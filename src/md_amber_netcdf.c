#include <md_amber_netcdf.h>

#include <md_system.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_common.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_platform.h>
#include <core/md_str.h>
#include <core/md_unit.h>

#include <netcdf.h>

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#define AMBER_NETCDF_PATH_CAP 4096

#if MD_PLATFORM_WINDOWS
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
static SRWLOCK netcdf_lock = SRWLOCK_INIT;
static void lock_netcdf(void) { AcquireSRWLockExclusive(&netcdf_lock); }
static void unlock_netcdf(void) { ReleaseSRWLockExclusive(&netcdf_lock); }
#elif MD_PLATFORM_UNIX
#include <pthread.h>
static pthread_mutex_t netcdf_lock = PTHREAD_MUTEX_INITIALIZER;
static void lock_netcdf(void) { pthread_mutex_lock(&netcdf_lock); }
static void unlock_netcdf(void) { pthread_mutex_unlock(&netcdf_lock); }
#else
#error "Amber NetCDF reader has no lock for this platform"
#endif

typedef struct amber_netcdf_source_t {
    int coordinate_var;
    bool has_fill_value;
    bool has_missing_value;
    float fill_value;
    float missing_value;
    size_t path_len;
    char path[AMBER_NETCDF_PATH_CAP];
} amber_netcdf_source_t;

typedef struct amber_netcdf_layout_t {
    size_t num_frames;
    size_t num_atoms;
    int coordinate_var;
    int time_var;
    int cell_lengths_var;
    int cell_angles_var;
    bool has_time_fill_value;
    bool has_time_missing_value;
    double time_fill_value;
    double time_missing_value;
    bool has_coordinate_fill_value;
    bool has_coordinate_missing_value;
    double coordinate_fill_value;
    double coordinate_missing_value;
    bool has_cell_lengths_fill_value;
    bool has_cell_lengths_missing_value;
    double cell_lengths_fill_value;
    double cell_lengths_missing_value;
    bool has_cell_angles_fill_value;
    bool has_cell_angles_missing_value;
    double cell_angles_fill_value;
    double cell_angles_missing_value;
    bool has_time;
    bool has_cell;
} amber_netcdf_layout_t;

static bool nc_ok(int status, const char* operation, str_t path) {
    if (status == NC_NOERR) {
        return true;
    }
    MD_LOG_ERROR("Amber NetCDF: %s failed for '" STR_FMT "': %s", operation, STR_ARG(path), nc_strerror(status));
    return false;
}

static bool nc_get_text_attribute(int ncid, int varid, const char* name, char* dst, size_t cap) {
    nc_type type;
    size_t len;
    if (nc_inq_att(ncid, varid, name, &type, &len) != NC_NOERR || type != NC_CHAR || len == 0 || len >= cap) {
        return false;
    }
    if (nc_get_att_text(ncid, varid, name, dst) != NC_NOERR) {
        return false;
    }
    dst[len] = '\0';
    return true;
}

static bool nc_variable_units_are(int ncid, int varid, const char* expected) {
    char units[64];
    return nc_get_text_attribute(ncid, varid, "units", units, sizeof(units)) && strcmp(units, expected) == 0;
}

static bool nc_variable_shape(int ncid, int varid, int rank, int dimids[]) {
    int actual_rank;
    if (nc_inq_var(ncid, varid, NULL, NULL, &actual_rank, dimids, NULL) != NC_NOERR) {
        return false;
    }
    return actual_rank == rank;
}

static bool nc_variable_is_float(int ncid, int varid) {
    nc_type type;
    return nc_inq_vartype(ncid, varid, &type) == NC_NOERR && (type == NC_FLOAT || type == NC_DOUBLE);
}

static bool nc_optional_numeric_attribute(int ncid, int varid, const char* name, bool* present, double* value) {
    nc_type type;
    size_t len;
    const int status = nc_inq_att(ncid, varid, name, &type, &len);
    if (status == NC_ENOTATT) {
        *present = false;
        return true;
    }
    if (status != NC_NOERR || len != 1 || type == NC_CHAR || type == NC_STRING) {
        return false;
    }
    *present = true;
    return nc_get_att_double(ncid, varid, name, value) == NC_NOERR && isfinite(*value);
}

static bool is_missing_double(double value, bool has_fill, double fill, bool has_missing, double missing) {
    return !isfinite(value) || value == NC_FILL_DOUBLE || value == (double)NC_FILL_FLOAT ||
        (has_fill && value == fill) || (has_missing && value == missing);
}

static bool nc_validate_spatial_labels(int ncid, int spatial_dim, str_t path) {
    int varid;
    const int status = nc_inq_varid(ncid, "spatial", &varid);
    if (status == NC_ENOTVAR) {
        return true;
    }
    int dims[NC_MAX_VAR_DIMS];
    nc_type type;
    if (status != NC_NOERR || nc_inq_vartype(ncid, varid, &type) != NC_NOERR || type != NC_CHAR ||
        !nc_variable_shape(ncid, varid, 1, dims) || dims[0] != spatial_dim) {
        MD_LOG_ERROR("Amber NetCDF: invalid spatial labels in '" STR_FMT "'", STR_ARG(path));
        return false;
    }
    char labels[3];
    if (nc_get_var_text(ncid, varid, labels) != NC_NOERR ||
        labels[0] != 'x' || labels[1] != 'y' || labels[2] != 'z') {
        MD_LOG_ERROR("Amber NetCDF: spatial labels must be x, y, z in '" STR_FMT "'", STR_ARG(path));
        return false;
    }
    return true;
}

static bool amber_netcdf_inspect(int ncid, str_t path, amber_netcdf_layout_t* layout) {
    int format;
    if (!nc_ok(nc_inq_format(ncid, &format), "format query", path)) {
        return false;
    }
    if (format != NC_FORMAT_CLASSIC && format != NC_FORMAT_64BIT_OFFSET) {
        MD_LOG_ERROR("Amber NetCDF: only classic NetCDF and 64-bit-offset NetCDF are supported in '" STR_FMT "'", STR_ARG(path));
        return false;
    }

    char conventions[64], convention_version[64];
    if (!nc_get_text_attribute(ncid, NC_GLOBAL, "Conventions", conventions, sizeof(conventions)) ||
        strcmp(conventions, "AMBER") != 0 ||
        !nc_get_text_attribute(ncid, NC_GLOBAL, "ConventionVersion", convention_version, sizeof(convention_version)) ||
        strcmp(convention_version, "1.0") != 0) {
        MD_LOG_ERROR("Amber NetCDF: '" STR_FMT "' does not declare the AMBER convention", STR_ARG(path));
        return false;
    }

    int frame_dim, atom_dim, spatial_dim;
    size_t frames, atoms, spatial;
    if (!nc_ok(nc_inq_dimid(ncid, "frame", &frame_dim), "frame dimension lookup", path) ||
        !nc_ok(nc_inq_dimid(ncid, "atom", &atom_dim), "atom dimension lookup", path) ||
        !nc_ok(nc_inq_dimid(ncid, "spatial", &spatial_dim), "spatial dimension lookup", path) ||
        !nc_ok(nc_inq_dimlen(ncid, frame_dim, &frames), "frame dimension query", path) ||
        !nc_ok(nc_inq_dimlen(ncid, atom_dim, &atoms), "atom dimension query", path) ||
        !nc_ok(nc_inq_dimlen(ncid, spatial_dim, &spatial), "spatial dimension query", path)) {
        return false;
    }
    if (!frames || !atoms || atoms > UINT32_MAX || frames > UINT32_MAX ||
        atoms > SIZE_MAX / 3 || frames > SIZE_MAX / (9 * sizeof(float)) || spatial != 3) {
        MD_LOG_ERROR("Amber NetCDF: invalid frame/atom/spatial dimensions in '" STR_FMT "'", STR_ARG(path));
        return false;
    }
    if (!nc_validate_spatial_labels(ncid, spatial_dim, path)) {
        return false;
    }

    int coordinate_var;
    int coordinate_dims[NC_MAX_VAR_DIMS];
    if (!nc_ok(nc_inq_varid(ncid, "coordinates", &coordinate_var), "coordinates variable lookup", path) ||
        !nc_variable_shape(ncid, coordinate_var, 3, coordinate_dims) ||
        !nc_variable_is_float(ncid, coordinate_var) ||
        coordinate_dims[0] != frame_dim || coordinate_dims[1] != atom_dim || coordinate_dims[2] != spatial_dim ||
        !nc_variable_units_are(ncid, coordinate_var, "angstrom")) {
        MD_LOG_ERROR("Amber NetCDF: coordinates must have shape {frame, atom, spatial} and units 'angstrom' in '" STR_FMT "'", STR_ARG(path));
        return false;
    }
    double coordinate_fill_value = 0.0, coordinate_missing_value = 0.0;
    bool has_coordinate_fill_value, has_coordinate_missing_value;
    if (!nc_optional_numeric_attribute(ncid, coordinate_var, "_FillValue", &has_coordinate_fill_value, &coordinate_fill_value) ||
        !nc_optional_numeric_attribute(ncid, coordinate_var, "missing_value", &has_coordinate_missing_value, &coordinate_missing_value)) {
        MD_LOG_ERROR("Amber NetCDF: invalid coordinate fill/missing-value metadata in '" STR_FMT "'", STR_ARG(path));
        return false;
    }

    int time_var = -1;
    bool has_time = nc_inq_varid(ncid, "time", &time_var) == NC_NOERR;
    if (has_time) {
        int dims[NC_MAX_VAR_DIMS];
        size_t length;
        if (!nc_variable_shape(ncid, time_var, 1, dims) || !nc_variable_is_float(ncid, time_var) || dims[0] != frame_dim ||
            nc_inq_dimlen(ncid, dims[0], &length) != NC_NOERR || length != frames ||
            !nc_variable_units_are(ncid, time_var, "picosecond")) {
            MD_LOG_ERROR("Amber NetCDF: time must have one value per frame and units 'picosecond' in '" STR_FMT "'", STR_ARG(path));
            return false;
        }
    }
    double time_fill_value = 0.0, time_missing_value = 0.0;
    bool has_time_fill_value = false, has_time_missing_value = false;
    if (has_time &&
        (!nc_optional_numeric_attribute(ncid, time_var, "_FillValue", &has_time_fill_value, &time_fill_value) ||
         !nc_optional_numeric_attribute(ncid, time_var, "missing_value", &has_time_missing_value, &time_missing_value))) {
        MD_LOG_ERROR("Amber NetCDF: invalid time fill/missing-value metadata in '" STR_FMT "'", STR_ARG(path));
        return false;
    }

    int lengths_var = -1, angles_var = -1;
    const bool has_lengths = nc_inq_varid(ncid, "cell_lengths", &lengths_var) == NC_NOERR;
    const bool has_angles = nc_inq_varid(ncid, "cell_angles", &angles_var) == NC_NOERR;
    if (has_lengths != has_angles) {
        MD_LOG_ERROR("Amber NetCDF: cell_lengths and cell_angles must either both be present or both be absent in '" STR_FMT "'", STR_ARG(path));
        return false;
    }
    bool has_lengths_fill_value = false, has_lengths_missing_value = false;
    bool has_angles_fill_value = false, has_angles_missing_value = false;
    double lengths_fill_value = 0.0, lengths_missing_value = 0.0;
    double angles_fill_value = 0.0, angles_missing_value = 0.0;
    if (has_lengths) {
        int cell_spatial_dim, cell_angular_dim, lengths_dims[NC_MAX_VAR_DIMS], angles_dims[NC_MAX_VAR_DIMS];
        if (nc_inq_dimid(ncid, "cell_spatial", &cell_spatial_dim) != NC_NOERR ||
            nc_inq_dimid(ncid, "cell_angular", &cell_angular_dim) != NC_NOERR ||
            !nc_variable_shape(ncid, lengths_var, 2, lengths_dims) ||
            !nc_variable_shape(ncid, angles_var, 2, angles_dims) ||
            !nc_variable_is_float(ncid, lengths_var) || !nc_variable_is_float(ncid, angles_var) ||
            lengths_dims[0] != frame_dim || angles_dims[0] != frame_dim ||
            lengths_dims[1] != cell_spatial_dim || angles_dims[1] != cell_angular_dim) {
            MD_LOG_ERROR("Amber NetCDF: cell variables have invalid dimensions in '" STR_FMT "'", STR_ARG(path));
            return false;
        }
        size_t cell_spatial, cell_angular;
        if (nc_inq_dimlen(ncid, cell_spatial_dim, &cell_spatial) != NC_NOERR ||
            nc_inq_dimlen(ncid, cell_angular_dim, &cell_angular) != NC_NOERR ||
            cell_spatial != 3 || cell_angular != 3 ||
            !nc_variable_units_are(ncid, lengths_var, "angstrom") ||
            !nc_variable_units_are(ncid, angles_var, "degree")) {
            MD_LOG_ERROR("Amber NetCDF: cell variables require 3 spatial/angle values with angstrom/degree units in '" STR_FMT "'", STR_ARG(path));
            return false;
        }
        if (!nc_optional_numeric_attribute(ncid, lengths_var, "_FillValue", &has_lengths_fill_value, &lengths_fill_value) ||
            !nc_optional_numeric_attribute(ncid, lengths_var, "missing_value", &has_lengths_missing_value, &lengths_missing_value) ||
            !nc_optional_numeric_attribute(ncid, angles_var, "_FillValue", &has_angles_fill_value, &angles_fill_value) ||
            !nc_optional_numeric_attribute(ncid, angles_var, "missing_value", &has_angles_missing_value, &angles_missing_value)) {
            MD_LOG_ERROR("Amber NetCDF: invalid cell fill/missing-value metadata in '" STR_FMT "'", STR_ARG(path));
            return false;
        }
    }

    *layout = (amber_netcdf_layout_t){
        .num_frames = frames,
        .num_atoms = atoms,
        .coordinate_var = coordinate_var,
        .time_var = time_var,
        .cell_lengths_var = lengths_var,
        .cell_angles_var = angles_var,
        .has_time_fill_value = has_time_fill_value,
        .has_time_missing_value = has_time_missing_value,
        .time_fill_value = time_fill_value,
        .time_missing_value = time_missing_value,
        .has_coordinate_fill_value = has_coordinate_fill_value,
        .has_coordinate_missing_value = has_coordinate_missing_value,
        .coordinate_fill_value = coordinate_fill_value,
        .coordinate_missing_value = coordinate_missing_value,
        .has_cell_lengths_fill_value = has_lengths_fill_value,
        .has_cell_lengths_missing_value = has_lengths_missing_value,
        .cell_lengths_fill_value = lengths_fill_value,
        .cell_lengths_missing_value = lengths_missing_value,
        .has_cell_angles_fill_value = has_angles_fill_value,
        .has_cell_angles_missing_value = has_angles_missing_value,
        .cell_angles_fill_value = angles_fill_value,
        .cell_angles_missing_value = angles_missing_value,
        .has_time = has_time,
        .has_cell = has_lengths,
    };
    return true;
}

bool md_amber_netcdf_file_is_trajectory(str_t filename) {
    char path[AMBER_NETCDF_PATH_CAP];
    if (str_empty(filename) || filename.len >= AMBER_NETCDF_PATH_CAP) {
        return false;
    }
    memcpy(path, filename.ptr, filename.len);
    path[filename.len] = '\0';

    lock_netcdf();
    int ncid = -1;
    amber_netcdf_layout_t layout = {0};
    const bool ok = nc_open(path, NC_NOWRITE, &ncid) == NC_NOERR &&
        amber_netcdf_inspect(ncid, filename, &layout);
    if (ncid >= 0) {
        nc_close(ncid);
    }
    unlock_netcdf();
    return ok;
}

static bool cell_from_lengths_angles(float out[9], const double lengths[3], const double angles[3]) {
    const double pi = 3.14159265358979323846;
    for (int i = 0; i < 3; ++i) {
        if (!isfinite(lengths[i]) || lengths[i] <= 0.0 ||
            !isfinite(angles[i]) || angles[i] <= 0.0 || angles[i] >= 180.0) {
            return false;
        }
    }

    const double alpha = angles[0] * pi / 180.0;
    const double beta  = angles[1] * pi / 180.0;
    const double gamma = angles[2] * pi / 180.0;
    const double sin_gamma = sin(gamma);
    if (fabs(sin_gamma) < 1.0e-12) {
        return false;
    }

    const double c_x = lengths[2] * cos(beta);
    const double c_y = lengths[2] * (cos(alpha) - cos(beta) * cos(gamma)) / sin_gamma;
    double c_z_sq = lengths[2] * lengths[2] - c_x * c_x - c_y * c_y;
    if (c_z_sq < -1.0e-8 * lengths[2] * lengths[2]) {
        return false;
    }
    if (c_z_sq < 0.0) {
        c_z_sq = 0.0;
    }
    out[0] = (float)lengths[0]; out[1] = 0.0f; out[2] = 0.0f;
    out[3] = (float)(lengths[1] * cos(gamma)); out[4] = (float)(lengths[1] * sin_gamma); out[5] = 0.0f;
    out[6] = (float)c_x; out[7] = (float)c_y; out[8] = (float)sqrt(c_z_sq);
    for (int i = 0; i < 9; ++i) {
        if (!isfinite(out[i])) {
            return false;
        }
    }
    return out[8] > 0.0f;
}

static size_t amber_netcdf_position_provider(void* dst, size_t cap, const md_attribute_t* attr,
                                             const md_attribute_slice_t* slice, void* user_data,
                                             md_attribute_io_t* io) {
    (void)io;
    const amber_netcdf_source_t* source = (const amber_netcdf_source_t*)user_data;
    if (!source || !slice || slice->num_idx == 0 || slice->num_idx > 2) {
        return 0;
    }
    const size_t frame = slice->idx[0];
    const size_t num_atoms = attr->format.shape[1];
    size_t first_atom = 0, atom_count = num_atoms;
    if (slice->num_idx == 2) {
        first_atom = slice->idx[1];
        atom_count = 1;
    }
    if (frame >= attr->format.shape[0] || first_atom >= num_atoms || cap != atom_count * 3) {
        return 0;
    }

    const str_t source_path = {source->path, source->path_len};
    int ncid = -1;
    lock_netcdf();
    int status = nc_open(source->path, NC_NOWRITE, &ncid);
    if (status == NC_NOERR) {
        const size_t start[3] = { frame, first_atom, 0 };
        const size_t count[3] = { 1, atom_count, 3 };
        status = nc_get_vara_float(ncid, source->coordinate_var, start, count, (float*)dst);
        const int close_status = nc_close(ncid);
        if (status == NC_NOERR) {
            status = close_status;
        }
    }
    char error[128] = {0};
    if (status != NC_NOERR) {
        snprintf(error, sizeof(error), "%s", nc_strerror(status));
    }
    unlock_netcdf();
    if (status != NC_NOERR) {
        MD_LOG_ERROR("Amber NetCDF: failed reading frame %zu from '" STR_FMT "': %s", frame,
                     STR_ARG(source_path), error);
        return 0;
    }
    const float default_fill = NC_FILL_FLOAT;
    const float double_fill = (float)NC_FILL_DOUBLE;
    const float custom_fill = (float)source->fill_value;
    const float custom_missing = (float)source->missing_value;
    for (size_t i = 0; i < cap; ++i) {
        const float value = ((float*)dst)[i];
        if (!isfinite(value) || value == default_fill || value == double_fill ||
            (source->has_fill_value && value == custom_fill) ||
            (source->has_missing_value && value == custom_missing)) {
            MD_LOG_ERROR("Amber NetCDF: non-finite coordinate in frame %zu of '" STR_FMT "'", frame,
                         STR_ARG(source_path));
            return 0;
        }
    }
    return cap;
}

bool md_amber_netcdf_system_publish_run(md_system_t* sys, str_t filename, str_t run, uint32_t flags) {
    ASSERT(sys);
    // NetCDF files already index frames; this reader does not create sidecar caches.
    (void)flags;
    if (str_empty(run) || str_empty(filename) || filename.len >= AMBER_NETCDF_PATH_CAP) {
        MD_LOG_ERROR("Amber NetCDF: invalid run name or path");
        return false;
    }
    char path_buf[AMBER_NETCDF_PATH_CAP];
    memcpy(path_buf, filename.ptr, filename.len);
    path_buf[filename.len] = '\0';
    const str_t path = {path_buf, filename.len};

    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    bool result = false;
    int ncid = -1;
    amber_netcdf_layout_t layout = {0};

    lock_netcdf();
    int status = nc_open(path_buf, NC_NOWRITE, &ncid);
    if (status != NC_NOERR) {
        nc_ok(status, "open", path);
    }
    const bool inspected = status == NC_NOERR && amber_netcdf_inspect(ncid, path, &layout);
    bool data_ok = inspected;
    double* times = NULL;
    float* cells = NULL;
    if (inspected) {
        times = md_alloc(arena, layout.num_frames * sizeof(double));
        if (!times) {
            MD_LOG_ERROR("Amber NetCDF: failed to allocate frame times");
            data_ok = false;
        } else if (layout.has_time) {
            status = nc_get_var_double(ncid, layout.time_var, times);
            if (status != NC_NOERR) {
                nc_ok(status, "time read", path);
                times = NULL;
                data_ok = false;
            }
        } else {
            for (size_t f = 0; f < layout.num_frames; ++f) {
                times[f] = (double)f;
            }
        }
        if (times && !layout.has_cell) {
            // Nonperiodic: publish an all-zero cell, which extraction reads as "no box"
            cells = md_alloc(arena, layout.num_frames * 9 * sizeof(float));
            if (!cells) {
                MD_LOG_ERROR("Amber NetCDF: failed to allocate cell data");
                data_ok = false;
            } else {
                memset(cells, 0, layout.num_frames * 9 * sizeof(float));
            }
        } else if (times) {
            double* lengths = md_alloc(arena, layout.num_frames * 3 * sizeof(double));
            double* angles = md_alloc(arena, layout.num_frames * 3 * sizeof(double));
            cells = md_alloc(arena, layout.num_frames * 9 * sizeof(float));
            if (!lengths || !angles || !cells) {
                MD_LOG_ERROR("Amber NetCDF: failed to allocate cell data");
                data_ok = false;
                cells = NULL;
            } else {
                status = nc_get_var_double(ncid, layout.cell_lengths_var, lengths);
                if (status == NC_NOERR) {
                    status = nc_get_var_double(ncid, layout.cell_angles_var, angles);
                }
                if (status != NC_NOERR) {
                    nc_ok(status, "cell read", path);
                    data_ok = false;
                    cells = NULL;
                } else {
                    for (size_t f = 0; f < layout.num_frames; ++f) {
                        bool missing = false;
                        for (size_t i = 0; i < 3; ++i) {
                            missing |= is_missing_double(lengths[f * 3 + i], layout.has_cell_lengths_fill_value,
                                                         layout.cell_lengths_fill_value, layout.has_cell_lengths_missing_value,
                                                         layout.cell_lengths_missing_value);
                            missing |= is_missing_double(angles[f * 3 + i], layout.has_cell_angles_fill_value,
                                                         layout.cell_angles_fill_value, layout.has_cell_angles_missing_value,
                                                         layout.cell_angles_missing_value);
                        }
                        if (missing) {
                            MD_LOG_ERROR("Amber NetCDF: missing periodic-cell value in frame %zu of '" STR_FMT "'", f, STR_ARG(path));
                            data_ok = false;
                            cells = NULL;
                            break;
                        }
                        if (!cell_from_lengths_angles(cells + f * 9, lengths + f * 3, angles + f * 3)) {
                            MD_LOG_ERROR("Amber NetCDF: invalid periodic cell in frame %zu of '" STR_FMT "'", f, STR_ARG(path));
                            data_ok = false;
                            cells = NULL;
                            break;
                        }
                    }
                }
            }
        }
    }
    if (ncid >= 0) {
        const int close_status = nc_close(ncid);
        if (close_status != NC_NOERR) {
            nc_ok(close_status, "close", path);
            data_ok = false;
        }
    }
    unlock_netcdf();

    if (!data_ok || !times || !cells) {
        goto done;
    }
    for (size_t f = 0; f < layout.num_frames; ++f) {
        const double time = times[f];
        if (is_missing_double(time, layout.has_time_fill_value, layout.time_fill_value,
                              layout.has_time_missing_value, layout.time_missing_value) ||
            (f && time < times[f - 1])) {
            MD_LOG_ERROR("Amber NetCDF: frame times must be finite and nondecreasing in '" STR_FMT "'", STR_ARG(path));
            goto done;
        }
    }

    if (!sys->attributes.alloc) {
        sys->attributes.alloc = sys->alloc;
    }
    if (!sys->attributes.alloc) {
        MD_LOG_ERROR("Amber NetCDF: the system attribute allocator is not set");
        goto done;
    }
    amber_netcdf_source_t* source = md_attributes_alloc_user_data(&sys->attributes, sizeof(*source));
    if (!source) {
        MD_LOG_ERROR("Amber NetCDF: failed to allocate coordinate source");
        goto done;
    }
    source->coordinate_var = layout.coordinate_var;
    source->has_fill_value = layout.has_coordinate_fill_value;
    source->has_missing_value = layout.has_coordinate_missing_value;
    source->fill_value = (float)layout.coordinate_fill_value;
    source->missing_value = (float)layout.coordinate_missing_value;
    source->path_len = path.len;
    memcpy(source->path, path.ptr, path.len);
    source->path[path.len] = '\0';
    const md_attribute_virtual_t position = {
        .provider = amber_netcdf_position_provider,
        .user_data = source,
        .user_data_size = sizeof(*source),
    };
    int64_t* offsets = md_alloc(arena, layout.num_frames * sizeof(int64_t));
    int64_t* sizes = md_alloc(arena, layout.num_frames * sizeof(int64_t));
    if (!offsets || !sizes) {
        md_free(sys->attributes.alloc, source, sizeof(*source));
        MD_LOG_ERROR("Amber NetCDF: failed to allocate source index");
        goto done;
    }
    for (size_t f = 0; f < layout.num_frames; ++f) {
        offsets[f] = -1;
        sizes[f] = 0;
    }
    const md_run_desc_t desc = {
        .num_frames = layout.num_frames,
        .num_atoms = layout.num_atoms,
        .time = times,
        .time_unit = layout.has_time ? md_unit_picosecond() : md_unit_none(),
        .unitcell = cells,
        .source_path = path,
        .source_offset = offsets,
        .source_size = sizes,
        .position_virt = &position,
    };
    result = md_run_publish(sys, run, &desc);
    if (!result) {
        md_free(sys->attributes.alloc, source, sizeof(*source));
    }

done:
    md_arena_allocator_destroy(arena);
    return result;
}
