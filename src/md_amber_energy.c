#include <md_amber_energy.h>

#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_parse.h>
#include <core/md_unit.h>
#include <md_system.h>

#include <math.h>
#include <string.h>

#define AMBER_MAX_TERMS_PER_RECORD 128
#define AMBER_MAX_TERM_NAME 128
#define AMBER_FILE_SNIFF_SIZE (1024 * 1024)

typedef struct amber_energy_value_t {
    str_t name;
    str_t slug;
    double value;
} amber_energy_value_t;

typedef struct amber_energy_term_t {
    str_t name;
    str_t slug;
    double* values;
    size_t present;
} amber_energy_term_t;

typedef struct amber_energy_segment_t {
    size_t first_row;
    size_t num_rows;
} amber_energy_segment_t;

typedef struct amber_energy_record_t {
    bool active;
    bool has_step;
    bool has_time;
    bool has_total;
    double step;
    double time;
    amber_energy_value_t values[AMBER_MAX_TERMS_PER_RECORD];
    size_t num_values;
} amber_energy_record_t;

typedef struct amber_energy_data_t {
    double* time;
    amber_energy_term_t* terms;
    amber_energy_segment_t* segments;
    size_t num_rows;
    double previous_step;
    bool has_previous_step;
    bool force_new_segment;
} amber_energy_data_t;

static str_t amber_join(char* buf, size_t cap, str_t a, str_t b);

static bool ascii_contains_ignore_case(const char* text, size_t text_len, const char* needle, size_t needle_len) {
    if (!text || !needle || needle_len == 0 || text_len < needle_len) return false;
    for (size_t i = 0; i <= text_len - needle_len; ++i) {
        size_t j = 0;
        while (j < needle_len && to_lower((unsigned char)text[i + j]) == to_lower((unsigned char)needle[j])) ++j;
        if (j == needle_len) return true;
    }
    return false;
}

typedef enum amber_energy_format_t {
    AMBER_ENERGY_FORMAT_NONE = 0,
    AMBER_ENERGY_FORMAT_MDOUT,
    AMBER_ENERGY_FORMAT_MDEN,
} amber_energy_format_t;

// An mden file opens with its "L0  Nsteps  time(ps)  Etot ..." column header; an mdout has
// NSTEP/TIME(PS)/Etot assignments somewhere in its first megabyte.
static amber_energy_format_t amber_energy_detect(str_t filename) {
    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_READ)) return AMBER_ENERGY_FORMAT_NONE;
    md_temp_scope_t temp = md_temp_begin();
    char* buf = md_temp_alloc(temp, AMBER_FILE_SNIFF_SIZE);
    if (!buf) {
        MD_LOG_ERROR("Amber energy: failed to allocate the file detection buffer");
        md_file_close(&file);
        md_temp_end(temp);
        return AMBER_ENERGY_FORMAT_NONE;
    }
    const size_t n = md_file_read(file, buf, AMBER_FILE_SNIFF_SIZE);
    md_file_close(&file);

    amber_energy_format_t format = AMBER_ENERGY_FORMAT_NONE;
    size_t i = 0;
    while (i < n && is_whitespace(buf[i])) ++i;
    if (i + 2 < n && buf[i] == 'L' && buf[i + 1] == '0' && is_whitespace(buf[i + 2])) {
        size_t end = i;
        while (end < n && buf[end] != '\n') ++end;
        if (ascii_contains_ignore_case(buf + i, end - i, "NSTEPS", 6) &&
            ascii_contains_ignore_case(buf + i, end - i, "TIME(PS)", 8) &&
            ascii_contains_ignore_case(buf + i, end - i, "ETOT", 4)) {
            format = AMBER_ENERGY_FORMAT_MDEN;
        }
    } else if (ascii_contains_ignore_case(buf, n, "NSTEP", 5) &&
               ascii_contains_ignore_case(buf, n, "TIME(PS)", 8) &&
               ascii_contains_ignore_case(buf, n, "ETOT", 4)) {
        format = AMBER_ENERGY_FORMAT_MDOUT;
    }
    md_temp_end(temp);
    return format;
}

bool md_amber_energy_file_is_mdout(str_t filename) {
    return amber_energy_detect(filename) == AMBER_ENERGY_FORMAT_MDOUT;
}

bool md_amber_energy_file_is_mden(str_t filename) {
    return amber_energy_detect(filename) == AMBER_ENERGY_FORMAT_MDEN;
}

bool md_amber_energy_file_is_energy(str_t filename) {
    return amber_energy_detect(filename) != AMBER_ENERGY_FORMAT_NONE;
}

static size_t amber_make_slug(char* dst, size_t cap, str_t name) {
    if (!cap) return 0;
    size_t len = 0;
    bool separator = false;
    for (size_t i = 0; i < name.len; ++i) {
        const unsigned char c = (unsigned char)name.ptr[i];
        if (is_alpha(c) || is_digit(c)) {
            if (separator && len) {
                if (len + 2 >= cap) break;
                dst[len++] = '_';
            } else if (len + 1 >= cap) {
                break;
            }
            dst[len++] = (char)to_lower(c);
            separator = false;
        } else {
            separator = true;
        }
    }
    while (len && dst[len - 1] == '_') --len;
    dst[len] = '\0';
    return len;
}

static bool amber_parse_number(str_t token, double* value) {
    token = str_trim(token);
    if (token.len && (token.ptr[token.len - 1] == ',' || token.ptr[token.len - 1] == ';')) {
        --token.len;
    }
    size_t n = md_parse_f64(value, token);
    return n == token.len && n > 0 && isfinite(*value);
}

static size_t amber_term_find(const amber_energy_data_t* data, str_t slug) {
    for (size_t i = 0; i < md_array_size(data->terms); ++i) {
        if (str_eq(data->terms[i].slug, slug)) return i;
    }
    return SIZE_MAX;
}

static bool amber_record_find(const amber_energy_record_t* record, str_t slug, double* value) {
    for (size_t i = 0; i < record->num_values; ++i) {
        if (str_eq(record->values[i].slug, slug)) {
            if (value) *value = record->values[i].value;
            return true;
        }
    }
    return false;
}

static bool amber_record_append_value(amber_energy_record_t* record, str_t name, double value, md_allocator_i* alloc) {
    char slug_buf[AMBER_MAX_TERM_NAME];
    const size_t slug_len = amber_make_slug(slug_buf, sizeof(slug_buf), name);
    if (!slug_len) return true;
    const str_t slug = {slug_buf, slug_len};

    if (str_eq(slug, STR_LIT("time_ps"))) {
        record->time = value;
        record->has_time = true;
        return true;
    }
    if (str_eq(slug, STR_LIT("nstep"))) {
        record->step = value;
        record->has_step = true;
    }
    if (str_eq(slug, STR_LIT("etot"))) record->has_total = true;

    if (record->num_values == AMBER_MAX_TERMS_PER_RECORD) {
        MD_LOG_ERROR("Amber mdout: a record contains more than %u terms", AMBER_MAX_TERMS_PER_RECORD);
        return false;
    }
    if (amber_record_find(record, slug, NULL)) {
        MD_LOG_ERROR("Amber mdout: duplicate term '" STR_FMT "' in an energy record", STR_ARG(name));
        return false;
    }
    record->values[record->num_values++] = (amber_energy_value_t) {
        .name = str_copy(name, alloc),
        .slug = str_copy(slug, alloc),
        .value = value,
    };
    return true;
}

static bool amber_record_finalize(amber_energy_data_t* data, amber_energy_record_t* record, md_allocator_i* alloc) {
    if (!record->active) {
        return true;
    }
    if (!record->has_step || !record->has_time || !record->has_total) {
        MD_LOG_ERROR("Amber mdout: incomplete energy record; NSTEP, TIME(PS), and Etot are required");
        return false;
    }
    if (data->num_rows == UINT32_MAX) {
        MD_LOG_ERROR("Amber mdout: too many energy records");
        return false;
    }

    const bool new_segment = data->num_rows == 0 || data->force_new_segment ||
        record->time < data->time[data->num_rows - 1] ||
        (data->has_previous_step && record->step < data->previous_step);
    if (new_segment) {
        if (data->num_rows) {
            MD_LOG_INFO("Amber energy WARNING: restart detected at record %zu (step %.0f, %.9g ps); "
                        "publishing each run segment separately without changing its timestamps",
                        data->num_rows, record->step, record->time);
        }
        md_array_push(data->segments, ((amber_energy_segment_t){data->num_rows, 0}), alloc);
    }

    for (size_t i = 0; i < record->num_values; ++i) {
        const amber_energy_value_t* value = &record->values[i];
        if (amber_term_find(data, value->slug) == SIZE_MAX) {
            amber_energy_term_t term = {
                .name = str_copy(value->name, alloc),
                .slug = str_copy(value->slug, alloc),
            };
            for (size_t r = 0; r < data->num_rows; ++r) {
                md_array_push(term.values, NAN, alloc);
            }
            md_array_push(data->terms, term, alloc);
        }
    }

    md_array_push(data->time, record->time, alloc);
    for (size_t i = 0; i < md_array_size(data->terms); ++i) {
        amber_energy_term_t* term = &data->terms[i];
        double value = NAN;
        if (amber_record_find(record, term->slug, &value)) {
            ++term->present;
        }
        md_array_push(term->values, value, alloc);
    }
    ++data->num_rows;
    ++data->segments[md_array_size(data->segments) - 1].num_rows;
    data->previous_step = record->step;
    data->has_previous_step = true;
    data->force_new_segment = false;
    record->active = false;
    record->has_step = false;
    record->has_time = false;
    record->has_total = false;
    record->step = 0.0;
    record->time = 0.0;
    record->num_values = 0;
    return true;
}

static bool amber_energy_validate_run_coverage(const md_system_t* sys, str_t run,
                                              const amber_energy_data_t* data,
                                              md_allocator_i* alloc) {
    char path_buf[512];
    const str_t run_time_path = amber_join(path_buf, sizeof(path_buf), run, STR_LIT("time"));
    const md_attribute_t* run_axis = md_attributes_find(&sys->attributes, run_time_path);
    if (!run_axis || md_attributes_axis(&sys->attributes, run_axis) != run_axis) {
        return true;
    }

    md_attributes_t probe = {.alloc = alloc};
    const size_t num_segments = md_array_size(data->segments);
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_attribute_id_t* energy_axes = md_temp_alloc_array(temp, md_attribute_id_t, num_segments);
    if (!energy_axes) {
        md_attributes_free(&probe);
        md_temp_end(temp);
        MD_LOG_ERROR("Amber energy: failed to allocate segment time axes");
        return false;
    }

    for (size_t s = 0; s < num_segments; ++s) {
        const amber_energy_segment_t* segment = &data->segments[s];
        char axis_path[64];
        md_snprintf(axis_path, sizeof(axis_path), "segment_%zu/time", s + 1);
        const md_attribute_id_t axis_id = md_attributes_create(&probe, &(md_attribute_desc_t){
            .path = str_from_cstr(axis_path),
            .format = {.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = {(uint32_t)segment->num_rows}},
            .flags = MD_ATTRIBUTE_FLAG_TEMPORAL,
            .unit = md_unit_picosecond(),
            .label = STR_LIT("Time"),
            .data = data->time + segment->first_row,
            .byte_size = segment->num_rows * sizeof(double),
        });
        energy_axes[s] = axis_id;
        if (axis_id == MD_ATTRIBUTE_INVALID) {
            md_attributes_free(&probe);
            md_temp_end(temp);
            MD_LOG_ERROR("Amber energy: failed to create time axis for segment %zu", s + 1);
            return false;
        }
    }

    for (size_t f = 0; f < run_axis->format.shape[0]; ++f) {
        bool mapped = false;
        for (size_t s = 0; s < num_segments; ++s) {
            size_t row;
            const md_attribute_t* axis = md_attributes_get(&probe, energy_axes[s]);
            if (axis && md_attribute_axis_map(&row, run_axis, f, axis)) {
                mapped = true;
                break;
            }
        }
        if (!mapped) {
            double time = 0.0;
            md_attribute_extract_f64(&time, 1, run_axis, md_attribute_slice_1((uint32_t)f), md_unit_picosecond());
            MD_LOG_ERROR("Amber energy: frame %zu of '" STR_FMT "' (%g ps) has no matching energy record; the file may belong to another run",
                         f, STR_ARG(run), time);
            md_attributes_free(&probe);
            md_temp_end(temp);
            return false;
        }
    }
    md_attributes_free(&probe);
    md_temp_end(temp);
    return true;
}

static bool amber_record_set(amber_energy_data_t* data, amber_energy_record_t* record, str_t name, double value, md_allocator_i* alloc) {
    char slug_buf[AMBER_MAX_TERM_NAME];
    const size_t slug_len = amber_make_slug(slug_buf, sizeof(slug_buf), name);
    const str_t slug = {slug_buf, slug_len};
    if (str_eq(slug, STR_LIT("nstep"))) {
        if (!amber_record_finalize(data, record, alloc)) return false;
        record->active = true;
    } else if (!record->active) {
        return true;
    }
    return amber_record_append_value(record, name, value, alloc);
}

static bool amber_append_key(char* key, size_t cap, size_t* len, str_t token) {
    if (!token.len) return true;
    if (*len && *len + 1 < cap) key[(*len)++] = ' ';
    if (*len + token.len >= cap) return false;
    MEMCPY(key + *len, token.ptr, token.len);
    *len += token.len;
    key[*len] = '\0';
    return true;
}

static bool amber_parse_line(amber_energy_data_t* data, amber_energy_record_t* record, str_t line, md_allocator_i* alloc) {
    char key[256] = {0};
    size_t key_len = 0;
    size_t i = 0;
    while (i < line.len) {
        while (i < line.len && is_whitespace(line.ptr[i])) ++i;
        if (i == line.len) break;
        const size_t beg = i;
        while (i < line.len && !is_whitespace(line.ptr[i])) ++i;
        str_t token = {line.ptr + beg, i - beg};
        size_t eq = token.len;
        for (size_t j = 0; j < token.len; ++j) {
            if (token.ptr[j] == '=') { eq = j; break; }
        }

        if (eq < token.len) {
            if (eq && !amber_append_key(key, sizeof(key), &key_len, (str_t){token.ptr, eq})) {
                key_len = 0;
                key[0] = '\0';
                continue;
            }
            str_t number = {token.ptr + eq + 1, token.len - eq - 1};
            if (!number.len) {
                while (i < line.len && is_whitespace(line.ptr[i])) ++i;
                const size_t num_beg = i;
                while (i < line.len && !is_whitespace(line.ptr[i])) ++i;
                number = (str_t){line.ptr + num_beg, i - num_beg};
            }
            double value;
            if (key_len && amber_parse_number(number, &value)) {
                bool has_alpha = false;
                for (size_t j = 0; j < key_len; ++j) has_alpha |= is_alpha((unsigned char)key[j]);
                if (has_alpha && !amber_record_set(data, record, (str_t){key, key_len}, value, alloc)) return false;
            }
            key_len = 0;
            key[0] = '\0';
        } else {
            double ignored;
            if (amber_parse_number(token, &ignored)) {
                key_len = 0;
                key[0] = '\0';
            } else if (!amber_append_key(key, sizeof(key), &key_len, token)) {
                key_len = 0;
                key[0] = '\0';
            }
        }
    }
    return true;
}

static md_unit_t amber_energy_unit(str_t slug) {
    if (str_eq(slug, STR_LIT("nstep"))) return md_unit_none();
    if (str_eq(slug, STR_LIT("temp_k")) || str_eq(slug, STR_LIT("temp")) ||
        str_eq(slug, STR_LIT("t_solute")) || str_eq(slug, STR_LIT("t_solv"))) {
        return md_unit_kelvin();
    }
    if (str_eq(slug, STR_LIT("press")) || str_eq(slug, STR_LIT("pressure")) ||
        str_eq(slug, STR_LIT("pres_x")) || str_eq(slug, STR_LIT("pres_y")) || str_eq(slug, STR_LIT("pres_z"))) {
        return md_unit_scl(md_unit_bar(), 1.01325);
    }
    if (str_eq(slug, STR_LIT("density"))) {
        return md_unit_scl(md_unit_div(md_unit_kilogram(), md_unit_pow(md_unit_meter(), 3)), 1000.0);
    }
    if (str_eq(slug, STR_LIT("volume"))) return md_unit_pow(md_unit_angstrom(), 3);
    if (str_eq(slug, STR_LIT("boxx")) || str_eq(slug, STR_LIT("boxy")) || str_eq(slug, STR_LIT("boxz"))) {
        return md_unit_angstrom();
    }
    // Dimensionless pressure scaling factors and dipole moments (Debye, which md_unit has no name for)
    if (str_eq(slug, STR_LIT("pres_scal_solu")) || str_eq(slug, STR_LIT("pres_scal_solv")) ||
        str_eq(slug, STR_LIT("av_permmoment")) || str_eq(slug, STR_LIT("av_indmoment")) || str_eq(slug, STR_LIT("av_totmoment"))) {
        return md_unit_none();
    }
    return md_unit_scl(md_unit_div(md_unit_calorie(), md_unit_mole()), 1000.0);
}

static str_t amber_join(char* buf, size_t cap, str_t a, str_t b) {
    const size_t len = a.len + 1 + b.len;
    if (len + 1 > cap) return (str_t){0};
    MEMCPY(buf, a.ptr, a.len);
    buf[a.len] = '/';
    MEMCPY(buf + a.len + 1, b.ptr, b.len);
    buf[len] = '\0';
    return (str_t){buf, len};
}

static bool amber_mdout_section_header(str_t line, str_t number, str_t name) {
    line = str_trim(line);
    if (!str_begins_with(line, number)) return false;
    return str_eq_ignore_case(str_trim(str_substr(line, number.len, SIZE_MAX)), name);
}

static bool amber_mdout_record_separator(str_t line) {
    line = str_trim(line);
    if (line.len < 5) return false;
    for (size_t i = 0; i < line.len; ++i) {
        if (line.ptr[i] != '-') return false;
    }
    return true;
}

static bool amber_parse_mdout(amber_energy_data_t* data, str_t filename, md_allocator_i* alloc) {
    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_READ)) {
        MD_LOG_ERROR("Amber mdout: failed to open '" STR_FMT "'", STR_ARG(filename));
        return false;
    }
    char buf[64 * 1024];
    md_buffered_reader_t reader = md_buffered_reader_from_file(buf, sizeof(buf), file);
    amber_energy_record_t record = {0};
    str_t line;
    bool ok = true;
    bool skipping_summary = false;
    bool after_summary = false;
    while (md_buffered_reader_extract_line(&line, &reader)) {
        if (ascii_contains_ignore_case(line.ptr, line.len, "A V E R A G E S", sizeof("A V E R A G E S") - 1) ||
            ascii_contains_ignore_case(line.ptr, line.len, "R M S  F L U C T U A T I O N S", sizeof("R M S  F L U C T U A T I O N S") - 1)) {
            if (!amber_record_finalize(data, &record, alloc)) { ok = false; break; }
            skipping_summary = true;
            after_summary = false;
            continue;
        }
        if (amber_mdout_section_header(line, STR_LIT("4."), STR_LIT("RESULTS"))) {
            if (!amber_record_finalize(data, &record, alloc)) { ok = false; break; }
            if (data->num_rows) data->force_new_segment = true;
            skipping_summary = false;
            after_summary = false;
            continue;
        }
        if (skipping_summary) {
            if (amber_mdout_section_header(line, STR_LIT("5."), STR_LIT("TIMINGS"))) {
                after_summary = true;
                continue;
            }
            const str_t trimmed = str_trim_beg(line);
            if (!after_summary || trimmed.len < 5 ||
                !str_eq_ignore_case(str_substr(trimmed, 0, 5), STR_LIT("NSTEP")) ||
                !ascii_contains_ignore_case(trimmed.ptr, trimmed.len, "TIME(PS)", 8)) {
                continue;
            }
            data->force_new_segment = data->num_rows != 0;
            skipping_summary = false;
        }
        if (record.active && record.has_total && amber_mdout_record_separator(line)) {
            if (!amber_record_finalize(data, &record, alloc)) { ok = false; break; }
            continue;
        }
        if (!amber_parse_line(data, &record, line, alloc)) {
            ok = false;
            break;
        }
    }
    if (ok) ok = amber_record_finalize(data, &record, alloc);
    md_file_close(&file);
    if (!ok) return false;
    if (data->num_rows == 0) {
        MD_LOG_ERROR("Amber mdout: no complete NSTEP/TIME(PS)/Etot records in '" STR_FMT "'", STR_ARG(filename));
        return false;
    }
    return true;
}

#define AMBER_MDEN_MAX_LINES 10
#define AMBER_MDEN_MAX_COLUMNS 32

// mden names its columns differently from mdout. The terms both files have are published under
// the mdout names, so attr("amber_energy/etot") or ".../eptot" reads the same quantity from either.
static str_t amber_mden_term_name(str_t name) {
    static const struct { const char* mden; const char* mdout; } map[] = {
        {"Nsteps", "NSTEP"},       {"time(ps)", "TIME(PS)"},  {"Etot", "Etot"},
        {"EKinetic", "EKtot"},     {"Temp", "TEMP(K)"},       {"Pressure", "PRESS"},
        {"volume", "VOLUME"},      {"EKComTot", "EKCMT"},     {"VIRIAL_tot", "VIRIAL"},
        {"E_pot", "EPtot"},        {"E_vdw", "VDWAALS"},      {"E_el", "EELEC"},
        {"E_hbon", "EHBOND"},      {"E_bon", "BOND"},         {"E_angle", "ANGLE"},
        {"E_dih", "DIHED"},        {"E_14vdw", "1-4 NB"},     {"E_14el", "1-4 EEL"},
        {"E_const", "RESTRAINT"},  {"E_pol", "EPOLZ"},        {"Density", "Density"},
        {"dV/dlambda", "DV/DL"},
    };
    for (size_t i = 0; i < ARRAY_SIZE(map); ++i) {
        if (str_eq_cstr(name, map[i].mden)) return str_from_cstr(map[i].mdout);
    }
    return name;
}

static size_t amber_split_tokens(str_t* tokens, size_t cap, str_t line) {
    size_t count = 0;
    size_t i = 0;
    while (i < line.len) {
        while (i < line.len && is_whitespace(line.ptr[i])) ++i;
        if (i == line.len) break;
        const size_t beg = i;
        while (i < line.len && !is_whitespace(line.ptr[i])) ++i;
        if (count == cap) return cap + 1;
        tokens[count++] = (str_t){line.ptr + beg, i - beg};
    }
    return count;
}

// An mden file is a header of "L<k>  name name ..." lines followed by records of the same lines
// with numbers in place of the names; record k-th line values belong to the k-th header line.
static bool amber_parse_mden(amber_energy_data_t* data, str_t filename, md_allocator_i* alloc) {
    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_READ)) {
        MD_LOG_ERROR("Amber mden: failed to open '" STR_FMT "'", STR_ARG(filename));
        return false;
    }
    char buf[64 * 1024];
    md_buffered_reader_t reader = md_buffered_reader_from_file(buf, sizeof(buf), file);

    str_t header[AMBER_MDEN_MAX_LINES][AMBER_MDEN_MAX_COLUMNS] = {0};
    size_t header_count[AMBER_MDEN_MAX_LINES] = {0};
    amber_energy_record_t record = {0};
    size_t line_number = 0;
    str_t line;
    bool ok = true;

    while (ok && md_buffered_reader_extract_line(&line, &reader)) {
        ++line_number;
        line = str_trim(line);
        if (!line.len) continue;

        str_t tokens[AMBER_MDEN_MAX_COLUMNS + 1];
        const size_t num_tokens = amber_split_tokens(tokens, ARRAY_SIZE(tokens), line);
        const str_t label = tokens[0];
        if (num_tokens < 2 || num_tokens > ARRAY_SIZE(tokens) || label.len != 2 || label.ptr[0] != 'L' || !is_digit(label.ptr[1])) {
            MD_LOG_ERROR("Amber mden: line %zu is not an 'L<k>' header or record line", line_number);
            ok = false;
            break;
        }
        const size_t k = (size_t)(label.ptr[1] - '0');
        const str_t* values = tokens + 1;
        const size_t num_values = num_tokens - 1;

        double value;
        if (!amber_parse_number(values[0], &value)) {
            // A header line. A new header (a second run appended to the file) closes the open record.
            if (!amber_record_finalize(data, &record, alloc)) { ok = false; break; }
            if (k == 0 && data->num_rows) data->force_new_segment = true;
            for (size_t c = 0; c < num_values; ++c) {
                header[k][c] = str_copy(values[c], alloc);
            }
            header_count[k] = num_values;
            continue;
        }

        if (k == 0) {
            if (!amber_record_finalize(data, &record, alloc)) { ok = false; break; }
            record.active = true;
        } else if (!record.active) {
            MD_LOG_ERROR("Amber mden: line %zu (L%zu) comes before the L0 line that starts its record", line_number, k);
            ok = false;
            break;
        }
        if (header_count[k] == 0) {
            MD_LOG_ERROR("Amber mden: line %zu (L%zu) has no column header", line_number, k);
            ok = false;
            break;
        }
        if (num_values != header_count[k]) {
            MD_LOG_ERROR("Amber mden: line %zu (L%zu) has %zu values where its header names %zu", line_number, k, num_values, header_count[k]);
            ok = false;
            break;
        }
        for (size_t c = 0; ok && c < num_values; ++c) {
            if (!amber_parse_number(values[c], &value)) {
                MD_LOG_ERROR("Amber mden: line %zu has an unreadable value '" STR_FMT "'", line_number, STR_ARG(values[c]));
                ok = false;
                break;
            }
            ok = amber_record_append_value(&record, amber_mden_term_name(header[k][c]), value, alloc);
        }
    }
    if (ok) ok = amber_record_finalize(data, &record, alloc);
    md_file_close(&file);
    if (!ok) return false;
    if (data->num_rows == 0) {
        MD_LOG_ERROR("Amber mden: no complete records in '" STR_FMT "'", STR_ARG(filename));
        return false;
    }
    return true;
}

static bool amber_parse_file(amber_energy_data_t* data, str_t filename, md_allocator_i* alloc) {
    switch (amber_energy_detect(filename)) {
    case AMBER_ENERGY_FORMAT_MDOUT: return amber_parse_mdout(data, filename, alloc);
    case AMBER_ENERGY_FORMAT_MDEN:  return amber_parse_mden(data, filename, alloc);
    default:
        MD_LOG_ERROR("Amber energy: '" STR_FMT "' is neither an Amber mdout nor an mden file", STR_ARG(filename));
        return false;
    }
}

bool md_amber_energy_system_supplement_from_file(md_system_t* sys, str_t filename, str_t run) {
    ASSERT(sys);
    if (str_empty(run)) {
        MD_LOG_ERROR("Amber energy: no run to publish the energies under");
        return false;
    }
    if (!sys->attributes.alloc) sys->attributes.alloc = sys->alloc;

    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);
    amber_energy_data_t data = {0};
    bool result = false;
    if (!amber_parse_file(&data, filename, alloc)) goto done;

    char group_buf[512];
    const str_t group = amber_join(group_buf, sizeof(group_buf), run, STR_LIT("amber_energy"));
    if (str_empty(group)) {
        MD_LOG_ERROR("Amber energy: run path '" STR_FMT "' is too long", STR_ARG(run));
        goto done;
    }

    if (!amber_energy_validate_run_coverage(sys, run, &data, alloc)) goto done;

    md_attributes_remove_prefix(&sys->attributes, group);
    const size_t num_segments = md_array_size(data.segments);
    bool published = true;
    if (num_segments > 1) {
        MD_LOG_INFO("Amber energy WARNING: preserving %zu restarted segments under '" STR_FMT
                    "/segment_1, /segment_2, ...; timestamps are unchanged",
                    num_segments, STR_ARG(group));
    }

    for (size_t s = 0; published && s < num_segments; ++s) {
        const amber_energy_segment_t* segment = &data.segments[s];
        char segment_group_buf[512];
        char time_path_buf[512];
        char term_path_buf[512];
        char segment_name[32];
        str_t segment_group = group;
        if (num_segments > 1) {
            md_snprintf(segment_name, sizeof(segment_name), "segment_%zu", s + 1);
            segment_group = amber_join(segment_group_buf, sizeof(segment_group_buf), group, str_from_cstr(segment_name));
        }
        const str_t time_path = amber_join(time_path_buf, sizeof(time_path_buf), segment_group, STR_LIT("time"));
        published = !str_empty(segment_group) && !str_empty(time_path) &&
            md_attributes_create(&sys->attributes, &(md_attribute_desc_t){
                .path = time_path,
                .format = {.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = {(uint32_t)segment->num_rows}},
                .flags = MD_ATTRIBUTE_FLAG_TEMPORAL,
                .unit = md_unit_picosecond(),
                .label = STR_LIT("Time"),
                .data = data.time + segment->first_row,
                .byte_size = segment->num_rows * sizeof(double),
            }) != MD_ATTRIBUTE_INVALID;

        md_array(str_t) taken = 0;
        md_array_push(taken, STR_LIT("time"), alloc);
        md_array_push(taken, STR_LIT("source"), alloc);
        for (size_t i = 0; published && i < md_array_size(data.terms); ++i) {
            amber_energy_term_t* term = &data.terms[i];
            bool complete = true;
            for (size_t r = 0; r < segment->num_rows; ++r) {
                complete &= isfinite(term->values[segment->first_row + r]);
            }
            if (!complete) {
                MD_LOG_INFO("Amber energy: omitting '" STR_FMT "' from segment %zu because it is not present in every record",
                            STR_ARG(term->name), s + 1);
                continue;
            }
            char slug_buf[AMBER_MAX_TERM_NAME];
            size_t slug_len = amber_make_slug(slug_buf, sizeof(slug_buf) - 8, term->name);
            if (!slug_len) slug_len = (size_t)md_snprintf(slug_buf, sizeof(slug_buf), "term%zu", i + 1);
            const size_t base_len = slug_len;
            for (int n = 2;; ++n) {
                bool clash = false;
                for (size_t k = 0; k < md_array_size(taken); ++k) {
                    if (str_eq(taken[k], (str_t){slug_buf, slug_len})) { clash = true; break; }
                }
                if (!clash) break;
                slug_len = base_len + (size_t)md_snprintf(slug_buf + base_len, sizeof(slug_buf) - base_len, "_%d", n);
            }
            md_array_push(taken, str_copy((str_t){slug_buf, slug_len}, alloc), alloc);
            const str_t path = amber_join(term_path_buf, sizeof(term_path_buf), segment_group, (str_t){slug_buf, slug_len});
            published = !str_empty(path) && md_attributes_create(&sys->attributes, &(md_attribute_desc_t){
                .path = path,
                .format = {.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = {(uint32_t)segment->num_rows}},
                .flags = MD_ATTRIBUTE_FLAG_TEMPORAL,
                .unit = amber_energy_unit(term->slug),
                .label = term->name,
                .data = term->values + segment->first_row,
                .byte_size = segment->num_rows * sizeof(double),
            }) != MD_ATTRIBUTE_INVALID;
        }
        md_array_free(taken, alloc);
    }

    if (published) {
        char canonical_buf[4096];
        char source_buf[512];
        const size_t path_len = md_path_write_canonical(canonical_buf, sizeof(canonical_buf), filename);
        const str_t source_path = amber_join(source_buf, sizeof(source_buf), group, STR_LIT("source"));
        str_t source = path_len ? str_copy((str_t){canonical_buf, path_len}, sys->attributes.alloc) : (str_t){0};
        published = path_len && !str_empty(source_path) && md_attributes_create(&sys->attributes, &(md_attribute_desc_t){
            .path = source_path,
            .format = {.type = MD_ATTRIBUTE_TYPE_STR, .components = 1, .rank = 0},
            .unit = md_unit_none(),
            .label = STR_LIT("Source"),
            .data = &source,
            .byte_size = sizeof(str_t),
        }) != MD_ATTRIBUTE_INVALID;
    }
    if (!published) {
        MD_LOG_ERROR("Amber energy: failed to publish energy series under '" STR_FMT "'", STR_ARG(group));
        md_attributes_remove_prefix(&sys->attributes, group);
        goto done;
    }
    result = true;

done:
    md_temp_end(temp);
    return result;
}
