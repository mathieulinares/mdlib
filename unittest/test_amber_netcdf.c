#include "utest.h"

#include <md_amber_netcdf.h>
#include <md_system.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_unit.h>

#include <netcdf.h>

#include <math.h>
#include <stdio.h>
#include <string.h>

#define AMBER_TEST_FILE "md_unittest_amber_netcdf.nc"
#define AMBER_TEST_RUN  STR_LIT("run/amber")

static bool amber_test_file_create(const char* path, bool with_time, bool with_cell, const char* convention,
                                   bool decreasing_time, bool offset64, const char* spatial_labels,
                                   bool fill_coordinate, bool double_coordinates) {
    int ncid, frame_dim, atom_dim, spatial_dim, cell_spatial_dim = -1, cell_angular_dim = -1;
    int coordinates_var, spatial_var, time_var = -1, lengths_var = -1, angles_var = -1;
    int coordinate_dims[3], time_dims[1], length_dims[2], angle_dims[2];
    if (nc_create(path, NC_CLOBBER | (offset64 ? NC_64BIT_OFFSET : 0), &ncid) != NC_NOERR) {
        return false;
    }
    if (nc_put_att_text(ncid, NC_GLOBAL, "Conventions", strlen(convention), convention) != NC_NOERR ||
        nc_put_att_text(ncid, NC_GLOBAL, "ConventionVersion", 3, "1.0") != NC_NOERR ||
        nc_def_dim(ncid, "frame", 2, &frame_dim) != NC_NOERR ||
        nc_def_dim(ncid, "atom", 2, &atom_dim) != NC_NOERR ||
        nc_def_dim(ncid, "spatial", 3, &spatial_dim) != NC_NOERR) {
        nc_close(ncid);
        return false;
    }
    if (nc_def_var(ncid, "spatial", NC_CHAR, 1, &spatial_dim, &spatial_var) != NC_NOERR) {
        nc_close(ncid);
        return false;
    }
    coordinate_dims[0] = frame_dim;
    coordinate_dims[1] = atom_dim;
    coordinate_dims[2] = spatial_dim;
    if (nc_def_var(ncid, "coordinates", double_coordinates ? NC_DOUBLE : NC_FLOAT, 3, coordinate_dims, &coordinates_var) != NC_NOERR ||
        nc_put_att_text(ncid, coordinates_var, "units", 8, "angstrom") != NC_NOERR) {
        nc_close(ncid);
        return false;
    }
    if (with_time) {
        time_dims[0] = frame_dim;
        if (nc_def_var(ncid, "time", NC_DOUBLE, 1, time_dims, &time_var) != NC_NOERR ||
            nc_put_att_text(ncid, time_var, "units", 10, "picosecond") != NC_NOERR) {
            nc_close(ncid);
            return false;
        }
    }
    if (with_cell) {
        if (nc_def_dim(ncid, "cell_spatial", 3, &cell_spatial_dim) != NC_NOERR ||
            nc_def_dim(ncid, "cell_angular", 3, &cell_angular_dim) != NC_NOERR) {
            nc_close(ncid);
            return false;
        }
        length_dims[0] = frame_dim;
        length_dims[1] = cell_spatial_dim;
        angle_dims[0] = frame_dim;
        angle_dims[1] = cell_angular_dim;
        if (nc_def_var(ncid, "cell_lengths", NC_DOUBLE, 2, length_dims, &lengths_var) != NC_NOERR ||
            nc_put_att_text(ncid, lengths_var, "units", 8, "angstrom") != NC_NOERR ||
            nc_def_var(ncid, "cell_angles", NC_DOUBLE, 2, angle_dims, &angles_var) != NC_NOERR ||
            nc_put_att_text(ncid, angles_var, "units", 6, "degree") != NC_NOERR) {
            nc_close(ncid);
            return false;
        }
    }
    if (nc_enddef(ncid) != NC_NOERR) {
        nc_close(ncid);
        return false;
    }
    float coords[12] = {
        1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f,
        7.0f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f,
    };
    double double_coords[12] = {
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
        7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
    };
    if (fill_coordinate) {
        if (double_coordinates) {
            double_coords[0] = NC_FILL_DOUBLE;
        } else {
            coords[0] = NC_FILL_FLOAT;
        }
    }
    const double times[2] = { 0.0, decreasing_time ? -1.0 : 2.5 };
    const double lengths[6] = { 10.0, 20.0, 30.0, 11.0, 21.0, 31.0 };
    const double angles[6] = { 90.0, 90.0, 90.0, 90.0, 90.0, 60.0 };
    bool ok = nc_put_var_text(ncid, spatial_var, spatial_labels) == NC_NOERR &&
              (double_coordinates ? nc_put_var_double(ncid, coordinates_var, double_coords)
                                  : nc_put_var_float(ncid, coordinates_var, coords)) == NC_NOERR;
    if (with_time) ok = ok && nc_put_var_double(ncid, time_var, times) == NC_NOERR;
    if (with_cell) {
        ok = ok && nc_put_var_double(ncid, lengths_var, lengths) == NC_NOERR &&
             nc_put_var_double(ncid, angles_var, angles) == NC_NOERR;
    }
    return nc_close(ncid) == NC_NOERR && ok;
}

UTEST(amber_netcdf, publishes_streamed_coordinates_time_and_triclinic_cells) {
    remove(AMBER_TEST_FILE);
    ASSERT_TRUE(amber_test_file_create(AMBER_TEST_FILE, true, true, "AMBER", false, false, "xyz", false, false));
    const str_t path = STR_LIT(AMBER_TEST_FILE);
    EXPECT_TRUE(md_amber_netcdf_file_is_trajectory(path));

    md_allocator_i* alloc = md_vm_arena_create(MEGABYTES(16));
    md_system_t sys = {.alloc = alloc};
    ASSERT_TRUE(md_amber_netcdf_system_publish_run(&sys, path, AMBER_TEST_RUN, MD_RUN_FLAG_DISABLE_CACHE_WRITE));
    const md_attribute_t* time = md_attributes_find(&sys.attributes, STR_LIT("run/amber/time"));
    const md_attribute_t* cell = md_attributes_find(&sys.attributes, STR_LIT("run/amber/unitcell"));
    const md_attribute_t* pos = md_attributes_find(&sys.attributes, STR_LIT("run/amber/atom/position"));
    ASSERT_TRUE(time && cell && pos);
    ASSERT_EQ(2u, time->format.shape[0]);
    EXPECT_TRUE(md_unit_equal(time->unit, md_unit_picosecond()));
    EXPECT_NEAR(2.5, ((const double*)time->data)[1], 1.0e-12);
    ASSERT_EQ(2u, pos->format.shape[0]);
    ASSERT_EQ(2u, pos->format.shape[1]);

    float frame[6];
    ASSERT_EQ(6u, md_attribute_extract_f32(frame, ARRAY_SIZE(frame), pos, md_attribute_slice_1(1), md_unit_angstrom()));
    EXPECT_NEAR(7.0, frame[0], 1.0e-6);
    EXPECT_NEAR(12.0, frame[5], 1.0e-6);
    float atom[3];
    ASSERT_EQ(3u, md_attribute_extract_f32(atom, ARRAY_SIZE(atom), pos, md_attribute_slice_2(0, 1), md_unit_angstrom()));
    EXPECT_NEAR(4.0, atom[0], 1.0e-6);
    EXPECT_NEAR(6.0, atom[2], 1.0e-6);

    float box[9];
    ASSERT_EQ(9u, md_attribute_extract_f32(box, ARRAY_SIZE(box), cell, md_attribute_slice_1(1), md_unit_angstrom()));
    EXPECT_NEAR(11.0, box[0], 1.0e-6);
    EXPECT_NEAR(10.5, box[3], 1.0e-6);
    EXPECT_NEAR(18.1865, box[4], 1.0e-4);
    EXPECT_NEAR(31.0, box[8], 1.0e-6);

    const str_t paths[] = { STR_LIT("atom/position"), STR_LIT("unitcell") };
    md_system_extract_t* ex = md_system_extract_begin(&sys, AMBER_TEST_RUN, paths, ARRAY_SIZE(paths), alloc);
    ASSERT_TRUE(ex != NULL);
    md_system_state_t state = {.alloc = alloc};
    ASSERT_TRUE(md_system_state_init(&state, 2));
    ASSERT_TRUE(md_system_extract_frame(ex, 1, &state));
    EXPECT_NEAR(7.0, state.xyz[0].x, 1.0e-6);
    EXPECT_TRUE(md_unitcell_is_triclinic(&state.unitcell));
    EXPECT_NEAR(11.0, state.unitcell.x, 1.0e-6);
    EXPECT_NEAR(10.5, state.unitcell.xy, 1.0e-6);
    EXPECT_NEAR(18.1865, state.unitcell.y, 1.0e-4);
    md_system_extract_end(ex);
    md_system_state_free(&state);

    md_system_free(&sys);
    md_vm_arena_destroy(alloc);
    remove(AMBER_TEST_FILE);
}

UTEST(amber_netcdf, nonperiodic_frames_extract_and_clear_previous_cell) {
    remove(AMBER_TEST_FILE);
    ASSERT_TRUE(amber_test_file_create(AMBER_TEST_FILE, true, false, "AMBER", false, false, "xyz", false, false));
    md_allocator_i* alloc = md_vm_arena_create(MEGABYTES(16));
    md_system_t sys = {.alloc = alloc};
    ASSERT_TRUE(md_amber_netcdf_system_publish_run(&sys, STR_LIT(AMBER_TEST_FILE), AMBER_TEST_RUN, 0));
    const md_attribute_t* time = md_attributes_find(&sys.attributes, STR_LIT("run/amber/time"));
    ASSERT_TRUE(time != NULL);
    EXPECT_TRUE(md_unit_equal(time->unit, md_unit_picosecond()));
    EXPECT_NEAR(2.5, ((const double*)time->data)[1], 1.0e-12);

    const str_t paths[] = { STR_LIT("atom/position"), STR_LIT("unitcell") };
    md_system_extract_t* ex = md_system_extract_begin(&sys, AMBER_TEST_RUN, paths, ARRAY_SIZE(paths), alloc);
    ASSERT_TRUE(ex != NULL);
    md_system_state_t state = {.alloc = alloc};
    ASSERT_TRUE(md_system_state_init(&state, 2));
    for (int64_t f = 1; f >= 0; --f) {
        state.unitcell = md_unitcell_from_extent(10.0, 20.0, 30.0);
        ASSERT_TRUE(md_system_extract_frame(ex, f, &state));
        EXPECT_EQ((double)f, state.frame);
        EXPECT_NEAR(1.0 + 6.0 * f, state.xyz[0].x, 1.0e-6);
        EXPECT_NEAR(6.0 + 6.0 * f, state.xyz[1].z, 1.0e-6);
        EXPECT_FALSE(md_system_state_has_unitcell(&state));
        EXPECT_EQ(0.0, state.unitcell.x);
        EXPECT_EQ(0.0, state.unitcell.y);
        EXPECT_EQ(0.0, state.unitcell.z);
    }
    md_system_extract_end(ex);
    md_system_state_free(&state);
    md_system_free(&sys);
    md_vm_arena_destroy(alloc);
    remove(AMBER_TEST_FILE);
}

UTEST(amber_netcdf, missing_time_and_cells_use_unitless_ordinals_and_no_box) {
    remove(AMBER_TEST_FILE);
    ASSERT_TRUE(amber_test_file_create(AMBER_TEST_FILE, false, false, "AMBER", false, false, "xyz", false, false));
    md_allocator_i* alloc = md_vm_arena_create(MEGABYTES(16));
    md_system_t sys = {.alloc = alloc};
    ASSERT_TRUE(md_amber_netcdf_system_publish_run(&sys, STR_LIT(AMBER_TEST_FILE), AMBER_TEST_RUN, 0));
    const md_attribute_t* time = md_attributes_find(&sys.attributes, STR_LIT("run/amber/time"));
    ASSERT_TRUE(time != NULL);
    EXPECT_TRUE(md_unit_equal(time->unit, md_unit_none()));
    const md_attribute_t* cell = md_attributes_find(&sys.attributes, STR_LIT("run/amber/unitcell"));
    ASSERT_TRUE(cell != NULL);
    float box[9];
    ASSERT_EQ(9u, md_attribute_extract_f32(box, ARRAY_SIZE(box), cell, md_attribute_slice_1(1), md_unit_angstrom()));
    for (int i = 0; i < 9; ++i) EXPECT_EQ(0.0f, box[i]);
    EXPECT_EQ(1.0, ((const double*)time->data)[1]);
    md_system_free(&sys);
    md_vm_arena_destroy(alloc);
    remove(AMBER_TEST_FILE);
}

UTEST(amber_netcdf, rejects_wrong_convention_and_decreasing_time) {
    remove(AMBER_TEST_FILE);
    ASSERT_TRUE(amber_test_file_create(AMBER_TEST_FILE, true, false, "OTHER", false, false, "xyz", false, false));
    EXPECT_FALSE(md_amber_netcdf_file_is_trajectory(STR_LIT(AMBER_TEST_FILE)));
    md_allocator_i* alloc = md_vm_arena_create(MEGABYTES(16));
    md_system_t sys = {.alloc = alloc};
    EXPECT_FALSE(md_amber_netcdf_system_publish_run(&sys, STR_LIT(AMBER_TEST_FILE), AMBER_TEST_RUN, 0));
    EXPECT_TRUE(md_attributes_find(&sys.attributes, STR_LIT("run/amber/time")) == NULL);
    md_system_free(&sys);

    ASSERT_TRUE(amber_test_file_create(AMBER_TEST_FILE, true, false, "AMBER", true, false, "xyz", false, false));
    sys = (md_system_t){.alloc = alloc};
    EXPECT_TRUE(md_amber_netcdf_file_is_trajectory(STR_LIT(AMBER_TEST_FILE)));
    EXPECT_FALSE(md_amber_netcdf_system_publish_run(&sys, STR_LIT(AMBER_TEST_FILE), AMBER_TEST_RUN, 0));
    EXPECT_TRUE(md_attributes_find(&sys.attributes, STR_LIT("run/amber/time")) == NULL);
    md_system_free(&sys);
    md_vm_arena_destroy(alloc);
    remove(AMBER_TEST_FILE);
}

UTEST(amber_netcdf, accepts_64_bit_offset_and_rejects_bad_axis_labels_and_fill_values) {
    remove(AMBER_TEST_FILE);
    ASSERT_TRUE(amber_test_file_create(AMBER_TEST_FILE, true, false, "AMBER", false, true, "xyz", false, true));
    EXPECT_TRUE(md_amber_netcdf_file_is_trajectory(STR_LIT(AMBER_TEST_FILE)));
    md_allocator_i* alloc = md_vm_arena_create(MEGABYTES(16));
    md_system_t sys = {.alloc = alloc};
    ASSERT_TRUE(md_amber_netcdf_system_publish_run(&sys, STR_LIT(AMBER_TEST_FILE), AMBER_TEST_RUN, 0));
    const md_attribute_t* pos = md_attributes_find(&sys.attributes, STR_LIT("run/amber/atom/position"));
    ASSERT_TRUE(pos != NULL);
    float frame[6];
    ASSERT_EQ(6u, md_attribute_extract_f32(frame, ARRAY_SIZE(frame), pos, md_attribute_slice_1(1), md_unit_angstrom()));
    EXPECT_NEAR(7.0, frame[0], 1.0e-6);
    md_system_free(&sys);

    ASSERT_TRUE(amber_test_file_create(AMBER_TEST_FILE, true, false, "AMBER", false, false, "zyx", false, false));
    EXPECT_FALSE(md_amber_netcdf_file_is_trajectory(STR_LIT(AMBER_TEST_FILE)));

    ASSERT_TRUE(amber_test_file_create(AMBER_TEST_FILE, true, false, "AMBER", false, false, "xyz", true, false));
    sys = (md_system_t){.alloc = alloc};
    ASSERT_TRUE(md_amber_netcdf_system_publish_run(&sys, STR_LIT(AMBER_TEST_FILE), AMBER_TEST_RUN, 0));
    pos = md_attributes_find(&sys.attributes, STR_LIT("run/amber/atom/position"));
    ASSERT_TRUE(pos != NULL);
    float xyz[6];
    EXPECT_EQ(0u, md_attribute_extract_f32(xyz, ARRAY_SIZE(xyz), pos, md_attribute_slice_1(0), md_unit_angstrom()));
    md_system_free(&sys);
    md_vm_arena_destroy(alloc);
    remove(AMBER_TEST_FILE);
}
