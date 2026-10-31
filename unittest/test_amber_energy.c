#include "utest.h"

#include <md_amber_energy.h>
#include <md_system.h>
#include <core/md_allocator.h>

#include <stdio.h>
#include <string.h>

#define AMBER_ENERGY_FILE(name) STR_LIT(MD_UNITTEST_SOURCE_DIR "/amber_energy_data/" name)

static md_attribute_id_t amber_energy_publish_run_axis(md_system_t* sys, const double* times, uint32_t count) {
    return md_attributes_create(&sys->attributes, &(md_attribute_desc_t){
        .path = STR_LIT("run/test/time"),
        .format = {.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = {count}},
        .flags = MD_ATTRIBUTE_FLAG_TEMPORAL,
        .unit = md_unit_picosecond(),
        .data = times,
        .byte_size = count * sizeof(double),
    });
}

UTEST(amber_energy, recognizes_and_publishes_mdout_samples) {
    const str_t path = AMBER_ENERGY_FILE("amber-mdout.out");
    ASSERT_TRUE(md_amber_energy_file_is_mdout(path));

    md_allocator_i* alloc = md_get_heap_allocator();
    const double run_time[] = {0.0, 1.0, 2.0};
    md_system_t sys = {.alloc = alloc};
    sys.attributes.alloc = alloc;
    ASSERT_NE(MD_ATTRIBUTE_INVALID, amber_energy_publish_run_axis(&sys, run_time, 3));
    ASSERT_TRUE(md_amber_energy_system_supplement_from_file(&sys, path, STR_LIT("run/test")));

    const md_attribute_t* time = md_attributes_find(&sys.attributes, STR_LIT("run/test/amber_energy/time"));
    ASSERT_TRUE(time != NULL);
    EXPECT_EQ(3u, time->format.shape[0]);
    EXPECT_TRUE(md_unit_equal(time->unit, md_unit_picosecond()));

    const md_attribute_t* etot = md_attributes_find(&sys.attributes, STR_LIT("run/test/amber_energy/etot"));
    ASSERT_TRUE(etot != NULL);
    EXPECT_TRUE(str_eq_cstr(etot->label, "Etot"));
    md_unit_t energy_unit = md_unit_scl(md_unit_div(md_unit_calorie(), md_unit_mole()), 1000.0);
    EXPECT_TRUE(md_unit_equal(etot->unit, energy_unit));
    double values[3];
    ASSERT_EQ(3u, md_attribute_extract_f64(values, 3, etot, md_attribute_slice_all(), energy_unit));
    EXPECT_NEAR(-10.0, values[0], 1.0e-12);
    EXPECT_NEAR(-9.5, values[1], 1.0e-12);
    EXPECT_NEAR(-9.0, values[2], 1.0e-12);
    EXPECT_TRUE(md_attributes_find(&sys.attributes, STR_LIT("run/test/amber_energy/ke_trans")) == NULL);

    const md_attribute_t* pair_nb = md_attributes_find(&sys.attributes, STR_LIT("run/test/amber_energy/1_4_nb"));
    ASSERT_TRUE(pair_nb != NULL);
    ASSERT_EQ(3u, md_attribute_extract_f64(values, 3, pair_nb, md_attribute_slice_all(), energy_unit));
    EXPECT_NEAR(0.6, values[2], 1.0e-12);
    EXPECT_TRUE(md_attributes_find(&sys.attributes, STR_LIT("run/test/amber_energy/restraint")) == NULL);
    EXPECT_TRUE(md_attributes_find(&sys.attributes, STR_LIT("run/test/amber_energy/source")) != NULL);
    md_attributes_free(&sys.attributes);
}

UTEST(amber_energy, refuses_energy_output_that_does_not_cover_the_run) {
    const str_t path = AMBER_ENERGY_FILE("amber-mdout.out");
    md_allocator_i* alloc = md_get_heap_allocator();
    const double run_time[] = {0.0, 1.0, 3.0};
    md_system_t sys = {.alloc = alloc};
    sys.attributes.alloc = alloc;
    ASSERT_NE(MD_ATTRIBUTE_INVALID, amber_energy_publish_run_axis(&sys, run_time, 3));
    EXPECT_FALSE(md_amber_energy_system_supplement_from_file(&sys, path, STR_LIT("run/test")));
    EXPECT_TRUE(md_attributes_find(&sys.attributes, STR_LIT("run/test/amber_energy/time")) == NULL);
    md_attributes_free(&sys.attributes);
}

UTEST(amber_energy, preserves_restarted_runs_as_separate_segments) {
    const str_t path = AMBER_ENERGY_FILE("amber-mdout-restarted.out");
    const double run_time[] = {0.0, 1.0};
    md_system_t sys = {.alloc = md_get_heap_allocator()};
    sys.attributes.alloc = sys.alloc;
    ASSERT_NE(MD_ATTRIBUTE_INVALID, amber_energy_publish_run_axis(&sys, run_time, 2));
    ASSERT_TRUE(md_amber_energy_system_supplement_from_file(&sys, path, STR_LIT("run/test")));

    const md_attribute_t* first_time = md_attributes_find(&sys.attributes, STR_LIT("run/test/amber_energy/segment_1/time"));
    const md_attribute_t* second_time = md_attributes_find(&sys.attributes, STR_LIT("run/test/amber_energy/segment_2/time"));
    const md_attribute_t* first_etot = md_attributes_find(&sys.attributes, STR_LIT("run/test/amber_energy/segment_1/etot"));
    const md_attribute_t* second_etot = md_attributes_find(&sys.attributes, STR_LIT("run/test/amber_energy/segment_2/etot"));
    ASSERT_TRUE(first_time && second_time && first_etot && second_etot);
    EXPECT_TRUE(md_attributes_find(&sys.attributes, STR_LIT("run/test/amber_energy/time")) == NULL);
    EXPECT_EQ(2u, first_time->format.shape[0]);
    EXPECT_EQ(2u, second_time->format.shape[0]);

    const md_unit_t energy_unit = md_unit_scl(md_unit_div(md_unit_calorie(), md_unit_mole()), 1000.0);
    double values[2];
    ASSERT_EQ(2u, md_attribute_extract_f64(values, 2, first_etot, md_attribute_slice_all(), energy_unit));
    EXPECT_NEAR(-10.0, values[0], 1.0e-12);
    EXPECT_NEAR(-9.5, values[1], 1.0e-12);
    ASSERT_EQ(2u, md_attribute_extract_f64(values, 2, second_etot, md_attribute_slice_all(), energy_unit));
    EXPECT_NEAR(-20.0, values[0], 1.0e-12);
    EXPECT_NEAR(-19.5, values[1], 1.0e-12);
    EXPECT_TRUE(md_attributes_find(&sys.attributes, STR_LIT("run/test/amber_energy/source")) != NULL);
    md_attributes_free(&sys.attributes);
}

static void check_concatenated_energy_file(int* utest_result, str_t path) {
    const double run_time[] = {0.0, 1.0, 2.0};
    md_system_t sys = {.alloc = md_get_heap_allocator()};
    sys.attributes.alloc = sys.alloc;
    ASSERT_NE(MD_ATTRIBUTE_INVALID, amber_energy_publish_run_axis(&sys, run_time, 3));
    ASSERT_TRUE(md_amber_energy_system_supplement_from_file(&sys, path, STR_LIT("run/test")));

    const md_attribute_t* first = md_attributes_find(&sys.attributes, STR_LIT("run/test/amber_energy/segment_1/etot"));
    const md_attribute_t* second = md_attributes_find(&sys.attributes, STR_LIT("run/test/amber_energy/segment_2/etot"));
    ASSERT_TRUE(first && second);
    EXPECT_EQ(2u, first->format.shape[0]);
    EXPECT_EQ(2u, second->format.shape[0]);
    const md_unit_t unit = md_unit_scl(md_unit_div(md_unit_calorie(), md_unit_mole()), 1000.0);
    double values[2];
    ASSERT_EQ(2u, md_attribute_extract_f64(values, 2, first, md_attribute_slice_all(), unit));
    EXPECT_NEAR(-10.0, values[0], 1.0e-12);
    EXPECT_NEAR(-9.5, values[1], 1.0e-12);
    ASSERT_EQ(2u, md_attribute_extract_f64(values, 2, second, md_attribute_slice_all(), unit));
    EXPECT_NEAR(-20.0, values[0], 1.0e-12);
    EXPECT_NEAR(-19.5, values[1], 1.0e-12);
    md_attributes_free(&sys.attributes);
}

UTEST(amber_energy, reads_concatenated_mdout_after_summary_sections) {
    check_concatenated_energy_file(utest_result, AMBER_ENERGY_FILE("amber-mdout-concatenated.out"));
}

UTEST(amber_energy, resumes_mdout_at_records_after_timings) {
    check_concatenated_energy_file(utest_result, AMBER_ENERGY_FILE("amber-mdout-after-timings.out"));
}

UTEST(amber_energy, separates_mden_runs_with_repeated_boundary_times) {
    check_concatenated_energy_file(utest_result, AMBER_ENERGY_FILE("amber-mden-concatenated.mden"));
}

UTEST(amber_energy, validates_coverage_across_many_segments) {
    const char* filename = "md_unittest_amber_many_segments.out";
    FILE* file = fopen(filename, "wb");
    ASSERT_TRUE(file != NULL);
    for (int i = 0; i < 20; ++i) {
        fprintf(file, "NSTEP = 0 TIME(PS) = %d Etot = -10\n", i * 10);
        fprintf(file, "NSTEP = 1 TIME(PS) = %d Etot = -9\n", i * 10 + 1);
    }
    ASSERT_EQ(0, fclose(file));

    double run_time[20];
    for (size_t i = 0; i < 20; ++i) run_time[i] = (double)(i * 10);
    md_system_t sys = {.alloc = md_get_heap_allocator()};
    sys.attributes.alloc = sys.alloc;
    ASSERT_NE(MD_ATTRIBUTE_INVALID, amber_energy_publish_run_axis(&sys, run_time, 20));
    ASSERT_TRUE(md_amber_energy_system_supplement_from_file(&sys, str_from_cstr(filename), STR_LIT("run/test")));
    for (size_t i = 0; i < 20; ++i) {
        char path[128];
        md_snprintf(path, sizeof(path), "run/test/amber_energy/segment_%zu/time", i + 1);
        const md_attribute_t* axis = md_attributes_find(&sys.attributes, str_from_cstr(path));
        ASSERT_TRUE(axis != NULL);
        double values[2];
        ASSERT_EQ(2u, md_attribute_extract_f64(values, 2, axis, md_attribute_slice_all(), md_unit_picosecond()));
        EXPECT_NEAR(run_time[i], values[0], 1.0e-12);
        EXPECT_NEAR(run_time[i] + 1.0, values[1], 1.0e-12);
    }
    md_attributes_free(&sys.attributes);
    EXPECT_EQ(0, remove(filename));
}

UTEST(amber_energy, truncates_long_term_names_safely) {
    const char* filename = "md_unittest_amber_long_term.out";
    char name[129];
    memset(name, 'A', 126);
    name[126] = ' ';
    name[127] = 'Z';
    name[128] = '\0';
    FILE* file = fopen(filename, "wb");
    ASSERT_TRUE(file != NULL);
    fprintf(file, "NSTEP = 0 TIME(PS) = 0 Etot = -1 %s = 12\n", name);
    ASSERT_EQ(0, fclose(file));

    md_system_t sys = {.alloc = md_get_heap_allocator()};
    ASSERT_TRUE(md_amber_energy_system_supplement_from_file(&sys, str_from_cstr(filename), STR_LIT("run/test")));
    char slug[120];
    memset(slug, 'a', sizeof(slug) - 1);
    slug[sizeof(slug) - 1] = '\0';
    char path[256];
    md_snprintf(path, sizeof(path), "run/test/amber_energy/%s", slug);
    EXPECT_TRUE(md_attributes_find(&sys.attributes, str_from_cstr(path)) != NULL);
    md_attributes_free(&sys.attributes);
    EXPECT_EQ(0, remove(filename));
}

UTEST(amber_energy, rejects_incomplete_records_without_reusing_previous_time) {
    const str_t path = AMBER_ENERGY_FILE("amber-mdout-incomplete.out");
    md_system_t sys = {.alloc = md_get_heap_allocator()};
    EXPECT_FALSE(md_amber_energy_system_supplement_from_file(&sys, path, STR_LIT("run/test")));
    EXPECT_TRUE(md_attributes_find(&sys.attributes, STR_LIT("run/test/amber_energy/time")) == NULL);
    md_attributes_free(&sys.attributes);
}

UTEST(amber_energy, recognizes_and_publishes_mden_records) {
    const str_t path = AMBER_ENERGY_FILE("amber-mden.mden");
    EXPECT_TRUE(md_amber_energy_file_is_mden(path));
    EXPECT_FALSE(md_amber_energy_file_is_mdout(path));
    EXPECT_TRUE(md_amber_energy_file_is_energy(path));
    EXPECT_FALSE(md_amber_energy_file_is_mden(AMBER_ENERGY_FILE("amber-mdout.out")));

    md_allocator_i* alloc = md_get_heap_allocator();
    const double run_time[] = {0.0, 1.0, 2.0};
    md_system_t sys = {.alloc = alloc};
    sys.attributes.alloc = alloc;
    ASSERT_NE(MD_ATTRIBUTE_INVALID, amber_energy_publish_run_axis(&sys, run_time, 3));
    ASSERT_TRUE(md_amber_energy_system_supplement_from_file(&sys, path, STR_LIT("run/test")));

    const md_attribute_t* time = md_attributes_find(&sys.attributes, STR_LIT("run/test/amber_energy/time"));
    ASSERT_TRUE(time != NULL);
    EXPECT_EQ(3u, time->format.shape[0]);

    md_unit_t energy_unit = md_unit_scl(md_unit_div(md_unit_calorie(), md_unit_mole()), 1000.0);
    double values[3];
    // Shared terms take their mdout names
    const md_attribute_t* etot = md_attributes_find(&sys.attributes, STR_LIT("run/test/amber_energy/etot"));
    ASSERT_TRUE(etot != NULL);
    ASSERT_EQ(3u, md_attribute_extract_f64(values, 3, etot, md_attribute_slice_all(), energy_unit));
    EXPECT_NEAR(-10.0, values[0], 1.0e-9);
    EXPECT_NEAR(-9.5, values[1], 1.0e-9);
    EXPECT_NEAR(-9.0, values[2], 1.0e-9);

    const md_attribute_t* eptot = md_attributes_find(&sys.attributes, STR_LIT("run/test/amber_energy/eptot"));
    ASSERT_TRUE(eptot != NULL);
    EXPECT_TRUE(str_eq_cstr(eptot->label, "EPtot"));
    ASSERT_EQ(3u, md_attribute_extract_f64(values, 3, eptot, md_attribute_slice_all(), energy_unit));
    EXPECT_NEAR(-14.0, values[2], 1.0e-9);

    const md_attribute_t* temp = md_attributes_find(&sys.attributes, STR_LIT("run/test/amber_energy/temp_k"));
    ASSERT_TRUE(temp != NULL);
    EXPECT_TRUE(md_unit_equal(temp->unit, md_unit_kelvin()));
    ASSERT_EQ(3u, md_attribute_extract_f64(values, 3, temp, md_attribute_slice_all(), md_unit_kelvin()));
    EXPECT_NEAR(302.0, values[2], 1.0e-9);

    const md_attribute_t* nb14 = md_attributes_find(&sys.attributes, STR_LIT("run/test/amber_energy/1_4_nb"));
    ASSERT_TRUE(nb14 != NULL);
    ASSERT_EQ(3u, md_attribute_extract_f64(values, 3, nb14, md_attribute_slice_all(), energy_unit));
    EXPECT_NEAR(0.45, values[1], 1.0e-9);

    // Terms without an mdout counterpart keep their mden names and get their own units
    const md_attribute_t* boxx = md_attributes_find(&sys.attributes, STR_LIT("run/test/amber_energy/boxx"));
    ASSERT_TRUE(boxx != NULL);
    EXPECT_TRUE(md_unit_equal(boxx->unit, md_unit_angstrom()));
    EXPECT_TRUE(md_attributes_find(&sys.attributes, STR_LIT("run/test/amber_energy/dv_dl")) != NULL);
    EXPECT_TRUE(md_attributes_find(&sys.attributes, STR_LIT("run/test/amber_energy/density")) != NULL);
    EXPECT_TRUE(md_attributes_find(&sys.attributes, STR_LIT("run/test/amber_energy/source")) != NULL);
    md_attributes_free(&sys.attributes);
}

UTEST(amber_energy, rejects_mden_record_lines_that_do_not_match_their_header) {
    const str_t path = AMBER_ENERGY_FILE("amber-mden-bad-columns.mden");
    ASSERT_TRUE(md_amber_energy_file_is_mden(path));
    md_system_t sys = {.alloc = md_get_heap_allocator()};
    EXPECT_FALSE(md_amber_energy_system_supplement_from_file(&sys, path, STR_LIT("run/test")));
    EXPECT_TRUE(md_attributes_find(&sys.attributes, STR_LIT("run/test/amber_energy/time")) == NULL);
    md_attributes_free(&sys.attributes);
}
