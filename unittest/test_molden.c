#include "utest.h"

#include <md_molden.h>
#include <md_system.h>
#include <core/md_allocator.h>
#include <core/md_str.h>

UTEST(molden, create_destroy) {
    md_allocator_i* alloc = md_get_heap_allocator();
    
    md_molden_t* molden = md_molden_create(alloc);
    ASSERT_TRUE(molden != NULL);
    
    md_molden_destroy(molden);
}

UTEST(molden, parse_water_molecule) {
    md_allocator_i* alloc = md_get_heap_allocator();
    
    md_molden_t* molden = md_molden_create(alloc);
    ASSERT_TRUE(molden != NULL);
    
    // Parse water molecule test file
    str_t filename = STR_LIT(MD_UNITTEST_DATA_DIR "/water.molden");
    bool success = md_molden_parse_file(molden, filename);
    ASSERT_TRUE(success);
    
    // Check atom count
    size_t num_atoms = md_molden_number_of_atoms(molden);
    EXPECT_EQ(3, num_atoms);
    
    // Check atomic numbers
    const uint8_t* atomic_numbers = md_molden_atomic_numbers(molden);
    ASSERT_TRUE(atomic_numbers != NULL);
    EXPECT_EQ(8, atomic_numbers[0]);  // Oxygen
    EXPECT_EQ(1, atomic_numbers[1]);  // Hydrogen
    EXPECT_EQ(1, atomic_numbers[2]);  // Hydrogen
    
    // Check coordinates
    const dvec3_t* coords = md_molden_atom_coordinates(molden);
    ASSERT_TRUE(coords != NULL);
    
    // Oxygen should be near origin
    EXPECT_NEAR(0.0, coords[0].x, 0.0001);
    EXPECT_NEAR(0.0, coords[0].y, 0.0001);
    EXPECT_NEAR(0.117176, coords[0].z, 0.0001);
    
    // First hydrogen
    EXPECT_NEAR(0.0, coords[1].x, 0.0001);
    EXPECT_NEAR(0.758602, coords[1].y, 0.0001);
    EXPECT_NEAR(-0.468706, coords[1].z, 0.0001);
    
    // Second hydrogen
    EXPECT_NEAR(0.0, coords[2].x, 0.0001);
    EXPECT_NEAR(-0.758602, coords[2].y, 0.0001);
    EXPECT_NEAR(-0.468706, coords[2].z, 0.0001);
    
    md_molden_destroy(molden);
}

UTEST(molden, system_init) {
    md_allocator_i* alloc = md_get_heap_allocator();
    
    md_molden_t* molden = md_molden_create(alloc);
    ASSERT_TRUE(molden != NULL);
    
    str_t filename = STR_LIT(MD_UNITTEST_DATA_DIR "/water.molden");
    bool success = md_molden_parse_file(molden, filename);
    ASSERT_TRUE(success);
    
    // Initialize system from molden data
    md_system_t molecule = {0};
    success = md_molden_system_init(&molecule, molden, alloc);
    ASSERT_TRUE(success);
    
    // Verify system data
    EXPECT_EQ(3, molecule.atom.count);
    ASSERT_TRUE(molecule.atom.x != NULL);
    ASSERT_TRUE(molecule.atom.y != NULL);
    ASSERT_TRUE(molecule.atom.z != NULL);
    ASSERT_TRUE(molecule.atom.type.z != NULL);
    
    // Check atomic numbers in system
    uint8_t z0 = molecule.atom.type.z[0];
    uint8_t z1 = molecule.atom.type.z[1];
    uint8_t z2 = molecule.atom.type.z[2];
    EXPECT_EQ(8, z0);
    EXPECT_EQ(1, z1);
    EXPECT_EQ(1, z2);
    
    // Clean up
    md_system_free(&molecule, alloc);
    md_molden_destroy(molden);
}

UTEST(molden, reset) {
    md_allocator_i* alloc = md_get_heap_allocator();
    
    md_molden_t* molden = md_molden_create(alloc);
    ASSERT_TRUE(molden != NULL);
    
    str_t filename = STR_LIT(MD_UNITTEST_DATA_DIR "/water.molden");
    bool success = md_molden_parse_file(molden, filename);
    ASSERT_TRUE(success);
    EXPECT_EQ(3, md_molden_number_of_atoms(molden));
    
    // Reset should clear all data
    md_molden_reset(molden);
    EXPECT_EQ(0, md_molden_number_of_atoms(molden));
    
    md_molden_destroy(molden);
}
