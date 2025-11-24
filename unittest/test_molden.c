#include "utest.h"

#include <md_molden.h>
#include <core/md_allocator.h>
#include <core/md_str.h>

UTEST(molden, molden_parse) {
	md_molden_t* molden = md_molden_create(md_get_heap_allocator());
    ASSERT_TRUE(molden != NULL);
    
    MD_LOG_INFO("Testing Molden file parsing...");
    bool result = md_molden_parse_file(molden, STR_LIT(MD_UNITTEST_DATA_DIR "/molden/h2o.molden"));
    if (!result) {
        MD_LOG_ERROR("Parse file failed!");
    }
    ASSERT_TRUE(result);

    MD_LOG_INFO("Parse successful, checking data...");
    
    // Check basic molecular properties
    size_t num_atoms = md_molden_number_of_atoms(molden);
    MD_LOG_INFO("Number of atoms: %zu", num_atoms);
    EXPECT_EQ(3, num_atoms);
    EXPECT_EQ(10, md_molden_number_of_alpha_electrons(molden));
    EXPECT_EQ(10, md_molden_number_of_beta_electrons(molden));
    EXPECT_EQ(1, md_molden_spin_multiplicity(molden));

    // Check atomic numbers
    const md_element_t* atomic_numbers = md_molden_atomic_numbers(molden);
    ASSERT_TRUE(atomic_numbers != NULL);
    EXPECT_EQ(8, atomic_numbers[0]);  // Oxygen
    EXPECT_EQ(1, atomic_numbers[1]);  // Hydrogen
    EXPECT_EQ(1, atomic_numbers[2]);  // Hydrogen

    // Check coordinates (in Ångström)
    const dvec3_t* coords = md_molden_atom_coordinates(molden);
    ASSERT_TRUE(coords != NULL);
    EXPECT_NEAR(0.0, coords[0].x, 1.0e-5);
    EXPECT_NEAR(0.0, coords[0].y, 1.0e-5);
    EXPECT_NEAR(0.119262, coords[0].z, 1.0e-5);

    // Check MO data
    size_t num_aos = md_molden_scf_number_of_atomic_orbitals(molden);
    size_t num_mos = md_molden_scf_number_of_molecular_orbitals(molden);
    EXPECT_EQ(24, num_aos);
    EXPECT_EQ(5, num_mos);

    // Check MO energies
    const double* energies = md_molden_scf_mo_energy(molden, MD_MOLDEN_MO_TYPE_ALPHA);
    ASSERT_TRUE(energies != NULL);
    EXPECT_NEAR(-20.2515, energies[0], 1.0e-3);
    EXPECT_NEAR(-1.2577, energies[1], 1.0e-3);
    EXPECT_NEAR(-0.5976, energies[2], 1.0e-3);

    // Check MO occupancies
    const double* occup = md_molden_scf_mo_occupancy(molden, MD_MOLDEN_MO_TYPE_ALPHA);
    ASSERT_TRUE(occup != NULL);
    EXPECT_NEAR(2.0, occup[0], 1.0e-5);
    EXPECT_NEAR(2.0, occup[1], 1.0e-5);

    // Check HOMO/LUMO indices
    size_t homo_idx = md_molden_scf_homo_idx(molden, MD_MOLDEN_MO_TYPE_ALPHA);
    size_t lumo_idx = md_molden_scf_lumo_idx(molden, MD_MOLDEN_MO_TYPE_ALPHA);
    EXPECT_EQ(4, homo_idx);
    EXPECT_EQ(0, lumo_idx);  // Will be 0 since all MOs are occupied

    md_molden_destroy(molden);
}

UTEST(molden, molden_gto_extraction) {
	md_molden_t* molden = md_molden_create(md_get_heap_allocator());
    ASSERT_TRUE(molden != NULL);
    
    bool result = md_molden_parse_file(molden, STR_LIT(MD_UNITTEST_DATA_DIR "/molden/h2o.molden"));
    ASSERT_TRUE(result);

    // Extract GTO data
    md_gto_data_t gto_data;
    MEMSET(&gto_data, 0, sizeof(md_gto_data_t));
    
    result = md_molden_scf_extract_gto_data(&gto_data, molden, 0.0, md_get_heap_allocator());
    ASSERT_TRUE(result);
    
    // Check that we have GTOs
    EXPECT_GT(gto_data.num_pgtos, 0);
    EXPECT_GT(gto_data.num_cgtos, 0);

    md_molden_destroy(molden);
}
