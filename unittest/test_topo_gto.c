#include "utest.h"

#include "qm_test_util.h"

#include <md_molden.h>
#include <md_topo.h>

#include <core/md_allocator.h>
#include <core/md_str.h>

#include <math.h>
#include <stdatomic.h>
#include <string.h>

// Certified critical points of a GTO density (md_topo_compute_extremum_graph_gto), driven exactly the way
// a consumer drives it: a QM file is read into a system, and the basis, the total density matrix and the
// geometry come back out of the attribute table. Water, cc-pVDZ, RHF (see test_data/molden/README.md):
// the topology is three nuclear attractors and two O-H bond critical points, nothing else.

typedef struct topo_gto_input_t {
    qm_test_t      t;
    md_gto_basis_t basis;
    double*        density;
    float          xyz[3 * 3];  // Bohr
} topo_gto_input_t;

static bool topo_gto_load_water(topo_gto_input_t* in) {
    memset(in, 0, sizeof(*in));
    qm_test_init(&in->t, MEGABYTES(16));
    if (!md_molden_system_init_from_file(&in->t.sys, &in->t.state, str_from_cstr(MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz.molden"))) return false;
    if (!qm_test_basis(&in->basis, &in->t)) return false;
    size_t dim = 0;
    in->density = qm_test_matrix(&in->t, STR_LIT("orbital/total/density"), &dim);
    if (!in->density || dim != md_gto_basis_num_ao(&in->basis)) return false;
    double xyz[9] = {0};
    const md_attribute_t* coord = qm_test_attr(&in->t, STR_LIT("qm/atom/coordinate"));   // Angstrom
    if (!coord || md_attribute_extract_f64(xyz, 9, coord, md_attribute_slice_all(), md_unit_none()) != 9) return false;
    for (int i = 0; i < 9; ++i) in->xyz[i] = (float)(xyz[i] * QM_TEST_ANGSTROM_TO_BOHR);
    return true;
}

static md_topo_gto_desc_t topo_gto_desc(const topo_gto_input_t* in) {
    return (md_topo_gto_desc_t){
        .basis = &in->basis,
        .atom_xyz = in->xyz,
        .density_matrix = in->density,
        .rho_min = 1.0e-4,
        .h_min = 1.0e-4,
        .trace_separatrices = true,
    };
}

static double topo_gto_dist(const md_topo_vert_t* v, const float* p) {
    const double dx = v->x - p[0], dy = v->y - p[1], dz = v->z - p[2];
    return sqrt(dx * dx + dy * dy + dz * dz);
}

UTEST(topo_gto, water_topology_is_certified) {
    topo_gto_input_t in;
    ASSERT_TRUE(topo_gto_load_water(&in));
    md_topo_gto_desc_t desc = topo_gto_desc(&in);
    volatile float progress = -1.0f;
    desc.progress = &progress;

    md_topo_extremum_graph_t graph = { .alloc = in.t.alloc };
    md_topo_gto_info_t info;
    ASSERT_TRUE(md_topo_compute_extremum_graph_gto(&graph, &info, &desc));

    // Every cube was decided: nothing near-degenerate here, so the result is proven complete.
    EXPECT_TRUE(info.complete);
    EXPECT_EQ(1.0f, progress);
    EXPECT_EQ(0u, info.num_unresolved_boxes);
    EXPECT_EQ(1, info.poincare_hopf);
    EXPECT_GT(info.domain_pad, 3.0);

    uint32_t count[MD_TOPO_NUM_TYPES] = {0};
    md_topo_count_vertex_types(count, &graph);
    EXPECT_EQ(3u, count[MD_TOPO_MAXIMUM]);
    EXPECT_EQ(2u, count[MD_TOPO_SPLIT_SADDLE]);
    EXPECT_EQ(0u, count[MD_TOPO_JOIN_SADDLE]);
    EXPECT_EQ(0u, count[MD_TOPO_MINIMUM]);

    // The vertices in double with the Hessian there: the float vertices are them rounded, the eigenvalues are
    // ascending with signs that give the type, the eigenvectors orthonormal. A copy keeps them.
    ASSERT_TRUE(graph.vertices_f64 != NULL);
    for (uint32_t i = 0; i < graph.num_vertices; ++i) {
        const md_topo_vert_f64_t* d = &graph.vertices_f64[i];
        EXPECT_EQ(graph.vertices[i].x, (float)d->x);
        EXPECT_EQ(graph.vertices[i].y, (float)d->y);
        EXPECT_EQ(graph.vertices[i].z, (float)d->z);
        EXPECT_EQ(graph.vertices[i].value, (float)d->value);
        EXPECT_LE(d->hessian_eigenvalues[0], d->hessian_eigenvalues[1]);
        EXPECT_LE(d->hessian_eigenvalues[1], d->hessian_eigenvalues[2]);
        int npos = 0;
        for (int k = 0; k < 3; ++k) npos += d->hessian_eigenvalues[k] > 0.0;
        EXPECT_EQ(graph.types[i] == MD_TOPO_MAXIMUM ? 0 : 1, npos);   // water: maxima and bond CPs only
        for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) {
            const double* u = d->hessian_eigenvectors[a];
            const double* w = d->hessian_eigenvectors[b];
            EXPECT_NEAR(a == b ? 1.0 : 0.0, u[0] * w[0] + u[1] * w[1] + u[2] * w[2], 1e-12);
        }
    }
    {
        md_topo_extremum_graph_t copy = { .alloc = in.t.alloc };
        md_topo_extremum_graph_copy(&copy, &graph);
        ASSERT_TRUE(copy.vertices_f64 != NULL);
        EXPECT_EQ(0, memcmp(copy.vertices_f64, graph.vertices_f64, sizeof(md_topo_vert_f64_t) * graph.num_vertices));
        md_topo_extremum_graph_free(&copy);
        EXPECT_TRUE(copy.vertices_f64 == NULL);
    }

    // Each maximum is a nuclear attractor. A GTO density has no cusp, so attractors are displaced from the
    // nuclei: negligibly for O, by several hundredths of a Bohr for H.
    int max_of_atom[3] = { -1, -1, -1 };
    for (uint32_t i = 0; i < graph.num_vertices; ++i) {
        if (graph.types[i] != MD_TOPO_MAXIMUM) continue;
        int best = -1; double bd = 1e9;
        for (int a = 0; a < 3; ++a) {
            const double d = topo_gto_dist(&graph.vertices[i], in.xyz + 3 * a);
            if (d < bd) { bd = d; best = a; }
        }
        ASSERT_GE(best, 0);
        EXPECT_EQ(-1, max_of_atom[best]);
        max_of_atom[best] = (int)i;
        EXPECT_LT(bd, best == 0 ? 1.0e-3 : 0.15);
    }

    // Each bond critical point lies between O and one H, and its two separatrices (the bond path) end in
    // exactly those two attractors.
    uint32_t bcp_seen = 0;
    for (uint32_t i = 0; i < graph.num_vertices; ++i) {
        if (graph.types[i] != MD_TOPO_SPLIT_SADDLE) continue;
        bcp_seen++;
        EXPECT_GT(graph.vertices[i].value, 0.2f);   // a covalent O-H bond, rho_b ~ 0.37
        int ends[2] = { -1, -1 }, n = 0;
        for (uint32_t e = 0; e < graph.num_edges; ++e) {
            if (graph.edges[e].from == i && n < 2) ends[n++] = (int)graph.edges[e].to;
        }
        ASSERT_EQ(2, n);
        const bool o_end = ends[0] == max_of_atom[0] || ends[1] == max_of_atom[0];
        const bool h_end = ends[0] == max_of_atom[1] || ends[1] == max_of_atom[1] || ends[0] == max_of_atom[2] || ends[1] == max_of_atom[2];
        EXPECT_TRUE(o_end);
        EXPECT_TRUE(h_end);
    }
    EXPECT_EQ(2u, bcp_seen);
    EXPECT_EQ(4u, graph.num_edges);

    md_topo_extremum_graph_free(&graph);
    qm_test_free(&in.t);
}

// Same input, same output, bit for bit, whatever the thread count: cubes are decided independently and
// their results merged in cube order, with no atomics and no ordering that depends on scheduling.
// The kinetic energy density at points: at the water BCPs rho agrees with the search's and the local virial theorem
// gives an O-H bond's shared-shell signature (H < 0, |V|/G > 1). For a density of one orbital, D = c c^T, G is the
// von Weizsaecker 1/8 |grad rho|^2 / rho exactly, here with grad rho from central differences of rho.
UTEST(topo_gto, kinetic_energy_density) {
    topo_gto_input_t in;
    ASSERT_TRUE(topo_gto_load_water(&in));
    md_topo_gto_desc_t desc = topo_gto_desc(&in);
    md_topo_extremum_graph_t graph = { .alloc = in.t.alloc };
    md_topo_gto_info_t info;
    ASSERT_TRUE(md_topo_compute_extremum_graph_gto(&graph, &info, &desc));
    ASSERT_TRUE(graph.vertices_f64 != NULL);
    int bcps = 0;
    for (uint32_t i = 0; i < graph.num_vertices; ++i) {
        if (graph.types[i] != MD_TOPO_SPLIT_SADDLE) continue;
        const md_topo_vert_f64_t* v = &graph.vertices_f64[i];
        const double x[3] = { v->x, v->y, v->z };
        double G = 0.0, rho = 0.0;
        ASSERT_TRUE(md_topo_gto_kinetic_energy_density(&G, &rho, x, 1, &desc));
        EXPECT_NEAR(v->value, rho, 1e-10 * v->value);
        const double lap = v->hessian_eigenvalues[0] + v->hessian_eigenvalues[1] + v->hessian_eigenvalues[2];
        const double V = 0.25 * lap - 2.0 * G;
        EXPECT_GT(G, 0.0);
        EXPECT_LT(G + V, 0.0);
        EXPECT_GT(-V / G, 1.0);
        ++bcps;
    }
    EXPECT_EQ(2, bcps);

    const size_t n = md_gto_basis_num_ao(&in.basis);
    double* c = (double*)md_alloc(in.t.alloc, sizeof(double) * n);
    double* D = (double*)md_alloc(in.t.alloc, sizeof(double) * n * n);
    for (size_t i = 0; i < n; ++i) c[i] = 0.3 * sin(1.7 * (double)i + 0.4);
    for (size_t i = 0; i < n; ++i) for (size_t j = 0; j < n; ++j) D[i * n + j] = c[i] * c[j];
    desc.density_matrix = D;
    const double h = 1e-4;
    const double p0[3][3] = { {0.3, -0.2, 0.5}, {1.1, 0.4, -0.3}, {-0.6, 0.9, 0.2} };
    for (int k = 0; k < 3; ++k) {
        double pts[7][3];
        for (int t = 0; t < 7; ++t) for (int a = 0; a < 3; ++a) pts[t][a] = p0[k][a] + in.xyz[a];
        for (int a = 0; a < 3; ++a) { pts[1 + 2 * a][a] += h; pts[2 + 2 * a][a] -= h; }
        double G[7], rho[7];
        ASSERT_TRUE(md_topo_gto_kinetic_energy_density(G, rho, &pts[0][0], 7, &desc));
        double g2 = 0.0;
        for (int a = 0; a < 3; ++a) {
            const double d = (rho[1 + 2 * a] - rho[2 + 2 * a]) / (2.0 * h);
            g2 += d * d;
        }
        ASSERT_GT(rho[0], 1e-6);
        EXPECT_NEAR(g2 / (8.0 * rho[0]), G[0], 1e-6 * fmax(G[0], 1e-3));
    }
    md_free(in.t.alloc, D, sizeof(double) * n * n);
    md_free(in.t.alloc, c, sizeof(double) * n);
    md_topo_extremum_graph_free(&graph);
    qm_test_free(&in.t);
}

UTEST(topo_gto, density_derivatives) {
    topo_gto_input_t in;
    ASSERT_TRUE(topo_gto_load_water(&in));
    md_topo_gto_desc_t desc = topo_gto_desc(&in);
    md_topo_extremum_graph_t graph = { .alloc = in.t.alloc };
    md_topo_gto_info_t info;
    ASSERT_TRUE(md_topo_compute_extremum_graph_gto(&graph, &info, &desc));
    ASSERT_TRUE(graph.vertices_f64 != NULL);
    // At the critical points: rho as found, no gradient, and the Hessian's trace is the Laplacian
    for (uint32_t i = 0; i < graph.num_vertices; ++i) {
        const md_topo_vert_f64_t* v = &graph.vertices_f64[i];
        const double x[3] = { v->x, v->y, v->z };
        double rho, g[3], h[6];
        ASSERT_TRUE(md_topo_gto_density_derivatives(&rho, g, h, x, 1, &desc));
        EXPECT_NEAR(v->value, rho, 1e-10 * v->value);
        const double lap = v->hessian_eigenvalues[0] + v->hessian_eigenvalues[1] + v->hessian_eigenvalues[2];
        EXPECT_NEAR(lap, h[0] + h[3] + h[5], 1e-6 * fabs(lap) + 1e-9);
        EXPECT_LT(sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]), 1e-6 * fmax(rho, 1.0));
    }
    // Elsewhere: the gradient and Hessian are the finite differences of rho and the gradient
    const double step = 1e-4;
    const double p0[3][3] = { {0.3, -0.2, 0.5}, {1.1, 0.4, -0.3}, {-0.6, 0.9, 1.2} };
    static const int sym[3][3] = { {0, 1, 2}, {1, 3, 4}, {2, 4, 5} };
    for (int k = 0; k < 3; ++k) {
        double pts[7][3];
        for (int t = 0; t < 7; ++t) for (int a = 0; a < 3; ++a) pts[t][a] = p0[k][a] + in.xyz[a];
        for (int a = 0; a < 3; ++a) { pts[1 + 2 * a][a] += step; pts[2 + 2 * a][a] -= step; }
        double rho[7], g[7][3], h[7][6];
        ASSERT_TRUE(md_topo_gto_density_derivatives(rho, &g[0][0], &h[0][0], &pts[0][0], 7, &desc));
        for (int a = 0; a < 3; ++a) {
            const double fd = (rho[1 + 2 * a] - rho[2 + 2 * a]) / (2.0 * step);
            EXPECT_NEAR(fd, g[0][a], 1e-6 * fmax(fabs(g[0][a]), 1e-2));
            for (int b = 0; b < 3; ++b) {
                const double fdh = (g[1 + 2 * a][b] - g[2 + 2 * a][b]) / (2.0 * step);
                EXPECT_NEAR(fdh, h[0][sym[a][b]], 1e-5 * fmax(fabs(h[0][sym[a][b]]), 1e-1));
            }
        }
    }
    // Only rho
    const double x[3] = { in.xyz[0] + 0.5, in.xyz[1], in.xyz[2] };
    double r0 = 0.0, r1 = 0.0;
    ASSERT_TRUE(md_topo_gto_density_derivatives(&r0, NULL, NULL, x, 1, &desc));
    ASSERT_TRUE(md_topo_gto_density_derivatives(&r1, NULL, NULL, x, 1, &desc));
    EXPECT_GT(r0, 0.0);
    EXPECT_EQ(r0, r1);
    md_topo_extremum_graph_free(&graph);
    qm_test_free(&in.t);
}

typedef struct topo_pauli_input_t {
    md_gto_shell_t shells[2];
    float alpha[2], coeff[2], xyz[6];
    double density[4];
    md_gto_basis_t basis;
    md_topo_gto_desc_t desc;
} topo_pauli_input_t;

static void topo_pauli_init(topo_pauli_input_t* in) {
    memset(in, 0, sizeof(*in));
    in->shells[0] = (md_gto_shell_t){ .atom_idx = 0, .primitive_offset = 0, .num_primitives = 1 };
    in->shells[1] = (md_gto_shell_t){ .atom_idx = 1, .primitive_offset = 1, .num_primitives = 1 };
    in->alpha[0] = in->alpha[1] = 0.5f;
    in->coeff[0] = in->coeff[1] = 1.0f;
    in->xyz[0] = -1.0f;
    in->xyz[3] = 1.0f;
    in->density[0] = in->density[3] = 1.0;
    in->basis = (md_gto_basis_t){
        .num_shells = 2, .num_primitives = 2, .shells = in->shells, .alpha = in->alpha, .coeff = in->coeff,
    };
    in->desc = (md_topo_gto_desc_t){
        .basis = &in->basis, .atom_xyz = in->xyz, .density_matrix = in->density, .num_threads = 1,
    };
}

static void topo_pauli_gaussians(const double x[3], double phi[2], double grad[2][3]) {
    for (int i = 0; i < 2; ++i) {
        const double r[3] = { x[0] - (i ? 1.0 : -1.0), x[1], x[2] };
        phi[i] = exp(-0.5 * (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]));
        for (int a = 0; a < 3; ++a) grad[i][a] = -r[a] * phi[i];
    }
}

UTEST(topo_gto, pauli_rank_one_gaussians) {
    topo_pauli_input_t in;
    topo_pauli_init(&in);
    const double c[2] = {0.8, 0.4};
    for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) in.density[2 * i + j] = c[i] * c[j];
    const double xyz[4][3] = {{0, 0, 0}, {0.4, 0.2, -0.3}, {-0.7, -0.5, 0.1}, {1.3, 0.1, 0.2}};
    double P[4], G[4], rho[4], grad[4][3], old_G[4], old_rho[4];
    ASSERT_TRUE(md_topo_gto_pauli_kinetic_energy_density(P, G, rho, &grad[0][0], &xyz[0][0], 4, &in.desc));
    ASSERT_TRUE(md_topo_gto_kinetic_energy_density(old_G, old_rho, &xyz[0][0], 4, &in.desc));
    for (int t = 0; t < 4; ++t) {
        double phi[2], dphi[2][3];
        topo_pauli_gaussians(xyz[t], phi, dphi);
        const double psi = c[0] * phi[0] + c[1] * phi[1];
        double tau = 0.0;
        EXPECT_NEAR(psi * psi, rho[t], 1e-14);
        for (int a = 0; a < 3; ++a) {
            const double dpsi = c[0] * dphi[0][a] + c[1] * dphi[1][a];
            tau += 0.5 * dpsi * dpsi;
            EXPECT_NEAR(2.0 * psi * dpsi, grad[t][a], 1e-14);
        }
        EXPECT_NEAR(tau, G[t], 1e-14);
        EXPECT_NEAR(0.0, P[t], 1e-14);
        EXPECT_EQ(old_G[t], G[t]);
        EXPECT_EQ(old_rho[t], rho[t]);
    }
}

UTEST(topo_gto, pauli_two_orbital_gaussians) {
    topo_pauli_input_t in;
    topo_pauli_init(&in);
    const double xyz[4][3] = {{0, 0, 0}, {0.4, 0.2, -0.3}, {-0.7, -0.5, 0.1}, {1.3, 0.1, 0.2}};
    double P[4], G[4], rho[4], grad[4][3], ref_rho[4], ref_grad[4][3], only_P[4];
    ASSERT_TRUE(md_topo_gto_pauli_kinetic_energy_density(P, G, rho, &grad[0][0], &xyz[0][0], 4, &in.desc));
    ASSERT_TRUE(md_topo_gto_density_derivatives(ref_rho, &ref_grad[0][0], NULL, &xyz[0][0], 4, &in.desc));
    ASSERT_TRUE(md_topo_gto_pauli_kinetic_energy_density(only_P, NULL, NULL, NULL, &xyz[0][0], 4, &in.desc));
    for (int t = 0; t < 4; ++t) {
        double phi[2], dphi[2][3];
        topo_pauli_gaussians(xyz[t], phi, dphi);
        const double a = phi[0] * phi[0], b = phi[1] * phi[1];
        // For two independently occupied real orbitals: t_P = ab |grad ln phi0 - grad ln phi1|^2 / (2(a+b)).
        EXPECT_NEAR(2.0 * a * b / (a + b), P[t], 1e-14);
        EXPECT_GT(P[t], 0.0);
        EXPECT_EQ(P[t], only_P[t]);
        EXPECT_EQ(rho[t], ref_rho[t]);
        double tau = 0.0;
        for (int k = 0; k < 3; ++k) {
            tau += 0.5 * (dphi[0][k] * dphi[0][k] + dphi[1][k] * dphi[1][k]);
            EXPECT_EQ(grad[t][k], ref_grad[t][k]);
            double xp[3] = {xyz[t][0], xyz[t][1], xyz[t][2]};
            double xm[3] = {xyz[t][0], xyz[t][1], xyz[t][2]};
            const double h = 1e-5;
            xp[k] += h;
            xm[k] -= h;
            double rp, rm;
            ASSERT_TRUE(md_topo_gto_density_derivatives(&rp, NULL, NULL, xp, 1, &in.desc));
            ASSERT_TRUE(md_topo_gto_density_derivatives(&rm, NULL, NULL, xm, 1, &in.desc));
            EXPECT_NEAR((rp - rm) / (2.0 * h), grad[t][k], 1e-9);
        }
        EXPECT_NEAR(tau, G[t], 1e-14);
    }
    for (int a = 0; a < 3; ++a) EXPECT_NEAR(0.0, grad[0][a], 1e-14);
    EXPECT_EQ(G[0], P[0]);
    EXPECT_GT(G[0], 0.0);
}

UTEST(topo_gto, pauli_threads_cancel_and_vacuum) {
    topo_pauli_input_t in;
    topo_pauli_init(&in);
    enum { N = 257 };
    double xyz[N][3], P[2][N], G[2][N], rho[2][N], grad[2][N][3];
    for (int t = 0; t < N; ++t) {
        xyz[t][0] = -2.0 + 4.0 * t / (N - 1);
        xyz[t][1] = 0.2;
        xyz[t][2] = -0.1;
    }
    for (int i = 0; i < 2; ++i) {
        in.desc.num_threads = i ? 4 : 1;
        ASSERT_TRUE(md_topo_gto_pauli_kinetic_energy_density(P[i], G[i], rho[i], &grad[i][0][0], &xyz[0][0], N, &in.desc));
    }
    EXPECT_EQ(0, memcmp(P[0], P[1], sizeof(P[0])));
    EXPECT_EQ(0, memcmp(G[0], G[1], sizeof(G[0])));
    EXPECT_EQ(0, memcmp(rho[0], rho[1], sizeof(rho[0])));
    EXPECT_EQ(0, memcmp(grad[0], grad[1], sizeof(grad[0])));
    const double vacuum[3] = {1e6, 0, 0};
    ASSERT_TRUE(md_topo_gto_pauli_kinetic_energy_density(P[0], G[0], rho[0], grad[0][0], vacuum, 1, &in.desc));
    EXPECT_TRUE(isnan(P[0][0]));
    EXPECT_EQ(0.0, G[0][0]);
    EXPECT_EQ(0.0, rho[0][0]);
    for (int a = 0; a < 3; ++a) EXPECT_EQ(0.0, grad[0][0][a]);

    volatile int32_t cancel = 1;
    volatile float progress = -1.0f;
    in.desc.cancel = &cancel;
    in.desc.progress = &progress;
    EXPECT_FALSE(md_topo_gto_pauli_kinetic_energy_density(P[0], G[0], rho[0], &grad[0][0][0], &xyz[0][0], N, &in.desc));
    EXPECT_EQ(-1.0f, progress);
    for (int t = 0; t < N; ++t) {
        EXPECT_TRUE(isnan(P[0][t]));
        EXPECT_TRUE(isnan(G[0][t]));
        EXPECT_TRUE(isnan(rho[0][t]));
        for (int a = 0; a < 3; ++a) EXPECT_TRUE(isnan(grad[0][t][a]));
    }
    cancel = 0;
    EXPECT_TRUE(md_topo_gto_pauli_kinetic_energy_density(P[0], NULL, NULL, NULL, NULL, 0, &in.desc));
    EXPECT_TRUE(isnan(P[0][0]));
    in.desc.density_matrix = NULL;
    EXPECT_FALSE(md_topo_gto_pauli_kinetic_energy_density(P[0], NULL, NULL, NULL, &xyz[0][0], 1, &in.desc));
    EXPECT_TRUE(isnan(P[0][0]));
}

typedef struct topo_cancel_state_t {
    atomic_int calls;
    int limit;
} topo_cancel_state_t;

static bool topo_cancel_callback(void* data) {
    topo_cancel_state_t* state = (topo_cancel_state_t*)data;
    return atomic_fetch_add_explicit(&state->calls, 1, memory_order_relaxed) >= state->limit;
}

UTEST(topo_gto, pauli_callback_cancellation) {
    topo_pauli_input_t in;
    topo_pauli_init(&in);
    enum { N = 64 };
    double xyz[N][3] = {{0}}, P[N], G[N], rho[N], grad[N][3], hess[N][6];
    topo_cancel_state_t state = { .calls = ATOMIC_VAR_INIT(0), .limit = 0 };
    in.desc.num_threads = 4;
    in.desc.cancel_callback = topo_cancel_callback;
    in.desc.cancel_user_data = &state;
    EXPECT_FALSE(md_topo_gto_pauli_kinetic_energy_density(P, G, rho, &grad[0][0], &xyz[0][0], N, &in.desc));
    for (int t = 0; t < N; ++t) {
        EXPECT_TRUE(isnan(P[t]));
        EXPECT_TRUE(isnan(G[t]));
        EXPECT_TRUE(isnan(rho[t]));
        for (int a = 0; a < 3; ++a) EXPECT_TRUE(isnan(grad[t][a]));
    }
    atomic_store(&state.calls, 0);
    state.limit = 20;
    EXPECT_FALSE(md_topo_gto_pauli_kinetic_energy_density(P, G, rho, &grad[0][0], &xyz[0][0], N, &in.desc));
    int filled = 0;
    for (int t = 0; t < N; ++t) {
        if (isfinite(P[t])) {
            ++filled;
            EXPECT_TRUE(isfinite(G[t]));
            EXPECT_TRUE(isfinite(rho[t]));
            for (int a = 0; a < 3; ++a) EXPECT_TRUE(isfinite(grad[t][a]));
        } else {
            EXPECT_TRUE(isnan(G[t]));
            EXPECT_TRUE(isnan(rho[t]));
        }
    }
    EXPECT_GT(filled, 0);
    EXPECT_LT(filled, N);
    atomic_store(&state.calls, 0);
    state.limit = 0;
    EXPECT_FALSE(md_topo_gto_density_derivatives(rho, &grad[0][0], &hess[0][0], &xyz[0][0], N, &in.desc));
    for (int t = 0; t < N; ++t) {
        EXPECT_TRUE(isnan(rho[t]));
        for (int a = 0; a < 3; ++a) EXPECT_TRUE(isnan(grad[t][a]));
        for (int a = 0; a < 6; ++a) EXPECT_TRUE(isnan(hess[t][a]));
    }
    in.desc.cancel_callback = NULL;
    ASSERT_TRUE(md_topo_gto_pauli_kinetic_energy_density(P, G, rho, &grad[0][0], &xyz[0][0], N, &in.desc));
    for (int t = 0; t < N; ++t) EXPECT_TRUE(isfinite(P[t]));
}

UTEST(topo_gto, deterministic_across_thread_counts) {
    topo_gto_input_t in;
    ASSERT_TRUE(topo_gto_load_water(&in));
    md_topo_gto_desc_t da = topo_gto_desc(&in);
    md_topo_gto_desc_t db = topo_gto_desc(&in);
    da.num_threads = 1;
    db.num_threads = 4;

    md_topo_extremum_graph_t a = { .alloc = in.t.alloc };
    md_topo_extremum_graph_t b = { .alloc = in.t.alloc };
    md_topo_gto_info_t ia, ib;
    ASSERT_TRUE(md_topo_compute_extremum_graph_gto(&a, &ia, &da));
    ASSERT_TRUE(md_topo_compute_extremum_graph_gto(&b, &ib, &db));
    EXPECT_EQ(1u, ia.num_threads);
    EXPECT_EQ(4u, ib.num_threads);
    ASSERT_EQ(a.num_vertices, b.num_vertices);
    ASSERT_EQ(a.num_edges, b.num_edges);
    EXPECT_EQ(ia.num_box_evals, ib.num_box_evals);
    EXPECT_EQ(0, memcmp(a.vertices, b.vertices, sizeof(md_topo_vert_t) * a.num_vertices));
    EXPECT_EQ(0, memcmp(a.types, b.types, sizeof(md_topo_critical_point_type_t) * a.num_vertices));
    EXPECT_EQ(0, memcmp(a.edges, b.edges, sizeof(md_topo_edge_t) * a.num_edges));

    md_topo_extremum_graph_free(&a);
    md_topo_extremum_graph_free(&b);
    qm_test_free(&in.t);
}

// A density threshold above the bond critical points removes them, and only them: what is reported is
// exactly the set of CPs with rho >= rho_min.
UTEST(topo_gto, rho_min_is_respected) {
    topo_gto_input_t in;
    ASSERT_TRUE(topo_gto_load_water(&in));
    md_topo_gto_desc_t desc = topo_gto_desc(&in);
    desc.rho_min = 0.5;
    desc.trace_separatrices = false;

    md_topo_extremum_graph_t graph = { .alloc = in.t.alloc };
    md_topo_gto_info_t info;
    ASSERT_TRUE(md_topo_compute_extremum_graph_gto(&graph, &info, &desc));
    uint32_t count[MD_TOPO_NUM_TYPES] = {0};
    md_topo_count_vertex_types(count, &graph);
    EXPECT_EQ(1u, count[MD_TOPO_MAXIMUM]);     // only the O attractor has rho > 0.5
    EXPECT_EQ(0u, count[MD_TOPO_SPLIT_SADDLE]);
    for (uint32_t i = 0; i < graph.num_vertices; ++i) EXPECT_GE(graph.vertices[i].value, 0.5f);

    md_topo_extremum_graph_free(&graph);
    qm_test_free(&in.t);
}

// A cancel flag that is already raised stops the search before any cube is decided: the call reports
// failure, says it was cancelled, and never claims completeness, nor reports the progress of one.
UTEST(topo_gto, cancel) {
    topo_gto_input_t in;
    ASSERT_TRUE(topo_gto_load_water(&in));
    md_topo_gto_desc_t desc = topo_gto_desc(&in);
    volatile int32_t cancel = 1;
    volatile float progress = -1.0f;
    desc.cancel = &cancel;
    desc.progress = &progress;

    md_topo_extremum_graph_t graph = { .alloc = in.t.alloc };
    md_topo_gto_info_t info;
    EXPECT_FALSE(md_topo_compute_extremum_graph_gto(&graph, &info, &desc));
    EXPECT_TRUE(info.cancelled);
    EXPECT_FALSE(info.complete);
    EXPECT_GE(progress, 0.0f);
    EXPECT_LT(progress, 0.1f);

    md_topo_extremum_graph_free(&graph);
    qm_test_free(&in.t);
}

static bool topo_gto_same_graph(const md_topo_extremum_graph_t* a, const md_topo_extremum_graph_t* b) {
    // Same CPs (matched by position: symmetry-equivalent CPs of equal density may be listed in another
    // order) and the same edges between them.
    if (a->num_vertices != b->num_vertices || a->num_edges != b->num_edges || a->num_vertices > 64) return false;
    uint32_t map[64];
    bool used[64] = { false };
    for (uint32_t i = 0; i < a->num_vertices; ++i) {
        double best = 1e9; uint32_t bj = 0;
        for (uint32_t j = 0; j < b->num_vertices; ++j) {
            if (a->types[i] != b->types[j]) continue;
            const double dx = a->vertices[i].x - b->vertices[j].x, dy = a->vertices[i].y - b->vertices[j].y, dz = a->vertices[i].z - b->vertices[j].z;
            const double d = sqrt(dx * dx + dy * dy + dz * dz);
            if (d < best) { best = d; bj = j; }
        }
        if (best > 1e-5 || used[bj]) return false;
        if (fabs(a->vertices[i].value - b->vertices[bj].value) > 1e-6 * (1.0 + fabs(a->vertices[i].value))) return false;
        used[bj] = true;
        map[i] = bj;
    }
    for (uint32_t e = 0; e < a->num_edges; ++e) {
        bool found = false;
        for (uint32_t f = 0; f < b->num_edges && !found; ++f) found = b->edges[f].from == map[a->edges[e].from] && b->edges[f].to == map[a->edges[e].to];
        if (!found) return false;
    }
    return true;
}

// k waters 'sep' Bohr apart along x, D block-diagonal: an input whose AOs mostly do not meet.
#define TOPO_GTO_CHAIN_MAX 8
typedef struct topo_gto_pair_t {
    md_gto_shell_t shells[256];
    md_gto_basis_t basis;
    float          xyz[9 * TOPO_GTO_CHAIN_MAX];
    double*        D;
    size_t         N;
} topo_gto_pair_t;

static bool topo_gto_make_chain(topo_gto_pair_t* p, const topo_gto_input_t* in, int k, float sep, md_allocator_i* alloc) {
    const uint32_t ns = in->basis.num_shells;
    const size_t n = md_gto_basis_num_ao(&in->basis);
    if (k > TOPO_GTO_CHAIN_MAX || k * ns > 256) return false;
    for (int c = 0; c < k; ++c) {
        for (uint32_t s = 0; s < ns; ++s) {
            p->shells[c * ns + s] = in->basis.shells[s];
            p->shells[c * ns + s].atom_idx += 3 * c;
        }
        for (int i = 0; i < 9; ++i) p->xyz[9 * c + i] = in->xyz[i] + (i % 3 == 0 ? sep * c : 0.0f);
    }
    p->basis = (md_gto_basis_t){ .num_shells = k * ns, .num_primitives = in->basis.num_primitives, .shells = p->shells, .alpha = in->basis.alpha, .coeff = in->basis.coeff };
    p->N = k * n;
    p->D = (double*)md_alloc(alloc, sizeof(double) * p->N * p->N);
    memset(p->D, 0, sizeof(double) * p->N * p->N);
    for (int c = 0; c < k; ++c) {
        for (size_t i = 0; i < n; ++i) for (size_t j = 0; j < n; ++j) p->D[(c * n + i) * p->N + c * n + j] = in->density[i * n + j];
    }
    return md_gto_basis_num_ao(&p->basis) == p->N;
}

static md_topo_gto_desc_t topo_gto_pair_desc(const topo_gto_input_t* in, const topo_gto_pair_t* p, double rho_min) {
    md_topo_gto_desc_t desc = topo_gto_desc(in);
    desc.basis = &p->basis;
    desc.atom_xyz = p->xyz;
    desc.density_matrix = p->D;
    desc.rho_min = rho_min;
    return desc;
}

// Screening by region: a cube keeps the shells within the local radius, which lies below the far one,
// and is charged the bound of those in between. Two waters beyond each other's far radius: every CP is
// the water's, in place or moved by the separation, and so are the edges.
UTEST(topo_gto, screening_by_region) {
    topo_gto_input_t in;
    ASSERT_TRUE(topo_gto_load_water(&in));
    md_topo_gto_desc_t dw = topo_gto_desc(&in);
    dw.rho_min = 1.0e-3;
    topo_gto_pair_t pair;
    ASSERT_TRUE(topo_gto_make_chain(&pair, &in, 2, 60.0f, in.t.alloc));
    const md_topo_gto_desc_t dp = topo_gto_pair_desc(&in, &pair, 1.0e-3);

    md_topo_extremum_graph_t w = { .alloc = in.t.alloc }, g = { .alloc = in.t.alloc };
    md_topo_gto_info_t iw, ip;
    ASSERT_TRUE(md_topo_compute_extremum_graph_gto(&w, &iw, &dw));
    ASSERT_TRUE(md_topo_compute_extremum_graph_gto(&g, &ip, &dp));
    EXPECT_GT(iw.screening_rung, 0u);
    EXPECT_LT(iw.screening_radius_local, iw.screening_radius_far);
    EXPECT_GT(iw.local_ao_bound, 0u);
    EXPECT_LE(iw.local_ao_bound, (uint32_t)md_gto_basis_num_ao(&in.basis));
    EXPECT_LT(ip.screening_radius_far, 60.0);
    EXPECT_TRUE(ip.complete);
    ASSERT_EQ(2 * w.num_vertices, g.num_vertices);
    EXPECT_EQ(2 * w.num_edges, g.num_edges);
    uint32_t matched = 0;
    for (uint32_t i = 0; i < g.num_vertices; ++i) {
        for (uint32_t j = 0; j < w.num_vertices; ++j) {
            if (g.types[i] != w.types[j]) continue;
            const float shift = g.vertices[i].x > 30.0f ? 60.0f : 0.0f;
            const double dx = g.vertices[i].x - shift - w.vertices[j].x, dy = g.vertices[i].y - w.vertices[j].y, dz = g.vertices[i].z - w.vertices[j].z;
            if (sqrt(dx * dx + dy * dy + dz * dz) < 1e-4) { matched++; break; }
        }
    }
    EXPECT_EQ(g.num_vertices, matched);

    md_topo_extremum_graph_free(&w);
    md_topo_extremum_graph_free(&g);
    md_free(in.t.alloc, pair.D, sizeof(double) * pair.N * pair.N);
    qm_test_free(&in.t);
}

#if MD_ENABLE_GPU
#include <core/md_gpu.h>

// The GPU sweep (fp32 with rounding margins; roots polished in double on the CPU) reaches the same
// certified topology as the CPU reference, and is itself deterministic. The separatrices are traced on the
// GPU in fp32 (trace_main): on water they end where the CPU's do.
UTEST(topo_gto, gpu_matches_cpu) {
    md_gpu_device_t dev = md_gpu_device_create(&(md_gpu_device_desc_t){ .label = "topo_gto unittest" });
    if (!dev) UTEST_SKIP("no GPU device");
    md_gpu_stream_t stream = md_gpu_stream_default(dev, MD_GPU_STREAM_COMPUTE);

    topo_gto_input_t in;
    ASSERT_TRUE(topo_gto_load_water(&in));
    md_topo_gto_desc_t desc = topo_gto_desc(&in);
    desc.rho_min = 1.0e-3;      // the same 5 CPs in a smaller domain: quick on software Vulkan too
    volatile float progress = -1.0f;
    desc.progress = &progress;

    md_topo_extremum_graph_t cpu = { .alloc = in.t.alloc }, gpu = { .alloc = in.t.alloc }, gpu2 = { .alloc = in.t.alloc };
    md_topo_gto_info_t ic, ig, ig2;
    ASSERT_TRUE(md_topo_compute_extremum_graph_gto(&cpu, &ic, &desc));
    ASSERT_TRUE(md_topo_compute_extremum_graph_gto_gpu(&gpu, &ig, &desc, stream));
    EXPECT_EQ(1.0f, progress);
    EXPECT_TRUE(ig.used_gpu);
    EXPECT_TRUE(ig.complete);
    EXPECT_GT(ig.num_gpu_box_evals, 0u);
    EXPECT_EQ(1, ig.poincare_hopf);
    EXPECT_TRUE(ig.separatrices_on_gpu);
    EXPECT_EQ(4u, ig.num_separatrix_paths);
    EXPECT_EQ(0u, ig.num_separatrix_cpu_paths);
    EXPECT_EQ(ig.num_separatrix_evals, ig.num_separatrix_evals_gpu);
    EXPECT_TRUE(topo_gto_same_graph(&cpu, &gpu));

    ASSERT_TRUE(md_topo_compute_extremum_graph_gto_gpu(&gpu2, &ig2, &desc, stream));
    EXPECT_EQ(ig.num_box_evals, ig2.num_box_evals);
    ASSERT_EQ(gpu.num_vertices, gpu2.num_vertices);
    ASSERT_EQ(gpu.num_edges, gpu2.num_edges);
    EXPECT_EQ(0, memcmp(gpu.vertices, gpu2.vertices, sizeof(md_topo_vert_t) * gpu.num_vertices));
    EXPECT_EQ(0, memcmp(gpu.edges, gpu2.edges, sizeof(md_topo_edge_t) * gpu.num_edges));

    md_topo_extremum_graph_free(&cpu);
    md_topo_extremum_graph_free(&gpu);
    md_topo_extremum_graph_free(&gpu2);
    qm_test_free(&in.t);
    md_topo_gpu_shutdown();      // releases the kernel before its device
    md_gpu_device_destroy(dev);
}

// D on the device as the 16 x 16 tiles of AO pairs that can meet in a local list, chosen when they take
// at most half the dense size: eight waters 60 Bohr apart in a row with a block-diagonal D, where the
// tiles that hold only pairs across them are dropped. The same graph as the CPU's. (rho_min 1e-2 keeps
// the domain small; every CP of water lies far above it.)
UTEST(topo_gto, gpu_density_tiles) {
    md_gpu_device_t dev = md_gpu_device_create(&(md_gpu_device_desc_t){ .label = "topo_gto unittest" });
    if (!dev) UTEST_SKIP("no GPU device");
    md_gpu_stream_t stream = md_gpu_stream_default(dev, MD_GPU_STREAM_COMPUTE);

    topo_gto_input_t in;
    ASSERT_TRUE(topo_gto_load_water(&in));
    topo_gto_pair_t pair;
    ASSERT_TRUE(topo_gto_make_chain(&pair, &in, 8, 60.0f, in.t.alloc));
    const md_topo_gto_desc_t desc = topo_gto_pair_desc(&in, &pair, 1.0e-2);
    const size_t N = pair.N;

    md_topo_extremum_graph_t c = { .alloc = in.t.alloc }, g = { .alloc = in.t.alloc };
    md_topo_gto_info_t ic, ig;
    ASSERT_TRUE(md_topo_compute_extremum_graph_gto(&c, &ic, &desc));
    ASSERT_TRUE(md_topo_compute_extremum_graph_gto_gpu(&g, &ig, &desc, stream));
    EXPECT_TRUE(ig.used_gpu);
    EXPECT_TRUE(ig.complete);
    EXPECT_LE(2 * ig.gpu_density_bytes, sizeof(float) * N * N);   // tiles, not the dense matrix
    uint32_t count[MD_TOPO_NUM_TYPES] = {0};
    md_topo_count_vertex_types(count, &g);
    EXPECT_EQ(24u, count[MD_TOPO_MAXIMUM]);
    EXPECT_EQ(16u, count[MD_TOPO_SPLIT_SADDLE]);
    EXPECT_TRUE(topo_gto_same_graph(&c, &g));

    md_topo_extremum_graph_free(&c);
    md_topo_extremum_graph_free(&g);
    md_free(in.t.alloc, pair.D, sizeof(double) * N * N);
    qm_test_free(&in.t);
    md_topo_gpu_shutdown();
    md_gpu_device_destroy(dev);
}
#endif
