#include <md_amber_prmtop.h>
#include <md_system.h>
#include <md_util.h>
#if MD_AMBER_NETCDF
#include <md_amber_netcdf.h>
#endif

#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_parse.h>

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

typedef enum amber_block_kind_t {
    POINTERS, ATOM_NAME, CHARGE, MASS, ATOMIC_NUMBER, AMBER_ATOM_TYPE,
    RESIDUE_LABEL, RESIDUE_POINTER, BONDS_H, BONDS_NO_H, ATOMS_PER_MOLECULE,
    SOLVENT_POINTERS, BOX_DIMENSIONS, BLOCK_COUNT
} amber_block_kind_t;

typedef struct amber_block_t {
    bool present;
    char format;
    size_t width;
    size_t repeat;
    str_t* values;
} amber_block_t;

static const char* block_names[BLOCK_COUNT] = {
    "POINTERS", "ATOM_NAME", "CHARGE", "MASS", "ATOMIC_NUMBER", "AMBER_ATOM_TYPE",
    "RESIDUE_LABEL", "RESIDUE_POINTER", "BONDS_INC_HYDROGEN", "BONDS_WITHOUT_HYDROGEN",
    "ATOMS_PER_MOLECULE", "SOLVENT_POINTERS", "BOX_DIMENSIONS",
};

static bool parse_integer(str_t token, int64_t* value) {
    token = str_trim(token);
    return token.len && md_parse_i64(value, token) == token.len;
}

static bool parse_real(str_t token, double* value) {
    token = str_trim(token);
    char buf[128];
    if (!token.len || token.len >= sizeof(buf)) return false;
    memcpy(buf, token.ptr, token.len);
    for (size_t i = 0; i < token.len; ++i) {
        if (buf[i] == 'D' || buf[i] == 'd') buf[i] = 'E';
    }
    return md_parse_f64(value, (str_t){buf, token.len}) == token.len && isfinite(*value);
}

bool md_amber_prmtop_file_is_topology(str_t filename) {
    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_READ)) return false;
    char buf[512];
    const size_t len = md_file_read(file, buf, sizeof(buf));
    md_file_close(&file);
    return str_begins_with(str_trim_beg((str_t){buf, len}), STR_LIT("%VERSION"));
}

static bool parse_format(amber_block_t* block, str_t line) {
    char buf[128];
    str_copy_to_char_buf(buf, sizeof(buf), line);
    int repeat = 0, width = 0, end = 0;
    char format = 0;
    if (sscanf(buf, "%%FORMAT( %d%c%d%n", &repeat, &format, &width, &end) != 3 ||
        repeat <= 0 || repeat > 1000 || width <= 0 || width > 128) return false;
    format = (char)to_lower((unsigned char)format);
    if (format != 'a' && format != 'i' && format != 'e' && format != 'd' && format != 'f') return false;
    str_t suffix = str_trim(str_substr(line, (size_t)end, SIZE_MAX));
    if (format != 'a' && format != 'i') {
        if (str_empty(suffix) || suffix.ptr[0] != '.') return false;
        int64_t precision;
        suffix = str_substr(suffix, 1, SIZE_MAX);
        const size_t n = md_parse_i64(&precision, suffix);
        if (!n || precision < 0 || precision > width) return false;
        suffix = str_trim(str_substr(suffix, n, SIZE_MAX));
    }
    if (!str_eq(suffix, STR_LIT(")"))) return false;
    block->repeat = (size_t)repeat;
    block->width = (size_t)width;
    block->format = format;
    return true;
}

static bool parse_topology(amber_block_t blocks[BLOCK_COUNT], str_t filename, md_allocator_i* alloc) {
    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_READ)) {
        MD_LOG_ERROR("Amber topology: cannot open '" STR_FMT "'", STR_ARG(filename));
        return false;
    }
    char buf[64 * 1024];
    md_buffered_reader_t reader = md_buffered_reader_from_file(buf, sizeof(buf), file);
    str_t line;
    amber_block_t* block = NULL;
    size_t line_number = 0;
    bool version = false, ok = true;
    while (ok && md_buffered_reader_extract_line(&line, &reader)) {
        ++line_number;
        if (str_begins_with(line, STR_LIT("%VERSION"))) {
            version = true;
        } else if (str_begins_with(line, STR_LIT("%FLAG"))) {
            if (!version) { ok = false; break; }
            if (block && !block->width) { ok = false; break; }
            const str_t name = str_trim(str_substr(line, 5, SIZE_MAX));
            if (str_eq_cstr(name, "CTITLE") || str_begins_with(name, STR_LIT("CHARMM_"))) {
                MD_LOG_ERROR("Amber topology: CHAMBER topologies are not supported");
                ok = false;
                break;
            }
            block = NULL;
            for (size_t i = 0; i < BLOCK_COUNT; ++i) {
                if (!str_eq_cstr(name, block_names[i])) continue;
                if (blocks[i].present) { ok = false; break; }
                block = &blocks[i];
                block->present = true;
                break;
            }
        } else if (str_begins_with(line, STR_LIT("%FORMAT"))) {
            if (block && (block->width || !parse_format(block, line))) ok = false;
        } else if (str_begins_with(line, STR_LIT("%COMMENT"))) {
            continue;
        } else if (block && line.len) {
            if (!block->width || line.len > block->repeat * block->width) { ok = false; break; }
            for (size_t i = 0; i < line.len; i += block->width) {
                str_t value = str_trim(str_substr(line, i, MIN(block->width, line.len - i)));
                if (str_empty(value)) continue;
                if (block->format != 'a') {
                    double number;
                    int64_t integer;
                    if (!(block->format == 'i' ? parse_integer(value, &integer) : parse_real(value, &number))) {
                        ok = false;
                        break;
                    }
                }
                md_array_push(block->values, str_copy(value, alloc), alloc);
            }
        }
    }
    md_file_close(&file);
    if (!ok || !version || (block && !block->width)) {
        MD_LOG_ERROR("Amber topology: invalid or duplicate block near line %zu in '" STR_FMT "'", line_number, STR_ARG(filename));
        return false;
    }
    return true;
}

static bool block_size(const amber_block_t* block, amber_block_kind_t kind, size_t count, bool required) {
    if (!block->present && !required) return true;
    const char expected = (kind == ATOM_NAME || kind == AMBER_ATOM_TYPE || kind == RESIDUE_LABEL) ? 'a' :
                          (kind == CHARGE || kind == MASS || kind == BOX_DIMENSIONS) ? 'e' : 'i';
    const bool correct_format = expected == 'e' ? (block->format == 'e' || block->format == 'd' || block->format == 'f') :
                                                 block->format == expected;
    if (!block->present || !correct_format || md_array_size(block->values) != count) {
        MD_LOG_ERROR("Amber topology: %s must contain %zu correctly formatted values (found %zu)",
                     block_names[kind], count, md_array_size(block->values));
        return false;
    }
    return true;
}

static int64_t integer_at(const amber_block_t* block, size_t i) {
    int64_t value = 0;
    parse_integer(block->values[i], &value);
    return value;
}

static double real_at(const amber_block_t* block, size_t i) {
    double value = 0;
    parse_real(block->values[i], &value);
    return value;
}

static bool read_ascii_restart(md_system_state_t* state, str_t filename, size_t num_atoms, bool periodic) {
    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_READ)) {
        MD_LOG_ERROR("Amber restart: cannot open '" STR_FMT "'", STR_ARG(filename));
        return false;
    }
    md_temp_scope_t temp = md_temp_begin_avoid(state->alloc);
    md_allocator_i* alloc = md_temp_allocator(temp);
    char buf[64 * 1024];
    md_buffered_reader_t reader = md_buffered_reader_from_file(buf, sizeof(buf), file);
    str_t line;
    int64_t count = 0;
    bool ok = md_buffered_reader_extract_line(&line, &reader) && md_buffered_reader_extract_line(&line, &reader);
    if (ok) {
        line = str_trim(line);
        const size_t n = md_parse_i64(&count, line);
        if (!n || count <= 0 || (uint64_t)count != num_atoms) ok = false;
        str_t time = str_trim(str_substr(line, n, SIZE_MAX));
        double value;
        if (!str_empty(time) && !parse_real(time, &value)) ok = false;
    }
    double* values = NULL;
    const size_t max_values = 6 * num_atoms + 6;
    while (ok && md_buffered_reader_extract_line(&line, &reader)) {
        if (str_empty(str_trim(line))) continue;
        if (line.len > 72) { ok = false; break; }
        for (size_t i = 0; i < line.len; i += 12) {
            str_t token = str_trim(str_substr(line, i, MIN((size_t)12, line.len - i)));
            double value;
            if (str_empty(token)) continue;
            if (!parse_real(token, &value) || fabs(value) > FLT_MAX || md_array_size(values) == max_values) {
                ok = false;
                break;
            }
            md_array_push(values, value, alloc);
        }
    }
    md_file_close(&file);
    const size_t size = md_array_size(values);
    const size_t box_size = periodic ? 6 : 0;
    // IFBOX resolves the two-atom ambiguity: six trailing numbers are either a box or velocities.
    if (!ok || (size != 3 * num_atoms + box_size && size != 6 * num_atoms + box_size)) {
        MD_LOG_ERROR("Amber restart: invalid atom count or coordinate/velocity/box layout in '" STR_FMT
                     "' (topology: %zu atoms, %s)", STR_ARG(filename), num_atoms, periodic ? "periodic" : "nonperiodic");
        md_temp_end(temp);
        return false;
    }
    md_unitcell_t cell = {0};
    if (periodic) {
        const double* box = values + size - 6;
        cell = md_unitcell_from_extent_and_angles(box[0], box[1], box[2], box[3], box[4], box[5]);
        if (box[0] <= 0 || box[1] <= 0 || box[2] <= 0 ||
            box[3] <= 0 || box[3] >= 180 || box[4] <= 0 || box[4] >= 180 || box[5] <= 0 || box[5] >= 180 ||
            !isfinite(cell.x) || !isfinite(cell.y) || !isfinite(cell.z) || cell.z <= 0) {
            MD_LOG_ERROR("Amber restart: invalid periodic box in '" STR_FMT "'", STR_ARG(filename));
            md_temp_end(temp);
            return false;
        }
    }
    ok = md_system_state_init(state, num_atoms);
    if (ok) {
        for (size_t i = 0; i < num_atoms; ++i) {
            state->xyz[i] = vec3_set((float)values[3*i], (float)values[3*i+1], (float)values[3*i+2]);
        }
        state->unitcell = cell;
    }
    md_temp_end(temp);
    return ok;
}

typedef struct amber_molecule_span_t {
    size_t atom;
    size_t num_atoms;
    size_t residue;
    size_t num_residues;
} amber_molecule_span_t;

static bool same_molecule(const md_system_t* sys, const float* charge,
                          amber_molecule_span_t a, amber_molecule_span_t b) {
    if (a.num_atoms != b.num_atoms || a.num_residues != b.num_residues ||
        memcmp(sys->atom.type_idx + a.atom, sys->atom.type_idx + b.atom, a.num_atoms * sizeof(*sys->atom.type_idx)) ||
        memcmp(charge + a.atom, charge + b.atom, a.num_atoms * sizeof(*charge))) return false;
    for (size_t r = 0; r < a.num_residues; ++r) {
        if (!md_label_eq(sys->residue.name[a.residue+r], sys->residue.name[b.residue+r]) ||
            sys->residue.atom_offset[a.residue+r] - a.atom != sys->residue.atom_offset[b.residue+r] - b.atom) return false;
    }
    if (sys->bond.conn.offset) {
        for (size_t i = 0; i < a.num_atoms; ++i) {
            const size_t count = md_bond_conn_count(&sys->bond, a.atom+i);
            if (count != md_bond_conn_count(&sys->bond, b.atom+i)) return false;
            for (size_t j = 0; j < count; ++j) {
                const md_atom_idx_t x = md_bond_conn_atom_idx(&sys->bond, sys->bond.conn.offset[a.atom+i], (uint32_t)j);
                const md_atom_idx_t y = md_bond_conn_atom_idx(&sys->bond, sys->bond.conn.offset[b.atom+i], (uint32_t)j);
                if ((int64_t)x - (int64_t)a.atom != (int64_t)y - (int64_t)b.atom) return false;
            }
        }
    }
    return true;
}

typedef struct amber_topology_info_t {
    size_t num_atoms;
    size_t num_residues;
    size_t num_molecules;
    size_t num_bonds;
    uint32_t* residue_offsets;
} amber_topology_info_t;

static bool validate_topology(amber_block_t blocks[BLOCK_COUNT], amber_topology_info_t* info, md_temp_scope_t temp) {
    const amber_block_t* pointers = &blocks[POINTERS];
    if (!pointers->present || pointers->format != 'i' || md_array_size(pointers->values) < 31) {
        MD_LOG_ERROR("Amber topology: POINTERS requires at least 31 integers");
        return false;
    }

    const int64_t natom = integer_at(pointers, 0), nres = integer_at(pointers, 11);
    const int64_t nbond_h = integer_at(pointers, 2), nbond_no_h = integer_at(pointers, 3);
    const int64_t ifbox = integer_at(pointers, 27);
    if (natom <= 0 || natom > INT32_MAX / 6 || nres <= 0 || nres > natom ||
        nbond_h < 0 || nbond_no_h < 0 || nbond_h > INT32_MAX / 3 || nbond_no_h > INT32_MAX / 3 ||
        nbond_h + nbond_no_h > INT32_MAX / 3 || ifbox < 0 || ifbox > 2) {
        MD_LOG_ERROR("Amber topology: invalid atom, residue, bond counts or IFBOX");
        return false;
    }

    const size_t na = (size_t)natom, nr = (size_t)nres;
    if (!block_size(&blocks[ATOM_NAME], ATOM_NAME, na, true) ||
        !block_size(&blocks[CHARGE], CHARGE, na, true) ||
        !block_size(&blocks[MASS], MASS, na, true) ||
        !block_size(&blocks[ATOMIC_NUMBER], ATOMIC_NUMBER, na, false) ||
        !block_size(&blocks[AMBER_ATOM_TYPE], AMBER_ATOM_TYPE, na, false) ||
        !block_size(&blocks[RESIDUE_LABEL], RESIDUE_LABEL, nr, true) ||
        !block_size(&blocks[RESIDUE_POINTER], RESIDUE_POINTER, nr, true) ||
        !block_size(&blocks[BONDS_H], BONDS_H, 3 * (size_t)nbond_h, nbond_h != 0) ||
        !block_size(&blocks[BONDS_NO_H], BONDS_NO_H, 3 * (size_t)nbond_no_h, nbond_no_h != 0) ||
        !block_size(&blocks[BOX_DIMENSIONS], BOX_DIMENSIONS, 4, false) ||
        !block_size(&blocks[SOLVENT_POINTERS], SOLVENT_POINTERS, 3, ifbox != 0)) return false;

    uint32_t* offsets = md_temp_alloc_array(temp, uint32_t, nr + 1);
    if (!offsets) {
        MD_LOG_ERROR("Amber topology: failed to allocate residue offsets");
        return false;
    }
    for (size_t r = 0; r < nr; ++r) {
        const int64_t offset = integer_at(&blocks[RESIDUE_POINTER], r);
        if (offset < 1 || offset > natom || (r == 0 && offset != 1) || (r && offset - 1 <= offsets[r-1])) {
            MD_LOG_ERROR("Amber topology: invalid RESIDUE_POINTER at residue %zu", r + 1);
            return false;
        }
        offsets[r] = (uint32_t)(offset - 1);
    }
    offsets[nr] = (uint32_t)na;

    if (ifbox && (!blocks[ATOMS_PER_MOLECULE].present || blocks[ATOMS_PER_MOLECULE].format != 'i')) {
        MD_LOG_ERROR("Amber topology: periodic topology requires ATOMS_PER_MOLECULE");
        return false;
    }
    size_t molecule_atoms = 0;
    const size_t nmol = md_array_size(blocks[ATOMS_PER_MOLECULE].values);
    if (blocks[ATOMS_PER_MOLECULE].present && blocks[ATOMS_PER_MOLECULE].format != 'i') goto invalid_molecule;
    for (size_t m = 0; m < nmol; ++m) {
        const int64_t count = integer_at(&blocks[ATOMS_PER_MOLECULE], m);
        if (count <= 0 || (uint64_t)count > na - molecule_atoms) goto invalid_molecule;
        molecule_atoms += (size_t)count;
        bool boundary = false;
        for (size_t r = 1; r <= nr; ++r) boundary |= offsets[r] == molecule_atoms;
        if (!boundary) goto invalid_molecule;
    }
    if (blocks[ATOMS_PER_MOLECULE].present && molecule_atoms != na) goto invalid_molecule;
    if (blocks[SOLVENT_POINTERS].present) {
        const int64_t solute = integer_at(&blocks[SOLVENT_POINTERS], 0);
        const int64_t molecules = integer_at(&blocks[SOLVENT_POINTERS], 1);
        const int64_t solvent = integer_at(&blocks[SOLVENT_POINTERS], 2);
        if (solute < 0 || solute > nres || molecules != (int64_t)nmol || solvent < 1 || solvent > molecules + 1) goto invalid_molecule;
    }
    if (ifbox && blocks[BOX_DIMENSIONS].present) {
        const amber_block_t* box = &blocks[BOX_DIMENSIONS];
        if (real_at(box, 0) <= 0 || real_at(box, 0) >= 180 ||
            real_at(box, 1) <= 0 || real_at(box, 2) <= 0 || real_at(box, 3) <= 0) {
            MD_LOG_ERROR("Amber topology: invalid BOX_DIMENSIONS");
            return false;
        }
    }
    for (size_t i = 0; i < na; ++i) {
        const double mass = real_at(&blocks[MASS], i);
        const double charge = real_at(&blocks[CHARGE], i) / 18.2223;
        const int64_t z = blocks[ATOMIC_NUMBER].present ? integer_at(&blocks[ATOMIC_NUMBER], i) : 0;
        if (mass < 0 || mass > FLT_MAX || fabs(charge) > FLT_MAX || z < 0 || z > 118) {
            MD_LOG_ERROR("Amber topology: invalid mass, charge or atomic number at atom %zu", i + 1);
            return false;
        }
    }
    for (size_t k = BONDS_H; k <= BONDS_NO_H; ++k) {
        const amber_block_t* block = &blocks[k];
        for (size_t b = 0; b < md_array_size(block->values); b += 3) {
            const int64_t a = integer_at(block, b), c = integer_at(block, b+1), parameter = integer_at(block, b+2);
            if (a < 0 || c < 0 || a % 3 || c % 3 || a / 3 >= natom || c / 3 >= natom || a == c ||
                parameter <= 0 || parameter > integer_at(pointers, 15)) {
                MD_LOG_ERROR("Amber topology: invalid bond in %s", block_names[k]);
                return false;
            }
        }
    }

    info->num_atoms = na;
    info->num_residues = nr;
    info->num_molecules = nmol;
    info->num_bonds = (size_t)(nbond_h + nbond_no_h);
    info->residue_offsets = offsets;
    return true;

invalid_molecule:
    MD_LOG_ERROR("Amber topology: invalid molecule or solvent pointers (molecules must end at residue boundaries)");
    return false;
}

bool md_amber_prmtop_system_init_from_file(md_system_t* sys, md_system_state_t* state,
                                         str_t filename, str_t coordinates) {
    ASSERT(sys && state);
    if (!sys->alloc || !state->alloc || str_empty(coordinates)) {
        MD_LOG_ERROR("Amber topology: set allocators and supply a matching rst7/inpcrd or NetCDF coordinate file");
        return false;
    }
    md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
    md_allocator_i* alloc = md_temp_allocator(temp);
    amber_block_t blocks[BLOCK_COUNT] = {0};
    bool result = false;
    amber_topology_info_t topology = {0};
    if (!parse_topology(blocks, filename, alloc) || !validate_topology(blocks, &topology, temp)) goto done;
    const size_t na = topology.num_atoms, nr = topology.num_residues, nmol = topology.num_molecules;
    const uint32_t* offsets = topology.residue_offsets;
    const int64_t ifbox = integer_at(&blocks[POINTERS], 27);
    md_system_state_t initial = {.alloc = alloc};
    bool netcdf = false;
    md_file_t coordinate_file = {0};
    if (md_file_open(&coordinate_file, coordinates, MD_FILE_READ)) {
        unsigned char magic[4] = {0};
        md_file_read(coordinate_file, magic, sizeof(magic));
        md_file_close(&coordinate_file);
        netcdf = !memcmp(magic, "CDF", 3) || !memcmp(magic, "\211HDF", 4);
    }
    if (netcdf) {
#if MD_AMBER_NETCDF
        if (!md_amber_netcdf_state_init_from_file(&initial, coordinates, na)) goto done;
        if (ifbox && !md_system_state_has_unitcell(&initial)) {
            MD_LOG_ERROR("Amber topology: periodic topology requires a cell in its NetCDF coordinates");
            goto done;
        }
#else
        MD_LOG_ERROR("Amber topology: NetCDF coordinates require MD_ENABLE_AMBER_NETCDF");
        goto done;
#endif
    } else if (!read_ascii_restart(&initial, coordinates, na, ifbox != 0)) goto done;

    md_system_reset(sys);
    if (!md_system_state_copy(state, &initial)) goto done;
    float* charge = md_temp_alloc_array(temp, float, na);
    float* mass = md_temp_alloc_array(temp, float, na);
    md_particle_type_add(&sys->atom.type, STR_LIT("Unk"), (str_t){0}, 0, 0, 0, 0, 0, sys->alloc);
    const size_t capacity = ROUND_UP(na, 16);
    md_array_resize(sys->atom.type_idx, capacity, sys->alloc);
    md_array_resize(sys->atom.flags, capacity, sys->alloc);
    memset(sys->atom.type_idx, 0, capacity * sizeof(*sys->atom.type_idx));
    memset(sys->atom.flags, 0, capacity * sizeof(*sys->atom.flags));
    size_t r = 0, inferred = 0;
    for (size_t i = 0; i < na; ++i) {
        while (i >= offsets[r+1]) ++r;
        mass[i] = (float)real_at(&blocks[MASS], i);
        charge[i] = (float)(real_at(&blocks[CHARGE], i) / 18.2223);
        md_atomic_number_t z = blocks[ATOMIC_NUMBER].present ? (md_atomic_number_t)integer_at(&blocks[ATOMIC_NUMBER], i) : 0;
        if (!z && mass[i] > 0) {
            z = md_atomic_number_infer_from_mass((float)mass[i]);
            if (!z) z = md_atomic_number_infer_from_label(blocks[ATOM_NAME].values[i], blocks[RESIDUE_LABEL].values[r], offsets[r+1] - offsets[r]);
            inferred += z != 0;
        }
        const str_t name = blocks[ATOM_NAME].values[i];
        const str_t ff = blocks[AMBER_ATOM_TYPE].present ? blocks[AMBER_ATOM_TYPE].values[i] : (str_t){0};
        const md_particle_type_flags_t flags = md_particle_type_flags_set_particle_kind(0,
            mass[i] == 0 && z == 0 ? MD_PARTICLE_VIRTUAL_SITE : MD_PARTICLE_ATOM);
        md_particle_type_idx_t type = 0;
        for (size_t t = 1; t < sys->atom.type.count; ++t) {
            if (str_eq(name, md_label_str(&sys->atom.type.name[t])) && str_eq(ff, md_label16_str(&sys->atom.type.ff_type[t])) &&
                z == sys->atom.type.z[t] && (float)mass[i] == sys->atom.type.mass[t] && flags == sys->atom.type.flags[t]) {
                type = (md_particle_type_idx_t)t;
                break;
            }
        }
        if (!type) type = md_particle_type_add(&sys->atom.type, name, ff, z, (float)mass[i],
                                              md_util_element_vdw_radius(z), md_util_element_cpk_color(z), flags, sys->alloc);
        sys->atom.type_idx[i] = type;
    }
    sys->atom.count = na;
    if (inferred) MD_LOG_INFO("Amber topology: inferred elements from masses/names for %zu atoms", inferred);
    for (size_t i = 0; i < nr; ++i) {
        md_array_push(sys->residue.name, make_label(blocks[RESIDUE_LABEL].values[i]), sys->alloc);
        md_array_push(sys->residue.seq_id, (md_sequence_id_t)(i + 1), sys->alloc);
        md_array_push(sys->residue.flags, MD_RESIDUE_FLAG_NONE, sys->alloc);
        md_array_push(sys->residue.atom_offset, offsets[i], sys->alloc);
    }
    md_array_push(sys->residue.atom_offset, offsets[nr], sys->alloc);
    sys->residue.count = nr;
    for (size_t k = BONDS_H; k <= BONDS_NO_H; ++k) {
        const amber_block_t* block = &blocks[k];
        for (size_t b = 0; b < md_array_size(block->values); b += 3) {
            const md_atom_pair_t pair = {{(md_atom_idx_t)(integer_at(block, b)/3), (md_atom_idx_t)(integer_at(block, b+1)/3)}};
            md_array_push(sys->bond.pairs, pair, sys->alloc);
            md_array_push(sys->bond.flags, md_bond_flags_set_origin(0, MD_BOND_ORIGIN_TOPOLOGY), sys->alloc);
            ++sys->bond.count;
        }
    }
    md_bond_build_connectivity(&sys->bond, na, sys->alloc);
    amber_molecule_span_t* types = NULL;
    size_t atom = 0, residue = 0;
    for (size_t m = 0; m < nmol; ++m) {
        const size_t count = (size_t)integer_at(&blocks[ATOMS_PER_MOLECULE], m);
        size_t end_residue = residue;
        while (end_residue < nr && offsets[end_residue] < atom + count) ++end_residue;
        const amber_molecule_span_t span = {atom, count, residue, end_residue - residue};
        size_t type = 0;
        while (type < md_array_size(types) && !same_molecule(sys, charge, types[type], span)) ++type;
        if (type == md_array_size(types)) {
            md_array_push(types, span, alloc);
            md_molecule_type_data_push(&sys->molecule.type, blocks[RESIDUE_LABEL].values[residue],
                                      MD_MOLECULE_TYPE_FLAG_NONE, sys->alloc);
        }
        md_array_push(sys->instance.id, (md_label_t){0}, sys->alloc);
        md_array_push(sys->instance.auth_id, (md_label_t){0}, sys->alloc);
        md_array_push(sys->instance.res_offset, (uint32_t)residue, sys->alloc);
        md_array_push(sys->instance.part_type_idx, -1, sys->alloc);
        md_array_push(sys->instance.mol_type_idx, (md_molecule_type_idx_t)type, sys->alloc);
        atom += count;
        residue = end_residue;
        ++sys->instance.count;
    }
    if (nmol) md_array_push(sys->instance.res_offset, (uint32_t)nr, sys->alloc);
    md_util_system_augment_particle_types(sys);
    md_util_system_infer_residue_flags(sys);
    md_attributes_publish_atom_column(&sys->attributes, STR_LIT("atom/charge"), md_unit_elementary_charge(), 1, charge, na);
    md_attributes_publish_atom_column(&sys->attributes, STR_LIT("atom/mass"), md_unit_dalton(), 1, mass, na);
    result = true;
    goto done;

done:
    md_temp_end(temp);
    return result;
}

static bool amber_name_matches(str_t structure, str_t topology) {
    if (str_eq_ignore_case(structure, topology)) return true;
    return structure.len >= 3 && topology.len > structure.len &&
           str_eq_ignore_case(structure, str_substr(topology, 0, structure.len));
}

bool md_amber_prmtop_system_supplement_from_file(md_system_t* sys, str_t filename) {
    ASSERT(sys);
    if (!sys->alloc || !sys->atom.count || !sys->atom.type_idx || !sys->atom.type.count ||
        !sys->atom.type.name || !sys->atom.type.ff_type || !sys->atom.type.z ||
        !sys->atom.type.mass || !sys->atom.type.flags || !sys->residue.name ||
        !sys->residue.atom_offset) {
        MD_LOG_ERROR("Amber topology: load a structure with atoms and residues before applying this topology");
        return false;
    }

    md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
    md_allocator_i* alloc = md_temp_allocator(temp);
    amber_block_t blocks[BLOCK_COUNT] = {0};
    amber_topology_info_t topology = {0};
    bool result = false;
    if (!parse_topology(blocks, filename, alloc) || !validate_topology(blocks, &topology, temp)) goto done;

    if (topology.num_atoms != sys->atom.count) {
        MD_LOG_ERROR("Amber topology: %zu atoms do not match the loaded structure's %zu atoms",
                     topology.num_atoms, sys->atom.count);
        goto done;
    }
    if (topology.num_residues != sys->residue.count) {
        MD_LOG_ERROR("Amber topology: %zu residues do not match the loaded structure's %zu residues",
                     topology.num_residues, sys->residue.count);
        goto done;
    }

    bool truncated_name = false;
    for (size_t r = 0; r < topology.num_residues; ++r) {
        const md_urange_t structure_range = md_system_residue_atom_range(sys, r);
        if (structure_range.beg != topology.residue_offsets[r] ||
            structure_range.end != topology.residue_offsets[r + 1] ||
            !amber_name_matches(md_system_residue_name(sys, r), blocks[RESIDUE_LABEL].values[r])) {
            MD_LOG_ERROR("Amber topology: residue %zu or its atom range does not match the loaded structure",
                         r + 1);
            goto done;
        }
        truncated_name |= !str_eq_ignore_case(md_system_residue_name(sys, r), blocks[RESIDUE_LABEL].values[r]);
    }
    for (size_t i = 0; i < topology.num_atoms; ++i) {
        const str_t structure_name = md_atom_name(&sys->atom, i);
        const str_t topology_name = blocks[ATOM_NAME].values[i];
        if (!amber_name_matches(structure_name, topology_name)) {
            MD_LOG_ERROR("Amber topology: atom %zu name '" STR_FMT "' does not match the loaded structure's '" STR_FMT "'",
                         i + 1, STR_ARG(topology_name), STR_ARG(structure_name));
            goto done;
        }
        truncated_name |= !str_eq_ignore_case(structure_name, topology_name);
    }
    if (truncated_name) {
        MD_LOG_INFO("Amber topology: some names are truncated or differ only by case; atom and residue order was checked");
    }

    float* charge = md_temp_alloc_array(temp, float, topology.num_atoms);
    float* mass = md_temp_alloc_array(temp, float, topology.num_atoms);
    md_atomic_number_t* elements = md_temp_alloc_array(temp, md_atomic_number_t, topology.num_atoms);
    if (!charge || !mass || !elements) {
        MD_LOG_ERROR("Amber topology: failed to allocate per-atom properties");
        goto done;
    }
    for (size_t i = 0; i < topology.num_atoms; ++i) {
        const md_particle_type_idx_t old_type = sys->atom.type_idx[i];
        if (old_type >= sys->atom.type.count) {
            MD_LOG_ERROR("Amber topology: loaded structure has an invalid atom type at atom %zu", i + 1);
            goto done;
        }
        const float atom_mass = (float)real_at(&blocks[MASS], i);
        md_atomic_number_t z = blocks[ATOMIC_NUMBER].present ? (md_atomic_number_t)integer_at(&blocks[ATOMIC_NUMBER], i) : 0;
        if (!z) z = sys->atom.type.z[old_type];
        if (!z && atom_mass > 0) {
            z = md_atomic_number_infer_from_mass(atom_mass);
            if (!z) {
                const size_t residue = (size_t)md_system_residue_find_by_atom_idx(sys, i);
                z = md_atomic_number_infer_from_label(blocks[ATOM_NAME].values[i],
                    blocks[RESIDUE_LABEL].values[residue],
                    topology.residue_offsets[residue + 1] - topology.residue_offsets[residue]);
            }
        }
        const md_atomic_number_t structure_z = sys->atom.type.z[old_type];
        if (structure_z && z && structure_z != z) {
            MD_LOG_ERROR("Amber topology: element mismatch at atom %zu (structure %u, topology %u)",
                         i + 1, (unsigned)structure_z, (unsigned)z);
            goto done;
        }
        elements[i] = z;
        mass[i] = atom_mass;
        charge[i] = (float)(real_at(&blocks[CHARGE], i) / 18.2223);
    }

    const size_t old_bond_count = sys->bond.count;
    const size_t num_noncovalent = sys->bond.noncovalent_count;
    size_t num_user_bonds = 0;
    for (size_t i = 0; i < old_bond_count; ++i) {
        num_user_bonds += md_bond_origin(sys->bond.flags[i]) == MD_BOND_ORIGIN_USER;
    }
    const size_t new_bond_count = topology.num_bonds + num_user_bonds;
    md_atom_pair_t* pairs = md_temp_alloc_array(temp, md_atom_pair_t,
                                                 new_bond_count + num_noncovalent + 1);
    md_bond_flags_t* flags = md_temp_alloc_array(temp, md_bond_flags_t,
                                                  new_bond_count + num_noncovalent + 1);
    if (!pairs || !flags) {
        MD_LOG_ERROR("Amber topology: failed to allocate topology bonds");
        goto done;
    }

    for (size_t i = 0; i < topology.num_atoms; ++i) {
        const md_particle_type_idx_t old_type = sys->atom.type_idx[i];
        const float atom_mass = mass[i];
        const md_atomic_number_t z = elements[i];
        const str_t atom_name = blocks[ATOM_NAME].values[i];
        const str_t ff_type = blocks[AMBER_ATOM_TYPE].present ? blocks[AMBER_ATOM_TYPE].values[i] : (str_t){0};
        const md_particle_type_flags_t flags = md_particle_type_flags_set_particle_kind(
            sys->atom.type.flags[old_type],
            atom_mass == 0 && z == 0 ? MD_PARTICLE_VIRTUAL_SITE : MD_PARTICLE_ATOM);
        md_particle_type_idx_t type = 0;
        for (size_t t = 1; t < sys->atom.type.count; ++t) {
            if (str_eq(atom_name, md_label_str(&sys->atom.type.name[t])) &&
                str_eq(ff_type, md_label16_str(&sys->atom.type.ff_type[t])) &&
                z == sys->atom.type.z[t] && atom_mass == sys->atom.type.mass[t] &&
                flags == sys->atom.type.flags[t]) {
                type = (md_particle_type_idx_t)t;
                break;
            }
        }
        if (!type) {
            type = md_particle_type_add(&sys->atom.type, atom_name, ff_type, z, atom_mass,
                md_util_element_vdw_radius(z), md_util_element_cpk_color(z), flags, sys->alloc);
        }
        sys->atom.type_idx[i] = type;
    }

    size_t bond_index = 0;
    for (size_t i = 0; i < old_bond_count; ++i) {
        if (md_bond_origin(sys->bond.flags[i]) != MD_BOND_ORIGIN_USER) continue;
        pairs[bond_index] = sys->bond.pairs[i];
        flags[bond_index++] = sys->bond.flags[i];
    }
    for (size_t k = BONDS_H; k <= BONDS_NO_H; ++k) {
        const amber_block_t* block = &blocks[k];
        for (size_t b = 0; b < md_array_size(block->values); b += 3) {
            pairs[bond_index] = (md_atom_pair_t){{
                (md_atom_idx_t)(integer_at(block, b) / 3),
                (md_atom_idx_t)(integer_at(block, b + 1) / 3),
            }};
            flags[bond_index++] = md_bond_flags_set_origin(MD_BOND_FLAG_NONE, MD_BOND_ORIGIN_TOPOLOGY);
        }
    }
    if (num_noncovalent) {
        memcpy(pairs + bond_index, sys->bond.pairs + old_bond_count, num_noncovalent * sizeof(*pairs));
        memcpy(flags + bond_index, sys->bond.flags + old_bond_count, num_noncovalent * sizeof(*flags));
    }

    md_bond_data_clear(&sys->bond);
    md_array_resize(sys->bond.pairs, new_bond_count + num_noncovalent, sys->alloc);
    md_array_resize(sys->bond.flags, new_bond_count + num_noncovalent, sys->alloc);
    if (new_bond_count) {
        memcpy(sys->bond.pairs, pairs, new_bond_count * sizeof(*pairs));
        memcpy(sys->bond.flags, flags, new_bond_count * sizeof(*flags));
    }
    if (num_noncovalent) {
        memcpy(sys->bond.pairs + new_bond_count, pairs + bond_index, num_noncovalent * sizeof(*pairs));
        memcpy(sys->bond.flags + new_bond_count, flags + bond_index, num_noncovalent * sizeof(*flags));
    }
    sys->bond.count = new_bond_count;
    sys->bond.noncovalent_count = num_noncovalent;
    md_bond_build_connectivity(&sys->bond, sys->atom.count, sys->alloc);
    md_util_system_infer_coordination(sys);
    md_util_system_infer_elastic_bonds(sys);

    if (!sys->attributes.alloc) sys->attributes.alloc = sys->alloc;
    const md_attribute_t* old_charge = md_attributes_find(&sys->attributes, STR_LIT("atom/charge"));
    const md_attribute_t* old_mass = md_attributes_find(&sys->attributes, STR_LIT("atom/mass"));
    const md_attribute_id_t old_charge_id = old_charge ? old_charge->id : MD_ATTRIBUTE_INVALID;
    const md_attribute_id_t old_mass_id = old_mass ? old_mass->id : MD_ATTRIBUTE_INVALID;
    if (old_charge_id != MD_ATTRIBUTE_INVALID) md_attributes_remove(&sys->attributes, old_charge_id);
    if (old_mass_id != MD_ATTRIBUTE_INVALID) md_attributes_remove(&sys->attributes, old_mass_id);
    md_attributes_publish_atom_column(&sys->attributes,
        STR_LIT("atom/charge"), md_unit_elementary_charge(), 1, charge, topology.num_atoms);
    md_attributes_publish_atom_column(&sys->attributes,
        STR_LIT("atom/mass"), md_unit_dalton(), 1, mass, topology.num_atoms);

    md_system_topology_changed(sys);
    md_util_system_infer(sys, &sys->reference,
        MD_UTIL_INFER_PART_BIT | MD_UTIL_INFER_MOLECULE_BIT | MD_UTIL_INFER_BACKBONE_BIT);
    md_util_system_infer_lipids(sys);
    md_util_system_infer_residue_types(sys);
    result = true;

done:
    md_temp_end(temp);
    return result;
}
