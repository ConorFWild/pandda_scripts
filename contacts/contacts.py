"""
Dependencies: fire, numpy, gemmi

"""

from pathlib import Path
import json

import fire
import numpy as np
import gemmi


class Constants:
    DEBUG = True


def get_out_file(out_path, pdb_path):
    """
    If you give a dir return the Path to dir/{protein_name}.json
    If you give a .json file, return the Path to file
    :param out_path:
    :param pdb_path:
    :return:
    """

    if out_path.is_dir():
        pdb_stem = pdb_path.stem
        out_file = out_path / f"{pdb_stem}.json"
        return out_file
    elif (out_path.suffix == ".json") and (out_path.parent.is_dir()):
        return out_path
    else:
        raise Exception("Out path must either be a .json file in an extant directory or a directory")


def get_structure(pdb_path):
    """
    Load a structure from a path
    :param pdb_path:
    :return:
    """
    structure = gemmi.read_structure(str(pdb_path))
    structure.setup_entities()
    return structure


def get_residues(structure, selection):
    """
    Get a dictionary of residues whose name are equal to the selection ("LIG" by default)
    :param structure:
    :param selection:
    :return:
    """
    residues = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.name == selection:
                    resid = f"{model.name}_{chain.name}_{residue.seqid.num}"
                    residues[resid] = residue

    return residues


def get_grid(structure):
    """
    Get the grid associated with a structure's resolution/unit cell/symmetry
    :param structure:
    :return:
    """
    cell = structure.cell

    res = structure.resolution

    spacegroup = gemmi.find_spacegroup_by_name(structure.spacegroup_hm)

    hkl = gemmi.make_miller_array(cell, spacegroup, res)

    mtz = gemmi.Mtz(with_base=True)
    mtz.spacegroup = spacegroup
    mtz.set_cell_for_all(cell)
    mtz.set_data(hkl)

    size = mtz.get_size_for_hkl()

    grid = gemmi.Int8Grid(*size)

    grid.spacegroup = spacegroup
    grid.set_unit_cell(cell)

    return grid


def copy_grid(grid):
    """
    Create a new grid with the same parameters, but zero everywhere
    :param grid:
    :return:
    """
    new_grid = gemmi.Int8Grid(grid.nu, grid.nv, grid.nw)

    new_grid.spacegroup = grid.spacegroup

    new_grid.set_unit_cell(grid.unit_cell)

    return new_grid


def double_grid(grid):
    """
    Return a grid with unit cell lengths and number of indexes doubled (effectively eight grids stacked in a cuboid)
    :param grid:
    :return:
    """
    new_grid = gemmi.Int8Grid(2 * grid.nu, 2 * grid.nv, 2 * grid.nw)

    new_grid.spacegroup = grid.spacegroup

    cell = grid.unit_cell
    new_cell = gemmi.UnitCell(2 * cell.a, 2 * cell.b, 2 * cell.c, cell.alpha, cell.beta, cell.gamma)

    new_grid.set_unit_cell(new_cell)

    return new_grid


def get_symmetry_mask(structure, grid, radius):
    """
    Basically looks for every symmetry atom and masks around it in a new grid

    :param structure:
    :param grid:
    :param radius:
    :return:
    """
    symmetry_mask = copy_grid(grid)

    symmetries = list(grid.spacegroup.operations())
    if Constants.DEBUG: print(f"Found: {len(symmetries)} symmetry operations")

    for model in structure:
        for chain in model:
            for residue in chain.get_polymer():
                for atom in residue:
                    pos = atom.pos
                    fractional = grid.unit_cell.fractionalize(pos)

                    for operation in symmetries[1:]:
                        sym_pos = operation.apply_to_xyz([fractional.x, fractional.y, fractional.z])

                        fractional_sym_pos = gemmi.Fractional(*sym_pos)

                        orthogonal_sym_pos = grid.unit_cell.orthogonalize(fractional_sym_pos)

                        symmetry_mask.set_points_around(orthogonal_sym_pos, radius=radius, value=1)

    return symmetry_mask


def get_cell_mask(structure, grid, radius: float):
    """
    Partitions protein atoms by the operation to put them into the unit cell, and generates
    masks for each of these
    :param structure:
    :param grid:
    :param radius:
    :return:
    """

    def partition_protein_atoms(structure, grid):
        partitioning = {}

        for model in structure:
            for chain in model:
                for residue in chain.get_polymer():
                    for atom in residue:
                        pos = atom.pos
                        fractional = grid.unit_cell.fractionalize(pos)
                        wrapped_fraction = fractional.wrap()
                        orthogonal_wrapped_fractional = grid.unit_cell.orthogonalize(wrapped_fraction)

                        x_trans = fractional.x // 1
                        y_trans = fractional.y // 1
                        z_trans = fractional.z // 1

                        block = (x_trans, y_trans, z_trans)

                        if block in partitioning:
                            partitioning[block].append(orthogonal_wrapped_fractional)

                        else:
                            partitioning[block] = []
                            partitioning[block].append(orthogonal_wrapped_fractional)

        return partitioning

    partitioning = partition_protein_atoms(structure, grid)

    partition_grids = {}
    for partitioning_key, partition_atoms in partitioning.items():

        # Setup grids
        doubled_grid = double_grid(grid)
        partition_grids[partitioning_key] = copy_grid(grid)

        # Mask in the doubled grid
        for atom_pos in partition_atoms:
            doubled_grid.set_points_around(atom_pos, radius=radius, value=1)

        # Pull out the unit cell component
        doubled_grid_array = np.array(doubled_grid, copy=False)
        grid_array = np.array(partition_grids[partitioning_key], copy=False)

        grid_size = (grid.nu, grid.nv, grid.nw)
        grid_array[:, :, :] = doubled_grid_array[:grid_size[0] - 1, :grid_size[1] - 1, :grid_size[2] - 1]

    # Combine the grids
    cell_mask = copy_grid(grid)
    cell_mask_array = np.array(cell_mask, copy=False)

    for partitioning_key, partition_grid in partition_grids.items():
        grid_array = np.array(partition_grid, copy=False)
        cell_mask_array = cell_mask_array + grid_array

    cell_mask_array[cell_mask_array <= 1] = 0
    cell_mask_array[cell_mask_array != 0] = 1

    return cell_mask


def get_protein_mask(structure, grid, radius):
    """
    Return a grid set to 1 at all points within radius angstroms of protein backbone atoms
    :param structure:
    :param grid:
    :param radius:
    :return:
    """

    protein_mask = copy_grid(grid)
    for model in structure:
        for chain in model:
            for residue in chain.get_polymer():
                for atom in residue:
                    pos = atom.pos
                    protein_mask.set_points_around(pos, radius=radius, value=1)

    return protein_mask


def combine_masks(symmetry_mask, cell_mask, protein_mask):
    """
    Return a mask which is 1 wherever the protein mask overlaps the cell image mask or the smmetry mask
    :param symmetry_mask: equals 1 where a symmetry atom is nearby
    :param cell_mask: equals 1 where a unit cell protein image is nearby
    :param protein_mask: equals 1 where a protein atom is nearby
    :return:
    """

    combined_mask = copy_grid(protein_mask)
    combined_array = np.array(combined_mask, copy=False)

    symmetry_array = np.array(symmetry_mask, copy=False)
    cell_array = np.array(cell_mask, copy=False)
    protein_array = np.array(protein_mask, copy=False)

    # Mask the overlap of points near the protein and points near the symmetry atoms
    combined_array[(protein_array == 1) & (symmetry_array == 1)] = 1
    # Mask the overlap of points near the protein and protein images from other cells
    combined_array[(protein_array == 1) & (cell_array == 1)] = 1

    return combined_mask


def get_contact_mask(structure, grid, radius):
    """
    Get the contact masks
    :param structure:
    :param grid:
    :param radius:
    :return:
    """
    symmetry_mask = get_symmetry_mask(structure, grid, radius)
    cell_mask = get_cell_mask(structure, grid, radius)
    protein_mask = get_protein_mask(structure, grid, radius)

    contact_mask = combine_masks(symmetry_mask, cell_mask, protein_mask)

    return contact_mask


def get_overlap(residue, contact_mask):
    """
    Return the fraction of residue atoms inside the contact mask
    :param residue:
    :param contact_mask:
    :return:
    """
    vals = []
    for atom in residue:
        pos = atom.pos
        val = contact_mask.interpolate_point(pos)
        vals.append(val)

    return sum(vals) / len(vals)


def get_overlap_dict(residues, contact_mask):
    """
    Get the dictionary of overlaps of selected residue with contacts
    :param structure:
    :param selection:
    :param radius:
    :return:
    """

    overlaps = {}
    for resid, residue in residues.items():
        overlap = get_overlap(residue, contact_mask)
        overlaps[resid] = overlap

    return overlaps


def write_out_file(pdb_path, out_file, selection, radius, overlap_dict):
    """
    Write a json containing all the overlaps and parameters
    :param pdb_path:
    :param out_file:
    :param selection:
    :param radius:
    :param overlap_dict:
    :return:
    """
    out_dict = {
        "pdb_path": str(pdb_path),
        "selection": str(selection),
        "radius": str(radius),
        "overlaps": {
            resid: overlap
            for resid, overlap
            in overlap_dict.items()
        }
    }

    with open(str(out_file), "w") as f:
        json.dump(out_dict, f)


def write_ccp4_mask(grid, file):
    ccp4 = gemmi.Ccp4Mask()
    ccp4.grid = grid
    ccp4.update_ccp4_header(0, True)
    ccp4.setup()
    ccp4.write_ccp4_map(str(file))


def get_contact_score(pdb_path, out_path=None, selection="LIG", radius=3.0, write_maps=True):
    """
    Get the score of the selection as a fraction inside of contact regions
    :param pdb_path:
    :param out_path:
    :param selection:
    :param radius:
    :return:
    """
    pdb_path = Path(pdb_path)
    if Constants.DEBUG: print(f"The supplied pdb path is: {pdb_path}")
    out_path = Path(out_path)
    if Constants.DEBUG: print(f"The supplied output  path is: {out_path}")

    out_file = get_out_file(out_path, pdb_path)
    if Constants.DEBUG: print(f"The json that output will be written to is: {out_file}")

    structure = get_structure(pdb_path)

    residues = get_residues(structure, selection)
    if Constants.DEBUG: print(f"Found {len(residues)} residue named {selection} to check for contacts")

    grid = get_grid(structure)
    if Constants.DEBUG: print(
        f"Found a grid for the structure with shape: {(grid.nu, grid.nv, grid.nw)}; spacegroup: {grid.spacegroup}; unit cell {grid.unit_cell}")

    symmetry_mask = get_symmetry_mask(structure, grid, radius)
    cell_mask = get_cell_mask(structure, grid, radius)
    protein_mask = get_protein_mask(structure, grid, radius)

    contact_mask = combine_masks(symmetry_mask, cell_mask, protein_mask)
    if write_maps:
        symmetry_mask_file = out_path.with_name("symmetry_mask.cpp4")
        cell_mask_file = out_path.with_name("cell_mask.ccp4")
        contact_mask_file = out_path.with_name("contact_mask.ccp4")
        print(
            (
                f"Writing ed maps to:\n"
                f"\tSymmetry mask: {symmetry_mask_file}\n"
                f"\tCell mask: {cell_mask_file}\n"
                f"\tContact mask: {contact_mask_file}\n"
            )
        )
        write_ccp4_mask(symmetry_mask, symmetry_mask_file)
        write_ccp4_mask(cell_mask, cell_mask_file)
        write_ccp4_mask(contact_mask, contact_mask_file)

    overlap_dict = get_overlap_dict(residues, contact_mask)

    write_out_file(pdb_path, out_file, selection, radius, overlap_dict)


if __name__ == "__main__":
    fire.Fire(get_contact_score)
