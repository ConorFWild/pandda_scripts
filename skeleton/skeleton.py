from pathlib import Path

import gemmi

import fire
import json

import numpy as np

from constants import Constants


def get_structure(structure_path: Path):
    structure: gemmi.Structure = gemmi.read_structure(str(structure_path))
    return structure


def get_xmap(event_map_file: Path):
    m = gemmi.read_ccp4_map(str(event_map_file))
    m.setup()

    grid_array = np.array(m.grid, copy=True)

    new_grid = gemmi.FloatGrid(*grid_array.shape)
    new_grid.spacegroup = m.grid.spacegroup  # gemmi.find_spacegroup_by_name("P 1")
    new_grid.set_unit_cell(m.grid.unit_cell)

    new_grid_array = np.array(new_grid, copy=False)
    new_grid_array[:, :, :] = grid_array[:, :, :]

    return new_grid


def select(structure, selection):
    sel = gemmi.parse_cid(selection)
    selected_atoms = []
    for model in sel.models(structure):
        for chain in sel.chains(model):
            for residue in sel.residues(chain):
                for atom in sel.atoms(residue):
                    selected_atoms.append(atom)

    return selected_atoms


def get_sample_positions(structure, atoms, xmap, radius=1.8):
    ns = gemmi.NeighborSearch(structure[0], structure.cell, radius).populate()

    positions = []
    for atom in atoms:
        pos = atom.pos
        positions.append(pos)
        ref = np.array([pos.x, pos.y, pos.z])
        for mark in ns.find_atoms(atom.pos,):
            moving = np.array([mark.x, mark.y, mark.z])

            if np.allclose(ref, moving):
                continue

            half_way = (ref + moving) / 2

            half_way_pos = gemmi.Position(half_way[0], half_way[1], half_way[2])

            positions.append(half_way_pos)

    fractional_positions = [xmap.unitcell.fractionalize(position) for position in positions]

    return fractional_positions


def sample(positions, xmap):
    samples = []
    for position in positions:
        sample = xmap.interpolate_value(position)
        samples.append(sample)

    return samples


def skeleton_score(structure, atoms, xmap, threshold=1.0):
    # Get positions of atoms and bonds in order to sample at them
    sample_positions = get_sample_positions(structure, atoms, xmap)
    print(f"Got {len(sample_positions)} samples")
    print(f"Got positions: {sample_positions}")

    # Get samples of electron density
    samples = sample(sample_positions, xmap)
    print(f"Got samples: {samples}")

    # Get the skeleton score from samples
    score = len([x for x in samples if x > threshold]) / len(samples)
    print(f"Got score: {score}")

    return score


def write_output(structure_path,
                 event_map_path,
                 selection,
                 score: float,
                 file: Path,
                 ):
    out_dict = {
        "structure_path": str(structure_path),
        "event_map_path": str(event_map_path),
        "selection": str(selection),
        "skeleton_score": score,
    }

    with open(str(file), "w") as f:
        json.dump(out_dict, f)


def skeleton(structure: str, event_map: str, out_dir: str, selection: str = "(LIG)"):
    # Type the input arguments
    structure_path = Path(structure)
    event_map_path = Path(event_map)
    out_dir = Path(out_dir)

    # Load the model
    structure = get_structure(structure_path)

    # Load the map
    xmap = get_xmap(event_map_path)

    # Perform selection
    atoms = select(structure, selection)

    # Calculate the skeleton score
    score = skeleton_score(structure, atoms, xmap)

    # Write output
    write_output(structure_path,
                 event_map_path,
                 selection,
                 score,
                 out_dir / Constants.SKELETON_SCORE_FILE,
                 )


if __name__ == "__main__":
    fire.Fire(skeleton)
