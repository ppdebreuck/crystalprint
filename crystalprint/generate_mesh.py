import argparse
import os

import numpy as np
import trimesh
from pymatgen.core import Structure
from trimesh.transformations import rotation_matrix

CYLINDER_SMOOTHNESS = 32
SPHERE_SMOOTHNESS = 3
Z_EPSILON = 0


# Parse command-line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generate 3D mesh from CIF file with bonds and atoms."
    )
    parser.add_argument("cif_path", type=str, help="Path to the input CIF file.")
    parser.add_argument(
        "--bond_cutoff",
        type=float,
        default=2.9,
        help="Bond cutoff distance. Default is 2.9.",
    )
    parser.add_argument(
        "--atom_radius",
        type=float,
        default=0.65,
        help="Radius of atoms. Default is 0.65.",
    )
    parser.add_argument(
        "--cylinder_diam",
        type=float,
        default=0.3,
        help="Diameter of the bond cylinder. Default is 0.3.",
    )
    parser.add_argument(
        "--out_dir",
        type=str,
        default=".",
        help="Directory for the output files.",
    )

    args = parser.parse_args()
    out_dir = os.path.abspath(args.out_dir)
    if not os.path.exists(out_dir):
        print(f"Warning: The directory '{out_dir}' does not exist.")

    args.out_dir = out_dir

    return args


# Function to add a bond with atoms
def add_bond(meshes, coords1, coords2, atom_radius, cylinder_diam):
    sphere1 = trimesh.primitives.Sphere(
        radius=atom_radius, center=coords1, subdivisions=SPHERE_SMOOTHNESS
    )
    sphere2 = trimesh.primitives.Sphere(
        radius=atom_radius, center=coords2, subdivisions=SPHERE_SMOOTHNESS
    )

    meshes.append(sphere1)
    meshes.append(sphere2)

    midpoint = (coords1 + coords2) / 2
    direction = coords1 - coords2
    height = np.linalg.norm(direction)
    direction = direction / height

    bond_cylinder = trimesh.primitives.Cylinder(
        radius=cylinder_diam, height=height, sections=CYLINDER_SMOOTHNESS
    )

    initial_direction = np.array([0, 0, 1])
    angle = np.arccos(np.dot(initial_direction, direction))
    axis = np.cross(initial_direction, direction)

    if axis.sum() == 0:
        print("---")
        print(initial_direction, direction)
        print("---")
        print(axis, angle)
        print("---")

    if (
        not np.isnan(angle)
        and angle != 0
        and not np.isnan(axis).any()
        and np.absolute(axis).sum() != 0
    ):
        print(axis, angle)
        axis = axis / np.linalg.norm(axis)
        rotation_mat = rotation_matrix(angle, axis)
        bond_cylinder.apply_transform(rotation_mat)

    bond_cylinder.apply_translation(midpoint)
    meshes.append(bond_cylinder)


# Function to check if the atom is within the unit cell
def in_unit_cell(frac_coord, supercell=2):
    return min(frac_coord) >= 0 and max(frac_coord) <= (1 / supercell)


# Main script execution
def main():
    args = parse_arguments()

    # Load CIF file and create structure
    cif_path = args.cif_path
    structure = Structure.from_file(cif_path)
    structure = structure.make_supercell(2)

    BOND_CUTOFF = args.bond_cutoff
    ATOM_RADIUS = args.atom_radius
    CYLINDER_DIAM = args.cylinder_diam

    SCALE_FACTOR = 1
    meshes = []
    unique_bonds = set()

    for i, site1 in enumerate(structure):
        neighbors = structure.get_neighbors(site1, BOND_CUTOFF)
        coords1 = site1.coords
        for neighbor in neighbors:
            coords2 = neighbor.coords
            if in_unit_cell(neighbor.frac_coords) and in_unit_cell(site1.frac_coords):
                bond_id = tuple((coords1 + coords2) / 2)
                if bond_id not in unique_bonds:
                    unique_bonds.add(bond_id)
                    add_bond(meshes, coords1, coords2, ATOM_RADIUS, CYLINDER_DIAM)

    # Combine meshes and apply transformations
    final_mesh = trimesh.util.concatenate(meshes)
    final_mesh.fill_holes()
    final_mesh.apply_scale(10 * SCALE_FACTOR)

    REPEATS = 2  # Number of repetitions for translation (can be adjusted)
    translated_meshes = []
    for i in range(REPEATS):
        for j in range(REPEATS):
            translated_mesh = final_mesh.copy()
            translation = np.array(
                [
                    i * structure.lattice.a * 10 * SCALE_FACTOR,
                    j * structure.lattice.b * 10 * SCALE_FACTOR,
                    0,
                ]
            )
            translated_mesh.apply_translation(translation)
            translated_meshes.append(translated_mesh)

    final_repeated_mesh = trimesh.util.concatenate(translated_meshes)

    # 3 mesh
    extended_translated_meshes = []
    for i in range(REPEATS):
        translated_mesh = final_repeated_mesh.copy()
        translation = np.array(
            [
                0,
                0,
                i * (structure.lattice.c / 2 + (2 * ATOM_RADIUS)) * 10 * SCALE_FACTOR
                + Z_EPSILON,
            ]
        )
        translated_mesh.apply_translation(translation)
        extended_translated_meshes.append(translated_mesh)
    # Combine all translated meshes
    final_extended_repeated_mesh = trimesh.util.concatenate(extended_translated_meshes)

    # saving to file
    out_path = os.path.join(
        args.out_dir, os.path.splitext(os.path.basename(cif_path))[0]
    )
    final_mesh.export(out_path + ".stl")
    final_repeated_mesh.export(out_path + "_repeated_2d.stl")
    final_extended_repeated_mesh.export(out_path + "_repeated_3d.stl")

    print(
        f"STL file with bonds and atoms has been created: {out_path}.stl and {out_path}_repeated.stl"
    )
    final_mesh.show()


# Ensure that this runs when invoked via `python -m mypackage`
if __name__ == "__main__":
    main()
