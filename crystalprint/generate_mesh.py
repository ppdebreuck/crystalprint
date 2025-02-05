import argparse
import os

import numpy as np
import trimesh
import yaml
from pymatgen.analysis.molecule_structure_comparator import CovalentRadius
from pymatgen.core import Structure
from pymatgen.io.xyz import XYZ
from trimesh.transformations import rotation_matrix

from .utils import structure_to_unit

CYLINDER_SMOOTHNESS = 32
SPHERE_SMOOTHNESS = 3
Z_EPSILON = 0


# Parse command-line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generate 3D mesh from CIF file with bonds and atoms."
    )
    parser.add_argument(
        "file_path", type=str, help="Path to the input XYZ or CIF file."
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=4,
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

    parser.add_argument(
        "--color",
        action="store_true",  # This makes it a flag
        default=False,  # Default value if the flag is not provided
        help="Enable color output",  # Help message for the flag
    )

    args = parser.parse_args()
    out_dir = os.path.abspath(args.out_dir)
    if not os.path.exists(out_dir):
        print(f"Warning: The directory '{out_dir}' does not exist.")

    args.out_dir = out_dir

    return args


with open("ElementColorSchemes.yaml") as fp:
    color_map = yaml.safe_load(fp)["Jmol"]
print(color_map)


# Function to add a bond with atoms
def add_bond(
    meshes,
    coords1,
    coords2,
    atom_radius,
    cylinder_diam,
    radius1,
    radius2,
    rank1="",
    rank2="",
    color=False,
):
    sphere1 = trimesh.primitives.Sphere(
        radius=radius1, center=coords1, subdivisions=SPHERE_SMOOTHNESS
    )
    sphere2 = trimesh.primitives.Sphere(
        radius=radius2, center=coords2, subdivisions=SPHERE_SMOOTHNESS
    )

    direction = coords2 - coords1
    height = np.linalg.norm(direction)
    direction = direction / height

    midpoint = (
        (coords2 - direction * radius2) - (coords1 + direction * radius1)
    ) * 0.5 + (coords1 + direction * radius1)

    height1 = np.linalg.norm(midpoint - coords1)
    height2 = np.linalg.norm(midpoint - coords2)

    shift1 = (midpoint + coords1) / 2
    shift2 = (midpoint + coords2) / 2

    cylinder_diam = min(radius1, radius2) / 3

    bond_cylinder1 = trimesh.primitives.Cylinder(
        radius=cylinder_diam, height=height1, sections=CYLINDER_SMOOTHNESS
    )
    bond_cylinder2 = trimesh.primitives.Cylinder(
        radius=cylinder_diam, height=height2, sections=CYLINDER_SMOOTHNESS
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
        bond_cylinder1.apply_transform(rotation_mat)
        bond_cylinder2.apply_transform(rotation_mat)

    bond_cylinder1.apply_translation(shift1)
    bond_cylinder2.apply_translation(shift2)

    if color:
        color1 = color_map[rank1]
        color2 = color_map[rank2]
        sphere1.visual.vertex_colors = np.tile(color1, (len(sphere1.vertices), 1))
        sphere2.visual.vertex_colors = np.tile(color2, (len(sphere2.vertices), 1))
        bond_cylinder1.visual.vertex_colors = np.tile(
            color1, (len(bond_cylinder1.vertices), 1)
        )
        bond_cylinder2.visual.vertex_colors = np.tile(
            color2, (len(bond_cylinder1.vertices), 1)
        )

    meshes.append((rank1, sphere1))
    meshes.append((rank2, sphere2))
    meshes.append((rank1, bond_cylinder1))
    meshes.append((rank2, bond_cylinder2))


# Main script execution
def main():
    args = parse_arguments()

    file_path = args.file_path

    if args.file_path.endswith(".cif"):
        unit = structure_to_unit(Structure.from_file(file_path))
    elif args.file_path.endswith(".xyz"):
        unit = XYZ.from_file(file_path).molecule
    else:
        raise ValueError("Unsupported file extension. Only .cif and .xyz are allowed.")

    BOND_CUTOFF = args.cutoff
    ATOM_RADIUS = args.atom_radius
    CYLINDER_DIAM = args.cylinder_diam

    SCALE_FACTOR = 1
    meshes = []
    unique_bonds = set()

    for i, site1 in enumerate(unit):
        neighbors = unit.get_neighbors(site1, BOND_CUTOFF)
        coords1 = site1.coords
        print(type(site1.species.elements[0]))
        # radius1 = site1.species.elements[0].atomic_radius_calculated * (3 / 4)
        radius1 = CovalentRadius.radius[site1.species.elements[0].symbol] * (4 / 5)
        for neighbor in neighbors:
            coords2 = neighbor.coords
            # radius2 = neighbor.species.elements[0].atomic_radius_calculated * (3 / 4)
            radius2 = CovalentRadius.radius[neighbor.species.elements[0].symbol] * (
                4 / 5
            )
            bond_id = tuple((coords1 + coords2) / 2)
            if bond_id not in unique_bonds:
                unique_bonds.add(bond_id)
                add_bond(
                    meshes,
                    coords1,
                    coords2,
                    ATOM_RADIUS,
                    CYLINDER_DIAM,
                    radius1,
                    radius2,
                    rank1=str(site1.species.elements[0].symbol),
                    rank2=str(neighbor.species.elements[0].symbol),
                    color=args.color,
                )

    # Combine meshes and apply transformations
    meshes = sorted(meshes, key=lambda x: x[0])
    meshes = [m[1] for m in meshes]
    final_mesh = trimesh.util.concatenate(meshes)
    final_mesh.fill_holes()
    final_mesh.apply_scale(10 * SCALE_FACTOR)

    # saving to file
    out_path = os.path.join(
        args.out_dir, os.path.splitext(os.path.basename(file_path))[0]
    )
    final_mesh.export(out_path + ".stl")
    if args.color:
        final_mesh.export(out_path + ".obj")
        # final_mesh.export(out_path + ".3mf")

    print("Done")
    final_mesh.show()


# Ensure that this runs when invoked via `python -m mypackage`
if __name__ == "__main__":
    main()
