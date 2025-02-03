from pymatgen.core import Molecule


# Function to check if the atom is within the unit cell
def in_unit_cell(frac_coord, supercell=2):
    return min(frac_coord) >= 0 and max(frac_coord) <= (
        (1 / supercell) + (0.01 / supercell)
    )


def structure_to_unit(structure):
    sites_to_keep = []
    structure = structure.make_supercell(2)
    for site in structure:
        if in_unit_cell(site.frac_coords):
            sites_to_keep.append(site)
    species = [s.specie for s in sites_to_keep]
    coords = [s.coords for s in sites_to_keep]
    unit = Molecule(species=species, coords=coords)

    return unit
