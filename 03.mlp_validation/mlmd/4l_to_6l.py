import numpy as np
from ase import io, build


atoms = io.read("../aimd_prl/coord.xyz")

# Pt slab
slab = build.fcc111("Pt", size=(6, 6, 6), a=3.976, vacuum=10.5, periodic=True)
cellpar = slab.cell.cellpar()
cellpar[-1] = 120.
slab.set_cell(cellpar)
slab.translate([0, 0, slab.cell[2, 2] / 2])
slab.wrap()

# water box
water = atoms[np.logical_not(atoms.symbols == "Pt")]
water.set_cell(slab.cell)
water.set_pbc(True)
water.center(axis=2)
# print(water)
# io.write("water.xyz", water)

new_atoms = water + slab
io.write("coord.xyz", new_atoms)
io.write("system.data", new_atoms, format="lammps-data", specorder=["O", "H", "Pt"])


