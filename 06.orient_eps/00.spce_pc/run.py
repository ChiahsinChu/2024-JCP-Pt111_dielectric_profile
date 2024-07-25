import os
from ase import io
import numpy as np

import MDAnalysis as mda
from MDAnalysis import transformations as trans

from intdielec.watanalysis.dielectric import InverseDielectricConstant as IDC


data_dir = os.environ.get("DATADIR")
if data_dir is None:
    raise ValueError("Please set env var DATADIR.")

dname = os.path.join(data_dir, "mlmd")

atoms = io.read(os.path.join(dname, "coord.xyz"))
topo = os.path.join(dname, "topo.data")
traj = os.path.join(dname, "dump.xtc")
u = mda.Universe(topo, traj, topology_format="DATA", format="XTC")
transform = trans.boxdimensions.set_dimensions(atoms.cell.cellpar())
u.trajectory.add_transformations(transform)
# print(u.dimensions)

ag = u.select_atoms("type 1 or type 2")
# print(ag)

metal_ids = np.where(atoms.symbols == "Pt")[0]
surf_ids = [metal_ids[-36:], metal_ids[:36]]
task = IDC(
    atomgroups=ag,
    bin_width=0.1,
    surf_ids=surf_ids,
    temperature=330,
    img_plane=0.935,
)
task.run()
task.save()
