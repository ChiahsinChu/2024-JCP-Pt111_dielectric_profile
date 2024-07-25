import os
from ase import io
import numpy as np

import MDAnalysis as mda


data_dir = os.environ.get("DATADIR")
if data_dir is None:
    raise ValueError("Please set env var DATADIR.")

dname = os.path.join(data_dir, "mlmd")

atoms = io.read(os.path.join(dname, "coord.xyz"))
metal_ids = np.where(atoms.symbols == "Pt")[0]
surf_ids = [metal_ids[-36:], metal_ids[:36]]

topo = os.path.join(dname, "topo.data")
traj = os.path.join(dname, "dump.xtc")
u = mda.Universe(topo, traj, topology_format="DATA", format="XTC")

z1 = 0
z2 = 0
for ts in u.trajectory:
    coord = ts.positions
    z1 += coord[surf_ids[0], 2].mean()
    z2 += coord[surf_ids[1], 2].mean()
z1 /= len(u.trajectory)
z2 /= len(u.trajectory)
print(z1, z2)
    