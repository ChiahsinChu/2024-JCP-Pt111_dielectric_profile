import glob
import os

from ase import io
import MDAnalysis as mda

from intdielec.workflow.ikkem import DualIterElecEps


data_dir = os.environ.get("DATADIR")
if data_dir is None:
    raise ValueError("Please set env var DATADIR.")

dname = os.path.join(data_dir, "mlmd")

# timestep in ps
dt = 100

atoms = io.read(os.path.join(dname, "coord.xyz"))
u = mda.Universe(
    os.path.join(dname, "coord.xyz"),
    os.path.join(dname, "dump.xtc"),
    topology_format="XYZ",
    format="XTC",
    in_memory=True,
    in_memory_step=2000,
)
for ii, ts in enumerate(u.trajectory):
    work_dir = "calc/{:0>10.3f}/pbc/".format(ii * dt)
    atoms.set_positions(ts.positions)
    os.makedirs(work_dir, exist_ok=True)
    io.write(os.path.join(work_dir, "coord.xyz"), atoms)

dnames = glob.glob("calc/*")
dnames.sort()
# print(dnames)
for dname in dnames:
    task = DualIterElecEps(work_dir=dname)
    task.workflow(configs="param.json")
