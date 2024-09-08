import os
import numpy as np

from ase import io
import MDAnalysis as mda
from MDAnalysis import transformations as trans

from toolbox.utils.utils import save_dict

from WatAnalysis.waterstructure import WaterStructure


if __name__ == "__main__":
    data_dir = os.environ.get("DATADIR")
    if data_dir is None:
        raise ValueError("Please set env var DATADIR.")

    dname = os.path.join(data_dir, "mlmd")
    topo = os.path.join(dname, "topo.data")
    traj = os.path.join(dname, "dump.xtc")
    u = mda.Universe(topo, traj, topology_format="DATA", format="XTC")
    atoms = io.read(os.path.join(dname, "coord.xyz"))

    workflow = []
    transform = trans.boxdimensions.set_dimensions(atoms.cell.cellpar())
    workflow.append(transform)
    ag = u.select_atoms("type 3")
    transform = mda.transformations.wrap(ag)
    workflow.append(transform)
    u.trajectory.add_transformations(*workflow)

    metal_ids = np.where(atoms.symbols == "Pt")[0]
    surf_ids = [metal_ids[-36:], metal_ids[:36]]
    task = WaterStructure(
        universe=u,
        surf_ids=surf_ids,
        oxygen_sel="type 1",
        hydrogen_sel="type 2",
        verbose=True,
    )
    task.run()

    out_dict = {
        "rho_water": task.results.rho_water,
        "geo_dipole_water": task.results.geo_dipole_water,
    }
    save_dict(out_dict, "output.pkl")
