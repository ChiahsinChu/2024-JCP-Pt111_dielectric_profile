from typing import Dict

import os
import numpy as np

from ase import io
import MDAnalysis as mda
from MDAnalysis import transformations as trans

from toolbox.utils.utils import save_dict

from WatAnalysis.waterstructure import WaterStructure, HBA
from WatAnalysis.waterdynamics import SP
from WatAnalysis.preprocess import make_selection
from WatAnalysis.hbonds.postprocess import count_by_time


def load_traj(dname):
    topo = os.path.join(dname, "topo.data")
    traj = os.path.join(dname, "dump.xtc")
    u = mda.Universe(topo, traj, topology_format="DATA", format="XTC")
    atoms = io.read(os.path.join(dname, "coord.xyz"))
    return u, atoms


def calc_water_structure(all_data: Dict):
    fname = "water_structure"

    kw = "mlmd"
    l_block = 20000
    n_block = 5
    for ii in range(1, n_block + 1):
        task, water_a_coverage, water_a_theta, water_b_coverage, water_b_theta = (
            _calc_water_structure(
                all_data[kw], start=ii * l_block, stop=(ii + 1) * l_block
            )
        )
        out_dict = {
            "rho_water": task.results.rho_water,
            "geo_dipole_water": task.results.geo_dipole_water,
            "water_a_coverage": water_a_coverage,
            "water_a_theta": water_a_theta,
            "water_b_coverage": water_b_coverage,
            "water_b_theta": water_b_theta,
        }
        save_dict(out_dict, "data/%s.%s.%03d.pkl" % (fname, kw, ii))


def _calc_water_structure(data: Dict, start=None, stop=None):
    metal_ids = np.where(data["atoms"].symbols == "Pt")[0]
    surf_ids = [metal_ids[-36:], metal_ids[:36]]
    task = WaterStructure(
        universe=data["universe"],
        surf_ids=surf_ids,
        oxygen_sel="type 1",
        hydrogen_sel="type 2",
        verbose=True,
    )
    task.run(start=start, stop=stop)

    water_a_coverage, water_a_theta = task.calc_sel_water([0.0, 2.7])
    water_b_coverage, water_b_theta = task.calc_sel_water([2.7, 4.5])
    return (
        task,
        water_a_coverage / 36 / 2,
        water_a_theta,
        water_b_coverage / 36 / 2,
        water_b_theta,
    )


def calc_water_sp(all_data: Dict):
    fname = "water_sp"

    kw = "mlmd"
    l_block = 20000
    n_block = 5
    for ii in range(1, n_block + 1):
        water_a_data, water_interface_data = _calc_water_sp(
            all_data[kw], start=ii * l_block, stop=(ii + 1) * l_block
        )
        out_dict = {
            "water_a": water_a_data,
            "water_a_and_b": water_interface_data,
        }
        save_dict(out_dict, "data/%s.%s.%03d.pkl" % (fname, kw, ii))


def _calc_water_sp(data: Dict, start=None, stop=None):
    metal_ids = np.where(data["atoms"].symbols == "Pt")[0]
    surf_ids = [metal_ids[-36:], metal_ids[:36]]

    task = SP(
        universe=data["universe"],
        surf_ids=surf_ids,
        sel_region=[0, 2.7],
        c_ag="type 1",
        verbose=True,
    )
    task.run(start=start, stop=stop, step=10, tau_max=20)
    water_a_data = np.concatenate([[task.tau_timeseries], [task.sp_timeseries]], axis=0)

    task = SP(
        universe=data["universe"],
        surf_ids=surf_ids,
        sel_region=[0, 4.5],
        c_ag="type 1",
        verbose=True,
    )
    task.run(start=start, stop=stop, step=10, tau_max=20)
    water_interface_data = np.concatenate(
        [[task.tau_timeseries], [task.sp_timeseries]], axis=0
    )
    return water_a_data, water_interface_data


def calc_water_hb(all_data: Dict):
    fname = "hb"

    # timestep in ps
    dt = 0.5 * 1e-3

    kw = "mlmd"
    l_block = 20000
    n_block = 5
    for ii in range(1, n_block + 1):
        wata_accpetor_hb, wata_donor_hb, watb_accpetor_hb, watb_donor_hb = (
            _calc_water_hb(all_data[kw], start=ii * l_block, stop=(ii + 1) * l_block)
        )

        np.save("data/%s.%s.%03d.wata_a.npy" % (fname, kw, ii), wata_accpetor_hb)
        output = count_by_time(wata_accpetor_hb, dt=dt)
        np.save("data/%s.%s.%03d.wata_a_count_by_time.npy" % (fname, kw, ii), output)

        np.save("data/%s.%s.%03d.wata_d.npy" % (fname, kw, ii), wata_donor_hb)
        output = count_by_time(wata_donor_hb, dt=dt)
        np.save("data/%s.%s.%03d.wata_d_count_by_time.npy" % (fname, kw, ii), output)

        np.save("data/%s.%s.%03d.watb_a.npy" % (fname, kw, ii), watb_accpetor_hb)
        output = count_by_time(watb_accpetor_hb, dt=dt)
        np.save("data/%s.%s.%03d.watb_a_count_by_time.npy" % (fname, kw, ii), output)

        np.save("data/%s.%s.%03d.watb_d.npy" % (fname, kw, ii), watb_donor_hb)
        output = count_by_time(watb_donor_hb, dt=dt)
        np.save("data/%s.%s.%03d.watb_d_count_by_time.npy" % (fname, kw, ii), output)


def _calc_water_hb(data: Dict, start=None, stop=None):
    metal_ids = np.where(data["atoms"].symbols == "Pt")[0]
    surf_ids = [metal_ids[-36:], metal_ids[:36]]

    # water A
    selection = make_selection(sel_region=[0, 2.7], surf_ids=surf_ids, c_ag="type 1")
    ## water A as acceptor
    hba = HBA(
        universe=data["universe"],
        donors_sel="type 1",
        hydrogens_sel="type 2",
        acceptors_sel=selection,
        d_h_cutoff=1.2,
        d_a_cutoff=3,
        d_h_a_angle_cutoff=150,
        update_acceptors=True,
    )
    hba.run(start=start, stop=stop)
    wata_accpetor_hb = hba.results.hbonds.copy()
    ## water A as donor
    hba = HBA(
        universe=data["universe"],
        donors_sel=selection,
        hydrogens_sel="type 1",
        acceptors_sel="type 2",
        d_h_cutoff=1.2,
        d_a_cutoff=3,
        d_h_a_angle_cutoff=150,
        update_donors=True,
    )
    hba.run(start=start, stop=stop)
    wata_donor_hb = hba.results.hbonds.copy()

    # water B
    selection = make_selection(sel_region=[2.7, 4.5], surf_ids=surf_ids, c_ag="type 1")
    ## water B as acceptor
    hba = HBA(
        universe=data["universe"],
        donors_sel="type 1",
        hydrogens_sel="type 2",
        acceptors_sel=selection,
        d_h_cutoff=1.2,
        d_a_cutoff=3,
        d_h_a_angle_cutoff=150,
        update_acceptors=True,
    )
    hba.run(start=start, stop=stop)
    watb_accpetor_hb = hba.results.hbonds.copy()
    ## water B as donor
    hba = HBA(
        universe=data["universe"],
        donors_sel=selection,
        hydrogens_sel="type 1",
        acceptors_sel="type 2",
        d_h_cutoff=1.2,
        d_a_cutoff=3,
        d_h_a_angle_cutoff=150,
        update_donors=True,
    )
    hba.run(start=start, stop=stop)
    watb_donor_hb = hba.results.hbonds.copy()
    return wata_accpetor_hb, wata_donor_hb, watb_accpetor_hb, watb_donor_hb


if __name__ == "__main__":
    u_mlmd, atoms_mlmd = load_traj("mlmd")
    transform = trans.boxdimensions.set_dimensions(atoms_mlmd.cell.cellpar())
    u_mlmd.trajectory.add_transformations(transform)

    all_data = {
        "mlmd": {"universe": u_mlmd, "atoms": atoms_mlmd},
    }

    calc_water_structure(all_data)
    calc_water_sp(all_data)
    calc_water_hb(all_data)
