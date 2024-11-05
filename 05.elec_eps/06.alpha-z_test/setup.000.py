"""
Build model: chemisorbed water @ Pt(111)
"""

import os

from ase import io, build
import numpy as np

from toolbox.io.cp2k import Cp2kInput


def build_model(d):
    assert d >= 1.0
    slab = build.fcc111(
        "Pt",
        size=(4, 4, 6),
        a=3.976,
        vacuum=15.0,
        orthogonal=True,
        periodic=True,
    )
    shift_vector = slab[-6].position.copy()
    shift_vector[2] = 0.0

    adsorbate = build.molecule("H2O")
    adsorbate.rotate(110, "y", center="COM")
    build.add_adsorbate(slab, adsorbate, d, "ontop")

    # shift water to the center of box
    coords = slab.get_positions()
    coords[-3:] += shift_vector
    slab.set_positions(coords)

    slab.set_pbc(True)
    return slab


def run(ii, fp_params, d0, delta_d, n):
    for jj in range(n):
        work_dir = "calc.000/task.%03d.%06d" % (ii, jj)
        atoms = build_model(d0 + jj * delta_d)

        cp2k_inp = Cp2kInput(
            atoms,
            hartree=True,
            eden=True,
            totden=True,
            wfn_restart=os.path.join(work_dir, "cp2k-RESTART.wfn"),
        )
        cp2k_inp.write(work_dir, fp_params=fp_params)


if __name__ == "__main__":
    # pzc
    fp_params = {
        "FORCE_EVAL": {
            "STRESS_TENSOR": "NONE",
            "DFT": {
                "LOCALIZE": {
                    "EPS_LOCALIZATION": 1e-02,
                    "RESTART": ".TRUE.",
                    "PRINT": {
                        "LOC_RESTART": {},
                        "WANNIER_CENTERS": {
                            "FILENAME": "=wannier.xyz",
                        },
                        "WANNIER_SPREADS": {
                            "FILENAME": "=wannier_spread.out",
                        },
                    },
                },
            }
        }
    }
    run(0, fp_params, 1.5, 0.2, 10)
    run(2, fp_params, 4.0, 0.5, 6)
    
    # add bias
    atoms = build_model(d=2.0)
    cell = atoms.get_cell()
    fp_params["FORCE_EVAL"]["DFT"]["POISSON"] = {
        "POISSON_SOLVER": "IMPLICIT",
        "IMPLICIT": {
            "BOUNDARY_CONDITIONS": "MIXED_PERIODIC",
            "DIRICHLET_BC": {
                "AA_PLANAR": [
                    {
                        "V_D": 0.0,
                        "PARALLEL_PLANE": "XY",
                        "X_XTNT": "0.0 %.4f" % cell[0, 0],
                        "Y_XTNT": "0.0 %.4f" % cell[1, 1],
                        "INTERCEPT": 5.0,
                        "PERIODIC_REGION": ".TRUE.",
                    },
                    {
                        "V_D": 5.0,
                        "PARALLEL_PLANE": "XY",
                        "X_XTNT": "0.0 %.4f" % cell[0, 0],
                        "Y_XTNT": "0.0 %.4f" % cell[1, 1],
                        "INTERCEPT": cell[2, 2] - 5.0,
                        "PERIODIC_REGION": ".TRUE.",
                    },
                ]
            },
        },
    }
    run(1, fp_params, 1.5, 0.2, 10)
    run(3, fp_params, 4.0, 0.5, 6)
