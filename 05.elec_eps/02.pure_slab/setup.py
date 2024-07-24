import glob
import numpy as np

from ase import io, build

from toolbox.io.cp2k import Cp2kInput


delta_vs = np.arange(-1.0, 1.1, 1.0)

atoms = build.fcc111("Pt", (6, 6, 6), 3.97605, periodic=True, vacuum=20.0)
cell_params = atoms.cell.cellpar()
cell_params[-1] = 120
atoms.set_cell(cell_params)
for ii, delta_v in enumerate(delta_vs):
    cell_params = atoms.cell.cellpar()
    fp_params = {
        "FORCE_EVAL": {
            "STRESS_TENSOR": "NONE",
            "DFT": {
                "POISSON": {
                    "POISSON_SOLVER": "IMPLICIT",
                    "IMPLICIT": {
                        "BOUNDARY_CONDITIONS": "MIXED_PERIODIC",
                        "DIRICHLET_BC": {
                            "AA_PLANAR": [
                                {
                                    "V_D": 0.0,
                                    "PARALLEL_PLANE": "XY",
                                    "X_XTNT": "0.0 %.4f" % atoms.cell[0][0],
                                    "Y_XTNT": "0.0 %.4f" % atoms.cell[1][1],
                                    "INTERCEPT": 5.0,
                                    "PERIODIC_REGION": ".TRUE.",
                                },
                                {
                                    "V_D": delta_v,
                                    "PARALLEL_PLANE": "XY",
                                    "X_XTNT": "0.0 %.4f" % atoms.cell[0][0],
                                    "Y_XTNT": "0.0 %.4f" % atoms.cell[1][1],
                                    "INTERCEPT": cell_params[2] - 5.0,
                                    "PERIODIC_REGION": ".TRUE.",
                                },
                            ]
                        },
                    },
                }
            },
            "PRINT": {"STRESS_TENSOR": {"_": "OFF"}},
        }
    }
    work_dir = "calc/ref.%03d" % ii
    # EXTENDED_FFT_LENGTHS
    cp2k_inp = Cp2kInput(
        atoms, hartree=True, totden=True, eden=True, extended_fft_lengths=True
    )
    cp2k_inp.write(work_dir, fp_params=fp_params)


fnames = glob.glob("../01.calc/calc/*/ref_lo/coord.xyz")
fnames.sort()
for ii, delta_v in enumerate(delta_vs):
    for jj, fname in enumerate(fnames):
        _atoms = io.read(fname)
        atoms = _atoms[_atoms.symbols == "Pt"]
        atoms.set_cell(cell_params)
        atoms.set_pbc(True)
        atoms.center()

        fp_params = {
            "FORCE_EVAL": {
                "STRESS_TENSOR": "NONE",
                "DFT": {
                    "POISSON": {
                        "POISSON_SOLVER": "IMPLICIT",
                        "IMPLICIT": {
                            "BOUNDARY_CONDITIONS": "MIXED_PERIODIC",
                            "DIRICHLET_BC": {
                                "AA_PLANAR": [
                                    {
                                        "V_D": 0.0,
                                        "PARALLEL_PLANE": "XY",
                                        "X_XTNT": "0.0 %.4f" % atoms.cell[0][0],
                                        "Y_XTNT": "0.0 %.4f" % atoms.cell[1][1],
                                        "INTERCEPT": 5.0,
                                        "PERIODIC_REGION": ".TRUE.",
                                    },
                                    {
                                        "V_D": delta_v,
                                        "PARALLEL_PLANE": "XY",
                                        "X_XTNT": "0.0 %.4f" % atoms.cell[0][0],
                                        "Y_XTNT": "0.0 %.4f" % atoms.cell[1][1],
                                        "INTERCEPT": cell_params[2] - 5.0,
                                        "PERIODIC_REGION": ".TRUE.",
                                    },
                                ]
                            },
                        },
                    }
                },
                "PRINT": {"STRESS_TENSOR": {"_": "OFF"}},
            }
        }

        work_dir = "calc/task.%03d.%06d" % (ii, jj)
        cp2k_inp = Cp2kInput(
            atoms,
            hartree=True,
            totden=True,
            eden=True,
            extended_fft_lengths=True,
            wfn_restart="calc/ref.%03d/cp2k-RESTART.wfn" % ii,
        )
        cp2k_inp.write(work_dir, fp_params=fp_params, save_dict=True)
