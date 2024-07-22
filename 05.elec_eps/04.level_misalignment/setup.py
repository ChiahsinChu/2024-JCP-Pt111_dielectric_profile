from toolbox.utils import *
from toolbox.utils.unit import *
from toolbox.io.cp2k import Cp2kInput


def make_fparam(atoms):
    fp_params = {"FORCE_EVAL": {"DFT": {"PRINT": {"PDOS": {"LDOS": []}}}}}
    oxygen_ids = np.where(atoms.symbols == "O")[0] + 1
    hydrogen_ids = np.where(atoms.symbols == "H")[0] + 1
    for ii in range(len(oxygen_ids)):
        fp_params["FORCE_EVAL"]["DFT"]["PRINT"]["PDOS"]["LDOS"].append(
            {
                "LIST": "%d %d %d"
                % (
                    oxygen_ids[ii],
                    hydrogen_ids[2 * ii],
                    hydrogen_ids[2 * ii + 1],
                )
            }
        )
    return fp_params


# PBC
atoms = io.read("pbc/coord.xyz")
fp_params = make_fparam(atoms)
task = Cp2kInput(
    atoms,
    hartree=True,
    eden=True,
    totden=True,
)
task.write(output_dir="pbc", fp_params=fp_params)


# dipole correction
atoms = io.read("dip_cor/coord.xyz")
fp_params = make_fparam(atoms)
task = Cp2kInput(
    atoms,
    hartree=True,
    eden=True,
    totden=True,
    dip_cor=True,
)
task.write(output_dir="dip_cor", fp_params=fp_params)
