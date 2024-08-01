# README

## to-do list

- [ ] add `requirements.txt` for python codes (`pipreqs myproj/`)

[info of publication]

## Before start...

1. Download the extra large files (e.g., trajectories) and set env var `DATADIR`

```bash
wget -c [url]
tar -zxvf large_files.tar.gz
cd large_files
export DATADIR=$(pwd)
```

2. Set env var `PUBDIR`

```bash
git clone https://github.com/ChiahsinChu/2024-JPCC-Pt111_dielectric_profile.git
cd 2024-JPCC-Pt111_dielectric_profile
export PUBDIR=$(pwd)
```

## `00.init_dataset`

## `01.dpgen`

Run workflow to explore the configuration space. Please adjust `param.json` according to the publication. 

- software: [`dpgen`](https://github.com/deepmodeling/dpgen)
- version: `0.11.0`

```bash
dpgen run param.json machine.json
```

## `02.mlp_train`

- software: [`deepmd-kit`](https://github.com/deepmodeling/deepmd-kit)
- version: `2.0.3`

```bash
dp train input.json
```

## `03.mlp_validation`

1. generate MLMD trajectory which will be compared with AIMD

   ```bash
   cd mlmd
   python make_config.py
   lmp_mpi -i input.lmp
   ```

2. Perform analyses

   ```bash
   mkdir data
   python test_aimd.py
   python test_mlmd.py
   ```

## `04.mlmd`

The MLMD simulations were performed with modified cp2k based on v9.0 (development version), which can be installed with the following commands [[ref]](https://wiki.cheng-group.net/wiki/software_installation/deepmd-kit/deepmd-kit_installation_191/#dp-cp2k):

```bash
git clone https://github.com/Cloudac7/cp2k.git -b deepmd_latest --recursive --depth=1
cd cp2k/tools/toolchain/
./install_cp2k_toolchain.sh --enable-cuda=no --with-deepmd=$deepmd_root --with-tfcc=$tensorflow_root --deepmd-mode=cuda --mpi-mode=no --with-libint=no --with-libxc=no --with-libxsmm=no
make -j 4 ARCH=local VERSION="ssmp sdbg"
```

The codes can also be found in `softwares/dp2.0-cp2k.tar.gz`.

The simulations can be run with the following commands:

```bash
module add deepmd/2.0
# only serial version is available!
cp2k.sopt input.inp
```

The output trajectory file is in `xyz` format, which can be converted to `xtc` format with the `xyz_to_xtc.py` script for further analysis.

## `05.elec_eps`

1. Perform calculation for water thickness test.

   ```bash
   cd 00.l_water_test
   cd 10A
   python run.py
   cd ../13A
   python run.py
   ```

2. Perform calculation for electronic dielectric constant.

   ```bash
   cd 01.eps_calc
   python run.py
   ```

3. Perform anaylsis for ``idealized" eps

   ```bash
   cd 03.water_polarizability
   python run.py
   ```

4. Perform calculation for water HOMO z-distribution

   ```bash
   cd 04.level_misalignment
   python setup.py
   cd pbc
   sbatch cp2k.slurm
   cd ../dip_cor
   sbatch cp2k.slurm
   ```

5. Calculate probability distribution of E_vac

   ```bash
   cd 05.efield_vac_distribution
   python run.py
   ```

6. Calculate maximal localized Wannier centers with different metal-water distances
   - `calc.000` for chemisorbed water
   - `calc.001` for physisorbed water
   - `calc.*/task.000.*` and `calc.*/task.002.*` for PBC
   - `calc.*/task.001.*` and `calc.*/task.003.*` for DBC

## `06.orient_eps`

- `00.spce`
   - Calculate inverse of orientational (ionic) dielectric constant via MLMD + SPCE charges (point charges).

- `01.dw_gaussian`
   - Calculate inverse of orientational (ionic) dielectric constant via MLMD + Wannier centroid from Deep Wannier (DW) model with spread from DFT calculation
   - DW model are downloaded from [aissquare](https://www.aissquare.com/datasets/detail?pageType=datasets&name=H2O-DPLR&id=17).

- `calc_z_surf.py`
   - Calculate surface average positions in MLMD

<!-- ## `07.manuscript` -->

<!--
## `trajs`

- `aimd.xyz`: AIMD trajectory of Pt(111)/water interface

## `softwares`

- `dp2.0-cp2k.tar.gz`: Deep Potential training code -->
