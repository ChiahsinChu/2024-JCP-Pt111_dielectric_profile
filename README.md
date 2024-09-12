# to-do list

- [ ] finish doc for `IntDielec`
- [ ] add `requirements.txt` for python codes (`pipreqs myproj/`)

# README

This repository contains the data and codes for the publication "Dielectric profile at the Pt(111)/water interface" (link), mainly including:

- training of a DP model (MLP) for Pt(111)/water interface
- validation of the DP model (via potential energy surface, thermodynamic properties, and dynamics properties)
- calculation of electronic and orientational dielectric profiles at the Pt(111)/water interface

Details for the files are shown in the following sections. The readers can reproduce the figures in the publication (including the main text and the supplementary information) with `pub_figs/make_figs.ipynb`.

## DOI

[placeholder]

## Title

Dielectric profile at the Pt(111)/water interface

## Authors

**Jia-Xin Zhu\***, **Jun Cheng\*** and **Katharina Doblhoff-Dier\***

## Contact e-mail

- jiaxinzhu@stu.xmu.edu.cn
- chengjun@xmu.edu.cn
- k.doblhoff-dier@lic.leidenuniv.nl

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

3. install required packages

```bash
git clone --recursive https://github.com/ChiahsinChu/IntDielec.git
cd IntDielec && pip install .
cd ..

git clone https://github.com/ChiahsinChu/mdanalysis.git -b jxzhu_dev
cd mdanalysis && pip install package/
cd ..

git clone https://github.com/ChiahsinChu/WatAnalysis.git
cd WatAnalysis && pip install .
cd ..
```

## `00.init_dataset`

## `01.dpgen`

Run workflow to explore the configuration space. Please adjust `param.json` according to the publication.

- software: [`dpgen`](https://github.com/deepmodeling/dpgen)
- version: `0.11.0`

```bash
dpgen run param.json machine.json
```

1. `param.json` is only applicable for the first 16 iterations. The users are expected to update the keyword `fp_task_max` as 100 and pad the `model_devi_jobs` for the following iterations to reproduce the routine shown in the publication.
2. The users are expected to update the `configs/POSCAR_*` to change the number of initial configurations for the exploration.

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

   Output files are in `data` directory. See `README.md` in the directory for details.

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

   Every directory in `calc` represents the calculation for a single snapshots. The usage of the workflow and the structures of the output data for a single snapshot is refered to [here](https://github.com/ChiahsinChu/IntDielec/blob/main/doc/elec_eps.md).

2. Perform calculation for electronic dielectric profiles at the Pt(111)/water interfaces.

   ```bash
   cd 01.eps_calc
   python run.py
   ```

   Every directory in `calc` represents the calculation for a single snapshots. The usage of the workflow and the structures of the output data for a single snapshot is refered to [here](https://github.com/ChiahsinChu/IntDielec/blob/main/doc/elec_eps.md).

3. Perform calculation for electronic dielectric profiles at the Pt(111)/vacuum interfaces.

   ```bash
   cd 02.pure_slab
   python setup.py
   python cp2k_calc.py
   ```

   Every directory in `calc` represents the calculation at a certain boundary condition. The post-processing of data for the dielectric profiles is finished by `data/run.py` and the output files are in `data` directory. The data structures are refered to [here](https://github.com/ChiahsinChu/IntDielec/blob/main/doc/elec_eps.md).

4. Perform anaylsis for ``idealized" eps

   ```bash
   cd 03.water_polarizability
   python run.py
   ```

   Every directory represents the analyses for a single snapshots. The prefixes of `lo` and `hi` refers to the lower and the upper interfaces in the snapshots.

   - `*_grid.npy`: grid in the z direction
   - `*_int_eps.npy`: interfacial dielectric profile
   - `*_wat_chi.npy`: water susceptibility defined as (interface dielectric profile) - (metal dielectric profiles)
   - `*_water_eden.npy`: valance electron density in water as defined in the publication
   - `*_wat_norm_chi.npy`: water susceptibility normalized by the valance electron density in water

5. Perform calculation for water HOMO z-distribution

   ```bash
   cd 04.level_misalignment
   python setup.py
   cd pbc
   sbatch cp2k.slurm
   cd ../dip_cor
   sbatch cp2k.slurm
   ```

   The post-processing of data is shown in `pub_figs/make_figs.ipynb`.

6. Calculate probability distribution of E_vac

   ```bash
   cd 05.efield_vac_distribution
   python run.py
   ```

7. Calculate displacements in the maximal localized Wannier centers (MLWCs) induced by E-fields (i.e., orbital polarizability) at different metal-water distances

   ```bash
   cd 06.alpha-z_test
   python setup.000.py
   python setup.001.py
   ```

   - `calc.000` for a water monomer at a Pt(111) slab with an ``O-down" configuration, which is similar to the chemisorbed water
   - `calc.001` for a water monomer at a Pt(111) slab with a ``H-down" configuration, which is similar to the physisorbed water
   - `calc.*/task.000.*` and `calc.*/task.002.*` for PBC
   - `calc.*/task.001.*` and `calc.*/task.003.*` for DBC

8. Calculate molecular orbitals for the maximal localized Wannier centers (MLWCs) cloest to the Pt(111) surface (`07.wannier_mo`)

   - `water_a` for a water monomer at a Pt(111) slab with an ``O-down" configuration, which is similar to the chemisorbed water
   - `water_b` for a water monomer at a Pt(111) slab with a ``H-down" configuration, which is similar to the physisorbed water
   - `water_a_ref` for a water monomer with the same coordinates as in `water_a` but without the Pt(111) surface
   - `water_b_ref` for a water monomer with the same coordinates as in `water_b` but without the Pt(111) surface

9. Water density profile for MLMD as reference

   ```python
   cd 08.water_density
   python run.py
   ```

   Read `output.pkl` get a dict, which contains keywords `rho_water` and `geo_dipole_water`. The data structures are identical with the output data in `03.mlp_validation/data/water_structure*.pkl`.

## `06.orient_eps`

- `00.spce`

  - Calculate inverse of orientational (ionic) dielectric constant via MLMD + SPCE charges (point charges).
  - The output data are stored in `inveps.pkl`, the structures of which are refered to [here](https://github.com/ChiahsinChu/IntDielec/blob/main/doc/orient_eps.md).

- `01.dw_gaussian`

  - Calculate inverse of orientational (ionic) dielectric constant via MLMD + Wannier centroid from Deep Wannier (DW) model with spread from DFT calculation.
  - The output data are stored in `inveps.pkl`, the structures of which are refered to [here](https://github.com/ChiahsinChu/IntDielec/blob/main/doc/orient_eps.md).
  - DW model are downloaded from [aissquare](https://www.aissquare.com/datasets/detail?pageType=datasets&name=H2O-DPLR&id=17).

- `calc_z_surf.py`
  - Calculate surface average positions in MLMD

<!-- ## `07.manuscript` -->

## `softwares`

- `dp2.0-cp2k.tar.gz`: CP2K for MLMD simulations
- `IntDielec.tar.gz`: IntDielec package
- `mdanalysis.tar.gz`: MDAnalysis package
- `WatAnalysis.tar.gz`: WatAnalysis package
