&FORCE_EVAL
  METHOD ${METHOD}
#  STRESS_TENSOR ANALYTICAL

  &MIXED
    GROUP_PARTITION 112 28
    MIXING_TYPE GENMIX
    &GENERIC
      # X: Energy force_eval 2
      # Y: Energy force_eval 3 (restraining water)
      MIXING_FUNCTION X+Y
      VARIABLES X Y
    &END
  &END MIXED

  &DFT
    BASIS_SET_FILE_NAME ./BASIS_MOLOPT
    POTENTIAL_FILE_NAME ./GTH_POTENTIALS 
#    CHARGE 1
#    UKS  T
#    WFN_RESTART_FILE_NAME ./pre-opt.wfn
    &MGRID
      CUTOFF 400
      NGRIDS 5
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-13
      EXTRAPOLATION ASPC 
      EXTRAPOLATION_ORDER 0
    &END QS
#    &SCF 
#       SCF_GUESS  RESTART 
#       EPS_SCF 3.0E-7 
#       MAX_SCF 500 
#       ADDED_MOS 500 
#       CHOLESKY INVERSE 
#       &SMEAR  ON 
#         METHOD FERMI_DIRAC 
#         ELECTRONIC_TEMPERATURE [K] 300 
#       &END SMEAR 
#       &DIAGONALIZATION 
#          ALGORITHM STANDARD 
#       &END DIAGONALIZATION 
#       &MIXING 
#           METHOD BROYDEN_MIXING 
#           ALPHA   0.3 
#           BETA    1.5 
#           NBROYDEN  14 
#       &END 
#       &PRINT 
#         &RESTART 
#           &EACH 
#               QS_SCF 50 
#           &END 
#           ADD_LAST NUMERIC 
#         &END 
#       &END 
#     &END SCF 

    &SCF
      SCF_GUESS RESTART 
      EPS_SCF 1.0E-7
      MAX_SCF 50
      MAX_SCF_HISTORY 5
      CHOLESKY INVERSE_DBCSR
      &OUTER_SCF
         EPS_SCF 1.0E-7
         MAX_SCF 20 
      &END
      &OT
         MINIMIZER DIIS 
         PRECOND_SOLVER INVERSE_UPDATE
         PRECONDITIONER FULL_SINGLE_INVERSE
         STEPSIZE 0.01
#         ENERGY_GAP 0.1
      &END
      &PRINT
        &RESTART
          &EACH
            MD 20 
          &END EACH
        &END RESTART
      &END PRINT
    &END SCF

    &XC
      &XC_FUNCTIONAL PBE 
      &END XC_FUNCTIONAL
#      &XC_GRID
#         XC_SMOOTH_RHO NN50
#         XC_DERIV NN50_SMOOTH
#      &END XC_GRID
       &VDW_POTENTIAL
         DISPERSION_FUNCTIONAL PAIR_POTENTIAL
         &PAIR_POTENTIAL
            TYPE DFTD3
            PARAMETER_FILE_NAME ./dftd3.dat
            CALCULATE_C9_TERM .TRUE.
            REFERENCE_C9_TERM .TRUE.
            R_CUTOFF [angstrom] 12.0
            REFERENCE_FUNCTIONAL PBE
         &END PAIR_POTENTIAL
      &END VDW_POTENTIAL
    &END XC
    &PRINT
#      &E_DENSITY_CUBE
#        ADD_LAST NUMERIC
#        &EACH
#          GEO_OPT 0
#        &END EACH
#      &END E_DENSITY_CUBE
#      &PDOS
#        COMPONENTS
#        NLUMO -1 
#        ADD_LAST  NUMERIC
#        &EACH
#          MD 50
#          GEO_OPT 0
#        &END EACH
#        COMMON_ITERATION_LEVELS 0
#        APPEND T
#      &END PDOS

      &MO_CUBES
        NHOMO -1
        NLUMO -1
        WRITE_CUBE F
        ADD_LAST NUMERIC
        &EACH
          MD 50
#          GEO_OPT 0
        &END EACH
      &END MO_CUBES

#      &V_HARTREE_CUBE ON
#         STRIDE 8 8 1
#         APPEND T
#         &EACH
#           MD 50
#          GEO_OPT 0
#         &END EACH
#        ADD_LAST NUMERIC
#      &END V_HARTREE_CUBE

   &END PRINT
  &END DFT

  &MM
    &FORCEFIELD
      DO_NONBONDED F
      &BOND
        ATOMS   H O
        KIND    CHARMM
        K       0.1
        R0      1.89
      &END BOND
      &CHARGE
        ATOM  H
        CHARGE 0
      &END CHARGE
      &CHARGE
        ATOM  O
        CHARGE 0
      &END CHARGE
      &CHARGE
        ATOM  F
        CHARGE 0
      &END CHARGE
      &CHARGE
        ATOM  Na
        CHARGE 0
      &END CHARGE
      &CHARGE
        ATOM  Pt
        CHARGE 0
      &END CHARGE
    &END FORCEFIELD

    &POISSON
      &EWALD
        EWALD_TYPE spme
        ALPHA .35
        GMAX 24
        O_SPLINE 6
      &END EWALD
    &END POISSON
  &END MM


  &SUBSYS
    &CELL
#      A  9.611  -5.549   0.000
#      B  0   11.098  0
#      C  0   0   21.796
      ABC [angstrom]  11.246   11.246   27.887     
      ALPHA_BETA_GAMMA  90 90 120
    &END CELL
    &COORD
       @INCLUDE ./coord.xyz
    &END COORD

    &TOPOLOGY
#      COORD_FILE_NAME ./coor
#      COORD_FILE_FORMAT XYZ
      CONN_FILE_NAME ./water-dump-1.psf
      CONNECTIVITY PSF
#      &DUMP_PSF
#      &END DUMP_PSF
#      &GENERATE
#         REORDER T
#      &END GENERATE
    &END TOPOLOGY

    &KIND O
      BASIS_SET DZVP-MOLOPT-SR-GTH 
      POTENTIAL GTH-PBE-q6 
    &END KIND
    &KIND Pt 
      BASIS_SET  DZVP-A5-Q10-323-MOL-T1-DERIVED_SET-1 
      POTENTIAL  GTH-PBE-q10
    &END KIND
    &KIND H
      BASIS_SET DZVP-MOLOPT-SR-GTH
      POTENTIAL GTH-PBE-q1
#      MASS 2.014
    &END KIND
  &END SUBSYS
&END FORCE_EVAL

