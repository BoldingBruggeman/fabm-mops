module mops_shared
   use fabm_types, only: rk
   real(rk), parameter :: vsafe = 1.0e-6_rk
   real(rk), parameter :: rcp = 117.0_rk       !redfield ratio C:P
   real(rk), parameter :: rnp = 16.0_rk        !redfield ratio N:P
   real(rk), parameter :: ro2ut = 167.0_rk !redfield -O2:P ratio
!   real(rk), parameter :: ro2ut = 151.13958_rk !redfield -O2:P ratio
! VS changed value of ro2ut for comparison with PETSC3.13 version in
! blogin/scratch/usr/shmvolki/models/TMM-MOPS/sim-single/MIT28/sim-01-02-2024
! (experiment OBS_NARR, Kriest et al., 2017)
! can I rather set this value in my "*.yaml" file, too ...?
   real(rk), parameter :: rhno3ut = 0.8_rk*ro2ut - rnp ! -HNO3:P ratio for denitrification
   real(rk), parameter :: bgc_dt = 1.0_rk      !max BGC timestep (d)
   real(rk), parameter :: convert_mol_to_mmol=1000.0_rk
   real(rk), parameter :: rho0=1024.5_rk
   real(rk), parameter :: permil=1.0_rk/rho0
   real(rk), parameter :: permeg=1.0e-6_rk
   real(rk), parameter :: alimit = 1.0d-3
   ! VS an aggregate variable for the total detritus production by plankton
   ! is to be used to calculate implicit CaCO3 divergences and their effect on DIC and Alk
   type (type_interior_standard_variable), parameter :: detritus_production_by_plankton = type_interior_standard_variable(name='detritus_production_by_plankton',units='mmol P/m3/d',aggregate_variable=.true.) 
end module
