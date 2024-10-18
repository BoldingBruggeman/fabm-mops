module mops_shared
   use fabm_types, only: rk
   real(rk), parameter :: vsafe = 1.0e-6_rk
   real(rk), parameter :: rcp = 117.0_rk       !redfield ratio C:P
   real(rk), parameter :: rnp = 16.0_rk        !redfield ratio N:P
   real(rk), parameter :: ro2ut = 151.13958_rk !redfield -O2:P ratio
   real(rk), parameter :: rhno3ut = 0.8_rk*ro2ut - rnp ! -HNO3:P ratio for denitrification
!   real(rk), parameter :: bgc_dt = 1.0_rk   !VS CORRECT FOR 1 STEP PER DAY
   real(rk), parameter :: bgc_dt = 0.5_rk   !VS CORRECT FOR 2 STEPS PER DAY
!   real(rk), parameter :: bgc_dt = 0.0625_rk   !VS CORRECT FOR 16 STEPS PER DAY
   real(rk), parameter :: convert_mol_to_mmol=1000.0_rk
   real(rk), parameter :: rho0=1024.5_rk
   real(rk), parameter :: permil=1.0_rk/rho0
   real(rk), parameter :: permeg=1.0e-6_rk
   real(rk), parameter :: alimit = 1.0d-3   
end module
