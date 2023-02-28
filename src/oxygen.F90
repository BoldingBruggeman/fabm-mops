#include "fabm_driver.h"

module mops_oxygen

   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_mops_oxygen
      type (type_dependency_id) :: id_bgc_theta, id_bgc_salt
      type (type_surface_dependency_id) :: id_bgc_atmosp, id_bgc_wind, id_bgc_seaice
      type (type_state_variable_id) :: id_oxy
   contains
      ! Model procedures
      procedure :: initialize
      procedure :: do_surface
   end type type_mops_oxygen

contains

   subroutine initialize(self, configunit)
      class (type_mops_oxygen), intent(inout), target :: self
      integer,                  intent(in)            :: configunit

      ! Register diagnostic variables
      call self%register_state_variable(self%id_oxy, 'c', 'mmol O2/m3', 'concentration')

      ! Register environmental dependencies
      call self%register_dependency(self%id_bgc_theta, standard_variables%temperature)
      call self%register_dependency(self%id_bgc_salt, standard_variables%practical_salinity)
      call self%register_dependency(self%id_bgc_atmosp, standard_variables%surface_air_pressure)
      call self%register_dependency(self%id_bgc_wind, standard_variables%wind_speed)
      call self%register_dependency(self%id_bgc_seaice, standard_variables%ice_area_fraction)

      self%dt = 86400.0_rk
   end subroutine

   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_mops_oxygen), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: surf_oxy, bgc_theta, bgc_salt, bgc_atmosp, bgc_wind, bgc_seaice, vgas660, o2gasex

      _SURFACE_LOOP_BEGIN_
         _GET_(self%id_oxy, surf_oxy)
         _GET_(self%id_bgc_theta, bgc_theta)
         _GET_(self%id_bgc_salt, bgc_salt)
         _GET_SURFACE_(self%id_bgc_atmosp, bgc_atmosp)
         _GET_SURFACE_(self%id_bgc_wind, bgc_wind)
         _GET_SURFACE_(self%id_bgc_seaice, bgc_seaice)

         bgc_atmosp = bgc_atmosp / 101325.0_rk   ! from Pa to atm

         vgas660=(0.337_rk*bgc_wind**2)*0.24_rk*(1.0_rk-bgc_seaice)
         CALL O2_SURFFORCING(vgas660,bgc_atmosp,bgc_theta,bgc_salt, &
             surf_oxy,o2gasex)
         _ADD_SURFACE_FLUX_(self%id_oxy, o2gasex)
      _SURFACE_LOOP_END_
   end subroutine do_surface

   SUBROUTINE O2_SURFFORCING(vgas660,atmosp,ttemp,stemp,soxy,o2ex)

! define Schmidt no. coefficients for O2
! based on Keeling et al [GBC, 12, 141, (1998)]
      real(rk), parameter :: sox1 = 1638.0_rk
      real(rk), parameter :: sox2 = -81.83_rk
      real(rk), parameter :: sox3 = 1.483_rk
      real(rk), parameter :: sox4 = -0.008004_rk

! coefficients for determining saturation O2
      real(rk), parameter :: oA0=  2.00907_rk
      real(rk), parameter :: oA1=  3.22014_rk
      real(rk), parameter :: oA2=  4.05010_rk
      real(rk), parameter :: oA3=  4.94457_rk
      real(rk), parameter :: oA4= -2.56847e-1_rk
      real(rk), parameter :: oA5=  3.88767_rk
      real(rk), parameter :: oB0= -6.24523e-3_rk
      real(rk), parameter :: oB1= -7.37614e-3_rk
      real(rk), parameter :: oB2= -1.03410e-2_rk
      real(rk), parameter :: oB3= -8.17083e-3_rk
      real(rk), parameter :: oC0= -4.88682e-7_rk

      real(rk), intent(in) :: vgas660,atmosp,ttemp,stemp,soxy
      real(rk), intent(out) :: o2ex

! local coefficients

      real(rk) :: SchmidtNoO2,aTT,aTK,aTS,aTS2,aTS3,aTS4,aTS5, &
            oCnew,o2s,O2sat,Kwexch

      SchmidtNoO2=sox1+sox2*ttemp+sox3*ttemp*ttemp &
         + sox4*ttemp*ttemp*ttemp

      KWexch = vgas660/sqrt(SchmidtNoO2/660.0_rk)

! Determine saturation O2
! using Garcia and Gordon (1992), L&O (mistake in original???)

      aTT  = 298.15_rk -ttemp
      aTK  = 273.15_rk +ttemp
      aTS  = log(aTT/aTK)
      aTS2 = aTS*aTS
      aTS3 = aTS2*aTS
      aTS4 = aTS3*aTS
      aTS5 = aTS4*aTS
      oCnew  = oA0+oA1*aTS+oA2*aTS2+oA3*aTS3+oA4*aTS4+oA5*aTS5 &
         + stemp*(oB0+oB1*aTS+oB2*aTS2+oB3*aTS3)+oC0*(stemp*stemp)
      o2s = EXP(oCnew)

! Convert from ml/l to mmol/m^3
! Note: o2 in mit is in mol/m3; I use mmol/m3, thus coonvert with 1d3
      O2sat = o2s/22391.6_rk * 1000.0_rk*1000.0_rk

! Determine flux, inc. correction for local atmos surface pressure
      o2ex = Kwexch*(atmosp*O2sat-soxy)

   END SUBROUTINE

end module mops_oxygen
