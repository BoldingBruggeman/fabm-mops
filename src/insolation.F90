#include "fabm_driver.h"

module mops_insolation

   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_mops_insolation
      type (type_horizontal_dependency_id)       :: id_YC
      type (type_global_dependency_id)           :: id_Time
      type (type_surface_diagnostic_variable_id) :: id_sfac, id_stau
   contains
      ! Model procedures
      procedure :: initialize
      procedure :: do_surface
   end type type_mops_insolation

contains

   subroutine initialize(self, configunit)
      class (type_mops_insolation), intent(inout), target :: self
      integer,                      intent(in)            :: configunit

      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_sfac, 'sfac', 'W m-2', 'daily mean downwelling shortwave')
      call self%register_diagnostic_variable(self%id_stau, 'stau', 'd', 'day length')

      ! Register environmental dependencies
      call self%register_dependency(self%id_YC, standard_variables%latitude)
      call self%register_dependency(self%id_Time, standard_variables%number_of_days_since_start_of_the_year)
   end subroutine

   ! find shortwave radiation as function of date and latitude
   ! based on paltridge and parson
   ! modified by SPK from MITGCM
   ! adapted by JB for FABM
   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_mops_insolation), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

! !INPUT PARAMETERS: ===================================================
! time                 :: current time
       real(rk) :: time,Yc
       real(rk) :: daysperyear = 365.2425_rk

! !LOCAL VARIABLES: ====================================================
       real(rk) ::  solar, albedo
       real(rk) ::  dayfrac, yday, delta
       real(rk) ::  lat, sun1, dayhrs
       real(rk) ::  cosz, frac, fluxi, fracmin, sfac
!EOP

!
      solar = 1360._rk   !solar constant
      albedo = 0.6_rk    !planetary albedo
!      par = 0.4_rk       !photosynthetically reactive frac

!
! find day (****NOTE for year starting in winter*****)
        _GET_GLOBAL_(self%id_Time, Time)
        dayfrac=mod(Time,daysperyear) &
                         /(daysperyear)          !fraction of year
        yday = 2.0_rk*3.1416_rk*dayfrac                         !convert to radians
        delta = (0.006918_rk - (0.399912_rk*cos(yday)) &     !cosine zenith angle
               +(0.070257_rk*sin(yday))             &    !(paltridge+platt)
               -(0.006758_rk*cos(2.0_rk*yday)) &
               +(0.000907_rk*sin(2.0_rk*yday)) &
               -(0.002697_rk*cos(3.0_rk*yday)) &
               +(0.001480_rk*sin(3.0_rk*yday)) )

       _SURFACE_LOOP_BEGIN_
! latitude in radians
         _GET_SURFACE_(self%id_YC, YC)
          lat=YC/180.0_rk*3.1416_rk
          sun1 = -sin(delta)/cos(delta) * sin(lat)/cos(lat)
          if (sun1.le.-0.999_rk) sun1=-0.999_rk
          if (sun1.ge. 0.999_rk) sun1= 0.999_rk
          dayhrs = abs(acos(sun1))
          cosz = ( sin(delta)*sin(lat)+ &             !average zenith angle
                 (cos(delta)*cos(lat)*sin(dayhrs)/dayhrs) )
          if (cosz.le.0.005_rk) cosz=0.005_rk
          frac = dayhrs/3.1416_rk               !fraction of daylight in day
! daily average photosynthetically active solar radiation just below surface
          fluxi = solar*(1.0_rk-albedo)*cosz*frac
!! VS nur kurz
!          print *, 'solar is', solar
!          print *, 'albedo is', albedo
!          print *, 'cosz is', cosz
          print *, 'Time is ', Time
          print *, 'dayperyear is ', daysperyear
          print *, 'dayfrac is ', dayfrac
          print *, 'frac is ', frac
!
! convert to sfac
          if (fluxi.gt.0.0_rk) sfac=fluxi
! very large for polar night
          if (fluxi.lt.0.00001_rk) sfac=0.00001_rk
          _SET_SURFACE_DIAGNOSTIC_(self%id_sfac, sfac)
! daylength; ensure that it lies between 0 and 1 (may be slightly
! out of this range in high latitudes)
          fracmin = MIN(frac,1.0_rk)
          _SET_SURFACE_DIAGNOSTIC_(self%id_stau, MAX(fracmin,0.0_rk))

      _SURFACE_LOOP_END_
   end subroutine do_surface

end module mops_insolation
