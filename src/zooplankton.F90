#include "fabm_driver.h"

module mops_zooplankton

   use fabm_types
   use mops_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_mops_zooplankton
      type (type_state_variable_id) :: id_c, id_phy, id_po4, id_din, id_oxy, id_det, id_dop, id_dic
      type (type_diagnostic_variable_id) :: id_graz

      real(rk) :: ACmuzoo, ACkphy, AClambdaz, AComniz, ACeff, graztodop, zlambda
   contains
      ! Model procedures
      procedure :: initialize
      procedure :: do
   end type

contains

   subroutine initialize(self, configunit)
      class (type_mops_zooplankton), intent(inout), target :: self
      integer,                       intent(in)            :: configunit

      call self%get_parameter(self%ACmuzoo, 'ACmuzoo', '1/d','max. grazing rate', default=1.893_rk)
      call self%get_parameter(self%ACkphy, 'ACkphy', 'mmol P/m3','half-saturation constant', default=SQRT(self%ACmuzoo/1.0_rk)/rnp)
      call self%get_parameter(self%ACeff, 'ACeff', '1','assimilation efficiency', default=0.75_rk)
      call self%get_parameter(self%graztodop, 'graztodop', '1','fraction of grazing that goes into DOP', default=0.0_rk)
      call self%get_parameter(self%AClambdaz, 'AClambdaz', '1/d','excretion', default=0.03_rk)
      call self%get_parameter(self%AComniz, 'AComniz', 'm3/(mmol P * day)','density dependent loss rate', default=4.548_rk)
      call self%get_parameter(self%zlambda, 'zlambda', '1/d','mortality', default=0.01_rk)

! VS nur kurz without minimum value to avoid clipping in TMM implementation
! (see Jorns mail on October 16, 2024)
      call self%register_state_variable(self%id_c, 'c', 'mmol P/m3', 'concentration') !, minimum=0.0_rk)

      call self%register_diagnostic_variable(self%id_graz, 'graz', 'mmol P/m3/d', 'grazing')

      call self%register_state_dependency(self%id_phy, 'phy', 'mmol P/m3', 'phytoplankton')
      call self%register_state_dependency(self%id_dop, 'dop', 'mmol P/m3', 'dissolved organic phosphorus')
      call self%register_state_dependency(self%id_det, 'det', 'mmol P/m3', 'detritus')
      call self%register_state_dependency(self%id_oxy, 'oxy', 'mmol O2/m3', 'oxygen')
      call self%register_state_dependency(self%id_din, 'din', 'mmol N/m3', 'dissolved inorganic nitrogen')
      call self%register_state_dependency(self%id_po4, 'pho', 'mmol P/m3', 'phosphate')
      call self%register_state_dependency(self%id_dic, 'dic', 'mmol C/m3', 'dissolved inorganic carbon')

      ! Register environmental dependencies
      call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_c)

      self%dt = 86400.0_rk
   end subroutine

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_mops_zooplankton), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: PHY, ZOO
      real(rk) :: graz0, graz, zooexu, zooloss
! VS nur kurz
      real(rk) :: term1, term2

      _LOOP_BEGIN_

      _GET_(self%id_phy, PHY)
      _GET_(self%id_c, ZOO)

       if(ZOO.gt.0.0_rk) then

         if(PHY.gt.0.0_rk) then

! Grazing of zooplankton, Holling III
           graz0=((self%ACmuzoo*(PHY*PHY))/(self%ACkphy*self%ACkphy+PHY*PHY))*ZOO
! Make sure not to graze more phytoplankton than available.
           graz = MIN(PHY,graz0*bgc_dt)/bgc_dt
         else !PHY < 0

           graz=0.0_rk

         endif !ZOO

! Zooplankton exudation
          zooexu = self%AClambdaz * ZOO

! Zooplankton mortality 
          zooloss = self%AComniz * ZOO * ZOO

       else !ZOO < 0

           graz   =0.0_rk
           zooexu = 0.0_rk
           zooloss = 0.0_rk

       endif !ZOO
!! VS nur kurz
!           print *, 'PHY in ZOO is ', PHY
!           print *, 'ZOO is ', ZOO
!           print *, 'ACkphy is ', self%ACkphy
!           print *, 'ACmuzoo is ', self%ACmuzoo
!           print *, 'graz0 is ', graz0 
!           print *, 'graz is ', graz
!           print *, ' '

       _SET_DIAGNOSTIC_(self%id_graz, graz)

! Collect all euphotic zone fluxes in these arrays.
! VS SETTING FLUXES TO ZERO (except for zooplankton, October 18, 2024)
        _ADD_SOURCE_(self%id_c, self%ACeff*graz-zooexu-zooloss)
        _ADD_SOURCE_(self%id_po4, zooexu)
        _ADD_SOURCE_(self%id_dop, self%graztodop*(1.0_rk-self%ACeff)*graz + self%graztodop*zooloss)
!        _ADD_SOURCE_(self%id_oxy, -zooexu*ro2ut)
        _ADD_SOURCE_(self%id_phy, -graz)
        _ADD_SOURCE_(self%id_det, (1.0_rk-self%graztodop)*(1.0_rk-self%ACeff)*graz + (1.0_rk-self%graztodop)*zooloss)
!! VS nur kurz
!        term1 = (1.0_rk-self%graztodop)*(1.0_rk-self%ACeff)*graz / 86400.0_rk 
!        print *, '(1.0_rk-self%graztodop)*(1.0_rk-self%ACeff)*graz / sec is', term1
!        term2 = (1.0_rk-self%graztodop)*zooloss / 86400.0_rk
!        print *, '(1.0_rk-self%graztodop)*zooloss / sec is', term2

        _ADD_SOURCE_(self%id_din, zooexu*rnp)
! VS nur kurz
!        print *, 'zooexu / sec is ', zooexu / 86400.0_rk

!        _ADD_SOURCE_(self%id_dic, zooexu*rcp)

         ZOO = MAX(ZOO - alimit*alimit, 0.0_rk)
         _ADD_SOURCE_(self%id_c, -self%zlambda*ZOO)
         _ADD_SOURCE_(self%id_dop, self%zlambda*ZOO)

      _LOOP_END_
   end subroutine do

end module mops_zooplankton
