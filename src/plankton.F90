#include "fabm_driver.h"

module mops_plankton

   use fabm_types
   use mops_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_mops_plankton
      type (type_dependency_id) :: id_bgc_theta, id_bgc_dz, id_ciz, id_att
      type (type_surface_dependency_id) :: id_bgc_tau
      type (type_state_variable_id) :: id_phy, id_zoo, id_po4, id_din, id_oxy, id_det, id_dop, id_dic
      type (type_diagnostic_variable_id) :: id_f1, id_f2, id_f6, id_chl

      real(rk) :: TempB, ACmuphy, ACik, ACkpo4, AClambda, AComni, plambda
      real(rk) :: ACmuzoo, ACkphy, AClambdaz, AComniz, ACeff, graztodop, zlambda
   contains
      ! Model procedures
      procedure :: initialize
      procedure :: do
   end type type_mops_plankton

   type (type_bulk_standard_variable), parameter :: total_chlorophyll = type_bulk_standard_variable(name='total_chlorophyll',units='mg/m^3',aggregate_variable=.true.)

contains

   subroutine initialize(self, configunit)
      class (type_mops_plankton), intent(inout), target :: self
      integer,                    intent(in)            :: configunit

      real(rk) :: ACkchl

      call self%get_parameter(self%TempB, 'TempB', 'degrees Celsius','ref. temperature for T-dependent growth', default=15.65_rk) 
      call self%get_parameter(self%ACmuphy, 'ACmuphy', '1/day','max. growth rate', default=0.6_rk) 
      call self%get_parameter(self%ACik, 'ACik', 'W/m2','light half-saturation constant', default=9.653_rk) 
      call self%get_parameter(self%ACkpo4, 'ACkpo4', 'mmol P/m3','half-saturation constant for PO4 uptake', default=0.4995_rk) 
      call self%get_parameter(ACkchl, 'ACkchl', '1/(m*mmol P/m3)','att. of Phy', default=0.03_rk*rnp) 
      call self%get_parameter(self%AClambda, 'AClambda', '1/day','exudation rate', default=0.03_rk) 
      call self%get_parameter(self%AComni, 'AComni', 'm3/(mmol P * day)','density dep. loss rate', default=0.0_rk) 
      call self%get_parameter(self%plambda, 'plambda', '1/d','phytoplankton mortality', default=0.01_rk) 

      call self%get_parameter(self%ACmuzoo, 'ACmuzoo', '1/d','max. grazing rate', default=1.893_rk)
      call self%get_parameter(self%ACkphy, 'ACkphy', 'mmol P/m3','zoo half-saturation constant', default=SQRT(self%ACmuzoo/1.0_rk)/rnp)
      call self%get_parameter(self%AClambdaz, 'AClambdaz', '1/d','zooplankton excretion', default=0.03_rk)
      call self%get_parameter(self%AComniz, 'AComniz', 'm3/(mmol P * day)','zooplankton mortality', default=4.548_rk)
      call self%get_parameter(self%ACeff, 'ACeff', '1','assimilation efficiency', default=0.75_rk)
      call self%get_parameter(self%graztodop, 'graztodop', '1','fraction grazing that goes into DOP', default=0.0_rk)
      call self%get_parameter(self%zlambda, 'zlambda', '1/d','zooplankton mortality', default=0.01_rk)

      call self%register_state_variable(self%id_phy, 'phy', 'mmol P/m3', 'phytoplankton', minimum=0.0_rk)
      call self%register_state_variable(self%id_zoo, 'zoo', 'mmol P/m3', 'zooplankton', minimum=0.0_rk)

      call self%register_diagnostic_variable(self%id_f1, 'f1', 'mmol P/m3/d', 'phytoplankton growth rate')
      call self%register_diagnostic_variable(self%id_f2, 'f2', 'mmol P/m3/d', 'zooplankton grazing')
      call self%register_diagnostic_variable(self%id_chl, 'chl', 'mg/m3/d', 'chlorophyll')

      call self%register_state_dependency(self%id_dop, 'dop', 'mmol P/m3', 'dissolved organic phosphorus')
      call self%register_state_dependency(self%id_det, 'det', 'mmol P/m3', 'detritus')
      call self%register_state_dependency(self%id_oxy, 'oxy', 'mmol O2/m3', 'oxygen')
      call self%register_state_dependency(self%id_din, 'din', 'mmol N/m3', 'dissolved inorganic nitrogen')
      call self%register_state_dependency(self%id_po4, 'po4', 'mmol P/m3', 'phosphate')
      call self%register_state_dependency(self%id_dic, 'dic', 'mmol C/m3', 'dissolved inorganic carbon')

      ! Register environmental dependencies
      call self%register_dependency(self%id_ciz, 'ciz', 'W m-2', 'PAR at top of the layer')
      call self%register_dependency(self%id_bgc_tau, 'tau', 'd', 'day length')
      call self%register_dependency(self%id_bgc_theta, standard_variables%temperature)
      call self%register_dependency(self%id_bgc_dz, standard_variables%cell_thickness)
      call self%register_dependency(self%id_att, standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux)

      call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, self%id_phy, scale_factor=ACkchl)
      call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_phy)
      call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_zoo)
      call self%add_to_aggregate_variable(total_chlorophyll, self%id_chl)

      self%dt = 86400.0_rk
   end subroutine

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_mops_plankton), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: bgc_theta, bgc_dz, ciz, bgc_tau, att, PO4, DIN, PHY, ZOO
      real(rk) :: tempscale, TACmuphy, TACik, atten, glbygd, flightlim, limnut, fnutlim, phygrow0, phygrow, phyexu, phyloss
      real(rk) :: graz0, graz, zooexu, zooloss
      real(rk) :: topo4

      _LOOP_BEGIN_

      _GET_(self%id_bgc_theta, bgc_theta)
      _GET_(self%id_bgc_dz, bgc_dz)
      _GET_(self%id_ciz, ciz)
      _GET_SURFACE_(self%id_bgc_tau, bgc_tau)
      _GET_(self%id_att, att)

      _GET_(self%id_po4, PO4)
      _GET_(self%id_din, DIN)
      _GET_(self%id_phy, PHY)
      _GET_(self%id_zoo, ZOO)

! temperature dependence of phytoplankton growth (Eppley)
! this affects the light-half-saturation constant via acik=acmuphy/alpha
       tempscale = EXP(bgc_theta/self%TempB)
       TACmuphy = self%ACmuphy*tempscale
       TACik = self%ACik*tempscale

! The light limitation function of phytoplankton.
! This function corresponds to Evans and Garcon, 1997.
! Note that the initial slope of the P-I curve, alpha, is ACMuPhy/ACIk
! flightlim thus gives the light limited growth rate, averaged over day 
! and layer, normalised by max. growth rate
       atten = att*bgc_dz !attenuation (dimensionless)
       glbygd = 2.0_rk*ciz/(TACik*bgc_tau)   ! 2 * G_L/G_D of EG97
       flightlim = bgc_tau/atten*(phi(glbygd)-phi(glbygd*exp(-atten)))

       if(PHY.gt.0.0_rk) then

         limnut = MIN(PO4,DIN/rnp)

         if(limnut.gt.vsafe) then

! The nutrient limitation of phytoplankton
           fnutlim = limnut/(self%ackpo4+limnut)

! The growth rate of phytoplankton: light*nutrient limitation.
           phygrow0 = TACmuphy*PHY*MIN(flightlim,fnutlim)

! Make sure not to take up more nutrients than available.
           phygrow = MIN(limnut,phygrow0*bgc_dt)/bgc_dt

         else !limnut < vsafe

           phygrow=0.0_rk

         endif !limnut

! The exudation of phytoplankton
         phyexu = self%AClambda * PHY

! Other losses of phytoplankton
         phyloss = self%AComni * PHY * PHY

       else !PHY < 0

         phygrow=0.0_rk
         phyexu =0.0_rk
         phyloss=0.0_rk

       endif !PHY

       if(ZOO.gt.0.0_rk) then

         if(PHY.gt.0.0_rk) then

! Grazing of zooplankton, Holling III
           graz0=self%ACmuzoo*PHY*PHY/(self%ACkphy*self%ACkphy+PHY*PHY)*ZOO

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

! Photosynthesis stored in this array for diagnostic purposes only.
       _SET_DIAGNOSTIC_(self%id_f1, phygrow)
       _SET_DIAGNOSTIC_(self%id_f2, graz)

! JB constant chlorophyll:carbon (IK pers comm 2023-03-07)
       _SET_DIAGNOSTIC_(self%id_chl, 50._rk * PHY)

! Collect all euphotic zone fluxes in these arrays.
        topo4 = -phygrow+zooexu
        _ADD_SOURCE_(self%id_po4, topo4)
        _ADD_SOURCE_(self%id_dop, self%graztodop*(1.0_rk-self%ACeff)*graz + self%graztodop*(phyexu+zooloss) + phyloss)
        _ADD_SOURCE_(self%id_oxy, (phygrow-zooexu)*ro2ut)
        _ADD_SOURCE_(self%id_phy, phygrow-graz-phyexu-phyloss)
        _ADD_SOURCE_(self%id_zoo, self%ACeff*graz-zooexu-zooloss)
        _ADD_SOURCE_(self%id_det, (1.0_rk-self%graztodop)*(1.0_rk-self%ACeff)*graz + (1.0_rk-self%graztodop)*(phyexu+zooloss))
        _ADD_SOURCE_(self%id_din, topo4*rnp)
        _ADD_SOURCE_(self%id_dic, topo4*rcp)

         PHY = MAX(PHY - alimit*alimit, 0.0_rk)
         ZOO = MAX(ZOO - alimit*alimit, 0.0_rk)
         _ADD_SOURCE_(self%id_phy, -self%plambda*PHY)
         _ADD_SOURCE_(self%id_zoo, -self%zlambda*ZOO)
         _ADD_SOURCE_(self%id_dop, self%plambda*PHY+self%zlambda*ZOO)

      _LOOP_END_
   end subroutine do

   elemental real(rk) FUNCTION phi(u)
      real(rk), intent(in) :: u
      
!      phi= u*(0.555588d0+0.004926d0*u)/(1.0d0+0.188721d0*u)

      if(u.gt.1.0e-6_rk) then
         phi= LOG(u+SQRT(1.0_rk+u*u))-(SQRT(1.0_rk+u*u)-1.0_rk)/u
      else
         phi=0.0_rk
      endif
   END FUNCTION

end module mops_plankton
