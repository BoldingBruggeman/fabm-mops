#include "fabm_driver.h"

module mops_nitrogen_fixation

   use fabm_types
   use mops_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_mops_nitrogen_fixation
      type (type_dependency_id) :: id_bgc_theta
      type (type_state_variable_id) :: id_po4, id_din
      type (type_diagnostic_variable_id) :: id_f6

      real(rk) :: tf0, tf1, tf2, tff, nfix
   contains
      ! Model procedures
      procedure :: initialize
      procedure :: do
   end type

contains

   subroutine initialize(self, configunit)
      class (type_mops_nitrogen_fixation), intent(inout), target :: self
      integer,                             intent(in)            :: configunit

! N2-Fixatioon
! Factors tf2, tf1 and tf0 are a polynomial (2nd order) 
! approximation to the functional relationship by Breitbarth et al. (2007),
! for temperature dependence of Trichodesmium growth, 
! their eq. (2), assuming that their powers relate to "e" (10^), and
! that the last but one sign has to be reversed.
! The relation will be scaled to their max. growth rate.
! Note that the second order approx. is basically similar to their
! function 2 for T-dependent nitrogen fixation multiplied by 4 
! (2 [N atoms per mole] * 12 [light hrs per day]/6 [C-atoms per N-atoms])
      call self%get_parameter(self%tf2, 'tf2', '1/(d degC^2)','quadratic coefficient for temperature dependence', default=-0.0042_rk)
      call self%get_parameter(self%tf1, 'tf1', '1/(d degC)','linear coefficient for temperature dependence', default=0.2253_rk)
      call self%get_parameter(self%tf0, 'tf0', '1/d','constant coefficient for temperature dependence', default=-2.7819_rk)
      call self%get_parameter(self%tff, 'tff', '1/d','normalization factor for temperature dependence', default=0.2395_rk)
      call self%get_parameter(self%nfix, 'nfix', 'mmol N/m3/d','max. fixation rate', default=0.002272073044_rk)

      call self%register_diagnostic_variable(self%id_f6, 'f6', 'mmol N/m3/d', 'nitrogen fixation')

      call self%register_state_dependency(self%id_din, 'din', 'mmol N/m3', 'dissolved inorganic nitrogen')
      call self%register_state_dependency(self%id_po4, 'po4', 'mmol P/m3', 'phosphate')

      ! Register environmental dependencies
      call self%register_dependency(self%id_bgc_theta, standard_variables%temperature)

      self%dt = 86400.0_rk
   end subroutine

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_mops_nitrogen_fixation), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: ttemp, PO4, DIN
      real(rk) :: nfixtfac, dinlim, nfixnfac, nfixation

      _LOOP_BEGIN_

         _GET_(self%id_bgc_theta, ttemp)
         _GET_(self%id_po4, PO4)
         _GET_(self%id_din, DIN)

! Relaxation of N:P to Redfield values (mimick cyanobacteria)
          if(PO4.gt.vsafe) then

            nfixtfac = MAX(0.0_rk,self%tf2*ttemp*ttemp + self%tf1*ttemp + self%tf0)/self%tff
            dinlim = MAX(0.0_rk,DIN)
            nfixnfac = MAX(0.0_rk, 1.0_rk-dinlim/(PO4*rnp))
            nfixation = nfixtfac*nfixnfac*self%nfix

          else

             nfixation = 0.0_rk

          endif

          _SET_DIAGNOSTIC_(self%id_f6, nfixation)
           _ADD_SOURCE_(self%id_din, nfixation)

      _LOOP_END_
   end subroutine do

end module mops_nitrogen_fixation
