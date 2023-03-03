#include "fabm_driver.h"

module mops_runoff

   use fabm_types
   use fabm_expressions
   use mops_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_mops_runoff
      type (type_bottom_dependency_id) :: id_burial
      type (type_horizontal_dependency_id) :: id_burial_am
      type (type_surface_diagnostic_variable_id) :: id_runoff
      type (type_state_variable_id) :: id_pho, id_din, id_dic
   contains
      procedure :: initialize
      procedure :: do_surface
   end type type_mops_runoff

contains

   subroutine initialize(self, configunit)
      class (type_mops_runoff), intent(inout), target :: self
      integer,                  intent(in)            :: configunit

      call self%register_dependency(self%id_burial, 'burial', 'mmol P/m2/d', 'burial')
      call self%register_dependency(self%id_burial_am, temporal_mean(self%id_burial, period=365*86400.0_rk, resolution=30*86400.0_rk))

      call self%register_diagnostic_variable(self%id_runoff, 'runoff', 'mmol P/m2/d', 'runoff')

      call self%register_state_dependency(self%id_din, 'din', 'mmol N/m3', 'dissolved inorganic nitrogen')
      call self%register_state_dependency(self%id_pho, 'pho', 'mmol P/m3', 'phosphate')
      call self%register_state_dependency(self%id_dic, 'dic', 'mmol C/m3', 'dissolved inorganic carbon')

      self%dt = 86400.0_rk
   end subroutine

   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_mops_runoff), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: burial_am

      _SURFACE_LOOP_BEGIN_
         _GET_BOTTOM_(self%id_burial_am, burial_am)
         _ADD_SURFACE_FLUX_(self%id_pho, burial_am)
         _ADD_SURFACE_FLUX_(self%id_din, burial_am*rnp)
         _ADD_SURFACE_FLUX_(self%id_dic, burial_am*rcp)
         _SET_SURFACE_DIAGNOSTIC_(self%id_runoff, burial_am)
      _SURFACE_LOOP_END_
   end subroutine

end module mops_runoff
