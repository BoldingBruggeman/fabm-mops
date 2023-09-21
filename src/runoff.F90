#include "fabm_driver.h"

module mops_runoff

   use fabm_types
   use fabm_expressions
   use mops_shared

   implicit none

   private

   type, extends(type_base_model) :: type_mops_runoff
      type (type_horizontal_dependency_id) :: id_source
      type (type_bottom_dependency_id) :: id_bottom_depth
      type (type_state_variable_id) :: id_pho, id_din, id_dic
   contains
      procedure :: initialize
   end type

   type, extends(type_mops_runoff), public :: type_mops_surface_runoff
   contains
      procedure :: do_surface
   end type

   type, extends(type_mops_runoff), public :: type_mops_interior_runoff
   contains
      procedure :: do
   end type

contains

   subroutine initialize(self, configunit)
      class (type_mops_runoff), intent(inout), target :: self
      integer,                  intent(in)            :: configunit

      call self%register_dependency(self%id_source, 'source', 'mmol P/m2/d', 'source')
      call self%register_dependency(self%id_bottom_depth, standard_variables%bottom_depth)

      call self%register_state_dependency(self%id_din, 'din', 'mmol N/m3', 'dissolved inorganic nitrogen')
      call self%register_state_dependency(self%id_pho, 'pho', 'mmol P/m3', 'phosphate')
      call self%register_state_dependency(self%id_dic, 'dic', 'mmol C/m3', 'dissolved inorganic carbon')

      self%dt = 86400.0_rk
   end subroutine

   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_mops_surface_runoff), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: source

      _SURFACE_LOOP_BEGIN_
         _GET_HORIZONTAL_(self%id_source, source)
         _ADD_SURFACE_FLUX_(self%id_pho, source)
         _ADD_SURFACE_FLUX_(self%id_din, source*rnp)
         _ADD_SURFACE_FLUX_(self%id_dic, source*rcp)
      _SURFACE_LOOP_END_
   end subroutine

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_mops_interior_runoff), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: source, bottom_depth

      _LOOP_BEGIN_
         _GET_HORIZONTAL_(self%id_source, source)
         _GET_BOTTOM_(self%id_bottom_depth, bottom_depth)
         source = source / bottom_depth
         _ADD_SOURCE_(self%id_pho, source)
         _ADD_SOURCE_(self%id_din, source*rnp)
         _ADD_SOURCE_(self%id_dic, source*rcp)
      _LOOP_END_
   end subroutine

end module mops_runoff
