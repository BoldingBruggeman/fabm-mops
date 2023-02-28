#include "fabm_driver.h"

module mops_tracer

   use fabm_types
   use mops_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_mops_tracer
      type (type_state_variable_id) :: id_c
   contains
      procedure :: initialize
   end type

contains

   subroutine initialize(self, configunit)
      class (type_mops_tracer), intent(inout), target :: self
      integer,                  intent(in)            :: configunit

      call self%register_state_variable(self%id_c, 'c', 'mmol P/m3', 'concentration')
      call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_c)
   end subroutine

end module mops_tracer
