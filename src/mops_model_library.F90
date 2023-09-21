module mops_model_library

   use fabm_types, only: type_base_model_factory, type_base_model

   use mops_carbon
   use mops_oxygen
   use mops_radiation
   use mops_insolation
   use mops_phytoplankton
   use mops_zooplankton
   use mops_nitrogen_fixation
   use mops_remineralization
   use mops_detritus
   use mops_tracer
   use mops_runoff

   implicit none

   private

   type, extends(type_base_model_factory) :: type_factory
   contains
      procedure :: create
   end type

   type (type_factory), save, target, public :: mops_model_factory

contains

   subroutine create(self,name,model)
      class (type_factory), intent(in) :: self
      character(*),         intent(in) :: name
      class (type_base_model), pointer :: model

      select case (name)
         case ('carbon');    allocate(type_mops_carbon::model)
         case ('oxygen');    allocate(type_mops_oxygen::model)
         case ('radiation'); allocate(type_mops_radiation::model)
         case ('insolation'); allocate(type_mops_insolation::model)
         case ('phytoplankton'); allocate(type_mops_phytoplankton::model)
         case ('zooplankton'); allocate(type_mops_zooplankton::model)
         case ('nitrogen_fixation'); allocate(type_mops_nitrogen_fixation::model)
         case ('remineralization'); allocate(type_mops_remineralization::model)
         case ('detritus'); allocate(type_mops_detritus::model)
         case ('tracer'); allocate(type_mops_tracer::model)
         case ('runoff'); allocate(type_mops_runoff::model)
         ! Add new models here
         case default
            call self%type_base_model_factory%create(name, model)
      end select
   end subroutine create

end module
