#include "fabm_driver.h"

module mops_radiation

   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_mops_radiation
      type (type_horizontal_dependency_id) :: id_bgc_swr, id_bgc_seaice
      type (type_dependency_id)            :: id_bgc_dz, id_att
      type (type_diagnostic_variable_id)   :: id_ciz, id_atten

      ! Parameters
      real(rk) :: parfrac, ACkw
   contains
      ! Model procedures
      procedure :: initialize
      procedure :: do_column
   end type type_mops_radiation

contains

   subroutine initialize(self, configunit)
      class (type_mops_radiation), intent(inout), target :: self
      integer,                     intent(in)            :: configunit

      call self%get_parameter(self%parfrac, 'parfrac', '1', 'PAR fraction of shortwave radiation', default=0.4_rk) 
      call self%get_parameter(self%ACkw, 'ACkw', '1/m','attenuation of water', default=0.04_rk) 

      call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, self%ACkw)

      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_ciz, 'ciz', 'W m-2', 'PAR at top of the layer', source=source_do_column)
      call self%register_diagnostic_variable(self%id_atten, 'atten', '-', 'attenuation across the layer', source=source_do_column)

      ! Register environmental dependencies
      call self%register_dependency(self%id_bgc_dz, standard_variables%cell_thickness)
      call self%register_dependency(self%id_bgc_swr, 'sfac', 'W m-2', 'net downwelling shortwave flux at water surface')
      call self%register_dependency(self%id_bgc_seaice, standard_variables%ice_area_fraction)
      call self%register_dependency(self%id_att, standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux)
end subroutine
   
   subroutine do_column(self, _ARGUMENTS_DO_COLUMN_)
      class (type_mops_radiation), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_

      real(rk) :: bgc_swr, bgc_seaice, bgc_dz, ciz, att, atten

      _GET_SURFACE_(self%id_bgc_swr,bgc_swr)
      _GET_SURFACE_(self%id_bgc_seaice,bgc_seaice)
      ciz = bgc_swr*(1.0_rk-bgc_seaice)*self%parfrac
      _DOWNWARD_LOOP_BEGIN_
         _SET_DIAGNOSTIC_(self%id_ciz,ciz)
         _GET_(self%id_bgc_dz,bgc_dz)     ! Layer height (m)
         _GET_(self%id_att,att)           ! Attenuation (1/m)
         atten = att*bgc_dz
         ciz = ciz * exp(-atten)
         _SET_DIAGNOSTIC_(self%id_atten,atten)
      _DOWNWARD_LOOP_END_
   end subroutine do_column

end module mops_radiation
