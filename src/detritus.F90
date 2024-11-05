#include "fabm_driver.h"

module mops_detritus

   use fabm_types
   use mops_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_mops_detritus
      type (type_dependency_id) :: id_bgc_z, id_bgc_dz
      type (type_bottom_dependency_id) :: id_bgc_z_bot
      type (type_state_variable_id) :: id_det
      type (type_bottom_diagnostic_variable_id) :: id_burial, id_detdiv

      real(rk) :: detlambda, detwb, detmartin
      real(rk) :: burdige_fac, burdige_exp
   contains
      ! Model procedures
      procedure :: initialize
!      procedure :: get_vertical_movement
      procedure :: do_column
      procedure :: do_bottom
   end type type_mops_detritus

contains

   subroutine initialize(self, configunit)
      class (type_mops_detritus), intent(inout), target :: self
      integer,                    intent(in)            :: configunit

      call self%get_parameter(self%detlambda, 'detlambda', '1/d','detritus remineralization rate', default=0.05_rk)
      call self%get_parameter(self%detwb, 'detwb', 'm/d','offset for detritus sinking', default=0.0_rk)
      call self%get_parameter(self%detmartin, 'detmartin', '-','exponent for Martin curve', default=0.8580_rk)
      call self%get_parameter(self%burdige_fac, 'burdige_fac', '-','factor for sediment burial (see Kriest and Oschlies, 2013)', default=1.6828_rk)
      call self%get_parameter(self%burdige_exp, 'burdige_exp', '-','exponent for sediment burial (see Kriest and Oschlies, 2013)', default=0.799_rk)

! VS nur kurz without minimum value to avoid clipping in TMM implementation
! (see Jorns mail on October 16, 2024)
      call self%register_state_variable(self%id_det, 'c', 'mmol P/m3', 'detritus')!, minimum=0.0_rk)
      call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_det)

      call self%register_diagnostic_variable(self%id_burial, 'burial', 'mmol P/m2/d', 'burial')
      call self%register_diagnostic_variable(self%id_detdiv, 'detdiv', 'mmol P/m3/d', 'divergence')

      ! Register environmental dependencies
      call self%register_dependency(self%id_bgc_z, standard_variables%depth)
      call self%register_dependency(self%id_bgc_dz, standard_variablles%cell_thickness)
      call self%register_dependency(self%id_bgc_z_bot, standard_variables%bottom_depth)

      self%dt = 86400.0_rk
   end subroutine

! VS replacing "get_vertical_movement" by "do_column"
!    in order to calculate detritus divergence like "PETSC-MOPS" does,
!    might want to adapt "PETSC-MOPS" later, instead
!   subroutine get_vertical_movement(self, _ARGUMENTS_DO_)
!      class (type_mops_detritus), intent(in) :: self
!      _DECLARE_ARGUMENTS_DO_
!
!      real(rk) :: detwa, bgc_z, wdet
!
!      detwa = self%detlambda/self%detmartin
!      _LOOP_BEGIN_
!         _GET_(self%id_bgc_z, bgc_z)
!         wdet = self%detwb + bgc_z*detwa
!         _ADD_VERTICAL_VELOCITY_(self%id_det, -wdet)
!      _LOOP_END_
!   end subroutine get_vertical_movement

   subroutine do_column(self, _ARGUMENTS_DO_COLUMN_)
      class (type_mops_radiation), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_

      real(rk) :: detwa, bgc_z, DET, wdet, fDET, flux_l

   end subroutine

   subroutine do_bottom(self, _ARGUMENTS_DO_BOTTOM_)
      class (type_mops_detritus), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: detwa, bgc_z, bgc_dz, DET, wdet, flux_u, flux_l

      detwa = self%detlambda/self%detmartin
      flux_u = 0.0_rk ! VS flux through upper box layer
      _DOWNWARD_LOOP_BEGIN_
         _GET_(self%id_bgc_z, bgc_z)
         _GET_(self%id_det, DET)
         DET = MAX(DET-alimit*alimit,0.0_rk)
         wdet = self%detwb + bgc_z*detwa
         fdet = wdet*DET
         flux_l = fdet ! VS flux thrugh lower box layer
         detdiv = () 
!         flux_l = MIN(1.0_rk,self%burdige_fac*fDET**self%burdige_exp)*fDET
         _SET_DIAGNOSTIC_(self%id_detdiv, detdiv)

!         ! VS nur kurz
!         print *, 'detwa is', detwa
!         print *, 'wdet is', wdet
!         print *, 'DET is', DET
!         print *, 'fDET is', fDET
!         print *, '-flux_l is', -flux_l
!         print *, '-flux_l / sec / layer_thickness is', -flux_l / 86400.0_rk / 50.0_rk
          flux_u = flux_l
      _DOWNWARD_LOOP_END_

   end subroutine

end module mops_detritus
