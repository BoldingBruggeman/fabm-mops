#include "fabm_driver.h"

module mops_detritus

   use fabm_types
   use mops_shared
   use fabm_builtin_depth_integral
   use fabm_expressions

   implicit none

   private

   type, extends(type_base_model), public :: type_mops_detritus
      type (type_state_variable_id) :: id_det
      type (type_state_variable_id) :: id_dic
      type (type_state_variable_id) :: id_alk
      type (type_dependency_id) :: id_bgc_z
      type (type_dependency_id) :: id_bgc_dz
      type (type_bottom_dependency_id) :: id_bgc_z_bot 
      type (type_dependency_id) :: id_det_prod ! VS for CaCO3 divergence
      type (type_horizontal_dependency_id) :: id_int_det_prod ! VS for CaCO3 divergence
      type (type_diagnostic_variable_id) :: id_fdiv_det ! VS detritus divergence
      type (type_diagnostic_variable_id) :: id_f9 ! VS CaCO3 production
      type (type_diagnostic_variable_id) :: id_fdiv_caco3 ! VS CaCO3 divergence
      type (type_dependency_id) :: id_fdiv_caco3_in ! an immediate dependency of the former
      type (type_bottom_diagnostic_variable_id) :: id_burial

      real(rk) :: detlambda, detwb, detmartin
      real(rk) :: burdige_fac, burdige_exp
      real(rk) :: frac_caco3, length_caco3 ! VS parameters for CaCO3 divergence
      integer :: file_unit, ierr ! VS only for short
   contains
      ! Model procedures
      procedure :: initialize
      procedure :: do_column ! VS to calculate divergence of detritus and CaCO3
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
      ! VS Parameter $\sigma_\mathrm{CaCO3}$ in Chien et al., 2022
      call self%get_parameter(self%frac_caco3, 'frac_caco3', 'mol CaCO3/mol C','calcite-to-organic-carbon ratio', default=0.032_rk)
      ! VS Parameter $l_\mathrm{CaCO3}$ in Cien et al., 2022
      call self%get_parameter(self%length_caco3, 'length_caco3', 'm','lenght scale for the e-fold function of dissolving CaCO3', default=4289.4_rk)

      ! VS detritus without minimum value to avoid clipping in TMM implementation
      ! (see Jorns mail on October 16, 2024)
      call self%register_state_variable(self%id_det, 'c', 'mmol P/m3', 'detritus')
      call self%register_state_variable(self%id_alk, 'alk', 'mmol/m3', 'alkalinity')
      call self%register_state_dependency(self%id_dic, 'dic', 'mmol C/m3', 'dissolved inorganic carbon')

      call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_det)

      call self%register_diagnostic_variable(self%id_fdiv_det, 'fdiv_det', 'mmol P/m3/d', 'divergence', source=source_do_column)
      call self%register_diagnostic_variable(self%id_burial, 'burial', 'mmol P/m2/d', 'burial')
      ! VS diagnostic variable f9
      call self%register_diagnostic_variable(self%id_f9, 'f9', &
         'mmol CaCO3/m3/d', 'production of CaCO3', source=source_do_column)
      ! VS diagnostic variable fdiv_caco3 (f8)
      call self%register_diagnostic_variable(self%id_fdiv_caco3, 'fdiv_caco3', &
         'mmol CaCO3/m3/d', 'divergence of CaCO3', source=source_do_column)
      ! VS remark: if a _SET_DIAGNOSTIC_ is done in another subroutine than do
      !    this must be flagged, e.g., as above by source=source_do_column

      ! Register environmental dependencies
      ! VS stuff needed for CaCO3 flux divergence in subroutine do_column
      call self%register_dependency(self%id_bgc_z, standard_variables%depth)
      call self%register_dependency(self%id_bgc_dz, standard_variables%cell_thickness)
      call self%register_dependency(self%id_det_prod, detritus_production_by_plankton)
      call self%register_expression_dependency(self%id_int_det_prod,vertical_integral(self%id_det_prod))
      ! VS we have an immediate dependency on diagnostic variable fdiv_caco3
      call self%register_dependency(self%id_fdiv_caco3_in, 'fdiv_caco3', 'mmol CaCO3/m3/d', 'divergence of CaCO3')

      call self%register_dependency(self%id_bgc_z_bot, standard_variables%bottom_depth)

      self%dt = 86400.0_rk

      ! VS borrowed the following commands from fabm-pisces,
      !    they appear to be necessary to make the _ADD_SOURCE_ commands
      !    work in subroutine do_column (default is subroutine do)
      self%id_alk%sms%link%target%source = source_do_column
      self%id_dic%sms%link%target%source = source_do_column
      self%id_det%sms%link%target%source = source_do_column
   end subroutine

   subroutine do_column(self, _ARGUMENTS_DO_COLUMN_)
      class (type_mops_detritus), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_

      real(rk) :: bgc_z, bgc_dz
      real(rk) :: detwa, DET, wdet
      real(rk) :: fdet_u, fdet_l, fdiv_det, fdet_u_bottom, fdiv_det_bottom, fdiv_caco3_bottom
      real(rk) :: fcaco3_u, fcaco3_l
      real(rk) :: det_prod ! produced detritus C in (euphotic) box [mol C/m3/d]
      real(rk) :: caco3_prod ! produced CaCO3 in (euphotic) box [mol CaCO3/m3/d]
      real(rk) :: int_det_prod ! total produced detritus C in water column [mol C/m2/d]
      real(rk) :: int_caco3_prod ! total produced CaCO3 in water column [mol CaCO3/m2/d]
      real(rk) :: fdiv_caco3 ! diagnostic value to be computed

      detwa = self%detlambda/self%detmartin
      fdet_u = 0.0_rk ! VS detritus flux through upper box layer
      _GET_HORIZONTAL_(self%id_int_det_prod, int_det_prod)
      int_caco3_prod = rcp * self%frac_caco3 * int_det_prod
      _DOWNWARD_LOOP_BEGIN_
         fdet_u_bottom = fdet_u
         _GET_(self%id_bgc_z, bgc_z)
         _GET_(self%id_bgc_dz, bgc_dz)
         _GET_(self%id_det, DET)
         _GET_(self%id_det_prod, det_prod )
         DET = MAX(DET-alimit*alimit,0.0_rk)
         wdet = self%detwb + bgc_z*detwa
         fdet_l = wdet*DET ! VS detritus flux through lower box layer
         fdiv_det = (fdet_u-fdet_l)/bgc_dz  ! divergence of detritus
         fdet_u = fdet_l ! VS outgoing flux is incoming flux of the box below
         fcaco3_u = exp( -( bgc_z - bgc_dz / 2._rk ) / self%length_caco3 ) ! flow portion through layer top
         fcaco3_l = exp( -( bgc_z + bgc_dz / 2._rk ) / self%length_caco3 ) ! flow portion through layer bottom 
         fdiv_caco3 = int_caco3_prod * ( fcaco3_u - fcaco3_l ) / bgc_dz ! CaCO3 flux divergence
         caco3_prod = rcp * self%frac_caco3 * det_prod ! CaCO3 portion of DET produced by plankton
         _SET_DIAGNOSTIC_(self%id_fdiv_det, fdiv_det)
         _SET_DIAGNOSTIC_( self%id_fdiv_caco3, fdiv_caco3 ) ! f8_out in original MOPS code
         _SET_DIAGNOSTIC_( self%id_f9, caco3_prod )
         _ADD_SOURCE_(self%id_det, fdiv_det)
         _ADD_SOURCE_(self%id_dic, fdiv_caco3 - caco3_prod )
         _ADD_SOURCE_(self%id_alk, 2._rk * ( fdiv_caco3 - caco3_prod ) )
      _DOWNWARD_LOOP_END_
      _MOVE_TO_BOTTOM_
      fdet_l = MIN(1.0_rk,self%burdige_fac*fdet_l**self%burdige_exp)*fdet_l ! VS flux through bottom layer
      fdiv_det_bottom = (fdet_u_bottom-fdet_l)/bgc_dz  ! VS divergence at bottom box
      ! VS at the seafloor, all incoming CaCO3 remains to be dissolved
      fdiv_caco3_bottom = int_caco3_prod * fcaco3_u / bgc_dz
      _SET_DIAGNOSTIC_(self%id_fdiv_det, fdiv_det_bottom)
      _SET_DIAGNOSTIC_( self%id_fdiv_caco3, fdiv_caco3_bottom )
      ! VS adding the following correction terms for the bottom layer:
      _ADD_SOURCE_(self%id_det, fdiv_det_bottom-fdiv_det) ! VS correcting the former one
      _ADD_SOURCE_(self%id_dic, fdiv_caco3_bottom-fdiv_caco3 )
      _ADD_SOURCE_(self%id_alk, 2._rk*fdiv_caco3_bottom-2._rk*fdiv_caco3 )
   end subroutine do_column

   subroutine do_bottom(self, _ARGUMENTS_DO_BOTTOM_)
      class (type_mops_detritus), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: detwa, bgc_z, DET, wdet, fdet, fdet_l

      detwa = self%detlambda/self%detmartin
      _BOTTOM_LOOP_BEGIN_
         _GET_BOTTOM_(self%id_bgc_z_bot, bgc_z)
         _GET_(self%id_det, DET)

         DET = MAX(DET-alimit*alimit,0.0_rk)

         wdet = self%detwb + bgc_z*detwa
         fdet = wdet * DET
         fdet_l = MIN(1.0_rk,self%burdige_fac*fdet**self%burdige_exp)*fdet

         _SET_BOTTOM_DIAGNOSTIC_(self%id_burial, fdet_l)
      _BOTTOM_LOOP_END_
   end subroutine

end module mops_detritus
