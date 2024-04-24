#include "fabm_driver.h"

module mops_detritus

   use fabm_types
   use mops_shared
   use fabm_builtin_depth_integral
   use fabm_expressions

   implicit none

   private

   type, extends(type_base_model), public :: type_mops_detritus
      type (type_dependency_id) :: id_bgc_z
      type (type_dependency_id) :: id_bgc_dz ! VS for CaCO3 divergence
      type (type_dependency_id) :: id_det_prod ! VS for CaCO3 divergence
      type (type_horizontal_dependency_id) :: id_int_det_prod ! VS for CaCO3 divergence
      type (type_bottom_dependency_id) :: id_bgc_z_bot
      type (type_state_variable_id) :: id_det, id_dic, id_alk
      type (type_diagnostic_variable_id) :: id_f8 ! VS CaCO3 production
      type (type_diagnostic_variable_id) :: id_fdiv_caco3 ! (f9) VS CaCO3 divergence
      type (type_dependency_id) :: id_fdiv_caco3_in ! an immediate dependency of the former
      type (type_bottom_diagnostic_variable_id) :: id_burial
      ! VS I am planning to use change rates id_det_prod
      !    and id_fdiv_caco3 (which is supposed to be calculated in "do_column")
      !    in procedure "do", using _ADD_TO_SOURCE_ for state variables DIC and Alk, later
      ! VS Can I use a diagnostic variable calculated in one type procedure (do_column)
      !    within another type procedure (do)?

      real(rk) :: detlambda, detwb, detmartin
      real(rk) :: burdige_fac, burdige_exp
      real(rk) :: frac_caco3, length_caco3 ! VS parameters for CaCO3 divergence
   contains
      ! Model procedures
      procedure :: initialize
      procedure :: get_vertical_movement
      procedure :: do_column ! VS to calculate CaCO3 divergence
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
      call self%get_parameter(self%frac_caco3, 'frac_caco3', 'mol CaCO3/mol C','calcite-to-organic-carbon ratio', default=0.32_rk)
      ! VS Parameter $l_\mathrm{CaCO3}$ in Cien et al., 2022
      call self%get_parameter(self%length_caco3, 'length_caco3', 'm','lenght scale for the e-fold function of dissolving CaCO3', default=4289.4_rk)

      call self%register_state_variable(self%id_det, 'c', 'mmol P/m3', 'detritus', minimum=0.0_rk)
      call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_det)

      call self%register_diagnostic_variable(self%id_burial, 'burial', 'mmol P/m2/d', 'burial')
      ! VS diagnostic variable f8
      call self%register_diagnostic_variable(self%id_f8, 'f8', 'mmol CaCO3/m3/d', 'production of CaCO3')
      ! VS diagnostic variable fdiv_caco3 (f9)
      call self%register_diagnostic_variable(self%id_fdiv_caco3, 'fdiv_caco3', 'mmol CaCO3/m3/d', 'divergence of CaCO3')

      ! Register environmental dependencies
      ! VS stuff needed for CaCO3 flux divergence in subroutine do_column
      call self%register_dependency(self%id_bgc_z, standard_variables%depth)
      call self%register_dependency(self%id_bgc_dz, standard_variables%cell_thickness)
      call self%register_dependency(self%id_det_prod, detritus_production_by_plankton)
      call self%register_expression_dependency(self%id_int_det_prod,vertical_integral(self%id_det_prod))
      ! direct dependency on corresponding diagnostic variable
      call self%register_dependency(self%id_fdiv_caco3_in, 'fdiv_caco3', 'mmol CaCO3/m3/d', 'divergence of CaCO3')

      call self%register_dependency(self%id_bgc_z_bot, standard_variables%bottom_depth)

      call self%register_state_variable(self%id_alk, 'alk', 'mmol Alk/m3', 'alkalinity', minimum=0.0_rk)

      self%dt = 86400.0_rk
   end subroutine

   subroutine get_vertical_movement(self, _ARGUMENTS_DO_)
      class (type_mops_detritus), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: detwa, bgc_z, wdet

      detwa = self%detlambda/self%detmartin
      _LOOP_BEGIN_
         _GET_(self%id_bgc_z, bgc_z)
         wdet = self%detwb + bgc_z*detwa
         _ADD_VERTICAL_VELOCITY_(self%id_det, -wdet)
      _LOOP_END_
   end subroutine get_vertical_movement

   subroutine do_column(self, _ARGUMENTS_DO_COLUMN_)
      class (type_mops_detritus), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_

      real(rk) :: z, dz, att_z, att_dz
      real(rk) :: int_det_prod ! total produced detritus C in water column (mol C/m2/d)
      real(rk) :: int_caco3_prod ! total produced CaCO3 in water column (mol CaCO3/m2/d)
      real(rk) :: fdiv_caco3 ! diagnostic value to be computed
      ! calculate CaCO3 divergence acc. to MOPS3.1 code (using "_GET_" for "space variate stuff"?)
      _GET_HORIZONTAL_(self%id_int_det_prod, int_det_prod) ! seems to be forbidden, here
      int_caco3_prod = rcp * self%frac_caco3 * int_det_prod
      _DOWNWARD_LOOP_BEGIN_
         _GET_(self%id_bgc_z, z)
         _GET_(self%id_bgc_dz, dz)
         att_z = exp( z / self%length_caco3 ) ! "attenuation" (of water columns produced CaCO3) at depth z
         att_dz = exp( dz / self%length_caco3 ) ! "attenuation" caused by a cell thickness dz
         ! VS Original MOPS3.1 code uses flux fractions fcaco3(k) of CaCO3 calculated as
         !    fcaco3(k) = exp( -zu(k)/length_caco3 ), with upper boundary depth zu(k) = z(k)-dz(k)/2.
         !    Since zu(k+1) is the same as the lower boundary depth z(k)+dz(k)/2 of layer k,
         !    we obtain (writing z=z(k), dz=dz(k), and L=length_caco3)
         !    fcaco3(k)-fcaco3(k+1) = exp((-z+dz/2)/L) - exp((-z-dz/2)/L)
         !                          = exp(-z/L) * (exp(dz/2/L)-exp(-dz/2/L))
         !                          = exp(-z/L) * (sqrt(exp(dz/L))-1/sqrt(exp(dz/L)))
         !                          = exp(-z/L) * (exp(dz/L)-1)/sqrt(exp(dz/L))
         !                          = 1 / att_z * ( att_dz - 1 ) / sqrt( att_dz )
         !    Now, fdiv_caco3(k) = int_caco3_prod * (fcaco3(k)-fcaco3(k+1)) / dz(k),
         !                       = int_caco3_prod * ( att_dz - 1 ) / sqrt( att_dz ) / att_z / dz
         fdiv_caco3 = int_caco3_prod * ( att_dz - 1._rk ) / sqrt( att_dz ) / att_z / dz
         _SET_DIAGNOSTIC_( self%id_fdiv_caco3, fdiv_caco3 )
      _DOWNWARD_LOOP_END_
      _MOVE_TO_BOTTOM_
      ! VS At the seafloor, all incoming CaCO3 remains to be dissolved,
      ! i.e., actually fdiv_caco3(k) = int_caco3_prod * fcaco3(k) / dz(k)
      !                              = int_caco3_prod * sqrt( att_dz ) / att_z / dz
      fdiv_caco3 = int_caco3_prod * sqrt( att_dz ) / att_z / dz
      _SET_DIAGNOSTIC_( self%id_fdiv_caco3, fdiv_caco3 )
   end subroutine do_column

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_mops_detritus), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: det_prod, caco3_prod, fdiv_caco3
      _LOOP_BEGIN_
         _GET_(self%id_det_prod, det_prod) ! detritus produced by plankton
         _GET_(self%id_fdiv_caco3_in, fdiv_caco3) ! implicite CaCO divergence
         caco3_prod = rcp * self%frac_caco3 * det_prod ! CaCO3 portion of detritus produced by plankton
         ! effects of CaCO3 processes on DIC and Alk:
         _ADD_SOURCE_(self%id_dic, fdiv_caco3-caco3_prod)
         _ADD_SOURCE_(self%id_alk, -2._rk*(fdiv_caco3-caco3_prod))
      _LOOP_END_
   end subroutine do

   subroutine do_bottom(self, _ARGUMENTS_DO_BOTTOM_)
      class (type_mops_detritus), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: detwa, bgc_z, DET, wdet, fDET, flux_l

      detwa = self%detlambda/self%detmartin
      _BOTTOM_LOOP_BEGIN_
         _GET_BOTTOM_(self%id_bgc_z_bot, bgc_z)
         _GET_(self%id_det, DET)
         wdet = self%detwb + bgc_z*detwa
         fDET = wdet*DET
         flux_l = MIN(1.0_rk,self%burdige_fac*fDET**self%burdige_exp)*fDET
         _ADD_BOTTOM_FLUX_(self%id_det, -flux_l)
         _SET_BOTTOM_DIAGNOSTIC_(self%id_burial, flux_l)
      _BOTTOM_LOOP_END_
   end subroutine

end module mops_detritus
