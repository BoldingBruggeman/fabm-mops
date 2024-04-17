#include "fabm_driver.h"

module mops_detritus

   use fabm_types
   use mops_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_mops_detritus
      type (type_dependency_id) :: id_bgc_z
      type (type_dependency_id) :: id_bgc_dz ! for CaCO3 divergence
      type (type_bottom_dependency_id) :: id_bgc_z_bot
      type (type_state_variable_id) :: id_det
      type (type_diagnostic_variable_id) :: id_fdiv_caco3 ! CaCO3 divergence
      type (type_bottom_diagnostic_variable_id) :: id_burial

      real(rk) :: detlambda, detwb, detmartin
      real(rk) :: burdige_fac, burdige_exp
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
      ! VS Parameter $l_\mathrm{CaCO3}$ in Cien et al., 2022
      call self%get_parameter(self%length_caco3, 'frac_caco3', 'mol CaCO3/mol C','calcite-to-organic-carbon ratio', default=0.32_rk)
      ! VS Parameter $\sigma_\mathrm{CaCO3}$ in Chien et al., 2022
      call self%get_parameter(self%length_caco3, 'length_caco3', 'm','lenght scale for the e-fold function of dissolving CaCO3', default=4289.4_rk)

      call self%register_state_variable(self%id_det, 'c', 'mmol P/m3', 'detritus', minimum=0.0_rk)
      call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_det)

      call self%register_diagnostic_variable(self%id_burial, 'burial', 'mmol P/m2/d', 'burial')
      ! VS diagnostic variable fdiv_caco3
      call self%register_diagnostic_variable(self%id_fdiv_caco3, 'fdiv_caco3', 'mmol CaCO3/m3/d', 'divergence of CaCO3', minimum=0.0_rk)

      ! Register environmental dependencies
      call self%register_dependency(self%id_bgc_z, standard_variables%depth)
      ! VS id_bgc_dz needed for CaCO3 flux divergence in subroutine do_column
      call self%register_dependency(self%id_bgc_dz, standard_variables%cell_thikness)
      call self%register_dependency(self%id_bgc_z_bot, standard_variables%bottom_depth)

      ! VS use depth integrated detritus production by plankton (mmol P/m2/d) in subroutine do_column
      call self%register_expression_dependency(self%id_intdet_prod,vertical_integral(self%id_det_prod))

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

   subroutine do_column(self, _ARGUMENTS_DO_)
      class (type_mops_detritus), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: z, dz, att_z, att_dz
      real(rk) :: int_det_prod ! total produced detritus C in water column (mol C/m2/d)
      real(rk) :: caco3_prod ! total produced CaCO3 in water column (mol CaCO3/m2/d)
      ! calculate CaCO3 divergence acc. to MOPS3.1 code
      _GET_(self%int_det_prod, int_det_prod) ! VS using "_GET_" for "space variate stuff"
      caco3_prod = self%rcp * self%frac_caco3 * int_det_prod
      _DOWNWARD_LOOP_BEGIN
         _GET_(self%id_bgc_z, z)
         _GET_(self%id_bgc_dz, dz)
         att_z = exp( z / self%length_caco3 ) ! "attenuation" (of water columns produced CaCO3) at depth z
         att_dz = exp( dz / self%length_caco3 ) ! "attenuation" caused by a cell thickness dz
         _SET_DIAGNOSTIC_( self%id_fdiv_caco3, caco3_prod * ( att_dz - 1._rk ) / sqrt( att_dz ) / att_z / dz )
         ! VS Original MOPS3.1 code uses flux fractions fcaco3(k) of CaCO3 calculated as
         !    fcaco3(k) = exp( -zu(k)/length_caco3 ), with upper boundary depth zu(k) = z(k)-dz(k)/2.
         !    Since zu(k+1) is the same as the lower boundary depth z(k)+dz(k)/2 of layer k,
         !    we obtain (writing z=z(k), dz=dz(k), L=length_caco3)
         !    fcaco3(k)-fcaco3(k+1) = exp((-z+dz/2)/L) - exp((-z-dz/2)/L)
         !                          = exp(-z/L) * (exp(dz/2/L)-exp(-dz/2/L))
         !                          = exp(-z/L) * (sqrt(exp(dz/L))-1/sqrt(exp(dz/L)))
         !                          = exp(-z/L) * (exp(dz/L)-1)/sqrt(exp(dz/L))
         !                          = 1 / att_z * ( att_dz - 1 ) / sqrt( att_dz )
         !    Now, fdiv_caco3 = caco3_prod * (fcaco3(k)-fcaco3(k+1)) / dz(k),
         !                    = caco3_prod * ( att_dz - 1 ) / sqrt( att_dz ) / att_z / dz
      _DOWNWARD_LOOP_END
 
   end subroutine

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
