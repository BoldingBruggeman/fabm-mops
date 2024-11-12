#include "fabm_driver.h"

module mops_remineralization

   use fabm_types
   use mops_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_mops_remineralization
      type (type_dependency_id) :: id_bgc_theta
      type (type_state_variable_id) :: id_dop, id_det, id_oxy, id_din, id_po4, id_dic
      type (type_diagnostic_variable_id) :: id_f4, id_f5, id_f7

      real(rk) :: dlambda, detlambda, subox, subdin, ACkbaco2, ACkbacdin
   contains
      ! Model procedures
      procedure :: initialize
      procedure :: do
   end type type_mops_remineralization

contains

   subroutine initialize(self, configunit)
      class (type_mops_remineralization), intent(inout), target :: self
      integer,                            intent(in)            :: configunit

      call self%get_parameter(self%dlambda, 'dlambda', '1/d','DOP remineralization rate', default=0.0005133333049196715389730860613828195004870736_rk)
      call self%get_parameter(self%detlambda, 'detlambda', '1/d','detritus remineralization rate', default=0.05_rk)
      call self%get_parameter(self%subox, 'subox', 'mmol/m3','minimum oxygen for oxic degradation', default=1.0_rk)
      call self%get_parameter(self%subdin, 'subdin', 'mmol/m3','minimum DIN for denitrification', default=16.0_rk)
      call self%get_parameter(self%ACkbaco2, 'ACkbaco2', 'mmol/m3','Half sat.-constant for oxic degradation (see Kriest and Oschlies, 2015)', default=1.145532_rk)
      call self%get_parameter(self%ACkbacdin, 'ACkbacdin', 'mmol/m3','Half sat.-constant for suboxic degradation (see Kriest and Oschlies, 2015)', default=23.083559_rk)

      call self%register_state_dependency(self%id_dop, 'dop', 'mmol P/m3', 'dissolved organic phosphorus')
      call self%register_state_dependency(self%id_det, 'det', 'mmol P/m3', 'detritus')
      call self%register_state_dependency(self%id_oxy, 'oxy', 'mmol O2/m3', 'oxygen')
      call self%register_state_dependency(self%id_din, 'din', 'mmol N/m3', 'dissolved inorganic nitrogen')
      call self%register_state_dependency(self%id_po4, 'pho', 'mmol P/m3', 'phosphate')
      call self%register_state_dependency(self%id_dic, 'dic', 'mmol C/m3', 'dissolved inorganic carbon')

      ! Register environmental dependencies
      call self%register_dependency(self%id_bgc_theta, standard_variables%temperature)

      call self%register_diagnostic_variable(self%id_f4, 'f4', 'mmol/m3/d', 'oxic remineralization')
      !call self%register_diagnostic_variable(self%id_f5, 'f5', 'mmol/m3/d', 'river input')
      call self%register_diagnostic_variable(self%id_f7, 'f7', 'mmol/m3/d', 'denitrification')

      self%dt = 86400.0_rk
   end subroutine

   subroutine do(self, _ARGUMENTS_DO_)
      class (type_mops_remineralization), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: DOP, DET, OXY, DIN
      real(rk) :: oxymm, o2req, o2usefrac, remindop, remindet
      real(rk) :: dinmm, dinreq, dinusefrac, denitdop, denitdet
      real(rk) :: topo4

      _LOOP_BEGIN_

      _GET_(self%id_dop, DOP)
      _GET_(self%id_det, DET)
      _GET_(self%id_oxy, OXY)
      _GET_(self%id_din, DIN)

      DOP = MAX(DOP - alimit*alimit, 0.0_rk)
      DET = MAX(DET - alimit*alimit, 0.0_rk)

! AEROBIC DECAY

! In contrast to the older (Kriest&Oschlies, 2013) version, this option:
! (1) does not degrade OM in the absence of O2, i.e. OM can accumulate 
! (2) uses a Michaelis-Menten Kinetic to slow down bacterial remineralisation under low O2
! (2) takes care not to use more O2 per timestep than available

! Michaelis-Menten limitation for oxic degradation: 

      OXY = MAX(OXY-self%subox,0.0_rk)
      oxymm = OXY*OXY/(OXY*OXY+self%ACkbaco2*self%ACkbaco2)

! O2 required for total remineralisation in a time step will then be:
      
      o2req = oxymm*(self%dlambda*DOP+self%detlambda*DET)*ro2ut*bgc_dt
      
! restrict remineralisation to amount of available oxygen
      
      if (o2req.gt.0.0_rk) then
         o2usefrac = MIN(OXY,o2req)/o2req
      else
         o2usefrac = 0.0_rk
      endif

      remindop = oxymm*self%dlambda*DOP*o2usefrac
      remindet = oxymm*self%detlambda*DET*o2usefrac

! ANAEROBIC DECAY INCL. ANAMMOX ETC.

      if(OXY.lt.36.0_rk) then

         DIN = MAX(DIN-self%subdin,0.0_rk)
         dinmm = DIN*DIN/(DIN*DIN+self%ACkbacdin*self%ACkbacdin)*(1.0_rk-oxymm)

! NO3 required for total remineralisation in a time step will then be:

         dinreq = dinmm*(self%dlambda*DOP+self%detlambda*DET)*rhno3ut*bgc_dt

! restrict remineralisation to amount of variable oxygen

         if (dinreq.gt.0.0_rk) then
            dinusefrac = MIN(DIN,dinreq)/dinreq
         else
            dinusefrac = 0.0_rk
         endif

! restrict anaerobic processes to regions with low oxygen concentration

         denitdop = dinmm*self%dlambda*DOP*dinusefrac
         denitdet = dinmm*self%detlambda*DET*dinusefrac

      else

         denitdop = 0.0_rk
         denitdet = 0.0_rk

      endif

      topo4 = remindop+remindet+denitdop+denitdet
      _ADD_SOURCE_(self%id_po4, topo4)
      _ADD_SOURCE_(self%id_dop, -remindop-denitdop)
      _ADD_SOURCE_(self%id_oxy, -(remindop+remindet)*ro2ut)
      _ADD_SOURCE_(self%id_det, -remindet-denitdet)

      _ADD_SOURCE_(self%id_din, +(remindop+remindet)*rnp-(denitdop+denitdet)*rhno3ut)
      _ADD_SOURCE_(self%id_dic, topo4*rcp)
      _SET_DIAGNOSTIC_(self%id_f4, remindop+remindet)
      _SET_DIAGNOSTIC_(self%id_f7, denitdop+denitdet)

!#ifdef RUNOFF
!      DO K=1,bgc_kloc
!        runoff(k) = bgc_globalrunoff * bgc_runoffvol(k)
!        _SET_DIAGNOSTIC_(self%id_f5, runoff(k)*bgc_dt)
!      ENDDO
!#else
!      runoff(1) = bgc_globalrunoff/bgc_dz(1)
!      _SET_DIAGNOSTIC_(self%id_f5, runoff(1)*bgc_dt)
!#endif

      _LOOP_END_
   end subroutine do

end module mops_remineralization
