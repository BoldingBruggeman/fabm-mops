#include "fabm_driver.h"

module mops_carbon

   use fabm_types
   use mops_shared

   implicit none

   private

   type, extends(type_base_model), public :: type_mops_carbon
      type (type_state_variable_id) :: id_dic
      type (type_dependency_id) :: id_pho, id_sil, id_bgc_salt, id_bgc_theta
      type (type_surface_dependency_id) :: id_bgc_wind, id_bgc_seaice, id_bgc_atmosp, id_pco2atm, id_surf_ph_in
      type (type_surface_diagnostic_variable_id) :: id_surf_ph, id_gasex

      ! Parameters
      real(rk) :: ocmip_alkfac
   contains
      ! Model procedures
      procedure :: initialize
      procedure :: do_surface
   end type type_mops_carbon

contains

   subroutine initialize(self, configunit)
      class (type_mops_carbon), intent(inout), target :: self
      integer,                  intent(in)            :: configunit

      call self%register_state_variable(self%id_dic, 'c', 'mmol C/m3', 'dissolved inorganic carbon')

      call self%get_parameter(self%ocmip_alkfac, 'ocmip_alkfac', 'meq/m3/PSU', 'alkalinity relative to salinity', default=2310.0_rk*1.0245_rk/34.88_rk)

      call self%register_diagnostic_variable(self%id_surf_ph, 'surf_ph', '-', 'surface pH', missing_value=8.0_rk)
      call self%register_diagnostic_variable(self%id_gasex, 'gasex', 'mmol C/m2/d', 'air-sea exchange of CO2')

      ! Register environmental dependencies
      call self%register_dependency(self%id_pho, 'pho', 'mmol P m-3', 'phosphate')
      call self%register_dependency(self%id_sil, 'sil', 'mmol Si m-3', 'silicate')
      call self%register_dependency(self%id_bgc_salt, standard_variables%practical_salinity)
      call self%register_dependency(self%id_bgc_theta, standard_variables%temperature)
      call self%register_dependency(self%id_bgc_wind, standard_variables%wind_speed)
      call self%register_dependency(self%id_bgc_seaice, standard_variables%ice_area_fraction)
      call self%register_dependency(self%id_bgc_atmosp, standard_variables%surface_air_pressure)
      call self%register_dependency(self%id_pco2atm, standard_variables%mole_fraction_of_carbon_dioxide_in_air)
      call self%register_dependency(self%id_surf_ph_in, 'surf_ph', '-', 'previous surface pH')

      self%dt = 86400.0_rk
   end subroutine

   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      class (type_mops_carbon), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: surf_dic, surf_pho, surf_sil, bgc_salt, bgc_theta, surf_alk, co2gasex
      real(rk) :: bgc_wind, bgc_seaice, bgc_atmosp, pco2atm, vgas660, surf_ph

      _SURFACE_LOOP_BEGIN_

         !CALL CAR_COEFFS(bgc_theta(1),bgc_salt(1),1)

   ! AIR-SEA GAS EXCHANGE OF CO2
         _GET_(self%id_dic, surf_dic)
         _GET_(self%id_pho, surf_pho)
         _GET_(self%id_sil, surf_sil)  ! normally constant Surface silicate from the OCMIP protocol
         _GET_(self%id_bgc_salt, bgc_salt)
         _GET_(self%id_bgc_theta, bgc_theta)
         _GET_SURFACE_(self%id_bgc_wind, bgc_wind)
         _GET_SURFACE_(self%id_bgc_seaice, bgc_seaice)
         _GET_SURFACE_(self%id_bgc_atmosp, bgc_atmosp)
         _GET_SURFACE_(self%id_pco2atm, pco2atm)
         _GET_SURFACE_(self%id_surf_ph_in, surf_ph)

   ! VS nur kurz? convert bgc_atmosp from Pa to bar
!         bgc_atmosp = bgc_atmosp / 101325.0_rk   ! from Pa to atm
         bgc_atmosp = bgc_atmosp / 100000.0_rk   ! from Pa to atm
   ! Surface total alkalinity from the OCMIP protocol
         surf_alk = self%ocmip_alkfac*bgc_salt

         vgas660=(0.337_rk*bgc_wind**2)*0.24_rk*(1.0_rk-bgc_seaice)
         CALL CO2_SURFFORCING(vgas660,bgc_atmosp, &
             surf_dic,surf_pho,surf_alk,surf_sil,bgc_theta,bgc_salt,pco2atm, &
             surf_ph,co2gasex)

! VS SETTING SURFACE FLUX TO ZERO
!         _ADD_SURFACE_FLUX_(self%id_dic, co2gasex)
         _SET_SURFACE_DIAGNOSTIC_(self%id_gasex, co2gasex)
         _SET_SURFACE_DIAGNOSTIC_(self%id_surf_ph, surf_ph)
      _SURFACE_LOOP_END_
   end subroutine do_surface


   elemental subroutine car_coeffs(t,s,bt,st,ft,ff,ak0,ak1,ak2,ak1p,ak2p,ak3p,aksi,akw,aks,akf,akb)
      real(rk), intent(in) :: t,s
      real(rk), intent(out) :: bt,st,ft,ff,ak0,ak1,ak2,ak1p,ak2p,ak3p,aksi,akw,aks,akf,akb

! local coefficients

      real(rk) :: tk,tk100,tk1002,invtk,dlogtk
      real(rk) :: is,is2,sqrtis,s2,sqrts,s15,scl

      tk = 273.15_rk + t
      tk100 = tk/100.0_rk
      tk1002=tk100*tk100
      invtk=1.0_rk/tk
      dlogtk=log(tk)

      is=19.924_rk*s/(1000.0_rk-1.005_rk*s)
      is2=is*is
      sqrtis=sqrt(is)
      s2=s*s
      sqrts=sqrt(s)
      s15=s**1.5_rk
      scl=s/1.80655_rk

! Calculate concentrations for borate, sulfate, and fluoride
! Uppstrom (1974), Morris & Riley (1966), Riley (1965)
      bt = 0.000232_rk * scl/10.811_rk
      st = 0.14_rk * scl/96.062_rk
      ft = 0.000067_rk * scl/18.9984_rk

! f = k0(1-pH2O)*correction term for non-ideality
! Weiss & Price (1980, Mar. Chem., 8, 347-359; Eq 13 with table 6 values)
      ff = exp(-162.8301_rk + 218.2968_rk/tk100  + &
             90.9241_rk*log(tk100) - 1.47696_rk*tk1002 + &
             s * (.025695_rk - .025225_rk*tk100 +  &
             0.0049867_rk*tk1002))


! K0 from Weiss 1974
      ak0 = exp(93.4517_rk/tk100 - 60.2409_rk +  &
           23.3585_rk * log(tk100) + &
           s * (0.023517_rk - 0.023656_rk*tk100 +  &
           0.0047036_rk*tk1002))


! k1 = [H][HCO3]/[H2CO3]
! k2 = [H][CO3]/[HCO3]     on hSWS
! Millero p.664 (1995) using Mehrbach et al. data on SEAWATER scale 
! (Original reference: Dickson and Millero, DSR, 1987)
      ak1=10**(-1._rk*(3670.7_rk*invtk -  &
             62.008_rk + 9.7944_rk*dlogtk - &
             0.0118_rk*s + 0.000116_rk*s2))
      ak2=10**(-1._rk*(1394.7_rk*invtk + 4.777_rk -  &
             0.0184_rk*s + 0.000118_rk*s2))


! k1p = [H][H2PO4]/[H3PO4] on hSWS
! Millero p.670 (1995)
      ak1p = exp(-4576.752_rk*invtk + 115.540_rk -  &
             18.453_rk*dlogtk +  &
   		    (-106.736_rk*invtk + 0.69171_rk)*sqrts + &
   		    (-0.65643_rk*invtk - 0.01844_rk)*s)


! k2p = [H][HPO4]/[H2PO4] on hSWS
! Millero p.670 (1995)
      ak2p = exp(-8814.715_rk*invtk + 172.1033_rk -  &
             27.927_rk*dlogtk + &
   		    (-160.340_rk*invtk + 1.3566_rk)*sqrts + &
   		    (0.37335_rk*invtk - 0.05778_rk)*s)


! k3p = [H][PO4]/[HPO4] on hSWS
! Millero p.670 (1995)
	   ak3p = exp(-3070.75_rk*invtk - 18.126_rk +  &
   		    (17.27039_rk*invtk + 2.81197_rk) * &
   		    sqrts + (-44.99486_rk*invtk - 0.09984_rk) * s)


! ksi = [H][SiO(OH)3]/[Si(OH)4] on hSWS
! Millero p.671 (1995) using data from Yao and Millero (1995)
! change to (mol/ kg soln)
      aksi = exp(-8904.2_rk*invtk + 117.400_rk -  &
             19.334_rk*dlogtk + &
   		    (-458.79_rk*invtk + 3.5913_rk) * sqrtis + &
   		    (188.74_rk*invtk - 1.5998_rk) * is + &
   		    (-12.1652_rk*invtk + 0.07871_rk) * is2 + &
   		    log(1._rk-0.001005_rk*s))


! kw = [H][OH] on hSWS
! Millero p.670 (1995) using composite data
      akw = exp(-13847.26_rk*invtk + 148.9802_rk -  &
             23.6521_rk*dlogtk + &
   		    (118.67_rk*invtk - 5.977_rk + 1.0495_rk * dlogtk) * &
   		    sqrts - 0.01615_rk * s)


! ks = [H][SO4]/[HSO4] on free H scale
! Dickson (1990, J. chem. Thermodynamics 22, 113)
! change to (mol/ kg soln)
      aks=exp(-4276.1_rk*invtk + 141.328_rk -  &
             23.093_rk*dlogtk + &
   		    (-13856._rk*invtk + 324.57_rk - 47.986_rk*dlogtk)*sqrtis + &
      		(35474._rk*invtk - 771.54_rk + 114.723_rk*dlogtk)*is - &
     		2698._rk*invtk*is**1.5_rk + 1776._rk*invtk*is2 + &
   		    log(1._rk - 0.001005_rk*s))


! kf = [H][F]/[HF] on free H scale
! Dickson and Riley (1979)
! change to (mol/ kg soln)
      akf=exp(1590.2_rk*invtk - 12.641_rk + 1.525_rk*sqrtis + &
   		    log(1._rk - 0.001005_rk*s)) 


! kb = [H][BO2]/[HBO2] on hSWS
! Dickson p.673 (1990)
! change from htotal to hSWS
      akb=exp( (-8966.90_rk - 2890.53_rk*sqrts - 77.942_rk*s + &
   		    1.728_rk*s15 - 0.0996_rk*s2)*invtk + &
     		(148.0248_rk + 137.1942_rk*sqrts + 1.62142_rk*s) + &
     		(-24.4344_rk - 25.085_rk*sqrts - 0.2474_rk*s) * &
     		dlogtk + 0.053105_rk*sqrts*tk + &
             log((1._rk+(st/aks)+(ft/akf))  &
             /(1._rk+(st/aks))) )
   end subroutine

   subroutine co2_surfforcing(vgas660,atmosp, &
                surf_dic,surf_pho,surf_alk,surf_sil,ttemp,stemp,pco2atm,surf_ph,co2ex)

      real(rk), intent(in) :: surf_dic,surf_pho,surf_alk,surf_sil
      real(rk), intent(in) :: vgas660,atmosp,ttemp,stemp,pco2atm
      real(rk), intent(inout) :: surf_ph
      real(rk), intent(out) :: co2ex

      real(rk), parameter :: scar1 = 2073.1_rk
      real(rk), parameter :: scar2 = 125.62_rk
      real(rk), parameter :: scar3 = 3.6276_rk
      real(rk), parameter :: scar4 = 0.043219_rk

! local coefficients

      real(rk) :: SchmidtNoCO2,kwexch,co2starair,co2star,co2sol
      real(rk) :: sdic,spho,ssil,salk   

      SchmidtNoCO2 = scar1 - scar2*ttemp + scar3*ttemp*ttemp &
     	 - scar4*ttemp*ttemp*ttemp
	
      KWexch = vgas660/sqrt(SchmidtNoCO2/660.0_rk)

! calculate co2star = pCO2 in the surface water 
! calculate co2sol = solubility of CO2 in water


      sdic = surf_dic/convert_mol_to_mmol
      spho = surf_pho/convert_mol_to_mmol
      ssil = surf_sil/convert_mol_to_mmol
      salk = surf_alk/convert_mol_to_mmol

      call co2_surface(ttemp,stemp,sdic,spho,ssil,salk,surf_ph,co2star,co2sol)
      
! sol is solubility of CO2 in mol/(m3*uatm)
! equilibrium [CO2]aq in mol/m^3 = sol*pCO2_atm*atmpres, where
! pCO2_atm = atmospheric mole fraction CO2 in dry air at 1 atm total pres (ppmv)
! atmpres= atmospheric pressure in atmospheres (1 atm==1013.25mbar)

      co2starair = co2sol*pco2atm*atmosp ! equilibrium CO2aq in mol/m^3
      co2ex=-KWexch*(co2star - co2starair)*convert_mol_to_mmol
!      !VS nur kurz
!      print *, 'KWexch is ', KWexch
!      print *, 'co2star is ', co2star
!      print *, 'co2starair is ', co2starair
!      print *, 'convert_mol_to_mmol ', convert_mol_to_mmol
!      print *, 'pco2atm is ', pco2atm
!      print *, 'atmops is ', atmosp
!      print *, 'co2ex is ', co2ex

   end subroutine


   subroutine co2_surface(t,s,sdic,spho,ssil,sta,sph,co2s,sol)
      real(rk), intent(in) :: t,s,sdic,spho,ssil,sta
      real(rk), intent(inout) :: sph
      real(rk), intent(out) :: co2s,sol

! local coefficients
      
      real(rk) :: pHlocal,btlocal,kblocal,k1local,k2local, &
            k1plocal,k2plocal,k3plocal,ksilocal,kwlocal,fflocal
      real(rk) :: pt,sit,ta,dic
      real(rk) :: phguess,hguess,bohg,stuff,h3po4g,h2po4g,hpo4g,po4g, &
            siooh3g,cag,gamm,hnew
      real(rk) :: stlocal, ftlocal, kslocal, kflocal, k0local

      call car_coeffs(t,s,btlocal,stlocal,ftlocal,fflocal,k0local,k1local,k2local,k1plocal,k2plocal,k3plocal,ksilocal,kwlocal,kslocal,kflocal,kblocal)
      pHlocal = sph

! change units from the input of mol/m^3 -> mol/kg:
! (1 mol/m^3)  x (1 m^3/1024.5 kg)
      pt=spho*permil
      sit=ssil*permil
      ta=sta*permil
      dic=sdic*permil

! first guess for ph, hydrogen ions, borate, phosphate, silicate
      phguess = phlocal
      hguess = 10.0_rk**(-phguess)
      bohg = btlocal*kblocal/(hguess+kblocal)
      stuff = hguess*hguess*hguess &
                + (k1plocal*hguess*hguess) &
                + (k1plocal*k2plocal*hguess) &
                + (k1plocal*k2plocal*k3plocal)
      h3po4g = (pt*hguess*hguess*hguess) / stuff
      h2po4g = (pt*k1plocal*hguess*hguess) / stuff
      hpo4g  = (pt*k1plocal*k2plocal*hguess) / stuff
      po4g   = (pt*k1plocal*k2plocal*k3plocal) / stuff
      siooh3g = sit*ksilocal / (ksilocal + hguess)

! estimate carbonate alkalinity
      cag = ta - bohg - (kwlocal/hguess) + hguess &
                - hpo4g - 2.0_rk*po4g + h3po4g &
                - siooh3g

! second guess of hydrogen ions
      gamm  = dic/cag
      stuff = (1.0_rk-gamm)*(1.0_rk-gamm)*k1local*k1local &
               - 4.0_rk*k1local*k2local*(1.0_rk-2.0_rk*gamm)
      hnew  = 0.5_rk*( (gamm-1.0_rk)*k1local + sqrt(stuff) )

! co2*
      co2s  = dic/ &
        (1.0_rk + (k1local/hnew) + (k1local*k2local/(hnew*hnew)))
! pH
      sph = -log10(hnew)

! co2* converted from mol/kg to mol/m3 
      co2s = co2s/permil
! VS nur kurz
!      print *, 'co2s is ', co2s

! fflocal is the solubility (computed in car_coeffs) in mol/(kg*atm)
! To convert to mol/(m^3*uatm), multiply ff by 1e-6*1024.5, i.e.

! solubility of CO2 in mol/(m3*uatm)
      sol=fflocal*permeg*rho0 

! equilibrium [CO2]aq in mol/m^3 = sol*pCO2_atm*atmpres, where
! pCO2_atm = atmospheric mole fraction CO2 in dry air at 1 atm total pres (ppmv)
! atmpres= atmospheric pressure in atmospheres (1 atm==1013.25mbar)

   end subroutine

end module mops_carbon
