module CLUBBmomentsMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module computes the surface higher order moments that the CLUBB
  ! scheme in CAM needs. These are computed accounting for land surface 
  ! heterogeneity, so all moments are first computed at the patch level 
  ! before scaling up to gridcell levels, according to Machulskaya and
  ! Mironov (2018), which is based on the original Andre et al. 
  ! 1978 equations. 
  ! 
  !
  ! !USES:
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use clm_varctl             , only : iulog
  use decompMod              , only : bounds_type
  !
  use CLUBBMomentsType       , only : clubbmoments_type
  use FrictionvelocityMod    , only : frictionvel_type
  use EnergyFluxType         , only : energyflux_type
  use atm2lndType            , only : atm2lnd_type
  use lnd2atmType            , only : lnd2atm_type
  use WaterDiagnosticBulkType , only : waterdiagnosticbulk_type
  use TemperatureType        , only : temperature_type
  use clm_varcon             , only : cpair,grav
  !
  use subgridAveMod          , only : p2g
  use GridcellType           , only : grc
  use LandunitType           , only : lun
  use ColumnType             , only : col
  use PatchType              , only : patch
  ! 
  ! use histFileMod            , only : hist_addfld1d
  
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CLUBBmoments            ! Calculate moments 

  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CLUBBmoments(bounds,clubbmoments_inst,frictionvel_inst,&
       energyflux_inst,atm2lnd_inst,lnd2atm_inst,&
       waterdiagnosticbulk_inst,temperature_inst)
    !
    ! !DESCRIPTION:
    !
    ! Compute moments that CLUBB uses at the surface:
    !    Start with wp2 (vertical variance) 
    !    Next up is thlp2 (temeprature variance).
    !
    ! ARGUMENTS:
    type(bounds_type)                      , intent(in)            :: bounds
    type(clubbmoments_type)                , intent(inout)         :: clubbmoments_inst
    type(frictionvel_type)                 , intent(in)            :: frictionvel_inst
    type(energyflux_type)                  , intent(in)            :: energyflux_inst
    type(atm2lnd_type)                     , intent(in)            :: atm2lnd_inst
    type(lnd2atm_type)                     , intent(in)            :: lnd2atm_inst
    type(waterdiagnosticbulk_type)         , intent(in)            :: waterdiagnosticbulk_inst
    type(temperature_type)                 , intent(in)            :: temperature_inst
    !  !LOCAL VARIABLES:
    integer  :: p                                        ! Patch indices
    integer  :: c,g                                      ! Column/gridcell indices
    real(r8) :: rhoair                                   ! Air density 
    real(r8) :: eflx_sh_pr_conversionPatch               ! Patch level eflx_sh correction term 
    real(r8) :: eflx_sh_il_conversionPatch               ! Patch level eflx_sh correction term
    real(r8) :: eflx_dynbal_grcPatch                     ! Patch level eflx_sh correction term
    real(r8) :: eflx_sh_patchCorrected                   ! Total patch eflx_sh, including correction terms
    real(r8) :: KH                                       ! Kinematic heat flux 

    real(r8) :: wp2(bounds%begp:bounds%endp)             ! Vertical velocityvariance (m2/s2)
    real(r8) :: thlp2(bounds%begp:bounds%endp)           ! Temperature variance (K2)
    real(r8) :: Tv(bounds%begp:bounds%endp)              ! Virtual temperature at the surface (column=patch level)
    !real(r8) :: wp2_weighted(bounds%begg:bounds%endg)    ! FINAL OUTPUT: Gridcell weighted mean of WP2 
    real(r8) :: thlp2_grid(bounds%begg:bounds%endg)      ! Gridcell weighted mean of THLP2
    real(r8) :: theta_grid(bounds%begg:bounds%endg)      ! Gridcell weighted mean of THETA (Tv)
    real(r8) :: theta_sqr(bounds%begp:bounds%endp)       ! Patch level difference in THETA from grid mean 
    real(r8) :: theta_sqr_grid(bounds%begg:bounds%endg)  ! Gridcell weighted mean of patch THETA differences 
    !real(r8) :: thlp2_weighted(bounds%begg:bounds%endg)  ! FINAL OUTPUT: Gridcell weighted mean  of THLP2 
   

    ! real(r8) :: sumwt(bounds%begg:bounds%endg)     ! sum of weights
    !------------------------------------------------------------------------

       do p = bounds%begp,bounds%endp
          if (patch%active(p)) then

             associate(                                                        &
                !wp2_weighted          => clubbmoments_inst%wp2_grid          ,  &      ! Output [real(r8) (:)   ]  vertical velocity variance [m**2/s**2]
                !thlp2_weighted        => clubbmoments_inst%thlp2_grid        ,  &      ! Output [real(r8) (:)   ]  temperature variance [K**2]  
                ustar                 => frictionvel_inst%ustar_patch        ,  &      ! Input: [real(r8) (:)   ]  friction velocity [m/s]  
                zeta                  => frictionvel_inst%zeta_patch         ,  &      ! Input: [real(r8) (:)   ]  dimensionless stability parameter 
                forc_rho              => atm2lnd_inst%forc_rho_downscaled_col,  &      ! Input: [real(r8) (:)   ]  density (kg/m**3)   
                q_ref2m               => waterdiagnosticbulk_inst%q_ref2m_patch,        &  ! Input: [real(r8) (:)   ]  2 m height surface specific humidity (kg/kg)
                t_ref2m               => temperature_inst%t_ref2m_patch,                &  ! Input: [real(r8) (:)   ]  2 m height surface air temperature (Kelvin)
                eflx_sh_tot           => energyflux_inst%eflx_sh_tot_patch   ,          &  ! Input:  [real(r8) (:)   ]  total sensible heat flux (W/m**2) [+ to atm]
                eflx_sh_pr_conversion => energyflux_inst%eflx_sh_precip_conversion_col, &  ! Input: [real(r8) (:) ] correction to eflx_sh
                eflx_sh_il_conversion => lnd2atm_inst%eflx_sh_ice_to_liq_col,           &  ! Input: [real(r8) (:) ] correction to eflx_sh  
                eflx_dynbal_grc       => energyflux_inst%eflx_dynbal_grc )                 ! Input: [real(r8) (:) ] correction eflx_sh 

               ! Compute the variance of vertical velocity for each patch 
               ! --------------------------------------------------------
               if (zeta(p) >= 0._r8) then  
                  wp2(p)=1.75_r8*(ustar(p)**2.0_r8)
               else
                  wp2(p)=(ustar(p)**2.0_r8)*(1.75_r8+(2.0_r8*(-zeta(p))**(2.0_r8/3.0_r8)))
               end if  

               !write(iulog,*)'MDF: Testing more - computed value of wp2: ',wp2(p)
               !write(iulog,*)'MDF: reading values... ustar= ',ustar(p)
               !write(iulog,*)'MDF: reading values... zeta=  ',zeta(p)

               ! Compute the variance of potential temperature for each patch
               ! ------------------------------------------------------------ 
               c = patch%column(p) 
               g = patch%gridcell(p)
    
               ! Use column or gridcell values for each patch 
               rhoair = forc_rho(c)    

               ! Correct total SHFLX as in lnd2atmMod
               eflx_sh_pr_conversionPatch = eflx_sh_pr_conversion(c) 
               eflx_sh_il_conversionPatch = eflx_sh_il_conversion(c)
               eflx_dynbal_grcPatch       = eflx_dynbal_grc(g) 
               eflx_sh_patchCorrected = eflx_sh_tot(p) + eflx_sh_pr_conversionPatch + eflx_sh_il_conversionPatch - eflx_dynbal_grcPatch
              
               ! Compute kinematic heat flux based on eflx_sh 
               KH = eflx_sh_patchCorrected / rhoair / cpair 

               !write(iulog,*)'MDF: patch type... ',patch%itype(p)        
               !write(iulog,*)'MDF: col type...   ',col%itype(c)
               !write(iulog,*)'MDF: patch weight: ',patch%wtgcell(p)

               if (zeta(p)>=0.0_r8) then 
                  thlp2(p) = (KH**2.0_r8/ustar(p)**2.0_r8)*4.0_r8
               else 
                  thlp2(p) = (KH**2.0_r8/ustar(p)**2.0_r8)*(4.0_r8*((1.0_r8 - 8.3_r8*zeta(p))**(-2.0_r8/3.0_r8)) )
               end if 

               ! Also compute virtual temperature, which we'll use as 'theta' 
               Tv(p)       = (1.0_r8 + 0.61_r8*q_ref2m(p))*t_ref2m(p)

               !write(iulog,*)'MDF: Testing thlp2 - computed value: ',thlp2(p)
               !write(iulog,*)'MDF:    value of KH:',KH
               !write(iulog,*)'MDF:    ustar:      ',ustar(p)
               !write(iulog,*)'MDF:    zeta:       ',zeta(p) 
               !write(iulog,*)'MDF: Testing theta - computed value: ', Tv(p)
               !write(iulog,*)'MDF:       qsfc / temp  ',q_ref2m(p),t_ref2m(p)
               !write(iulog,*)'MDF: rhoair:            ',rhoair 
               !write(iulog,*)'MDF: eflux_sh_total:    ',eflx_sh_patchCorrected
             end associate
         end if 
       end do      

    ! Get weighted means
    ! -----------------
    call p2g(bounds, &
         wp2(bounds%begp:bounds%endp), & 
         clubbmoments_inst%wp2_grid (bounds%begg:bounds%endg), & 
         p2c_scale_type='unity', c2l_scale_type= 'unity',l2g_scale_type='unity')    

!    call p2g(bounds, &
!         wp2(bounds%begp:bounds%endp), &
!         wp2_weighted(bounds%begg:bounds%endg), &
!         p2c_scale_type='unity', c2l_scale_type= 'unity',l2g_scale_type='unity')

    call  p2g(bounds, & 
          thlp2(bounds%begp:bounds%endp), &
          thlp2_grid(bounds%begg:bounds%endg), &
          p2c_scale_type='unity', c2l_scale_type='unity',l2g_scale_type='unity')

    call  p2g(bounds, &
          Tv(bounds%begp:bounds%endp), &
          theta_grid(bounds%begg:bounds%endg), &
          p2c_scale_type='unity', c2l_scale_type='unity',l2g_scale_type='unity')

    ! Some moments need a second term to represent SGS heterogeneity
    !   (The difference in patch values from gridcell mean) 
    do p = bounds%begp,bounds%endp
          if (patch%active(p)) then
             g = patch%gridcell(p)
             theta_sqr(p) = (Tv(p) - theta_grid(g))**2.0_r8
             !write(iulog,*)'MDF: Patch-level theta diff from grid mean theta: ',theta_sqr(p)
          end if
    end do 

    call  p2g(bounds, &
          theta_sqr(bounds%begp:bounds%endp), &
          theta_sqr_grid(bounds%begg:bounds%endg), &
          p2c_scale_type='unity', c2l_scale_type='unity',l2g_scale_type='unity')

    clubbmoments_inst%thlp2_grid(bounds%begg:bounds%endg) = thlp2_grid(bounds%begg:bounds%endg) + theta_sqr_grid(bounds%begg:bounds%endg)

    write(iulog,*)'MDF: Final value of wp2: ',clubbmoments_inst%wp2_grid
    !write(iulog,*)'MDF: Grid-mean thlp2: ',thlp2_grid
    !write(iulog,*)'MDF: Grid-mean theta: ',theta_grid
    !write(iulog,*)'MDF: Grid-mean theta diff: ',theta_sqr_grid 
    write(iulog,*)'MDF: Final value of thlp2: ',clubbmoments_inst%thlp2_grid

    ! Save files in way that can be written out to history file 
    ! ---------------------------------------------------------
    !   Not entirely sure how to do this...  



  end subroutine CLUBBmoments

end module CLUBBmomentsMod
