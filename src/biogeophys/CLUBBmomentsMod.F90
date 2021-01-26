module CLUBBmomentsMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module computes the surface higher order moments that the CLUBB
  ! scheme in CAM needs. These are computed accounting for land surface 
  ! heterogeneity, so all moments are first computed at the patch level 
  ! before scaling up to gridcell levels. 
  ! 
  !
  ! !USES:
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use clm_varctl             , only : iulog
  use decompMod              , only : bounds_type
  !
  use FrictionvelocityMod    , only : frictionvel_type
  use EnergyFluxType         , only : energyflux_type
  use clm_varcon             , only : cpair,grav
  !
  use subgridAveMod          , only : p2g
  use GridcellType           , only : grc
  use LandunitType           , only : lun
  use ColumnType             , only : col
  use PatchType              , only : patch

  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CLUBBmoments            ! Calculate moments 

  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CLUBBmoments(bounds,frictionvel_inst,energyflux_inst)
    !
    ! !DESCRIPTION:
    !
    ! Compute moments that CLUBB uses at the surface:
    !    Start with wp2 (vertical variance) 
    !    Next up is thlp2 (temeprature variance).
    !
    ! ARGUMENTS:
    type(bounds_type)                      , intent(in)            :: bounds
    type(frictionvel_type)                 , intent(inout)         :: frictionvel_inst
    type(energyflux_type)                  , intent(inout)         :: energyflux_inst
    !  !LOCAL VARIABLES:
    integer  :: p                                        ! Patch indices
    integer  :: c                                        ! Column indices
    real(r8) :: cpair                                    ! heat capacity dry airat const pres (j/kg/k)
    real(r8) :: rhoair                                   ! Air density 
    real(r8) :: wp2(bounds%begp:bounds%endp)             ! Vertical velocityvariance (m2/s2)
    real(r8) :: thlp2(bounds%begp:bounds%endp)           ! Temperature variance (K2)
    real(r8) :: wp2_grid(bounds%begg:bounds%endg)        ! Gridcell weighted mean of WP2 

    ! real(r8) :: sumwt(bounds%begg:bounds%endg)     ! sum of weights
    !------------------------------------------------------------------------


   !  sumwt(bounds%begg : bounds%endg) = 0._r8
    do p = bounds%begp,bounds%endp
       if (patch%active(p)) then

          associate(                                                      &
             ustar                  => frictionvel_inst%ustar_patch  , & ! Output: [real(r8) (:)   ]  friction velocity [m/s]  
             zeta                   => frictionvel_inst%zeta_patch ) ! Output: [real(r8) (:)   ]  dimensionless stability parameter 

            ! Compute the variance of vertical velocity for each patch 
            if (zeta(p) >= 0._r8) then  
               wp2(p)=1.75_r8*(ustar(p)**2.0_r8)
            else
               wp2(p)=(ustar(p)**2.0_r8)*(1.75_r8+(2.0_r8*(-zeta(p))**(2.0_r8/3.0_r8)))
            end if  

            ! Compute the variance of potential temperature for each patch 
            if (zeta(p)>=0.0_r8) then 
               thlp2 = 
            else 

            end if 

            !write(iulog,*)'MDF: Testing more - computed value of wp2: ',wp2(p)
            !write(iulog,*)'MDF: reading values... ustar= ',ustar(p)
            !write(iulog,*)'MDF: reading values... zeta=  ',zeta(p)
            !write(iulog,*)'MDF: Weight of patch: ',patch%wtgcell(p)
            !write(iulog,*)'MDF: Done with timestep for CLUBBmomentsMod...'
          end associate
      end if 
    end do      

    ! Get weighted mean  
    call p2g(bounds, &
         wp2(bounds%begp:bounds%endp), & 
         wp2_grid(bounds%begg:bounds%endg), & 
         p2c_scale_type='unity', c2l_scale_type= 'unity',l2g_scale_type='unity')    

    write(iulog,*)'MDF: Grid-mean wp2: ',wp2_grid


  end subroutine CLUBBmoments

end module CLUBBmomentsMod
