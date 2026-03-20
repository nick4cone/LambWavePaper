module ic_baro_dry_jw06
  !-----------------------------------------------------------------------
  !
  ! Purpose: Set idealized initial conditions for an isothermal acid test with a cosine mountain
  !
  !
  !-----------------------------------------------------------------------
  use cam_logfile,         only: iulog
  use shr_kind_mod,        only: r8 => shr_kind_r8
  use cam_abortutils,      only: endrun
  use spmd_utils,          only: masterproc
  use shr_sys_mod,         only: shr_sys_flush

  use physconst, only : rair, cpair, gravit, rearth, pi, omega
  use hycoef,    only : hyai, hybi, hyam, hybm, ps0

  implicit none
  private

  !=======================================================================
  !    Test case parameters: Acid test
  !=======================================================================
  real(r8), parameter, private ::             &
       T0                     = 275._r8,      & ! isothermal T in K
       p00                    = 1.e5_r8,      & ! reference surface pressure in Pa
       p_half                 = 1._r8/30._r8, & ! halfwidth of a Gaussian pressure perturbation without the Earth's radius=1 (a/30 with a=1) 
       pressure_amplitude     = 5000._r8,     & ! amplitude of the surface pressure perturbation 5000 Pa 
       mountain_amplitude     = 2000._r8,     & ! mountain amplitude of 2000 m
       mountain_radius        = pi/30._r8,    & ! width of the cosine mountain 
     ! center_longitude       = 180._r8,      & ! longitudinal center position in degrees, 180E
       center_longitude       = 270._r8,      & ! longitudinal center position in degrees, 270E
       center_latitude        = 0._r8           ! latitudinal center position in degrees, 0N (equator)

  real(r8), parameter :: deg2rad = pi/180._r8   ! conversion to radians

  ! Public interface
  public :: bc_dry_jw06_set_ic

contains

  subroutine bc_dry_jw06_set_ic(vcoord, latvals, lonvals, U, V, T, PS, PHIS, &
                                Q, m_cnst, mask, verbose)
    use dyn_tests_utils, only: vc_moist_pressure, vc_dry_pressure, vc_height
    use constituents,    only: cnst_name
    use const_init,      only: cnst_init_default

    !-----------------------------------------------------------------------
    !
    ! Purpose: Set initial values for the dry acid test
    !
    !-----------------------------------------------------------------------

    ! Dummy arguments
    integer, intent(in)               :: vcoord
    real(r8),           intent(in)    :: latvals(:) ! lat in degrees (ncol)
    real(r8),           intent(in)    :: lonvals(:) ! lon in degrees (ncol)
                                                    ! z_k for vccord 1)
    real(r8), optional, intent(inout) :: U(:,:)     ! zonal velocity
    real(r8), optional, intent(inout) :: V(:,:)     ! meridional velocity
    real(r8), optional, intent(inout) :: T(:,:)     ! temperature
    real(r8), optional, intent(inout) :: PS(:)      ! surface pressure
    real(r8), optional, intent(out)   :: PHIS(:)    ! surface geopotential
    real(r8), optional, intent(inout) :: Q(:,:,:)   ! tracer (ncol, lev, m)
    integer,  optional, intent(in)    :: m_cnst(:)  ! tracer indices (reqd. if Q)
    logical,  optional, intent(in)    :: mask(:)    ! Only init where .true.
    logical,  optional, intent(in)    :: verbose    ! For internal use
    ! Local variables
    logical, allocatable              :: mask_use(:)
    logical                           :: verbose_use
    logical                           :: lu,lv,lt,lq,l3d_vars
    integer                           :: i, k, m
    integer                           :: ncol
    integer                           :: nlev
    integer                           :: ncnst
    character(len=*), parameter       :: subname = 'BC_DRY_JW06_SET_IC'
    real(r8)                          :: r(size(latvals))         ! great circle distance (unit circle)
    real(r8)                          :: surface_pressure(size(latvals))
    real(r8)                          :: surface_height(size(latvals))
    real(r8)                          :: center_lon, center_lat
    real(r8)                          :: scale_height
    logical                           :: gaussian

    ! flag that distinguishes between the two perturbations
    ! gaussian          = .true.               ! flag that either selects the Gaussian pressure perturbation (.true.) or the cosine mountain (.false.)
    gaussian          = .false.                ! flag that either selects the Gaussian pressure perturbation (.true.) or the cosine mountain (.false.)

    scale_height = rair*T0/gravit ! scale height

    !*******************************
    ! Center position of the perturbation, convert to radians
    !*******************************
    center_lon = center_longitude * deg2rad    ! in radians
    center_lat = center_latitude  * deg2rad    ! in radians

    allocate(mask_use(size(latvals)))
    if (present(mask)) then
      if (size(mask_use) /= size(mask)) then
        call endrun(subname//': input, mask, is wrong size')
      end if
      mask_use = mask
    else
      mask_use = .true.
    end if

    if (present(verbose)) then
      verbose_use = verbose
    else
      verbose_use = .true.
    end if

    ncol = size(latvals, 1)
    nlev = -1

    !
    ! We do not yet handle height-based vertical coordinates
    if (vcoord == vc_height) then
      call endrun(subname//':  height-based vertical coordinate not currently supported')
    end if

    !*******************************
    !
    ! Initialize the surface height (topography) and surface pressure
    !
    !*******************************
    !
    if (gaussian) then
    !  Gaussian perturbation of the surface pressure without topography 
       where(mask_use)
    !    great circle distance without the Earth's radius (unit circle)
         r(:) = acos( sin(center_lat)*sin(latvals(:)) + cos(center_lat)*cos(latvals(:))*cos(lonvals(:)-center_lon))
         surface_height(:) = 0._r8
    !    surface pressure with an overlaid Gaussian pressure perturbation
         surface_pressure(:) = p00 + pressure_amplitude * exp(- (r(:)/p_half)**2._r8 )
       end where
    else 
    ! use a cosine mountain perturbation 
       do i = 1, size(latvals,1)
         if (mask_use(i)) then
    !       great circle distance without the Earth's radius (unit circle)
            r(i) = acos( sin(center_lat)*sin(latvals(i)) + cos(center_lat)*cos(latvals(i))*cos(lonvals(i)-center_lon))
            if (r(i) .lt. mountain_radius) then
              surface_height(i) = mountain_amplitude*0.5_r8*(1+cos(pi*r(i)/mountain_radius))
            else
              surface_height(i) = 0._r8
            endif
    !       surface pressure (in hydrostatic balance) for isothermal conditions 
            surface_pressure(i) = p00 * exp(-surface_height(i)/scale_height)
         endif
       end do
    endif

    !*******************************
    !
    ! Initialize PHIS
    !
    !*******************************
    !
     if (present(PHIS)) then
      where(mask_use)
      ! surface geopotential: mountain height times gravity 
        PHIS(:) = gravit*surface_height(:)
      end where
      if(masterproc .and. verbose_use) then
        write(iulog,*) '          PHIS initialized by "',subname,'"'
      end if
     end if
    !
    !*******************************
    !
    ! initialize surface pressure
    !
    !*******************************
    !
     if (present(PS)) then
      where(mask_use)
        PS(:) = surface_pressure(:)
      end where

      if(masterproc .and. verbose_use) then
        write(iulog,*) '          PS initialized by "',subname,'"'
      end if
     end if
    !
    !
    !*******************************
    !
    ! Initialize 3D vars
    !
    !
    !*******************************
    !
    lu = present(U)
    lv = present(V)
    lT = present(T)
    lq = present(Q)
    l3d_vars = lu .or. lv .or. lt .or.lq
    nlev = -1
    if (l3d_vars) then
      if (lu) nlev = size(U, 2)
      if (lv) nlev = size(V, 2)
      if (lt) nlev = size(T, 2)
      if (lq) nlev = size(Q, 2)

      if (lu) then
        do k = 1, nlev
          where(mask_use)
             U(:,k) = 0.0_r8
          end where
        end do
        if(masterproc.and. verbose_use) then
          write(iulog,*) '          U initialized by "',subname,'"'
        end if
      end if
      if (lv) then
        do k = 1, nlev
          where(mask_use)
            V(:,k) = 0.0_r8
          end where
        end do
        if(masterproc.and. verbose_use) then
          write(iulog,*) '          V initialized by "',subname,'"'
        end if
      end if
      if (lt) then
       do k = 1, nlev
         where(mask_use)
            T(:,k)   = T0 
         end where
       enddo
        if(masterproc.and. verbose_use) then
          write(iulog,*) '          T initialized by "',subname,'"'
        end if
      end if
      if (lq) then
        do k = 1, nlev
          where(mask_use)
            Q(:,k,1) = 0.0_r8
          end where
        end do
        if(masterproc.and. verbose_use) then
          write(iulog,*) '         ', trim(cnst_name(m_cnst(1))), ' initialized by "',subname,'"'
        end if
      end if
    end if

    if (lq) then
      ncnst = size(m_cnst, 1)
      if ((vcoord == vc_moist_pressure) .or. (vcoord == vc_dry_pressure)) then
        do m = 2, ncnst
          call cnst_init_default(m_cnst(m), latvals, lonvals, Q(:,:,m_cnst(m)),&
               mask=mask_use, verbose=verbose_use, notfound=.false.)
        end do
      end if
    end if

    deallocate(mask_use)

  end subroutine bc_dry_jw06_set_ic

end module ic_baro_dry_jw06
