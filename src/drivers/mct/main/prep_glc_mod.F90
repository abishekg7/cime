module prep_glc_mod

  use shr_kind_mod    , only: r8 => SHR_KIND_R8 
  use shr_kind_mod    , only: cs => SHR_KIND_CS
  use shr_kind_mod    , only: cl => SHR_KIND_CL
  use shr_sys_mod     , only: shr_sys_abort, shr_sys_flush
  use seq_comm_mct    , only: num_inst_glc, num_inst_lnd, num_inst_frc, &
                              num_inst_ocn
  use seq_comm_mct    , only: CPLID, GLCID, logunit
  use seq_comm_mct    , only: seq_comm_getData=>seq_comm_setptrs 
  use seq_infodata_mod, only: seq_infodata_type, seq_infodata_getdata  
  use seq_map_type_mod 
  use seq_map_mod
  use seq_flds_mod
  use t_drv_timers_mod
  use mct_mod
  use perf_mod
  use component_type_mod, only: component_get_x2c_cx, component_get_c2x_cx
  use component_type_mod, only: glc, lnd
  implicit none
  save
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: prep_glc_init
  public :: prep_glc_mrg_lnd

  public :: prep_glc_accum_lnd
  public :: prep_glc_accum_ocn
  public :: prep_glc_accum_avg

  public :: prep_glc_calc_l2x_gx
  public :: prep_glc_calc_o2x_gx

  public :: prep_glc_get_l2x_gx
  public :: prep_glc_get_l2gacc_lx
  public :: prep_glc_get_l2gacc_lx_cnt
  public :: prep_glc_get_mapper_SFl2g

  public :: prep_glc_get_mapper_So2g
  public :: prep_glc_get_mapper_Fo2g

  public :: prep_glc_calculate_subshelf_boundary_fluxes

  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------

  private :: prep_glc_merge

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  ! mappers
  type(seq_map), pointer :: mapper_SFl2g

  ! attribute vectors 
  type(mct_aVect), pointer :: l2x_gx(:) ! Lnd export, glc grid, cpl pes - allocated in driver
  type(mct_aVect), pointer :: o2x_gx(:) ! Ocn export, glc grid, cpl pes - allocated in driver

  ! accumulation variables
  
  type(mct_aVect), pointer :: x2gacc_gx(:) ! Glc export, glc grid, cpl pes - allocated in driver
  integer        , target :: x2gacc_gx_cnt ! x2gacc_gx: number of time samples accumulated 

  type(mct_aVect), pointer :: l2gacc_lx(:) ! Lnd export, lnd grid, cpl pes - allocated in driver
  integer        , target :: l2gacc_lx_cnt ! l2gacc_lx: number of time samples accumulated

  ! other module variables
  integer :: mpicom_CPLID  ! MPI cpl communicator
  !================================================================================================

contains

  !================================================================================================

  subroutine prep_glc_init(infodata, lnd_c2_glc, ocn_c2_glcshelf)

    !---------------------------------------------------------------
    ! Description
    ! Initialize module attribute vectors and mapping variables
    !
    ! Arguments
    type (seq_infodata_type) , intent(inout) :: infodata
    logical                  , intent(in)    :: lnd_c2_glc ! .true.  => lnd to glc coupling on
    logical                  , intent(in)    :: ocn_c2_glcshelf ! .true.  => ocn to glc coupling on
    !
    ! Local Variables
    integer                          :: eli, egi
    integer                          :: lsize_l
    integer                          :: lsize_g

    logical                          :: esmf_map_flag ! .true. => use esmf for mapping
    logical                          :: iamroot_CPLID ! .true. => CPLID masterproc
    character(CL)                    :: lnd_gnam      ! lnd grid
    character(CL)                    :: glc_gnam      ! glc grid
    character(CL)                    :: ocn_gnam      ! ocn grid

    type(mct_avect), pointer         :: l2x_lx
    type(mct_avect), pointer         :: x2g_gx
    type(mct_avect), pointer         :: o2x_ox

    character(*), parameter          :: subname = '(prep_glc_init)'
    character(*), parameter          :: F00 = "('"//subname//" : ', 4A )"
    !---------------------------------------------------------------

    call seq_infodata_getData(infodata , &
         esmf_map_flag=esmf_map_flag   , &
         glc_present=glc_present       , &
         lnd_gnam=lnd_gnam             , &
         glc_gnam=glc_gnam             , &
         ocn_gnam=ocn_gnam)

    allocate(mapper_SFl2g)

    if (glc_present) then

       call seq_comm_getData(CPLID, &
            mpicom=mpicom_CPLID, iamroot=iamroot_CPLID)

       l2x_lx => component_get_c2x_cx(lnd(1))
       lsize_l = mct_aVect_lsize(l2x_lx)
       
       x2g_gx => component_get_x2c_cx(glc(1))
       lsize_g = mct_aVect_lsize(x2g_gx)
       
       allocate(l2x_gx(num_inst_lnd))
       allocate(l2gacc_lx(num_inst_lnd))
       do eli = 1,num_inst_lnd
          call mct_aVect_initSharedFields(l2x_lx, x2g_gx, l2x_gx(eli) ,lsize=lsize_g)
          call mct_aVect_zero(l2x_gx(eli))
          
          call mct_aVect_initSharedFields(l2x_lx, x2g_gx, l2gacc_lx(eli), lsize=lsize_l)
          call mct_aVect_zero(l2gacc_lx(eli))
       enddo
       l2gacc_lx_cnt = 0
       
       if (iamroot_CPLID) then
          write(logunit,*) ' '
          write(logunit,F00) 'Initializing mapper_SFl2g'
       end if
       call seq_map_init_rearrolap(mapper_SFl2g, lnd(1), glc(1), 'mapper_SFl2g')
       call shr_sys_flush(logunit)

    end if

    if (glc_present .and. ocn_c2_glcshelf) then

       call seq_comm_getData(CPLID, &
            mpicom=mpicom_CPLID, iamroot=iamroot_CPLID)

       o2x_ox => component_get_c2x_cx(ocn(1))
       lsize_o = mct_aVect_lsize(o2x_ox)

       x2g_gx => component_get_x2c_cx(glc(1))
       lsize_g = mct_aVect_lsize(x2g_gx)

       allocate(o2x_gx(num_inst_ocn))
       do eoi = 1,num_inst_ocn
          call mct_aVect_init(o2x_gx(eoi), rList=seq_flds_o2x_fields, lsize=lsize_g)
          call mct_aVect_zero(o2x_gx(eoi))
       enddo

       allocate(x2gacc_gx(num_inst_glc))
       do egi = 1,num_inst_glc
          call mct_aVect_init(x2gacc_gx(egi), x2g_gx, lsize_g)
          call mct_aVect_zero(x2gacc_gx(egi))
       end do

       x2gacc_gx_cnt = 0
       samegrid_go = .true.
       if (trim(ocn_gnam) /= trim(glc_gnam)) samegrid_go = .false.
       if (iamroot_CPLID) then
          write(logunit,*) ' '
          write(logunit,F00) 'Initializing mapper_So2g'
       end if
       call seq_map_init_rcfile(mapper_So2g, ocn(1), glc(1), &
       'seq_maps.rc','ocn2glc_smapname:','ocn2glc_smaptype:',samegrid_go, &
       'mapper_So2g initialization',esmf_map_flag)
       if (iamroot_CPLID) then
          write(logunit,*) ' '
          write(logunit,F00) 'Initializing mapper_Fo2g'
       end if
       call seq_map_init_rcfile(mapper_Fo2g, ocn(1), glc(1), &
       'seq_maps.rc','ocn2glc_fmapname:','ocn2glc_fmaptype:',samegrid_go, &
       'mapper_Fo2g initialization',esmf_map_flag)

       !Initialize module-level arrays associated with compute_melt_fluxes
       allocate(oceanTemperature(lsize_g))
       allocate(oceanSalinity(lsize_g))
       allocate(oceanHeatTransferVelocity(lsize_g))
       allocate(oceanSaltTransferVelocity(lsize_g))
       allocate(interfacePressure(lsize_g))
       allocate(iceTemperature(lsize_g))
       allocate(iceTemperatureDistance(lsize_g))
       allocate(iceFloatingMask(lsize_g))
       allocate(outInterfaceSalinity(lsize_g))
       allocate(outInterfaceTemperature(lsize_g))
       allocate(outFreshwaterFlux(lsize_g))
       allocate(outOceanHeatFlux(lsize_g))
       allocate(outIceHeatFlux(lsize_g))
       ! TODO: Can we allocate these only while used or are we worried about performance hit?
       ! TODO: add deallocates!

       call shr_sys_flush(logunit)

    end if


  end subroutine prep_glc_init

  !================================================================================================

  subroutine prep_glc_accum(timer)

    !---------------------------------------------------------------
    ! Description
    ! Accumulate glc inputs from lnd
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eli
    type(mct_avect), pointer :: l2x_lx

    character(*), parameter :: subname = '(prep_glc_accum_lnd)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eli = 1,num_inst_lnd
       l2x_lx => component_get_c2x_cx(lnd(eli))
       if (l2gacc_lx_cnt == 0) then
          call mct_avect_copy(l2x_lx, l2gacc_lx(eli))
       else
          call mct_avect_accum(l2x_lx, l2gacc_lx(eli))
       endif
    end do
    l2gacc_lx_cnt = l2gacc_lx_cnt + 1
    call t_drvstopf  (trim(timer))

  end subroutine prep_glc_accum_lnd

  !================================================================================================

  subroutine prep_glc_accum_ocn(timer)

    !---------------------------------------------------------------
    ! Description
    ! Accumulate glc inputs from ocn
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: egi
    type(mct_avect), pointer :: x2g_gx

    character(*), parameter :: subname = '(prep_glc_accum_ocn)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do egi = 1,num_inst_glc
       x2g_gx => component_get_x2c_cx(glc(egi))
       if (x2gacc_gx_cnt == 0) then
          call mct_avect_copy(x2g_gx, x2gacc_gx(egi))
       else
          call mct_avect_accum(x2g_gx, x2gacc_gx(egi))
       endif
    end do
    x2gacc_gx_cnt = x2gacc_gx_cnt + 1
    call t_drvstopf  (trim(timer))

  end subroutine prep_glc_accum_ocn

  !================================================================================================


  subroutine prep_glc_accum_avg(timer)

    !---------------------------------------------------------------
    ! Description
    ! Finalize accumulation of glc inputs
    ! Note: There could be separate accum_avg routines for forcing coming
    ! from each component (LND and OCN), but they can be combined here
    ! by taking advantage of l2gacc_lx_cnt and x2gacc_gx_cnt variables
    ! that will only be greater than 0 if corresponding coupling is enabled.
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eli, egi
    type(mct_avect), pointer :: x2g_gx

    character(*), parameter :: subname = '(prep_glc_accum_avg)'
    !---------------------------------------------------------------

    ! Accumulation for LND
    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    if (l2gacc_lx_cnt > 1) then
       do eli = 1,num_inst_lnd
          call mct_avect_avg(l2gacc_lx(eli), l2gacc_lx_cnt)
       end do
    end if
    l2gacc_lx_cnt = 0

    ! Accumulation for OCN
    if (x2gacc_gx_cnt > 1) then
       do egi = 1,num_inst_glc
          ! temporary formation of average
          call mct_avect_avg(x2gacc_gx(egi), x2gacc_gx_cnt)

          ! ***NOTE***THE FOLLOWING ACTUALLY MODIFIES x2g_gx
          x2g_gx => component_get_x2c_cx(glc(egi))
          call mct_avect_copy(x2gacc_gx(egi), x2g_gx)
       enddo
    end if
    x2gacc_gx_cnt = 0

    call t_drvstopf  (trim(timer))
    
  end subroutine prep_glc_accum_avg

  !================================================================================================
  
  subroutine prep_glc_mrg(infodata, timer_mrg) 

    !---------------------------------------------------------------
    ! Description
    ! Merge glc inputs
    !
    ! Arguments
    type(seq_infodata_type) , intent(in)    :: infodata
    character(len=*)        , intent(in)    :: timer_mrg
    !
    ! Local Variables
    integer :: egi, eli
    type(mct_avect), pointer :: x2g_gx
    character(*), parameter  :: subname = '(prep_glc_mrg_lnd)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer_mrg),barrier=mpicom_CPLID)
    do egi = 1,num_inst_glc
       ! Use fortran mod to address ensembles in merge
       eli = mod((egi-1),num_inst_lnd) + 1
       x2g_gx => component_get_x2c_cx(glc(egi)) 
       call prep_glc_merge(l2x_gx(eli), x2g_gx)
    enddo
    call t_drvstopf  (trim(timer_mrg))

  end subroutine prep_glc_mrg_lnd

  !================================================================================================
  subroutine prep_glc_merge( s2x_g, x2g_g )

    !----------------------------------------------------------------------- 
    ! Arguments
    type(mct_aVect), intent(inout)  :: s2x_g  ! input
    type(mct_aVect), intent(inout)  :: x2g_g  ! output
    !----------------------------------------------------------------------- 

    integer       :: nflds,i,i1,o1
    logical       :: iamroot
    logical, save :: first_time = .true.
    character(CL),allocatable :: mrgstr(:)   ! temporary string
    character(CL) :: field   ! string converted to char
    type(mct_aVect_sharedindices),save :: s2x_sharedindices
    character(*), parameter   :: subname = '(prep_glc_merge) '

    !----------------------------------------------------------------------- 

    call seq_comm_getdata(CPLID, iamroot=iamroot)

    if (first_time) then
       nflds = mct_aVect_nRattr(x2g_g)

       allocate(mrgstr(nflds))
       do i = 1,nflds
          field = mct_aVect_getRList2c(i, x2g_g)
          mrgstr(i) = subname//'x2g%'//trim(field)//' ='
       enddo

       call mct_aVect_setSharedIndices(s2x_g, x2g_g, s2x_SharedIndices)

       !--- document copy operations ---
       do i=1,s2x_SharedIndices%shared_real%num_indices
          i1=s2x_SharedIndices%shared_real%aVindices1(i)
          o1=s2x_SharedIndices%shared_real%aVindices2(i)
          field = mct_aVect_getRList2c(i1, s2x_g)
          mrgstr(o1) = trim(mrgstr(o1))//' = s2x%'//trim(field)
       enddo
    endif

    ! Create input glc state directly from land snow output state
    call mct_aVect_copy(aVin=s2x_g, aVout=x2g_g, vector=mct_usevector, sharedIndices=s2x_SharedIndices)

    if (first_time) then
       if (iamroot) then
          write(logunit,'(A)') subname//' Summary:'
          do i = 1,nflds
             write(logunit,'(A)') trim(mrgstr(i))
          enddo
       endif
       deallocate(mrgstr)
    endif

    first_time = .false.

  end subroutine prep_glc_merge_lnd_forcing


  subroutine prep_glc_calc_o2x_gx(timer)
    !---------------------------------------------------------------
    ! Description
    ! Create o2x_gx

    ! Arguments
    character(len=*), intent(in) :: timer

    character(*), parameter :: subname = '(prep_glc_calc_o2x_gx)'
    ! Local Variables
    integer eoi
    type(mct_avect), pointer :: o2x_ox

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eoi = 1,num_inst_ocn
      o2x_ox => component_get_c2x_cx(ocn(eoi))
      call seq_map_map(mapper_So2g, o2x_ox, o2x_gx(eoi), &
                       fldlist=seq_flds_x2g_states_from_ocn,norm=.true.)
    enddo

    call t_drvstopf  (trim(timer))
  end subroutine prep_glc_calc_o2x_gx

  !================================================================================================


  !================================================================================================

  subroutine prep_glc_calc_l2x_gx(timer)
    !---------------------------------------------------------------
    ! Description
    ! Create l2x_gx (note that l2x_gx is a local module variable)
    ! Also l2x_gx is really the accumulated l2xacc_lx mapped to l2x_gx
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eli
    character(*), parameter :: subname = '(prep_glc_calc_l2x_gx)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eli = 1,num_inst_lnd
       call seq_map_map(mapper_SFl2g, l2gacc_lx(eli), l2x_gx(eli), norm=.true.)
    enddo
    call t_drvstopf  (trim(timer))
  end subroutine prep_glc_calc_l2x_gx

  !================================================================================================

  function prep_glc_get_l2x_gx()
    type(mct_aVect), pointer :: prep_glc_get_l2x_gx(:)
    prep_glc_get_l2x_gx => l2x_gx(:)   
  end function prep_glc_get_l2x_gx

  function prep_glc_get_l2gacc_lx()
    type(mct_aVect), pointer :: prep_glc_get_l2gacc_lx(:)
    prep_glc_get_l2gacc_lx => l2gacc_lx(:)   
  end function prep_glc_get_l2gacc_lx

  function prep_glc_get_l2gacc_lx_cnt()
    integer, pointer :: prep_glc_get_l2gacc_lx_cnt
    prep_glc_get_l2gacc_lx_cnt => l2gacc_lx_cnt
  end function prep_glc_get_l2gacc_lx_cnt

  function prep_glc_get_mapper_SFl2g()
    type(seq_map), pointer :: prep_glc_get_mapper_SFl2g
    prep_glc_get_mapper_SFl2g => mapper_SFl2g  
  end function prep_glc_get_mapper_SFl2g

  function prep_glc_get_mapper_So2g()
    type(seq_map), pointer :: prep_glc_get_mapper_So2g
    prep_glc_get_mapper_So2g=> mapper_So2g
  end function prep_glc_get_mapper_So2g

  function prep_glc_get_mapper_Fo2g()
    type(seq_map), pointer :: prep_glc_get_mapper_Fo2g
    prep_glc_get_mapper_Fo2g=> mapper_Fo2g
  end function prep_glc_get_mapper_Fo2g

!***********************************************************************
!
!  routine compute_melt_fluxes
!
!> \brief   Computes ocean and ice melt fluxes, etc.
!> \author  Xylar Asay-Davis
!> \date    3/27/2015
!>  This routine computes melt fluxes (melt rate, temperature fluxes
!>  into the ice and the ocean, and salt flux) as well as the interface
!>  temperature and salinity.  This routine expects an ice temperature
!>  in the bottom layer of ice and ocean temperature and salinity in
!>  the top ocean layer as well as the pressure at the ice/ocean interface.
!>
!>  The ocean heat and salt transfer velocities are determined based on
!>  observations of turbulent mixing rates in the under-ice boundary layer.
!>  They should be the product of the friction velocity and a (possibly
!>  spatially variable) non-dimenional transfer coefficient.
!>
!>  The iceTemperatureDistance is the distance between the location
!>  where the iceTemperature is supplied and the ice-ocean interface,
!>  used to compute a temperature gradient.  The ice thermal conductivity,
!>  SHR_CONST_KAPPA_LAND_ICE, is zero for the freezing solution from Holland and Jenkins
!>  (1999) in which the ice is purely insulating.
!
!-----------------------------------------------------------------------

  subroutine compute_melt_fluxes( &
       oceanTemperature, &
       oceanSalinity, &
       oceanHeatTransferVelocity, &
       oceanSaltTransferVelocity, &
       interfacePressure, &
       iceTemperature, &
       iceTemperatureDistance, &
       iceFloatingMask, &
       outInterfaceSalinity, &
       outInterfaceTemperature, &
       outFreshwaterFlux, &
       outOceanHeatFlux, &
       outIceHeatFlux, &
       gsize)

    use shr_const_mod,     only: SHR_CONST_CPICE,  &
                                 SHR_CONST_CPSW,   &
                                 SHR_CONST_LATICE, &
                                 SHR_CONST_RHOICE, &
                                 SHR_CONST_RHOSW,  &
                                 SHR_CONST_DTF_DP, &
                                 SHR_CONST_DTF_DS, &
                                 SHR_CONST_DTF_DPDS, &
                                 SHR_CONST_TF0,    &
                                 SHR_CONST_KAPPA_LAND_ICE

    !-----------------------------------------------------------------
    !
    ! input variables
    !
    !-----------------------------------------------------------------

    real (kind=r8), dimension(:), intent(in) :: &
         oceanTemperature, &          !< Input: ocean temperature in top layer
         oceanSalinity, &             !< Input: ocean salinity in top layer
         oceanHeatTransferVelocity, & !< Input: ocean heat transfer velocity
         oceanSaltTransferVelocity, & !< Input: ocean salt transfer velocity
         interfacePressure, &         !< Input: pressure at the ice-ocean interface
         iceTemperature, &            !< Input: ice temperature in bottom layer
         iceTemperatureDistance       !< Input: distance to ice temperature from ice-ocean interface
    integer, dimension(:), intent(in) :: &
         iceFloatingMask              !< Input: mask of cells that contain floating ice

    integer, intent(in) :: gsize !< Input: number of values in each array

    !-----------------------------------------------------------------
    !
    ! output variables
    !
    !-----------------------------------------------------------------

    real (kind=r8), dimension(:), intent(out) :: &
         outInterfaceSalinity, &    !< Output: ocean salinity at the interface
         outInterfaceTemperature, & !< Output: ice/ocean temperature at the interface
         outFreshwaterFlux, &   !< Output: ocean thickness flux (melt rate)
         outOceanHeatFlux, & !< Output: the temperature flux into the ocean
         outIceHeatFlux      !< Output: the temperature flux into the ice

    !-----------------------------------------------------------------
    !
    ! local variables
    !
    !-----------------------------------------------------------------

    real (kind=r8) :: T0, transferVelocityRatio, Tlatent, nu, a, b, c, eta, &
                         iceHeatFluxCoeff, iceDeltaT, dTf_dS
    integer :: n
    character(*), parameter :: subname = '(compute_melt_fluxes)'

    real (kind=r8), parameter :: minInterfaceSalinity = 0.001_r8

    real (kind=r8), parameter :: referencePressure = 0.0_r8 ! Using reference pressure of 0

    real (kind=r8) :: pressureOffset

    Tlatent = SHR_CONST_LATICE/SHR_CONST_CPSW
    do n = 1, gsize
       if (iceFloatingMask(n) == 0) cycle ! Only calculate on floating cells

       if (oceanHeatTransferVelocity(n) == 0.0_r8) then
          write(logunit,*) 'compute_melt_fluxes ERROR: oceanHeatTransferVelocity value of 0 causes divide by 0 at index ', n
          call shr_sys_abort('compute_melt_fluxes ERROR: oceanHeatTransferVelocity value of 0 causes divide by 0')
       end if

       iceHeatFluxCoeff = SHR_CONST_RHOICE*SHR_CONST_CPICE*SHR_CONST_KAPPA_LAND_ICE/iceTemperatureDistance(n)
       nu = iceHeatFluxCoeff/(SHR_CONST_RHOSW*SHR_CONST_CPSW*oceanHeatTransferVelocity(n))
       pressureOffset = max(interfacePressure(n) - referencePressure, 0.0_r8)
       T0 = SHR_CONST_TF0 + SHR_CONST_DTF_DP * pressureOffset
            !Note: These two terms for T0 are not needed because we are evaluating at salinity=0:
            !+ SHR_CONST_DTF_DS * oceanSalinity(n) + SHR_CONST_DTF_DPDS * pressureOffset * oceanSalinity(n)
       iceDeltaT = T0 - iceTemperature(n)
       dTf_dS = SHR_CONST_DTF_DS + SHR_CONST_DTF_DPDS * pressureOffset

       transferVelocityRatio = oceanSaltTransferVelocity(n)/oceanHeatTransferVelocity(n)

       a = -1.0_r8 * dTf_dS * (1.0_r8 + nu)
       b = transferVelocityRatio*Tlatent - nu*iceDeltaT + oceanTemperature(n) - T0
       c = -transferVelocityRatio*Tlatent*max(oceanSalinity(n), 0.0_r8)
       ! a is non-negative; c is strictly non-positive so we never get imaginary roots.
       ! Since a can be zero, we need a solution of the quadratic equation for 1/Si instead of Si.
       ! Following: https://people.csail.mit.edu/bkph/articles/Quadratics.pdf
       ! Since a and -c are are non-negative, the term in the square root is also always >= |b|.
       ! In all reasonable cases, b will be strictly positive, since transferVelocityRatio*Tlatent ~ 2 C,
       ! T0 ~ -1.8 C and oceanTemperature should never be able to get below about -3 C
       ! As long as either b or both a and c are greater than zero, the strictly non-negative root is
       outInterfaceSalinity(n) = max(-(2.0_r8*c)/(b + sqrt(b**2 - 4.0_r8*a*c)), minInterfaceSalinity)

       outInterfaceTemperature(n) = dTf_dS*outInterfaceSalinity(n)+T0

       outFreshwaterFlux(n) = SHR_CONST_RHOSW*oceanSaltTransferVelocity(n) &
            * (oceanSalinity(n)/outInterfaceSalinity(n) - 1.0_r8)

       ! According to Jenkins et al. (2001), the temperature fluxes into the ocean are:
       !   1. the advection of meltwater into the top layer (or removal for freezing)
       !   2. the turbulent transfer of heat across the boundary layer, based on the termal driving
       outOceanHeatFlux(n) = SHR_CONST_CPSW*(outFreshwaterFlux(n)*outInterfaceTemperature(n) &
            - SHR_CONST_RHOSW*oceanHeatTransferVelocity(n)*(oceanTemperature(n)-outInterfaceTemperature(n)))

       ! the temperature fluxes into the ice are:
       !   1. the advection of ice at the interface temperature out of the domain due to melting
       !      (or in due to freezing)
       !   2. the diffusion (if any) of heat into the ice, based on temperature difference between
       !      the reference point in the ice (either the surface or the middle of the bottom layer)
       !      and the interface
       outIceHeatFlux(n) = -SHR_CONST_CPICE*outFreshwaterFlux(n)*outInterfaceTemperature(n)

       outIceHeatFlux(n) = outIceHeatFlux(n) &
            - iceHeatFluxCoeff*(iceTemperature(n) - outInterfaceTemperature(n))

    end do

  !--------------------------------------------------------------------
  end subroutine compute_melt_fluxes

end module prep_glc_mod
