module source_sub

!!$  use precision, only: dp, si, i8b
  use precision, only: dp, si, li
  use my_mpi
  use file_admin, only: logf, results_dir
  use cgsconstants, only: m_p
  use sizes, only: mesh
  use astroconstants, only: M_SOLAR
  use cosmology_parameters, only: Omega_B, Omega0, h
  use nbody, only: id_str, M_grid, dir_src, zred_array
  use material, only: xh, dens_ND, dens_ND_prev
  use grid, only: x,y,z, vol_cMpc3
  use c2ray_parameters, only: phot_per_atom, fstar, &
       lifetime, S_star_nominal, StillNeutral

  implicit none

  integer, public :: NumAGrid !< Number of Active Grids having subgrid sources
  integer, public, parameter :: MHflag = 2
!!$  integer, public, parameter :: MHflag = 0

  integer,dimension(mesh(1),mesh(2),mesh(3)) :: AGflag

  real(kind=dp) :: densNDcrit, densNDcrit_prev !< critical nondimensional density over which cell can be minihalo populated

  integer,dimension(:,:),allocatable     :: sub_srcpos !< mesh position of subgrid sources
  real(kind=dp),dimension(:,:),allocatable :: rsub_srcpos !< grid position of subgrid sources
  real(kind=dp),dimension(:),allocatable :: ssM_msun !< array of mass of subgrid sources (one source) in solar mass unit. Similar to sM00_msun.
  real(kind=dp),dimension(:),allocatable :: subsrcMass !< masses of subgrid sources. Similar to srcMass.
  real(kind=dp),dimension(:),allocatable :: subNormFlux !< normalized ionizing flux of subgrid sources
  integer,dimension(:),allocatable       :: subSrcSeries !< a randomized list of sources
  
  ! Definitions for reading precalculated # density of minihalos
  integer(kind=li),  public :: Nzdata, Ndelta1data  !< LGnMH data size integer header
  real   (kind=dp),  public :: zmin, zmax, LGdelta1min, LGdelta1max  !< LGnMH data domain real*8 header
  real   (kind=dp), public :: dzred, dLGdelta1  !< Needed to search LGnMH data
  real   (kind=dp), dimension(:,:), allocatable, public :: LGnMH_Mpc3 !< log10 of number of MHs per comoving (Mpc)^3, at given redshfit and delta+1, where delta is overdensity.
  real   (kind=dp), parameter, public :: M_PIIIstar_msun = 300d0 !< mass of PopIII star in solar mass unit.
  real   (kind=dp), public :: tot_subsrcM_msun !< total STELLAR mass in all active (survived after suppression) minihalos
  character(len=100), parameter :: MHfit_file="zred_halodelta1_nMHMpc3_Planck"


contains

  ! =======================================================================

subroutine AGrid_properties(nz,AGlifetime,jLWgrid)

   ! Input routine: establish the subgrid source properties
   ! Author: Kyungjin Ahn
   ! Update: 5-Sep-2009, 15-Jul-2010

   use perm  ! for random permutation of sources

   integer,      intent(in) :: nz
   real(kind=dp),intent(in) :: AGlifetime ! time step
   real(kind=dp),dimension(mesh(1),mesh(2),mesh(3)),intent(in) :: jLWgrid

   real(kind=dp)            :: zred_prev,zred_now
   real(kind=dp)            :: jLWcrit_now, jLWcrit_min
   real(kind=dp)            :: diff_subsrcMsun

   character(len=512) :: sourcelistfile_sub
   character(len=6)   :: z_str

   integer :: ns
   integer :: i, j, k

#ifdef MPI
   integer :: mympierror
#endif
  ! CRITICAL: Add debug BEFORE any IF statement or array access
   write(6,*) 'DEBUG: AGrid_properties called, rank=', rank, ' nz=', nz
   write(6,*) 'DEBUG: AGlifetime=', AGlifetime
   write(6,*) 'DEBUG: mesh=', mesh(1), mesh(2), mesh(3)
   
   ! Check if we can access jLWgrid at all
   write(6,*) 'DEBUG: About to access jLWgrid(1,1,1)...'
   write(6,*) 'DEBUG: jLWgrid(1,1,1)=', jLWgrid(1,1,1)
   write(6,*) 'DEBUG: jLWgrid access successful'

   if (rank == 0) then
      write(6,*) 'DEBUG: Entering AGrid_properties, nz=', nz
      write(6,*) 'DEBUG: Checking array allocations...'
      write(6,*) 'DEBUG: allocated(dens_ND)=', allocated(dens_ND)
      write(6,*) 'DEBUG: allocated(dens_ND_prev)=', allocated(dens_ND_prev)
      write(6,*) 'DEBUG: allocated(LGnMH_Mpc3)=', allocated(LGnMH_Mpc3)
      write(6,*) 'DEBUG: allocated(xh)=', allocated(xh)
      if (allocated(dens_ND)) write(6,*) 'DEBUG: dens_ND(1,1,1)=', dens_ND(1,1,1)
   endif

   if (rank == 0) write(6,*) 'DEBUG: Entering AGrid_properties, nz=', nz

   ! Ask for input
   if (rank == 0) then

      if (nz > 1) zred_prev = zred_array(nz-1)
      zred_now  = zred_array(nz)
      write(z_str,'(f6.3)') zred_now
      sourcelistfile_sub=trim(adjustl(results_dir))//&
            trim(adjustl(z_str))//"-"//&
            trim(adjustl(id_str))//"_SUBsources.dat"

      if (nz > 1) call get_denscrit(zred_prev, densNDcrit_prev)
      call get_denscrit(zred_now,  densNDcrit)
      write(logf,*) 'DEBUG: zred, densNDcrit', zred_now, densNDcrit
      write(6,*) 'DEBUG: zred, densNDcrit', zred_now, densNDcrit

      jLWcrit_now = jLWcrit(zred_now)
      jLWcrit_min = 0.1d0 * jLWcrit_now 

      ! Initialize flags
      NumAGrid = 0
      AGflag   = 0

      if (rank == 0) write(6,*) 'DEBUG: Starting AGflag loop'
      do k=1,mesh(3)
         do j=1,mesh(2)
            do i=1,mesh(1)
               if (nz > 1) then
                  diff_subsrcMsun = &
                        subsrcM_msun(zred_now,dens_ND(i,j,k),densNDcrit) &
                        - subsrcM_msun(zred_prev,dens_ND_prev(i,j,k),densNDcrit_prev)
               else
                  diff_subsrcMsun = subsrcM_msun(zred_now,dens_ND(i,j,k),densNDcrit)
               endif

               ! DEBUG: print first cell values
               if (i==1 .and. j==1 .and. k==1) then
                  write(6,*) 'DEBUG: First cell check, diff_subsrcMsun=', diff_subsrcMsun, &
                              ' jLWgrid=', jLWgrid(1,1,1), ' xh=', xh(1,1,1,1)
               endif

               if (xh(i,j,k,1) < StillNeutral .and. jLWgrid(i,j,k) < jLWcrit_now &
                     .and. diff_subsrcMsun > 0d0) then
                  AGflag(i,j,k) = 1
                  NumAGrid      = NumAGrid + 1
               endif
            enddo
         enddo
      enddo
      if (rank == 0) write(6,*) 'DEBUG: Finished AGflag loop, NumAGrid=', NumAGrid

      ! Allocate subgrid arrays
      if (allocated(ssM_msun)) deallocate(ssM_msun)
      if (allocated(subsrcMass)) deallocate(subsrcMass)
      if (allocated(subNormFlux)) deallocate(subNormFlux)
      if (allocated(subSrcSeries)) deallocate(subSrcSeries)
      if (allocated(sub_srcpos)) deallocate(sub_srcpos)
      if (allocated(rsub_srcpos)) deallocate(rsub_srcpos)

      if (rank == 0) write(6,*) 'DEBUG: Allocating arrays, NumAGrid=', NumAGrid
      allocate(ssM_msun(NumAGrid))
      allocate(subsrcMass(NumAGrid))
      allocate(subNormFlux(NumAGrid))
      allocate(subSrcSeries(NumAGrid))
      allocate(sub_srcpos(3,NumAGrid))
      allocate(rsub_srcpos(3,NumAGrid))
      if (rank == 0) write(6,*) 'DEBUG: Allocation successful'

      ssM_msun         = 0d0
      NumAGrid         = 0
      tot_subsrcM_msun = 0d0

      ! Second loop: fill subgrid masses
      if (rank == 0) write(6,*) 'DEBUG: Starting subgrid mass loop'
      do k=1,mesh(3)
         do j=1,mesh(2)
            do i=1,mesh(1)
               if (nz > 1) then
                  diff_subsrcMsun = &
                        subsrcM_msun(zred_now,dens_ND(i,j,k),densNDcrit) &
                        - subsrcM_msun(zred_prev,dens_ND_prev(i,j,k),densNDcrit_prev)
               else
                  diff_subsrcMsun = subsrcM_msun(zred_now,dens_ND(i,j,k),densNDcrit)
               endif

               if (AGflag(i,j,k) == 1) then
                  NumAGrid = NumAGrid + 1
                  sub_srcpos(1, NumAGrid) = i
                  sub_srcpos(2, NumAGrid) = j
                  sub_srcpos(3, NumAGrid) = k
                  rsub_srcpos(1, NumAGrid) = x(i)
                  rsub_srcpos(2, NumAGrid) = y(j)
                  rsub_srcpos(3, NumAGrid) = z(k)

                  if (jLWgrid(i,j,k) <= jLWcrit_min) then
                     ssM_msun(NumAGrid) = diff_subsrcMsun
                  else
                     ssM_msun(NumAGrid) = diff_subsrcMsun * &
                                          (jLWcrit_now - jLWgrid(i,j,k)) / &
                                          (jLWcrit_now - jLWcrit_min)
                  endif
               endif
            enddo
         enddo
      enddo
      if (rank == 0) write(6,*) 'DEBUG: Finished subgrid mass loop, tot NumAGrid=', NumAGrid

      tot_subsrcM_msun = sum(ssM_msun)

      ! Convert to subsrcMass
      if (MHflag == 1) then
         subsrcMass = ssM_msun * (M_SOLAR/M_grid) * phot_per_atom(3)
      elseif (MHflag == 2) then
         subsrcMass = ssM_msun * (M_SOLAR/M_grid) * phot_per_atom(3)/fstar(3)
      endif

      ! Save subgrid source list
      open(unit=49,file=sourcelistfile_sub,status='unknown')
      write(49,*) NumAGrid
      write(49,*) MHflag
      do ns=1,NumAGrid
         write(49,333) sub_srcpos(1,ns),sub_srcpos(2,ns),sub_srcpos(3,ns), &
                        subsrcMass(ns), ssM_msun(ns)
333       format(3I5,2e13.4)
      enddo
      close(49)

      write(logf,*) "Number of Active Grids: ", NumAGrid

      ! Compute subNormFlux
      if (MHflag == 1) then
         subNormFlux = subsrcMass * M_grid * Omega_B/(Omega0*m_p) / S_star_nominal / AGlifetime
      elseif (MHflag == 2) then
         subNormFlux = subsrcMass * M_grid / m_p / S_star_nominal / AGlifetime
      endif

      write(6,*) 'DEBUG: Finished subNormFlux computation'

   endif ! end rank==0

#ifdef MPI
   ! Distribute to other nodes
   call MPI_BCAST(NumAGrid,1,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
   if (rank /= 0) then
      if (allocated(ssM_msun)) deallocate(ssM_msun)
      if (allocated(subsrcMass)) deallocate(subsrcMass)
      if (allocated(subNormFlux)) deallocate(subNormFlux)
      if (allocated(subSrcSeries)) deallocate(subSrcSeries)
      if (allocated(sub_srcpos)) deallocate(sub_srcpos)
      if (allocated(rsub_srcpos)) deallocate(rsub_srcpos)

      allocate(ssM_msun(NumAGrid))
      allocate(subsrcMass(NumAGrid))
      allocate(subNormFlux(NumAGrid))
      allocate(subSrcSeries(NumAGrid))
      allocate(sub_srcpos(3,NumAGrid))
      allocate(rsub_srcpos(3,NumAGrid))
   endif
   call MPI_BCAST(ssM_msun, NumAGrid, MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
   call MPI_BCAST(subsrcMass, NumAGrid, MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
   call MPI_BCAST(subNormFlux, NumAGrid, MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
   call MPI_BCAST(sub_srcpos, 3*NumAGrid, MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
   call MPI_BCAST(rsub_srcpos, 3*NumAGrid, MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
#endif

   if (rank == 0) then
      write(logf,*) 'Subgrid Source lifetime=', AGlifetime/3.1536e13
      write(logf,*) 'Subgrid Total flux= ',sum(subNormFlux)

      ! Prepare random series
      do ns=1,NumAGrid
         subSrcSeries(ns)=ns
      enddo
      call permi(NumAGrid, subSrcSeries)
      write(6,*) 'DEBUG: Finished random permutation'
   endif

#ifdef MPI
   call MPI_BCAST(subSrcSeries,NumAGrid,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
#endif

   if (rank == 0) write(6,*) 'DEBUG: Exiting AGrid_properties'

   return
end subroutine AGrid_properties


  ! =======================================================================

  subroutine get_denscrit(zred, densNDcrit)

    ! Find out the critical dimensionless density, over which cells
    ! can populate minihalos. This is set by matching
    ! the global number density of all minihalos to high-res, small-box
    ! simulation result which resolves minihalos.

    real(kind=dp), intent(in) :: zred

    ! For total number of minihalos in the small box
    integer       :: Ntable, itable
    real(kind=dp) :: NMH_smallbox, NMH_simbox
    real(kind=dp) :: size_smallbox     ! In h^-1 cMpc 
    real(kind=dp) :: volratio    ! (small box volume / simulation box volume)
    real(kind=dp), allocatable :: ztable(:), numMHtable(:)
    real(kind=dp) :: dz_table

    ! For iteration
    integer, parameter :: Nmaxiter = 14
    integer            :: iiter, i, j, k
    real(kind=dp)      :: LGdensND_L, LGdensND_R, LGdensND_N
    real(kind=dp)      :: frac

    ! For fitting
    integer            :: idx_zred, idx_LGdelta1
    real(kind=dp)      :: LGdelta1

    ! what we want
    real(kind=dp), intent(out) :: densNDcrit


    if (MHflag == 1) then
       write(6,*)  'Do something!!!!!!!!!!!!!!!!!!!!!!!'
       stop
    elseif (MHflag == 2) then
       ! Read and find the total number of minihalos in small box
       open(unit=55, file="z_numMH_6.3Mpc_full", status="old")
       read(55,*) size_smallbox  ! in unit of h^-1 cMpc
       read(55,*) Ntable

       allocate(ztable(Ntable))
       allocate(numMHtable(Ntable))

       do itable = 1, Ntable
          read(55,*) ztable(itable), numMHtable(itable)
       enddo
       close(55)

       volratio = (size_smallbox/h)**3 &
            / (vol_cMpc3*real(mesh(1)*mesh(2)*mesh(3),dp))

       ! Interpolate by finding the nearest neighbor index. Data table given 
       ! should be in the ascending order in z, and has to be in very fine, 
       ! equal size bin. This is the sheer total number of minihalos in the
       ! small box, NOT the number density of minihalos.
       dz_table     = (ztable(Ntable)-ztable(1))/real(Ntable-1, dp)
       NMH_smallbox = numMHtable(int((zred-ztable(1))/dz_table) + 1)
       
       deallocate(ztable)
       deallocate(numMHtable)

       ! Now iterate until we find the right number for the bigbox
       idx_zred     = int((zred - zmin) /dzred) + 1
       ! Very poor extrapolation. Just limit your simulation
       ! within given redshift range.
       if (idx_zred     <= 0         ) idx_zred = 1
       if (idx_zred     > Nzdata     ) idx_zred = Nzdata


       ! bisecting(?) interpolation
       LGdensND_L = -1d0 ! starting L value, should be small enough.
       LGdensND_R = 1d0  ! starting R value, should be large enough.
       iiter      = 0
       do
          iiter = iiter + 1
          LGdensND_N = (LGdensND_R + LGdensND_L)*0.5d0

          NMH_simbox = 0d0
          do k=1,mesh(3)
             do j=1,mesh(2)
                do i=1,mesh(1)
                   if (log10(dens_ND(i,j,k)) >= LGdensND_N) then
                      idx_LGdelta1 = int((log10(dens_ND(i,j,k)) - LGdelta1min)&
                           /dLGdelta1) + 1
                      ! Extrapolation, very poor, but only if needed. 
                      if (idx_LGdelta1 <= 0         ) idx_LGdelta1=1
                      if (idx_LGdelta1 > Ndelta1data) idx_LGdelta1=Ndelta1data

                      NMH_simbox = NMH_simbox &
                           + 10d0**LGnMH_Mpc3(idx_zred,idx_LGdelta1)
                   endif
                enddo
             enddo
          enddo
          ! Normalize to the simulation box to get the total number of
          ! minihalos in the simulation box.
          !!!!!!!!! Careful
          ! Because fitting function uses h=0.7, we need to scale further with used h.
          ! SHOULD BE FIXED LATER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!
          NMH_simbox = NMH_simbox * vol_cMpc3 * (h/0.7d0)**3
          
          ! fractional difference
          frac = (NMH_simbox*volratio - NMH_smallbox)/NMH_smallbox

          ! convergence is 1%, chosen just arbitrarily.
          if (abs(frac) <= 0.01d0 .or. iiter == Nmaxiter) exit
          if (NMH_simbox*volratio < NMH_smallbox) then
             LGdensND_R = LGdensND_N
          else
             LGdensND_L = LGdensND_N
          endif
       enddo
       densNDcrit = 10d0**LGdensND_N
    endif   ! end of MHflag=2 case

  end subroutine get_denscrit

  ! =======================================================================

  function subsrcM_msun(zred, densND, densNDcrit)
    
    ! Get subgrid src mass (those not resolved in N-body) in solar mass unit
    ! Cell characteristic given by zred & dimensionless density

    integer                   :: idx_zred, idx_LGdelta1

    real(kind=dp)             :: subsrcM_msun
    real(kind=dp), intent(in) :: zred
    real(kind=si), intent(in) :: densND
    real(kind=dp), intent(in) :: densNDcrit
    real(kind=dp)             :: LGdelta1
    real(kind=dp)             :: N_MH


    if (MHflag == 1) then
       write(6,*)  'Do something!!!!!!!!!!!!!!!!!!!!!!!'
       subsrcM_msun = 0d0
    elseif (MHflag == 2) then
       LGdelta1 = log10(dble(densND))
       if (densND <= 0d0) LGdelta1 = -5  ! for numerical safety.

       idx_zred     = int((zred     - zmin       ) /dzred    ) + 1
       ! Very poor extrapolation. Just limit your simulation
       ! within given redshift range.
       if (idx_zred     <= 0         ) idx_zred = 1
       if (idx_zred     > Nzdata     ) idx_zred = Nzdata

       idx_LGdelta1 = int((LGdelta1 - LGdelta1min) /dLGdelta1) + 1
       ! Extrapolation. LGdelta1 covers -1 to 1, which is not too bad.
       if (idx_LGdelta1 <= 0         ) idx_LGdelta1 = 1
       if (idx_LGdelta1 > Ndelta1data) idx_LGdelta1 = Ndelta1data

       ! densNDcrit is CHOSEN to match the total number of minihalos in
       ! small-box, high-res simulation. So N_MH can be sometimes less than 1,
       ! and more generally fractional numbers.
       if (densND >= densNDcrit) then
          ! number of minihalos contained in the cell, with volume vol_cMpc3 
          ! in comoving Mpc^3 unit.
          ! number of minihalos on the cell with given volume
          ! !!!!!!!! Careful
          ! Because fitting function uses h=0.7, we need to scale further with used h.
          ! SHOULD BE FIXED LATER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! !!!!!!!!
          N_MH  = 10d0**LGnMH_Mpc3(idx_zred,idx_LGdelta1) * vol_cMpc3 &
               * (h/0.7d0)**3

          ! Total mass of stellar baryon in the cell.
          subsrcM_msun = N_MH * M_PIIIstar_msun
       else
          subsrcM_msun = 0d0
       endif
    endif

  end function subsrcM_msun

  ! =======================================================================

  subroutine read_LGnMH_Mpc3

    use nbody, only: binaryaccess
#ifdef MPI
    integer :: mympierror
#endif

    if (rank == 0) then
!!$       open(unit=20, file=trim(adjustl(MHfit_file)), form="binary", status="old")
       open(unit=20, file=trim(adjustl(MHfit_file)), form="unformatted", access=binaryaccess, status="old")
       read(20) Nzdata, Ndelta1data  
       read(20) zmin, zmax, LGdelta1min, LGdelta1max 
       if (rank==0) write(6,*)  '.................'
       if (rank==0) write(6,*)  nzdata, ndelta1data, zmin, zmax, LGdelta1min, LGdelta1max
       
       if (allocated(LGnMH_Mpc3)) deallocate(LGnMH_Mpc3)
       allocate(LGnMH_Mpc3(Nzdata, Ndelta1data))

       read(20) LGnMH_Mpc3
       close(20)

       dzred     = (zmax       -zmin       ) /real(Nzdata     -1, 8)
       dLGdelta1 = (LGdelta1max-LGdelta1min) /real(Ndelta1data-1, 8)
    endif

#ifdef MPI
    call MPI_BCAST(Nzdata,     1,                  MPI_INTEGER,         0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(Ndelta1data,1,                  MPI_INTEGER,         0,MPI_COMM_NEW,mympierror)
    if (rank /= 0) then
       if (allocated(LGnMH_Mpc3)) deallocate(LGnMH_Mpc3)
       allocate(LGnMH_Mpc3(Nzdata, Ndelta1data))
    endif
    call MPI_BCAST(zmin,       1,                  MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(zmax,       1,                  MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(LGdelta1min,1,                  MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(LGdelta1max,1,                  MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(LGnMH_Mpc3, Nzdata*Ndelta1data, MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
#endif

  end subroutine read_LGnMH_Mpc3

  ! =======================================================================

  function jLWcrit(zred)

    ! Get critical jLW at zred

    real(kind=dp)             :: jLWcrit
    real(kind=dp), intent(in) :: zred

    jLWcrit = 0.05d0 * 1d-21

  end function jLWcrit


end module source_sub
