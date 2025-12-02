!>
!! \brief Main program for C2Ray-3Dm (3D, multiple sources)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 22-May-2008
!<
Program C2Ray

  ! Authors: Garrelt Mellema, Ilian Iliev

  ! Date: 22-May-2008 (06-Mar-2008 (30-Jan-2008 (8-Dec-2007 (23-Sep-2006))

  ! Goal:
  ! Cosmological reionization simulation using precomputed density fields
  ! and source lists.
  
  ! Does not include hydrodynamics
  ! Assumes constant time step

  ! Needs following `modules'
  ! c2ray_parameters : all the tunable parameters for the code
  ! my_mpi : sets up the MPI parallel environment
  ! output_module : output routines
  ! grid : sets up the grid
  ! radiation : radiation tools
  ! nbody : interface to N-body output
  ! cosmology : cosmological utilities
  ! material : material properties
  ! times : time and time step utilities
  ! sourceprops : source properties
  ! evolve : evolve grid in time

  use precision, only: dp, si
  use file_admin, only: stdinput, logf, file_input, flag_for_file_input
  use c2ray_parameters, only: cosmological, type_of_clumping
  use astroconstants, only: YEAR
  use my_mpi !, only: mpi_setup, mpi_end, rank
  use output_module, only: setup_output,output,close_down
  use grid, only: grid_ini
  use sizes, only: mesh
  use radiation, only: rad_ini
  use nbody, only: nbody_type, nbody_ini, NumZred, zred_array, snap, &
       NumZred_coarse, zred_array_coarse
  use cosmology, only: cosmology_init, redshift_evol, cosmo_evol, &
       time2zred, zred2time, zred
  use material, only: mat_ini, xfrac_ini, dens_ini, set_clumping
#ifdef MH
  use material, only: dens_ND, dens_ND_prev, nz0_coarse, xh
#endif
  use times, only: time_ini, set_timesteps
  use sourceprops, only: source_properties_ini, source_properties, NumSrc
#ifdef MH
  use source_sub, only: AGrid_properties, MHflag, read_LGnMH_Mpc3, NumAGrid, tot_subsrcM_msun
  use jLWgreen, only: get_HcOm
  use getjLW, only: get_jLW, jLW, read_jLW
  use mergesrc, only: merge_allsrc
#endif
  use evolve, only: evolve3D

#ifdef XLF
  USE XLFUTILITY, only: iargc, getarg, flush => flush_
#endif

  implicit none

#ifdef PGI
  include 'lib3f.h' ! for iargc, getargc
#endif

  ! Start and end time for CPU report
  real :: cputime1 !< Start time for CPU report
  real :: cputime2 !< End time for CPU report
  real(kind=dp) :: cpu_seconds=0.0
  integer :: cpu_hours=0
  integer :: cpu_minutes=0

#ifdef MH
  real(kind=dp) :: zred_coarse !< coarsely spaced redshifts, only used when minihalo redshifts(zred) are finer than more massive halo redshifts(zred_coarse)
#endif

  ! Wall clock time variables
  integer :: cntr1 !< Start time wall clock
  integer :: cntr2 !< End time wall clock
  integer :: countspersec !< counts per second (for wall clock time)
  real(kind=dp) :: clock_seconds=0.0
  integer :: clock_hours=0
  integer :: clock_minutes=0

  integer :: restart=0 !< restart flag
  integer :: iter_restart=0 !< restart from iteration flag
  integer :: nz !< loop counter for loop over redshift list
#ifdef MH
  integer :: nz_coarse !< loop counter for loop over coarse redshift list
#endif
  integer :: nz0 !< index of starting redshift from redshift list
  integer :: ierror !< error flag
  integer :: photcons_flag=0 !< photon conservation flag, non-zero if photon conservation is violated. This stops the simulation

  ! end_time - end time of the simulation (s)
  ! dt - time step (s)
  ! sim_time - actual time (s)
  ! end_time - end time of calculation (s)
  ! output_time - time interval between outputs (s)
  ! next_output_time - time of next output (s)
  real(kind=dp) :: end_time !< end time of the simulation (s)
  real(kind=dp) :: sim_time !< actual time (s)
  real(kind=dp) :: output_time !< time interval between outputs (s)
  real(kind=dp) :: next_output_time !< time of next output (s)
  real(kind=dp) :: dt !< calculated time step
  real(kind=dp) :: actual_dt !< actual time step (s)
  real(kind=dp) :: zred_interm !< intermediate redshift (for restart)
  real(kind=dp) :: interm_zred !< calculated intermediate redshift (for restart)
  real(kind=dp) :: est_mem_gb

#ifdef MPI
  integer :: mympierror
#endif

  ! Input file
  character(len=512) :: inputfile !< name of input file
  character(len=1) :: answer !< y or n answer


#ifdef MH
  ! set the constant HcOm for LW calculation
  call get_HcOm
#endif

  ! Initialize cpu timer
  call cpu_time(cputime1)

  ! Initialize wall clock times
  call system_clock(cntr1)

  ! Set up MPI structure
  call mpi_setup()

  if (rank == 0) print*,"This is rank 0"
   if (rank == 1) print*,"This is rank 1"
   if (rank == 2) print*,"This is rank 2"

  ! Set up input stream (either standard input or from file given
  ! by first argument)
  if (iargc() > 0) then
     call getarg(1,inputfile)
     if (rank == 0) then
        write(logf,*) "reading input from ",trim(adjustl(inputfile))
        open(unit=stdinput,file=inputfile)
        call flag_for_file_input(.true.)
     endif
  endif

  ! Initialize output
  call setup_output ()
  if (rank == 0) call flush(logf)

  ! Initialize grid
  call grid_ini ()
  ! Initialize photo-ionization calculation
  call rad_ini ( )

  ! Initialize the material properties
  call mat_ini (restart, nz0, ierror)
 print*,"Pre Nbody"
  ! Find the redshifts we are dealing with
  call nbody_ini ()
print*,"Post Nbody"
#ifdef MH
  ! Read in precalculated minihalo data
  if (MHflag == 1) then
     write(6,*) 'do something!!!!!!!!!!'
  elseif (MHflag == 2) then
     call read_LGnMH_Mpc3
  endif
#endif
  ! Initialize the source model
  call source_properties_ini ()
  ! Initialize time step parameters
  call time_ini ()
  if (rank == 0) call flush(logf)

  ! Set time to zero
  sim_time=0.0

  ! Initialize cosmology
  call cosmology_init(zred_array(nz0),sim_time)

  ! If a restart, inquire whether to restart from iteration
  if (restart /= 0) then
     if (rank == 0) then
        if (.not.file_input) &
             write(*,"(A,$)") "Restart from iteration dump (y/n)? : "
        read(stdinput,*) answer
        write(logf,*) "restart answer: ",answer
        if (answer == "y" .or. answer == "Y") then
           ! Set flag, this is passed to evolve3d
           iter_restart=1
           write(logf,*) "Restarting from iteration dump."
        endif
     endif
#ifdef MPI
     call MPI_BCAST(iter_restart,1,MPI_INTEGER,0,MPI_COMM_NEW, &
          mympierror)
#endif
  endif

  ! If a restart, read ionization fractions from file
  if (restart == 1) call xfrac_ini(zred_array(nz0))
  if (restart == 2) then
     if (rank == 0) then
        if (.not.file_input) &
             write(*,"(A,$)") "At which redshift to restart x_frac?: "
        read(stdinput,*) zred_interm
     endif
#ifdef MPI
     call MPI_BCAST(zred_interm,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
          mympierror)
#endif
     call xfrac_ini(zred_interm)
  end if
  
#ifdef MH
  ! Assume zero initial LW intensity.
  jLW        = 0d0
  print*, "mesh(1),mesh(2),mesh(3)", mesh(1),mesh(2),mesh(3)
  print *, "Size in X-direction = ", size(jLW, 1)
  print *, "Size in Y-direction = ", size(jLW, 2)
  print *, "Size in Z-direction = ", size(jLW, 3)
  ! Read in LW intensity if restart.
  if (restart /= 0) call read_jLW(zred_array(nz0))
#endif

#ifdef MH
  if (zred_array_coarse(nz0_coarse) < zred_array(nz0)) then
     if (rank == 0) write(logf,*) "No, the starting zred_coarse should be equal to or larger than zred(fine)"
     if (rank == 0) write(6,*)    "No, the starting zred_coarse should be equal to or larger than zred(fine)"
     stop
  endif
  ! nz_coarse initialization
  nz_coarse = nz0_coarse
#endif
  ! Loop over redshifts
  do nz=nz0,NumZred-1

     zred        = zred_array(nz)
#ifdef MH
     zred_coarse = zred_array_coarse(nz_coarse)
#endif
     if (rank == 0) then
        write(logf,*) "-------------------------------------"
        write(logf,*) "Doing redshift: ",zred," to ", &
             zred_array(nz+1)
        write(logf,*) "-------------------------------------"
     endif
     if (rank == 0) then 
        write(6,*) "-------------------------------------"
        write(6,*) "Doing redshift: ",zred," to ", &
             zred_array(nz+1)
        write(6,*) "-------------------------------------"
     endif
     ! Initialize time parameters
     call set_timesteps(zred,zred_array(nz+1), &
          end_time,dt,output_time)
     if (rank == 0) write(logf,*) &
          "This is time ",sim_time/YEAR," to ",end_time/YEAR
     if (rank == 0 .and. nbody_type == "LG") &
          write(logf,*) "This is snapshot ", snap(nz)
         
     ! Initialize source position
#ifdef MPILOG     
     write(logf,*) 'Calling source_properties'
#endif 
     if (rank ==0) print *, 'calling source'
#ifdef MH
     ! When finer redshifts are used for MHs, use a coarse redshift
     ! to read in more-massive-than-MH sources and keep using it for
     ! every fine redshift (therefore end_time-sim_time corresponds to
     ! fine time spacing, which also requires change of fgamma) until
     ! new coarse redshift overlaps with fine redshifts. This way 
     ! source suppression of intermediate mass halos(10^8 < M/Msun < 10^9) 
     ! occurs only every coarse redshift.
     if (zred_coarse == zred) then
        call source_properties(zred_coarse,nz_coarse,end_time-sim_time,restart)
        if (rank==0) print *, 'source props read. numsrc:', numsrc
        nz_coarse     = nz_coarse + 1
     elseif (zred_coarse < zred) then
        ! Need to reassign NumSrc, SrcSeries, srcpos, NormFlux back to
        ! those at old coarse redshift, because these quantities are renewed by 
        ! merge_allsrc every fine redshift.
        ! Need ionized fraction at right redshift
        if (nz_coarse-1 == 1) then
           xh(:,:,:,0) = 1_dp
           xh(:,:,:,1) = 0_dp
        else
           call xfrac_ini(zred_array_coarse(nz_coarse-1))
        endif
        call source_properties(zred_array_coarse(nz_coarse-1),nz_coarse-1,end_time-sim_time,restart)
        ! Revert to the current-redshift ionized fraction for MH syppression.
        call xfrac_ini(zred) 
     ! When restart and when redshifts mismatch in the first place:
     elseif (restart > 0 .and. zred_coarse > zred_array(nz0)) then
        ! Need ionized fraction at right redshift
        call xfrac_ini(zred_coarse)  
        call source_properties(zred_coarse,nz0_coarse,end_time-sim_time,restart)
        ! Revert to the current-redshift ionized fraction for MH syppression.
        call xfrac_ini(zred) 
        nz_coarse = nz_coarse + 1  
        !Now becomes normal run
     endif
#else
     call source_properties(zred,nz,end_time-sim_time,restart)
#endif

     if (rank ==0) print *, 'calling dens_ini', ' zred, nz, nz0', zred, nz, nz0
     ! Initialize density field
     call dens_ini(zred,nz)
     if (rank ==0) print *, 'finishing dens_ini'
     if (type_of_clumping == 5) call set_clumping(zred)

#ifdef MH
     NumAGrid = 0
     ! Based upon LW intensity, xfrac and density, get subgrid sources.
     if (MHflag == 1 .OR. MHflag == 2) then
        if (rank ==0) write(6,*) 'before agrid_prop, jlw(1,1,1)', jlw(1,1,1)
        ! Before calling AGrid_properties, assign dens_ND_prev
        if (restart == 0) then
           if (nz > nz0) then !dens_ND assigned before
              dens_ND_prev = dens_ND
           elseif (nz > 1) then
              if (rank == 0) write(logf,*) "You are neglecting previous history. If not a restart, please start from the first redshift."
              stop
           else !The very first start.
              dens_ND_prev = 0_si
           endif
        else
           if (nz > nz0) then
              dens_ND_prev = dens_ND
           elseif (nz > 1) then !When restart just commenced
              call dens_ini(zred_array(nz-1),nz-1)
              dens_ND_prev = dens_ND
              call dens_ini(zred_array(nz)  ,nz  )
           else
              if (rank == 0) write(logf,*) "Are you sure? You said restart, but how can it start from the first redshift? Stopping!"
              stop
           endif
        endif
        call AGrid_properties(nz,end_time-sim_time,jLW)
        
     endif
     if (rank == 0) then
        write(logf,*) 'Number of active subgrids:', NumAGrid
        write(6,*)  'Number of active subgrids:', NumAGrid
        write(logf,*) 'Total active stellar mass in subgrids, in solar mass:', tot_subsrcM_msun
        write(6,*) 'Total active stellar mass in subgrids, in solar mass:', tot_subsrcM_msun
     endif

     ! Get LW intensity jLW at nz+1; Note carefully you should call
     ! get_jLW with "nz" for the jLW at "nz+1" !!! 
     call get_jLW(nz0, nz)
     if (rank ==0) write(6,*) 'after get_jLW, jlw(1,1,1)', jlw(1,1,1)

     ! Now merge subgrid source and real source lists into one.
     ! New NumSrc, srcpos, SrcSeries, NormFlux created.
     if (MHflag == 1 .OR. MHflag == 2) then
        call merge_allsrc
     endif

     if (rank == 0) write(6,*) 'sources merged'
#endif

     ! Set time if restart at intermediate time
     ! Set next output time
     ! Note: this assumes that there are ALWAYS two time steps
     ! between redshifts. Should be made more general.
     if (restart >= 2) then
        interm_zred=time2zred(zred2time(zred)+dt)
        if (abs(interm_zred-zred_interm) < 0.001) then
           sim_time=zred2time(zred)+dt
        else
           sim_time=zred2time(zred_interm)
        endif
        next_output_time=end_time
     else
        next_output_time=sim_time+output_time
     endif
     ! Reset restart flag now that everything has been dealt with
     restart=0

     ! Loop until end time is reached
     do
        ! Make sure you produce output at the correct time
        actual_dt=min(next_output_time-sim_time,dt)
        
        ! Report time and time step
        if (rank == 0) write(logf,"(A,2(1pe10.3,x),A)") "Time, dt:", &
             sim_time/YEAR,actual_dt/YEAR," (years)"

        ! For cosmological simulations evolve proper quantities
        if (cosmological) then
           call redshift_evol(sim_time+0.5*actual_dt)
           call cosmo_evol()
        endif
        ! Do not call in case of position dependent clumping,
        ! the clumping grid should have been initialized above
        if (type_of_clumping /= 5) call set_clumping(zred)

        ! Take one time step
        if (NumSrc > 0) call evolve3D(actual_dt,iter_restart)

        ! Reset flag for restart from iteration 
        ! (evolve3D is the last routine affected by this)
        iter_restart=0

        ! Update time
        sim_time=sim_time+actual_dt
        if (rank==0) then
        write(6,*) '............'
        write(6,*) 'sim_time, output_time', sim_time/(1e6*year), next_output_time/(1e6*year)
        write(6,*) '.............'
        write(6,*) 'output redshift', time2zred(sim_time)
        write(6,*) '.............'
        write(6,*) 'actual_dt', actual_dt/(1e6*year)
        endif
        ! Write output
        if (rank==0)        write(6,*) 'before output'
        if (abs(sim_time-next_output_time) <= 1e-6*sim_time) then
           if (NumSrc > 0) call output(time2zred(sim_time),sim_time,actual_dt,&
                photcons_flag)
           next_output_time=next_output_time+output_time
           if (photcons_flag /= 0 .and. rank == 0) &
                write(logf,*) &
                "Exiting because of photon conservation violation"
           if (photcons_flag /= 0) exit ! photon conservation violated
        endif
        if (rank==0)        write(6,*) 'after output'
        ! end time for this redshift interval reached
        if (abs(sim_time-end_time) <= 1e-6*end_time) exit
        if (rank == 0) call flush(logf)
     enddo
     if (rank==0) write(6,*) 'inner loop finished'
     if (rank==0) write(6,*) 'photcons_flag:', photcons_flag
     ! Get out: photon conservation violated
     if (photcons_flag /= 0 .and. rank == 0) &
          write(logf,*) "Exiting because of photon conservation violation"
     if (photcons_flag /= 0) exit ! photon conservation violated

     ! Scale to the current redshift
     if (cosmological) then
        call redshift_evol(sim_time)
        call cosmo_evol()
     endif

     ! Find out intermediate CPU time (to avoid overflowing the counter)
     call cpu_time(cputime2)
     cpu_seconds=cpu_seconds+real(cputime2-cputime1)
     cputime1=cputime2
     cpu_minutes = cpu_minutes + int(cpu_seconds) / 60
     cpu_seconds = MOD ( cpu_seconds , 60.0 )
     cpu_hours = cpu_hours + cpu_minutes / 60
     cpu_minutes = MOD ( cpu_minutes , 60 )

     call system_clock(cntr2,countspersec)
     clock_seconds=clock_seconds+real(cntr2-cntr1)/real(countspersec)
     cntr1=cntr2
     clock_minutes = clock_minutes + int(clock_seconds) / 60
     clock_seconds = MOD ( clock_seconds , 60.0 )
     clock_hours = clock_hours + clock_minutes / 60
     clock_minutes = MOD ( clock_minutes , 60 )

  enddo

  ! Write final output

  if (photcons_flag == 0) call output(zred,sim_time,actual_dt,photcons_flag)

  call close_down ()
  
  ! Find out CPU time
  call cpu_time(cputime2)
  cpu_seconds=cpu_seconds+real(cputime2-cputime1,dp)
  cpu_minutes = cpu_minutes + int(cpu_seconds) / 60
  cpu_seconds = MOD ( cpu_seconds , 60.0 )
  cpu_hours = cpu_hours + cpu_minutes / 60
  cpu_minutes = MOD ( cpu_minutes , 60 )

  ! Find out wall clock time
  call system_clock(cntr2,countspersec)
  clock_seconds=clock_seconds+real(cntr2-cntr1,dp)/real(countspersec,dp)
  clock_minutes = clock_minutes + int(clock_seconds) / 60
  clock_seconds = MOD ( clock_seconds , 60.0 )
  clock_hours = clock_hours + clock_minutes / 60
  clock_minutes = MOD ( clock_minutes , 60 )
  
  if (rank == 0) then
     write(logf,*) "CPU time: ",cpu_hours,' hours',cpu_minutes,' minutes', &
          cpu_seconds,' seconds.'
     write(logf,*) "Wall clock time: ",clock_hours,' hours', &
          clock_minutes,' minutes',clock_seconds,' seconds.'
  endif

  ! End the run
  call mpi_end ()

end Program C2Ray
