module mergesrc
  use precision,   only: dp
  use my_mpi
  use file_admin,  only: logf
  use sourceprops, only: NumSrc,      SrcSeries,     srcpos,    NormFlux
  use source_sub,  only: NumAGrid, subSrcSeries, sub_srcpos, subNormFlux


contains
! ======================================================================
  subroutine merge_allsrc

    ! temporary, randomized list of sources
    integer,dimension(:),  allocatable     :: TempSrcSeries
    ! temporary src position
    integer,dimension(:,:),allocatable     :: Tempsrcpos
    ! temporary normalized ionizing flux of sources
    real(kind=dp),dimension(:),allocatable :: TempNormFlux
    ! temporary number of sources
    integer                                :: TempNumSrc

#ifdef MPI
    integer :: mympierror
#endif

    if (NumAGrid > 0 .and. NumSrc > 0) then

       if (rank == 0) then
          write(logf,*) 'BOTH SOURCES (larger & smaller than 10^8 Msun) exist'
          write(logf,*) Numsrc, ' sources have mass > 10^8 Msun'
          write(logf,*) NumAGrid, ' grids(!) have mass < 10^8 Msun'
          write(6,*)  'BOTH SUBSOURCES exist', Numsrc+NumAGrid
          write(6,*)   Numsrc, ' sources have mass > 10^8 Msun'
          write(6,*)   NumAGrid, ' grids(!) have mass < 10^8 Msun'
       endif

       if (allocated(TempSrcSeries)) deallocate(TempSrcSeries)
       if (allocated(Tempsrcpos   )) deallocate(Tempsrcpos   )
       if (allocated(TempNormFlux )) deallocate(TempNormFlux )

       allocate(TempSrcSeries(NumSrc + NumAGrid))
       allocate(Tempsrcpos(3, NumSrc + NumAGrid))
       allocate(TempNormFlux (NumSrc + NumAGrid))

       TempSrcSeries(   1:NumSrc) = SrcSeries   (   1:NumSrc)
       Tempsrcpos   (:, 1:NumSrc) = srcpos      (:, 1:NumSrc)
       TempNormFlux (   1:NumSrc) = NormFlux    (   1:NumSrc)

       TempSrcSeries(   NumSrc+1:NumSrc+NumAGrid) = subSrcSeries(   1:NumAGrid) + NumSrc
       Tempsrcpos   (:, NumSrc+1:NumSrc+NumAGrid) = sub_srcpos  (:, 1:NumAGrid)
       TempNormFlux (   NumSrc+1:NumSrc+NumAGrid) = subNormFlux (   1:NumAGrid)

       ! NumSrc is now counting all the sources.
       TempNumSrc = NumSrc + NumAGrid

       CALL MPI_BARRIER(MPI_COMM_NEW,mympierror)
       NumSrc     = TempNumSrc

       deallocate(SrcSeries)
       deallocate(srcpos)
       deallocate(NormFlux)
       
       allocate(SrcSeries(NumSrc))
       allocate(srcpos(3, NumSrc))
       allocate(NormFlux (NumSrc))
       
       SrcSeries = TempSrcSeries
       srcpos    = Tempsrcpos
       NormFlux  = TempNormFlux
       
       deallocate (TempSrcSeries)
       deallocate (Tempsrcpos)
       deallocate (TempNormFlux)
       
       CALL MPI_BARRIER(MPI_COMM_NEW,mympierror)

    elseif (NumAGrid > 0 .and. NumSrc == 0) then
       if (rank == 0) then
          write(logf,*) 'ONLY SUBGRID SOURCES (smaller than 10^8 Msun) exist'
          write(logf,*) NumAGrid, ' grids(!) have mass < 10^8 Msun'
          write(6,*)  'ONLY SUBSOURCES exist, whose number is ', NumAGrid
       endif

       if (allocated(SrcSeries)) deallocate(SrcSeries)
       if (allocated(srcpos   )) deallocate(srcpos   )
       if (allocated(NormFlux )) deallocate(NormFlux )
       
       NumSrc = NumAGrid
       
       allocate(SrcSeries(NumSrc))
       allocate(srcpos(3, NumSrc))
       allocate(NormFlux (NumSrc))
       

       SrcSeries = subSrcSeries
       srcpos    = sub_srcpos  
       NormFlux  = subNormFlux 
    endif

  end subroutine merge_allsrc

end module mergesrc
