!***********************************************!
!THE PROGRAM BEGINS.
!***********************************************!
program fes_calc
implicit none
integer i, iter, ncv, umbr_n, nbin1, nbin2
real*8, allocatable :: grid0(:,:), v(:), prob(:,:),biased_prob(:,:,:),grid(:,:)
real*8 :: kt, toler, dummy
logical :: parent
integer :: i_umbr, i_s1, i_s2
character (len=50) cvfile,outputfile,f1,f3
real*8 cnvg
real*8, allocatable :: umbr_mean(:), umbr_k(:),a(:)
integer, allocatable :: nmax(:)
real*8 :: dummy_1,dummy_2,dummy_3
integer :: rank, gleng1_min,gleng1_max,gleng2,ngrid

CALL MPI_Start
CALL Set_Parent(parent)


if(parent)then
  open (unit=1, file='input',status='old')
  !read number of cv, number of umbrella,kt in energy unit
  read(1,*)ncv, umbr_n, kt
end if

!broadcast ncv
CALL IBcast(ncv,1)
CALL IBcast(umbr_n,1)
CALL RBcast(kt,1)

!allocate grid0
allocate(grid(3,ncv))
allocate(grid0(3,ncv))
allocate(a(umbr_n))
allocate(umbr_mean(umbr_n))
allocate(umbr_k(umbr_n))
allocate(nmax(umbr_n))
a=1.d0

if (parent) then
 !read grid_min, grid_max, grid_bin
 do i=1,ncv
    read(1,*)grid0(1:3,i)
    write(*,'(i10,3f16.6)')i,grid0(1:3,i)
 end do
end if

if (parent) then
if (ncv.eq.2) then
nbin1=nint((grid0(2,1)-grid0(1,1))/grid0(3,1))+1 
nbin2=nint((grid0(2,2)-grid0(1,2))/grid0(3,2))+1
else if (ncv.ge.3)then
STOP '3 or more CVs not implemented'
end if
end if

!write(*,*) nbin1,nbin2 

CALL RBcast(grid0,3*ncv)
CALL IBcast(nbin1,1)
CALL IBcast(nbin2,1)

allocate(biased_prob(nbin1,nbin2,umbr_n))
allocate(prob(nbin1,nbin2))

if (parent) then
  open( unit =2, file= 'whaminput', status = 'old') 
   read(2,*)toler
   do i_umbr=1, umbr_n 
   write(*,*) 'umbrella simulation #', i_umbr
   !reads probability file name
   read(2,'(a)') cvfile
   write(*,*)'CVFILE',cvfile
   !reads force constant, r0, max points in that umbrella
   read(2,*) umbr_mean(i_umbr),umbr_k(i_umbr),nmax(i_umbr) 
       umbr_k(i_umbr)=umbr_k(i_umbr)*627.51 !converted in kcal
       umbr_mean(i_umbr)=umbr_mean(i_umbr)/0.529 !converted in au
   !reads the probability written in prob file for each umbrella
      f1 = cvfile 
      open( unit=3, file=f1, status='old' )
      do i_s1=1,nbin1 !US
      do i_s2=1,nbin2 !MTD
      read(3,*)dummy,dummy,dummy
      biased_prob(i_s1,i_s2,i_umbr)=dexp(-dummy/kt)
      end do
      end do
  enddo
end if

call RBcast(toler,1)
call RBcast(umbr_mean,umbr_n)
call RBcast(umbr_k,umbr_n)
call IBcast(nmax,umbr_n)
call RBcast(biased_prob,nbin1*nbin2*umbr_n)

!PERFORMS WHAM.
 if (parent)  write(*,*) 'wham begins'
 CALL DistributeGrids(ncv,grid0,grid,rank,gleng1_min,gleng1_max)

 gleng2=nint((grid(2,2)-grid(1,2))/grid(3,2))+1

 write(*,*)'new_grid', gleng1_min, gleng1_max,gleng2, rank, grid(1,1)
  iter=0
  scf_loop : do
       iter=iter+1
       call wham_scf(a,umbr_k,umbr_mean,nmax,grid0,biased_prob,umbr_n,nbin1,nbin2,ncv,kt,prob,cnvg,gleng1_min,gleng1_max,gleng2)
       if (parent) write(*,*) 'iteration #', iter, 'convergence =',cnvg
       if (mod(iter,100) .eq. 0 ) then
       !prints the free energy.
       call print_pmf(ncv,nbin1,nbin2,grid0,kt,prob,parent)
       endif
       if(cnvg.lt.toler)then
       write(*,*)'** convergence achieved **'
       exit scf_loop
       end if
   end do scf_loop

    !prints the free energy.
    call print_pmf(ncv,nbin1,nbin2,grid0,kt,prob,parent)
call MPI_Stop
end program 

!***********************************************!


 subroutine print_pmf(ncv,nbin1,nbin2,grid0,kt,prob,parent)
!prints 
 implicit none
 integer i_s1, i_s2
 integer nbin1, nbin2, ncv
 real*8 :: prob(nbin1,nbin2), grid0(3,ncv)
 real*8 s1, s2, kt
 real*8, allocatable :: prob1(:,:)
 logical parent
 character (len=50) f2
 
 allocate(prob1(nbin1,nbin2))
 CALL GlobSumR(prob,prob1,nbin1*nbin2)
 if (parent)then
 f2= 'free_energy'
 open( unit =7 , file = f2, status =  'unknown' )
 s1=grid0(1,1)
 do i_s1=1,nbin1 !US cv
 s2=grid0(1,2)
  do i_s2=1,nbin2 !MTD cv
    write(7,*) s1, s2, -kt*dlog( prob1(i_s1,i_s2) )
    s2= s2 + grid0(3,2)
  enddo
    s1= s1 + grid0(3,1)
    write(7,*)
 enddo
 write(*,*) 'free energy written in ',f2
 close(7)
 end if
 endsubroutine
!**************************************************************************************************!


 subroutine wham_scf(a,umbr_k,umbr_mean,nmax,grid0,biased_prob,umbr_n,nbin1,nbin2,ncv,kt,prob,cnvg,gleng1_min,gleng1_max,gleng2)
 !performs wham scf.
 implicit none
 integer i_s1, i_s2, i_umbr, nbin1, nbin2, umbr_n,ncv
 integer  :: nmax(umbr_n) 
 real*8 :: umbr_k(umbr_n), umbr_mean(umbr_n)
 real*8 :: num, den, dummy_v, dummy_s1, avg, del_s1, kt, dummy, cnvg
 real*8, allocatable :: dummy_a(:) 
 real*8 ::  prob(nbin1,nbin2),a(umbr_n)
 real*8 ::  grid0(3,ncv), biased_prob(nbin1,nbin2,umbr_n),grid(3,ncv)
 integer :: rank, gleng1_min,gleng1_max,gleng2,ngrid
 real*8,allocatable :: dummy_a1(:)

 allocate(dummy_a(umbr_n))
 allocate (dummy_a1(umbr_n) )

 dummy_a = 0.0d0

 !calculates probability at each grid_point.
 do i_s1 =gleng1_min,gleng1_max !over US cv
 do i_s2 =1,gleng2 !over MTD cv

 !sets the numerator and denominator of prob to 0.
     num = 0.0d0
     den = 0.0d0
     dummy_s1 = grid0(1,1)+dfloat(i_s1-1)*grid0(3,1)
 !calculates probability.
   do i_umbr=1,umbr_n
     del_s1 = dummy_s1 - umbr_mean(i_umbr)
     dummy_v = dexp(-(0.50d0*umbr_k(i_umbr)*del_s1*del_s1/kt))
     num = num + dfloat(nmax(i_umbr))* biased_prob(i_s1,i_s2,i_umbr)
     den = den + dfloat(nmax(i_umbr))*a(i_umbr)*dummy_v
   enddo
 prob(i_s1,i_s2) = num/den
 if(prob(i_s1,i_s2).ne.prob(i_s1,i_s2)) prob(i_s1,i_s2) = 1.0e-15 !remove NaN
 if(prob(i_s1,i_s2)+1.eq.prob(i_s1,i_s2)) prob(i_s1,i_s2) = 1.0e-15 !remove infinity

 !CALCULATES A's.
   do i_umbr=1,umbr_n
      del_s1 = dummy_s1 - umbr_mean(i_umbr)
      dummy_v = dexp(-(0.50d0*umbr_k(i_umbr)*del_s1*del_s1/kt))
      dummy_a(i_umbr) = dummy_a(i_umbr) +grid0(3,1)*grid0(3,2)*dummy_v*prob(i_s1,i_s2)
   enddo
 enddo !end of MTD cv loop
 enddo !end of US cv loop

   CALL GlobSumR(dummy_a,dummy_a1,umbr_n)
   do i_umbr=1,umbr_n
   dummy_a(i_umbr)=dummy_a1(i_umbr)
   end do

!FINDS CONVERGENCE AND UPDATES A's.
 avg = 0.0d0
 cnvg = 0.0d0
 do i_umbr=1,umbr_n
 dummy_a(i_umbr) = 1.0d0/dummy_a(i_umbr)
 cnvg = cnvg + dabs(dlog(dummy_a(i_umbr))-dlog(a(i_umbr)))
 a(i_umbr) = dummy_a(i_umbr)
 enddo
 cnvg = kt*cnvg
 end subroutine
!**************************************************************************************************!



SUBROUTINE DistributeGrids(ncv,grid0,grid,rank,gleng1_min,gleng1_max)
!Distribute X grid over processors by mapping grid0 to grid 
IMPLICIT NONE
INTEGER :: ncv,rank
REAL*8 :: grid0(3,ncv), grid(3,ncv)
!
INTEGER :: i,ncpu,icpu,ngrids,ngrids_m,ngrids_y,ngrids_z,ngrids_o
INTEGER :: gleng1_min, gleng1_max

CALL Get_ncpu(ncpu)
CALL Get_cpuid(icpu)

write(*,*)'NCPUS',ncpu
DO i=1,ncv
  grid(1:3,i)=grid0(1:3,i)
END DO
rank=0

IF(ncv.EQ.1)THEN
 ngrids_y=1
 ngrids_z=1
ELSE IF(ncv.EQ.2)THEN
 ngrids_y=(nint((grid0(2,2)-grid0(1,2))/grid0(3,2))+1)
 ngrids_z=1
ELSE IF(ncv.EQ.3)THEN
 ngrids_y=(nint((grid0(2,2)-grid0(1,2))/grid0(3,2))+1)
 ngrids_z=(nint((grid0(2,3)-grid0(1,3))/grid0(3,3))+1)
END IF

if(ncpu.eq.1) then 
gleng1_min=1
gleng1_max=nint((grid0(2,1)-grid0(1,1))/grid0(3,1))+1
end if

!Distribute X grids
if(icpu.eq.0)WRITE(*,'(3A12,3A16)') 'CPU','CV', 'GRID SIZE', 'GRID MIN', 'GRID MAX', 'GRID BIN'
CALL Sync_procs
IF(ncpu.GT.1)THEN
  ngrids=nint((grid0(2,1)-grid0(1,1))/grid0(3,1))+1
  write(*,*)'NEW_GRID',ngrids,icpu,ncpu
  ngrids_o=ngrids
  ngrids=ngrids/ncpu
  IF(icpu.eq.ncpu-1)THEN
    ngrids_m=ngrids+mod(ngrids_o,ncpu)
    grid(1,1)=grid(1,1)+DFLOAT(icpu*ngrids)*grid0(3,1)
    grid(2,1)=grid(1,1)+DFLOAT(ngrids_m-1)*grid0(3,1)
  ELSE
    ngrids_m=ngrids
    grid(1,1)=grid(1,1)+DFLOAT(icpu*ngrids)*grid0(3,1)
    grid(2,1)=grid(1,1)+DFLOAT(ngrids-1)*grid0(3,1)
  END IF
  CALL Sync_procs
!  WRITE(*,'(3I12,3F16.6)') icpu, 1, ngrids_m, grid(1,1), grid(2,1), grid(3,1)
  WRITE(*,*) 'ALL_INFO',icpu, 1, ngrids_m, grid(1,1), grid(2,1), grid(3,1)
  rank=ngrids_z*ngrids_y*ngrids*icpu
  gleng1_min=ngrids*icpu+1
  gleng1_max=ngrids*(icpu+1)
  if(icpu.eq.ncpu-1) gleng1_max=ngrids*(icpu+1)+mod(ngrids_o,ncpu)
END IF 
END 
!****************************************************************************************!



!****************************************************************************************!

SUBROUTINE MPI_Start()
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER :: i_err
!#if defined (_PARALLEL)
  call MPI_INIT(i_err)
!#endif
END
!****************************************************************************************!

SUBROUTINE MPI_Stop()
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER :: i_err
!#if defined (_PARALLEL)
call MPI_FINALIZE(i_err)
!#endif
END
!****************************************************************************************!

SUBROUTINE get_ncpu(ncpu)
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER :: ncpu, i_err
!ncpu=1
!#if defined (_PARALLEL)
call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,i_err)
!#endif
END
!****************************************************************************************!

SUBROUTINE get_cpuid(icpu)
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER :: icpu, i_err
!icpu=0
!#if defined (_PARALLEL)
call MPI_COMM_RANK(MPI_COMM_WORLD,icpu,i_err)
!#endif
END
!****************************************************************************************!

SUBROUTINE IBcast(myint,leng)
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER :: leng, myint(*), i_err
!#if defined (_PARALLEL)
CALL MPI_BCAST(myint,leng,MPI_INTEGER,0,MPI_COMM_WORLD,i_err)
!#endif
END
!****************************************************************************************!

SUBROUTINE RBcast(myreal,leng)
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER :: myreal(*)
INTEGER :: leng, i_err
!#if defined (_PARALLEL)
CALL MPI_BCAST(myreal,leng,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,i_err)
!#endif
END
!****************************************************************************************!

SUBROUTINE Sync_Procs
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER i_err
!#if defined (_PARALLEL)
call MPI_Barrier(MPI_COMM_WORLD,i_err)
!#endif
END
!****************************************************************************************!

SUBROUTINE Set_Parent(parent)
IMPLICIT NONE
INCLUDE 'mpif.h'
LOGICAL :: parent
INTEGER :: icpu, i_err
parent=.false.
!icpu=0
!#if defined (_PARALLEL)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,icpu,i_err)
!#endif
IF(icpu.eq.0)parent=.true.
END
!****************************************************************************************!

SUBROUTINE GlobSumR(myreal_in,myreal_out,leng)
IMPLICIT NONE
!#if defined (_PARALLEL)
INCLUDE 'mpif.h'
!#endif
REAL*8 :: myreal_in(*), myreal_out(*)
INTEGER :: leng,i_err
!#if defined (_PARALLEL)
CALL MPI_Allreduce(myreal_in,myreal_out,leng,MPI_DOUBLE_PRECISION, &
                   MPI_SUM,MPI_COMM_WORLD,i_err)
!#endif
END
!****************************************************************************************!

SUBROUTINE GlobSumI(myint_in,myint_out,leng)
IMPLICIT NONE
!#if defined (_PARALLEL)
INCLUDE 'mpif.h'
!#endif
INTEGER :: myint_in(*), myint_out(*)
INTEGER :: leng,i_err
!#if defined (_PARALLEL)
call MPI_Allreduce(myint_in,myint_out,leng,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,i_err)
!#endif
END
!****************************************************************************************!
