!---------------------------------------------------------------------!
!Written by Shalini Awasthi (ashalini@iitk.ac.in)
!---------------------------------------------------------------------!
      PROGRAM WSMTD_rw_2D
      IMPLICIT NONE
      REAL*8 gridmin1, gridmax1, griddif1, dummy2,v,
     &       gridmin2, gridmax2, griddif2,num,
     &       prob,den,alpha,fes,fes1,grid,Ro,
     &       cv1,cv2,ht,kt,ktb,bias_fact,ct,vbias,hill,width,
     &       diff_s2,ds2,ss,hh,dum,s1,s2,dummy11
      ALLOCATABLE cv1(:),cv2(:),ht(:),vbias(:),ct(:),hill(:),width(:),
     &            prob(:,:),fes(:,:),fes1(:),grid(:)
      INTEGER mtd_steps,md_steps,dummy1,i,j,index1,index2,k,
     &        t_min,t_max,
     &        nbin1, nbin2,
     &        w_cv,w_hill,
     &        i_mtd, i_md, i_s2, i_s1,mtd_max,
     &        narg
      LOGICAL pmf,inpgrid
      CHARACTER*120 :: arg 
      REAL*8, PARAMETER :: kb=1.9872041E-3 !kcal K-1 mol-1
      REAL*8, PARAMETER :: au_to_kcal = 627.51 

!      character(len=50) f1

      OPEN(11,FILE='colvar_mtd',STATUS='unknown')
      OPEN(12,FILE='parvar_mtd',STATUS='unknown')
      OPEN(13,FILE='CONSTRAINT',STATUS='unknown')
      OPEN(14,FILE='cvmdck_mtd',STATUS='unknown')
!
      CALL get_steps(11,mtd_steps)
      CALL get_steps(14,md_steps)
!
!      CALL check_files(11,w_hill)
!      CALL check_files(14,w_cv)
!      CALL check_files(13,i)
!      IF(i.ne.w_cv)THEN
!        print *, 
!     &    '!!Print Freq. of cvmdck_mtd & CONSTRAINTS are not equal!!'
!        print *, '    cvmdck_mtd Print Frq. =',w_cv
!        print *, '    CONSTRAINTS Print Frq. =',i
!        STOP 
!      END IF 

      kt=300.D0
      bias_fact=1500.D0
      t_min=1
      t_max=md_steps
      pmf=.FALSE.
      inpgrid=.false.
      t_max=md_steps
      narg = IARGC()
      DO i=1,narg
        CALL GETARG(i,arg)
        IF(INDEX(arg,'-T').NE.0)THEN
           CALL GETARG(i+1,arg)
           READ(arg,*)kt
        ELSE IF(INDEX(arg,'-dT').NE.0)THEN
           CALL GETARG(i+1,arg)
           READ(arg,*)bias_fact
        ELSE IF(INDEX(arg,'-tmin').NE.0)THEN
           CALL GETARG(i+1,arg)
           READ(arg,*)t_min
        ELSE IF(INDEX(arg,'-tmax').NE.0)THEN
           CALL GETARG(i+1,arg)
           READ(arg,*)t_max
           IF(t_max.gt.md_steps)STOP '!!ERROR: t_max > total MD steps'
        ELSE IF(INDEX(arg,'-pmf').NE.0)THEN
           pmf=.true.
        ELSE IF(INDEX(arg,'-grid').NE.0)THEN
            CALL GETARG(i+1,arg)
            READ(arg,*)gridmin1
            CALL GETARG(i+2,arg)
            READ(arg,*)gridmax1
            CALL GETARG(i+3,arg)
            READ(arg,*)griddif1
            CALL GETARG(i+4,arg)
            READ(arg,*)gridmin2
            CALL GETARG(i+5,arg)
            READ(arg,*)gridmax2
            CALL GETARG(i+6,arg)
            READ(arg,*)griddif2
            inpgrid=.true.
        ELSE IF(INDEX(arg,'-pfrqMD').NE.0)THEN
            CALL GETARG(i+1,arg)
            READ(arg,*)w_cv
        ELSE IF(INDEX(arg,'-dtMTD').NE.0)THEN
            CALL GETARG(i+1,arg)
            READ(arg,*)w_hill
        END IF
      END DO

      WRITE(*,'(A,I10)')'No: of MTD steps        =',mtd_steps
      WRITE(*,'(A,I10)')'No: of MD  steps        =',md_steps
      WRITE(*,'(A,F9.2)')'T (K)                   =',kt
      WRITE(*,'(A,F9.2)')'Bias Factor (K)         =',bias_fact
      WRITE(*,'(A,I10)')'Print Freq. cvmdck_mtd  =',w_cv
      WRITE(*,'(A,I10)')'Freq. of Hill Update    =',w_hill
      WRITE(*,'(A,I10)')'Reweigtht: Min step     =',t_min
      WRITE(*,'(A,I10)')'Reweigtht: Max step     =',t_max
      IF(pmf)WRITE(*,'(A)')'PMF will be written in PMF.dat'
!
      bias_fact=(kt+bias_fact)/kt
      kt=kb*kt
      ktb=kt*bias_fact

!Make CV file containing CV1(t), CV2(t) from CONSTRAINT and cvmdck_mtd
      ALLOCATE(cv1(md_steps),cv2(md_steps),vbias(md_steps))
      ALLOCATE(ht(mtd_steps),ct(mtd_steps),hill(mtd_steps),
     &         width(mtd_steps))

      OPEN(16,file='cv.dat',status='unknown')
      DO i_md=1,md_steps
        READ(13,*)
        READ(13,*)
        READ(13,*)
        READ(13,*)
        READ(13,*)
        READ(13,*)
        READ(13,*)
        READ(13,*)
        READ(13,*)dummy1,dummy1,Ro,cv1(i_md)
        READ(14,*)dummy1,dummy1,dummy11,cv2(i_md)
        cv1(i_md)=cv1(i_md)+Ro
        WRITE(16,'(I10,2F16.6)')i_md,cv1(i_md),cv2(i_md)
      END DO
      WRITE(*,'(A)')'CV values written in cv.dat'
!
      IF(.NOT.inpgrid) 
     &  CALL get_gridmin_max(16,gridmin1,gridmax1,griddif1,
     &                          gridmin2,gridmax2,griddif2)
      

      alpha=bias_fact/(bias_fact-1.D0)

      nbin1 = NINT((gridmax1-gridmin1)/griddif1)+1
      nbin2 = NINT((gridmax2-gridmin2)/griddif2)+1
      WRITE(*,'(7X,4A10)')'GRIDMIN','GRIDMAX','GRIDBIN','GRIDSIZE'
      WRITE(*,'(A10,3F8.4,I10)')'US  COORD:',
     &          gridmin1,gridmax1,griddif1,nbin1
      WRITE(*,'(A10,3F8.4,I10)')'MTD COORD:',
     &          gridmin2,gridmax2,griddif2,nbin2

      ALLOCATE(prob(nbin1,nbin2),fes(nbin1,nbin2))
      ALLOCATE(fes1(nbin2),grid(nbin2))

      DO i_mtd=1,mtd_steps
        READ(11,*) dummy1,hill(i_mtd),dummy2
        READ(12,*) dummy1,dummy2,width(i_mtd),ht(i_mtd)
        ht(i_mtd)=ht(i_mtd)*au_to_kcal
      END DO

      CLOSE(11 12 13 14 16)


!calculate c(t)
      DO i_s2=1,nbin2 !grid over cv2 on which MTD is being done
        grid(i_s2)=gridmin2+DFLOAT(i_s2-1)*griddif2
      END DO

      OPEN(21,FILE='ct.dat',STATUS='unknown')
      DO i_mtd=1,mtd_steps
        ds2=width(i_mtd)*width(i_mtd)
        ss=hill(i_mtd)
        hh=ht(i_mtd)
        num=0.D0
        den=0.D0
        DO i_s2=1,nbin2
           diff_s2=grid(i_s2)-ss
           diff_s2=diff_s2*diff_s2*0.5D0
           fes1(i_s2)=fes1(i_s2)-hh*DEXP(-diff_s2/ds2)
           num=num+DEXP(-fes1(i_s2)/kt)
           den=den+DEXP(-fes1(i_s2)/ktb)
        END DO
        ct(i_mtd)=kt*DLOG(num/den)
        WRITE(21,'(I10,F16.8)')i_mtd,ct(i_mtd)
      END DO
      CLOSE(21)
      WRITE(*,'(A)')'CV values written in cv.dat'

!calculate v(s,t)
      DO i_md=1,md_steps
         mtd_max=(i_md*w_cv/w_hill)
         ss=cv2(i_md)
         dum=0.d0
         DO i_mtd=1,mtd_max
           ds2=width(i_mtd)*width(i_mtd)
           hh=ht(i_mtd)/alpha
           diff_s2=ss-hill(i_mtd)
           diff_s2=diff_s2*diff_s2*0.5D0
           dum=dum+hh*DEXP(-diff_s2/ds2)
         END DO
         vbias(i_md)=dum
      END DO

!calculate prob (unbiased from MTD potential)       
      den=0.d0
      prob=0.d0
      DO i_md=1,md_steps
        IF((i_md.GT.t_min).AND.(i_md.LT.t_max))THEN
          index1 = nint((cv1(i_md)-gridmin1)/griddif1) +1
          index2 = nint((cv2(i_md)-gridmin2)/griddif2) +1
          i_mtd=i_md*w_cv/w_hill
          dum=vbias(i_md) - ct(i_mtd)
          prob(index1,index2)=prob(index1,index2)+DEXP(dum/kt)
          den=den+DEXP(dum/kt)
        END IF
      END DO
      dum=den*griddif1*griddif2
      den=1.D0/dum

      OPEN(2,FILE='Pu.dat',STATUS='unknown')
      DO i_s1=1,nbin1
        s1=DFLOAT(i_s1-1)*griddif1+gridmin1
        DO i_s2=1,nbin2
          s2=DFLOAT(i_s2-1)*griddif2+gridmin2
          prob(i_s1,i_s2)=prob(i_s1,i_s2)*den
          WRITE(2,'(3E16.8)')s1,s2,prob(i_s1,i_s2)
        END DO
        WRITE(2,*)
      END DO
      WRITE(*,'(A)')'Unbiased distribution written in Pu.dat'

!writes PMF
      IF(PMF)THEN 
        OPEN(3,FILE='PMF.dat',STATUS='unknown')
        DO i_s1=1,nbin1
           s1=DFLOAT(i_s1-1)*griddif1+gridmin1
           DO i_s2=1,nbin2
             s2=DFLOAT(i_s2-1)*griddif2+gridmin2
             prob(i_s1,i_s2)=MAX(prob(i_s1,i_s2),10E-16)
             dum=-kt*DLOG(prob(i_s1,i_s2))!calculating fes
             WRITE(3,'(3E16.8)')s1,s2,dum
           END DO
           WRITE(3,*)
        END DO
      WRITE(*,'(A)')'PMF written in PMF.dat'
!      CLOSE(3)
      END IF
      CLOSE(2 3)

      END PROGRAM WSMTD_rw_2D
!---------------------------------------------------------------------!


!---------------------------------------------------------------------!
      SUBROUTINE get_steps(iunit,nsteps)
      IMPLICIT NONE
      INTEGER iunit, nsteps
      INTEGER ios
      nsteps=0
      REWIND(iunit)
      Read_Loop: DO
         READ(iunit,*,IOSTAT=ios)
         IF(ios.ne.0)EXIT Read_Loop
         nsteps=nsteps+1
      END DO Read_Loop 
      REWIND(iunit)
      END 
!---------------------------------------------------------------------!
      SUBROUTINE check_files(iunit,dt)
      IMPLICIT NONE
      INTEGER iunit, dt
      INTEGER ii, jj,i,ios
      dt=0
      i=2
      REWIND(iunit)
      READ(iunit,*)ii
      READ(iunit,*)jj
      dt=jj-ii
      ii=jj
      RLoop: DO 
        i=i+1
        READ(iunit,*,IOSTAT=ios)jj
        IF(ios.ne.0)EXIT RLoop
        IF(jj.ne.ii+dt)THEN
           print *, '!!ERROR: Steps are not at constant stride!!'
           print *, '!!       Unit No:',iunit,'!!'
           print *, '!!       Line No:',i,'!!'
           print *, '!! Expected stride =', dt,'!!'
           print *, '!! Actual stride =', jj-ii,'!!'
           STOP
        END IF
        ii=jj
      END DO RLoop
      REWIND(iunit)
      END 
!---------------------------------------------------------------------!
      SUBROUTINE get_gridmin_max(iunit,gridmin1,gridmax1,griddif1,
     &                           gridmin2,gridmax2,griddif2)
      IMPLICIT NONE 
      INTEGER :: iunit 
      REAL*8  :: gridmin1, gridmax1, griddif1, 
     &           gridmin2, gridmax2, griddif2
      INTEGER :: ii, ios
      REAL*8  :: cv1, cv2
      INTEGER, PARAMETER :: Def_Grid_Size=101
      REWIND(iunit)
      READ(iunit,*,IOSTAT=ios)ii,cv1,cv2
      if(ios.ne.0)stop 'ERROR reading CV.dat'
      gridmin1=cv1
      gridmax1=cv1
      gridmin2=cv2
      gridmax2=cv2
      RLoop: DO 
        READ(iunit,*,IOSTAT=ios)ii,cv1,cv2
        if(ios.ne.0)EXIT RLoop
        gridmin1=MIN(gridmin1,cv1)
        gridmin2=MIN(gridmin2,cv2)
        gridmax1=MAX(gridmax1,cv1)
        gridmax2=MAX(gridmax2,cv2)
      END DO RLoop
      griddif1=(gridmax1-gridmin1)/DFLOAT(Def_Grid_Size)
      griddif2=(gridmax2-gridmin2)/DFLOAT(Def_Grid_Size)
      END
!---------------------------------------------------------------------!

