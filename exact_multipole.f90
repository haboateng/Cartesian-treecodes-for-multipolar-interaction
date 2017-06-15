      PROGRAM TREEDRIVER
      IMPLICIT NONE

      INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(12)

! runtime parameters

      INTEGER :: numpars
      REAL(KIND=r8) :: t1,t2,rmsE,temp 

! arrays for coordinates, charge, force calculations 

      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:,:) :: mpole
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: x,y,z ! particles
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: denergy

! variables for potential energy computation

      REAL(KIND=r8) :: tpeng,dpeng

! variables needed for f90 DATE_AND_TIME intrinsic

      INTEGER,DIMENSION(8) :: time1,time2 
      CHARACTER (LEN=8)  :: datec
      CHARACTER (LEN=10) :: timec
      CHARACTER (LEN=5)  :: zonec
      CHARACTER (LEN=20) :: sampin1,sampin2,sampout
      REAL(KIND=r8)      :: totaltime,timedirect

! variables for error calculations

      REAL(KIND=r8) :: inferr,relinferr

! local variables

      INTEGER :: i,j,err,mpdim

! EXECUTABLE STATEMENTS
      WRITE(6,*) 'Enter the name of input file 1'
      READ(5,*) sampin1
      WRITE(6,*) 'Enter the name of output file'
      READ(5,*) sampout
      WRITE(6,*) 'Enter NUMPARS : '
      READ(5,*)  numpars

      mpdim=35
      ALLOCATE(x(numpars),y(numpars),z(numpars),mpole(numpars,1:mpdim),STAT=err)
      IF (err .NE. 0) THEN
          WRITE(6,*) 'Error allocating arrays for x, y, z and mpole! '
          STOP
      END IF


      ALLOCATE(denergy(numpars),STAT=err)
      IF (err .NE. 0) THEN
          WRITE(6,*) 'Error allocating tenergy or denergy! '
          STOP
      END IF

      OPEN(unit=82,file=sampin1,status='old',action='read')
 
 
      DO i=1,numpars
         READ(82,*) x(i),y(i),z(i)!,qS(i)
      END DO
     
      CLOSE(82)

      CALL RANDOM_SEED()
      CALL RANDOM_NUMBER(mpole); mpole = -1.0_r8 + 2.0_r8*mpole

      OPEN(unit=85,file=sampout,status='unknown',action='write',&
           position='append')

      CALL DATE_AND_TIME(datec,timec,zonec,time1)

      CALL DIRECT_ENG(x,y,z,mpole,numpars,mpdim,denergy,dpeng)

      CALL DATE_AND_TIME(datec,timec,zonec,time2)
      CALL TTIME(time1,time2,timedirect)


      WRITE(85,13)timedirect,dpeng
      DO j=1,numpars
         WRITE(85,14)j,denergy(j)
      END DO

      CLOSE(unit=85)

      DEALLOCATE(x,y,z,mpole,denergy)
 13   FORMAT(2X,E24.16,2X,E24.16)
 14   FORMAT(I12,2X,E24.16)
      STOP
      END PROGRAM TREEDRIVER

!!!!!!!!!!!!!!


      SUBROUTINE DIRECT_ENG(x,y,z,mpole,numpars,mpdim,denergy,dpeng) 

      IMPLICIT NONE

      INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(12)
      INTEGER,INTENT(IN) :: numpars,mpdim
      REAL(KIND=r8),DIMENSION(numpars),INTENT(IN) :: x,y,z
      REAL(KIND=r8),DIMENSION(numpars,1:mpdim),INTENT(IN) :: mpole
      REAL(KIND=r8),DIMENSION(numpars),INTENT(INOUT) :: denergy
      REAL(KIND=r8),INTENT(INOUT) :: dpeng

! local variables 
   
      INTEGER :: i,j
      REAL(KIND=r8) :: tx,ty,tz,xi,yi,zi,teng
      REAL(KIND=r8) :: jmp(1:mpdim),rr,rr2,rr3,rr5,rr7,frr2,tfrr2,srr2
      REAL(KIND=r8) :: tx2,ty2,tz2,txy,txz,tyz,txyz
      REAL(KIND=r8) :: imp1,imp2,imp3,imp4,imp5,imp6,imp7,imp8,imp9,imp10
      REAL(KIND=r8) :: imp11,imp12,imp13,imp14,imp15,imp16,imp17,imp18,imp19
      REAL(KIND=r8) :: imp20,imp21,imp22,imp23,imp24,imp25,imp26,imp27,imp28
      REAL(KIND=r8) :: imp29,imp30,imp31,imp32,imp33,imp34,imp35

      dpeng=0.0_r8; denergy=0.0_r8
      DO i=1,numpars-1
         xi=x(i)
         yi=y(i)
         zi=z(i)
         imp1=mpole(i,1);imp2=mpole(i,2);imp3=mpole(i,3);imp4=mpole(i,4);imp5=mpole(i,5);imp6=mpole(i,6)
         imp7=mpole(i,7);imp8=mpole(i,8);imp9=mpole(i,9);imp10=mpole(i,10);imp11=mpole(i,11);imp12=mpole(i,12)
         imp13=mpole(i,13);imp14=mpole(i,14);imp15=mpole(i,15);imp16=mpole(i,16);imp17=mpole(i,17);imp18=mpole(i,18)
         imp19=mpole(i,19);imp20=mpole(i,20);imp21=mpole(i,21);imp22=mpole(i,22);imp23=mpole(i,23);imp24=mpole(i,24)
         imp25=mpole(i,25);imp26=mpole(i,26);imp27=mpole(i,27);imp28=mpole(i,28);imp29=mpole(i,29);imp30=mpole(i,30)
         imp31=mpole(i,31);imp32=mpole(i,32);imp33=mpole(i,33);imp34=mpole(i,34);imp35=mpole(i,35)

         teng=0.0_r8
         DO j=i+1,numpars
            tx=xi-x(j)
            ty=yi-y(j)
            tz=zi-z(j)
            jmp=mpole(j,:)
            tx2=tx*tx; ty2=ty*ty; tz2=tz*tz; txy=tx*ty; txz=tx*tz; tyz=ty*tz; txyz=tx*ty*tz
            rr2=1.0_r8/(tx*tx+ty*ty+tz*tz); rr=SQRT(rr2); frr2=5.0_r8*rr2; srr2=7.0_r8*rr2; tfrr2=35.0_r8*rr2
            rr3=rr2*rr; rr5=3.0_r8*rr2*rr3; rr7=5.0_r8*rr5*rr2

            ! Differentiate with respect to y (The L_j operator)
 
            teng=teng+jmp(1)*rr+rr3*(jmp(2)*tx+jmp(3)*ty+jmp(4)*tz)-rr3*(jmp(5)+jmp(8)+jmp(10))+&
                 rr5*(tx*(jmp(5)*tx+jmp(6)*ty+jmp(7)*tz)+ty*(jmp(8)*ty+jmp(9)*tz)+tz*jmp(10)*tz-&
                 jmp(11)*tx*(3.0_r8-tx2*frr2)-jmp(17)*ty*(3.0_r8-ty2*frr2)-jmp(20)*tz*(3.0_r8-tz2*frr2)-&
                 (jmp(12)*ty+jmp(13)*tz)*(1.0_r8-tx2*frr2)-(jmp(14)*tx+jmp(18)*tz)*(1.0_r8-ty2*frr2)-&
                 (jmp(16)*tx+jmp(19)*ty)*(1.0_r8-tz2*frr2))+jmp(15)*txyz*rr7+&
                 rr5*(jmp(21)*(3.0_r8+tx2*rr2*(tfrr2*tx2-30.0_r8))+jmp(24)*(1.0_r8-frr2*(tx2+ty2-srr2*tx2*ty2))+&
                 jmp(26)*(1.0_r8-frr2*(tx2+tz2-srr2*tx2*tz2))+jmp(31)*(3.0_r8+ty2*rr2*(tfrr2*ty2-30.0_r8))+&
                 jmp(33)*(1.0_r8-frr2*(ty2+tz2-srr2*ty2*tz2))+jmp(35)*(3.0_r8+tz2*rr2*(tfrr2*tz2-30.0_r8)))+&
                 rr7*(jmp(22)*txy*(srr2*tx2-3.0_r8)+jmp(23)*txz*(srr2*tx2-3.0_r8)+jmp(25)*tyz*(srr2*tx2-1.0_r8)+&
                 jmp(27)*txy*(srr2*ty2-3.0_r8)+jmp(28)*txz*(srr2*ty2-1.0_r8)+jmp(29)*txy*(srr2*tz2-1.0_r8)+&
                 jmp(30)*txz*(srr2*tz2-3.0_r8)+jmp(32)*tyz*(srr2*ty2-3.0_r8)+jmp(34)*tyz*(srr2*tz2-3.0_r8))

            ! Differentiate with respect to x (The L_i operator)

            denergy(j)=denergy(j)+imp1*rr-rr3*(imp2*tx+imp3*ty+imp4*tz)-rr3*(imp5+imp8+imp10)+&
                 rr5*(tx*(imp5*tx+imp6*ty+imp7*tz)+ty*(imp8*ty+imp9*tz)+tz*imp10*tz+&
                 imp11*tx*(3.0_r8-tx2*frr2)+imp17*ty*(3.0_r8-ty2*frr2)+imp20*tz*(3.0_r8-tz2*frr2)+&
                 (imp12*ty+imp13*tz)*(1.0_r8-tx2*frr2)+(imp14*tx+imp18*tz)*(1.0_r8-ty2*frr2)+&
                 (imp16*tx+imp19*ty)*(1.0_r8-tz2*frr2))-imp15*txyz*rr7+&
                 rr5*(imp21*(3.0_r8+tx2*rr2*(tfrr2*tx2-30.0_r8))+imp24*(1.0_r8-frr2*(tx2+ty2-srr2*tx2*ty2))+&
                 imp26*(1.0_r8-frr2*(tx2+tz2-srr2*tx2*tz2))+imp31*(3.0_r8+ty2*rr2*(tfrr2*ty2-30.0_r8))+&
                 imp33*(1.0_r8-frr2*(ty2+tz2-srr2*ty2*tz2))+imp35*(3.0_r8+tz2*rr2*(tfrr2*tz2-30.0_r8)))+&
                 rr7*(imp22*txy*(srr2*tx2-3.0_r8)+imp23*txz*(srr2*tx2-3.0_r8)+imp25*tyz*(srr2*tx2-1.0_r8)+&
                 imp27*txy*(srr2*ty2-3.0_r8)+imp28*txz*(srr2*ty2-1.0_r8)+imp29*txy*(srr2*tz2-1.0_r8)+&
                 imp30*txz*(srr2*tz2-3.0_r8)+imp32*tyz*(srr2*ty2-3.0_r8)+imp34*tyz*(srr2*tz2-3.0_r8))

         END DO
         denergy(i)=denergy(i)+teng
      END DO
      
      dpeng=SUM(denergy)
      
      END SUBROUTINE DIRECT_ENG  

!!!!!!!!!!!!!!!!!!
      SUBROUTINE TTIME(timebeg,timeend,totaltime)
      IMPLICIT NONE
!
! TTIME computes the time difference in seconds between
! the timestamps TIMEBEG and TIMEEND returned by the 
! f90 intrinsic DATE_AND_TIME
!
      INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(12)
      INTEGER,DIMENSION(8),INTENT(INOUT) :: timebeg,timeend
      REAL(KIND=r8),INTENT(OUT) :: totaltime

! TIMEEND is modifed by borrowing in case each of its fields
! are not .GE. to the corresponding field in TIMEBEG (up to
! and including days) 

      IF (timeend(8) .LT. timebeg(8)) THEN
          timeend(8)=timeend(8)+1000
          timeend(7)=timeend(7)-1
      END IF
      IF (timeend(7) .LT. timebeg(7)) THEN
          timeend(7)=timeend(7)+60
          timeend(6)=timeend(6)-1
      END IF
      IF (timeend(6) .LT. timebeg(6)) THEN
          timeend(6)=timeend(6)+60
          timeend(5)=timeend(5)-1
      END IF
      IF (timeend(5) .LT. timebeg(5)) THEN
          timeend(5)=timeend(5)+24
          timeend(3)=timeend(3)-1
      END IF

      totaltime=  REAL(timeend(8)-timebeg(8),KIND=r8) +          &
            1000.0_r8*( REAL(timeend(7)-timebeg(7),KIND=r8) +    &
              60.0_r8*( REAL(timeend(6)-timebeg(6),KIND=r8) +    &
              60.0_r8*( REAL(timeend(5)-timebeg(5),KIND=r8) +    &
              24.0_r8*( REAL(timeend(3)-timebeg(3),KIND=r8)))))
      totaltime=totaltime/1000.0_r8


      RETURN
      END SUBROUTINE TTIME
   
