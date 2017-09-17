      PROGRAM EXACT_SUM
      IMPLICIT NONE

      INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(12)

! runtime parameters

      INTEGER :: numpars
      REAL(KIND=r8) :: t1,t2,rmsE,temp 

! arrays for coordinates, charge, force calculations 

      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:,:) :: mpole
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: x,y,z ! particles
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: denergy,fxx,fyy,fzz
      INTEGER,      ALLOCATABLE,DIMENSION(:) :: ord

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

      REAL(KIND=r8) :: randi,inferr,relinferr

! local variables

      INTEGER :: i,j,err,m,mpdim

! EXECUTABLE STATEMENTS
      WRITE(6,*) 'Enter the name of input file 1'
      READ(5,*) sampin1
      WRITE(6,*) 'Enter the name of output file'
      READ(5,*) sampout
      WRITE(6,*) 'Enter NUMPARS : '
      READ(5,*)  numpars

      mpdim=10
      ALLOCATE(x(numpars),y(numpars),z(numpars),&
               ord(numpars),mpole(numpars,1:mpdim),STAT=err)
      IF (err .NE. 0) THEN
          WRITE(6,*) 'Error allocating arrays for x, y, z, ord and mpole! '
          STOP
      END IF


      ALLOCATE(denergy(numpars),fxx(numpars),fzz(numpars),fyy(numpars),STAT=err)
      IF (err .NE. 0) THEN
          WRITE(6,*) 'Error allocating denergy or forces! '
          STOP
      END IF

      OPEN(unit=82,file=sampin1,status='old',action='read')
 
 
      mpole=0.0_r8
      CALL RANDOM_SEED()
      DO i=1,numpars
         READ(82,*) x(i),y(i),z(i)!,qS(i)
         CALL RANDOM_NUMBER(randi); randi = 3.0_r8*randi
         ord(i)=FLOOR(randi)
         m=(ord(i)+1)*(ord(i)+2)*(ord(i)+3)/6
         CALL RANDOM_NUMBER(mpole(i,1:m)); mpole(i,1:m)=-1.0_r8+2.0_r8*mpole(i,1:m) 
      END DO
     
      CLOSE(82)

      
!      CALL RANDOM_NUMBER(mpole); mpole = -1.0_r8 + 2.0_r8*mpole

      OPEN(unit=85,file=sampout,status='unknown',action='write',&
           position='append')

      CALL DATE_AND_TIME(datec,timec,zonec,time1)

      CALL DIRECT_ENG_FRC(x,y,z,ord,mpole,numpars,mpdim,denergy,&
                      fxx,fyy,fzz,dpeng)

      CALL DATE_AND_TIME(datec,timec,zonec,time2)
      CALL TTIME(time1,time2,timedirect)


      WRITE(85,13)timedirect,dpeng
      DO j=1,numpars
         WRITE(85,14)j,denergy(j),fxx(j),fyy(j),fzz(j)
      END DO

      CLOSE(unit=85)

      DEALLOCATE(x,y,z,ord,mpole,denergy,fxx,fyy,fzz)
 13   FORMAT(2X,E24.16,2X,E24.16)
 14   FORMAT(I12,4(2X,E24.16))
      STOP
      END PROGRAM EXACT_SUM 

!!!!!!!!!!!!!!


      SUBROUTINE DIRECT_ENG_FRC(x,y,z,ord,mpole,numpars,mpdim,denergy,&
                            fxx,fyy,fzz,dpeng) 

      IMPLICIT NONE

      INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(12)
      INTEGER,INTENT(IN) :: numpars,mpdim
      REAL(KIND=r8),DIMENSION(numpars),INTENT(IN) :: x,y,z
      REAL(KIND=r8),DIMENSION(numpars,1:mpdim),INTENT(IN) :: mpole
      REAL(KIND=r8),DIMENSION(numpars),INTENT(INOUT) :: denergy,fxx,fyy,fzz
      REAL(KIND=r8),INTENT(INOUT) :: dpeng
      INTEGER,DIMENSION(numpars),INTENT(IN) :: ord 
! local variables 
   
      INTEGER :: i

      dpeng=0.0_r8; denergy=0.0_r8; fxx=0.0_r8; fyy=0.0_r8; fzz=0.0_r8

      DO i=1,numpars-1
         
         SELECT CASE (ord(i))
           CASE (0) 
             
             CALL COMPUTE_ORD_ZERO(x,y,z,ord,mpole,i,numpars,mpdim,denergy,&
                            fxx,fyy,fzz,dpeng) 

           CASE (1)
           
             CALL COMPUTE_ORD_ONE(x,y,z,ord,mpole,i,numpars,mpdim,denergy,&
                            fxx,fyy,fzz,dpeng)

           CASE (2)
          
             CALL COMPUTE_ORD_TWO(x,y,z,ord,mpole,i,numpars,mpdim,denergy,&
                            fxx,fyy,fzz,dpeng)

           CASE DEFAULT
             WRITE(*,*)'The order ',ord(i),' is out of range'
             STOP
        END SELECT           

      END DO
      
      dpeng=SUM(denergy)
      
      END SUBROUTINE DIRECT_ENG_FRC  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE COMPUTE_ORD_ZERO(x,y,z,ord,mpole,id,numpars,mpdim,denergy,&
                            fxx,fyy,fzz,dpeng) 

      IMPLICIT NONE

      INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(12)
      INTEGER,INTENT(IN) :: numpars,mpdim,id
      REAL(KIND=r8),DIMENSION(numpars),INTENT(IN) :: x,y,z
      REAL(KIND=r8),DIMENSION(numpars,1:mpdim),INTENT(IN) :: mpole
      REAL(KIND=r8),DIMENSION(numpars),INTENT(INOUT) :: denergy,fxx,fyy,fzz
      REAL(KIND=r8),INTENT(INOUT) :: dpeng
      INTEGER,DIMENSION(numpars),INTENT(IN) :: ord

! local variables 
   
      INTEGER :: i,j,mord
      REAL(KIND=r8) :: xi,yi,zi,teng,peng
      REAL(KIND=r8) :: rr,rr2,rr3,rr5,rr7,frr2,tfrr2,srr2
      REAL(KIND=r8) :: xy,xz,yz,xyz,fx,fy,fz,fix,fiy,fiz,imp1
      REAL(KIND=r8) :: jmp1,jmp2,jmp3,jmp4,jmp5,jmp6,jmp7,jmp8,jmp9,jmp10
      REAL(KIND=r8) :: b0,b1,b2,b3,b4,b5,rsq
      REAL(KIND=r8) :: dpx,dpy,dpz,dpxx,dpyy,dpzz,dpxy,dpxz,     &
                       dpyz,dpxxx,dpyyy,dpzzz,dpxxy,dpxxz,       &
                       dpxyy,dpxzz,dpxyz,dpyyz,dpyzz,            &
                       dpxzzz,dpxxxx,dpyyyy,dpzzzz,              &
                       dpxxxy,dpxxxz,dpxxyy,dpxxzz,dpxxyz,       &
                       dpxyyy,dpyzzz,dpyyzz,dpxyyz,dpxyzz,       &
                       dpyyyz,dpxxxxx,dpyyyyy,dpzzzzz,dpxxxxy,   &
                       dpxxxxz,dpyyyyz,dpxzzzz,dpyzzzz,dpxyyyy,  &
                       dpxxxyy,dpxxxzz,dpyyyzz,dpyyzzz,dpxxyyy,  &
                       dpxxzzz,dpxyyzz,dpxxyyz,dpxxyzz,dpxyyyz,  &
                       dpxyzzz,dpxxxyz,                          &
                       thrb2,thrb3,thrb4,sixb4,tenb4,fifb3,      &
                       edq,eqq,ijmp,bijmp,xx,yy,zz,xx2,yy2,      &
                       zz2,tm1,tm2,tm3,tm4,tm5,tm6,tm7,tm8,      &
                       tm9,tm10,tm11,tm12,tm13,tm14,tm15
      
      i=id
      imp1=mpole(i,1)
      xi=x(i); yi=y(i); zi=z(i)
      teng=0.0_r8; fix=0.0_r8; fiy=0.0_r8; fiz=0.0_r8

         DO j=i+1,numpars 
            mord=ord(j)
            peng=0.0_r8; fx=0.0_r8; fy=0.0_r8; fz=0.0_r8

            xx=xi-x(j)
            yy=yi-y(j)
            zz=zi-z(j)

            xx2=xx*xx; yy2=yy*yy; zz2=zz*zz; xy=xx*yy; xz=xx*zz; yz=yy*zz; xyz=xx*yy*zz
            rsq=1.0_r8/(xx2+yy2+zz2); rr=SQRT(rsq); ! frr2=5.0_r8*rsq; srr2=7.0_r8*rsq; tfrr2=35.0_r8*rsq
            rr3=rsq*rr; rr5=3.0_r8*rsq*rr3;! rr7=5.0_r8*rr5*rsq

! compute recurrence terms
 
            b0 = rr 
            b1 = b0*rsq
            b2 = 3.0_r8*b1 * rsq
            b3 = 5.0_r8*b2 * rsq
            b4 = 7.0_r8*b3 * rsq
            b5 = 9.0_r8*b4 * rsq

! charge-charge interaction

            jmp1 = mpole(j,1)

            ijmp= imp1*jmp1

   ! potential on i due to j (The L_j operator)

            teng=teng+jmp1*rr

   ! potential on j due to i (The L_i operator. In this case it's just a multiplication by the charge)
 
            denergy(j)=denergy(j)+imp1*rr
   ! forces

            bijmp = ijmp*b1
            fx  = bijmp*xx
            fy  = bijmp*yy
            fz  = bijmp*zz

! charge-dipole interactions

            IF (mord .GE. 1) THEN

               jmp2=mpole(j,2); jmp3=mpole(j,3); jmp4=mpole(j,4)


   ! potential on i due to j

               teng=teng+rr3*(jmp2*xx+jmp3*yy+jmp4*zz)

   ! forces

               tm1 = -imp1*jmp2
               tm2 = -imp1*jmp3
               tm3 = -imp1*jmp4

               dpxx= xx2*b2-b1
               dpyy= yy2*b2-b1
               dpzz= zz2*b2-b1
               dpxy= xx*yy*b2
               dpxz= xx*zz*b2
               dpyz= yy*zz*b2


               fx = fx - (dpxx*tm1 + dpxy*tm2 + dpxz*tm3)
               fy = fy - (dpxy*tm1 + dpyy*tm2 + dpyz*tm3)
               fz = fz - (dpxz*tm1 + dpyz*tm2 + dpzz*tm3)

            END IF

! charge-quadrupole 

            IF (mord .EQ. 2) THEN

               jmp5=mpole(j,5); jmp6=mpole(j,6); jmp7=mpole(j,7)
               jmp8=mpole(j,8); jmp9=mpole(j,9); jmp10=mpole(j,10)

   ! potential on i due to j

               teng=teng-rr3*(jmp5+jmp8+jmp10)+&
               rr5*(xx*(jmp5*xx+jmp6*yy+jmp7*zz)+yy*(jmp8*yy+jmp9*zz)+zz*jmp10*zz)

               dpxxx = xx*(3.0_r8*b2-xx2*b3)
               dpyyy = yy*(3.0_r8*b2-yy2*b3)
               dpzzz = zz*(3.0_r8*b2-zz2*b3)
               dpxxy = yy*(b2-xx2*b3)
               dpxxz = zz*(b2-xx2*b3)
               dpxyy = xx*(b2-yy2*b3)
               dpyyz = zz*(b2-yy2*b3)
               dpxzz = xx*(b2-zz2*b3)
               dpyzz = yy*(b2-zz2*b3)
               dpxyz = -xx*yy*zz*b3

               tm1 = imp1*jmp5; tm2=imp1*jmp6
               tm3 = imp1*jmp7; tm4=imp1*jmp8
               tm5 = imp1*jmp9; tm6=imp1*jmp10


              fx = fx - (tm1*dpxxx+tm2*dpxxy+tm3*dpxxz+tm4*dpxyy+tm5*dpxyz+tm6*dpxzz)
              fy = fy - (tm1*dpxxy+tm2*dpxyy+tm3*dpxyz+tm4*dpyyy+tm5*dpyyz+tm6*dpyzz)
              fz = fz - (tm1*dpxxz+tm2*dpxyz+tm3*dpxzz+tm4*dpyyz+tm5*dpyzz+tm6*dpzzz)

            END IF

            fix = fix + fx
            fiy = fiy + fy
            fiz = fiz + fz
 
            fxx(j) = fxx(j) - fx
            fyy(j) = fyy(j) - fy
            fzz(j) = fzz(j) - fz

         END DO
 
        denergy(i)=denergy(i)+teng
        fxx(i) = fxx(i)+fix
        fyy(i) = fyy(i)+fiy
        fzz(i) = fzz(i)+fiz

        END SUBROUTINE COMPUTE_ORD_ZERO
!!!!!!!!!!!!!!!!!!
      SUBROUTINE COMPUTE_ORD_ONE(x,y,z,ord,mpole,id,numpars,mpdim,denergy,&
                            fxx,fyy,fzz,dpeng) 

      IMPLICIT NONE

      INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(12)
      INTEGER,INTENT(IN) :: numpars,mpdim,id
      REAL(KIND=r8),DIMENSION(numpars),INTENT(IN) :: x,y,z
      REAL(KIND=r8),DIMENSION(numpars,1:mpdim),INTENT(IN) :: mpole
      REAL(KIND=r8),DIMENSION(numpars),INTENT(INOUT) :: denergy,fxx,fyy,fzz
      REAL(KIND=r8),INTENT(INOUT) :: dpeng
      INTEGER,DIMENSION(numpars),INTENT(IN) :: ord

! local variables 
   
      INTEGER :: i,j,mord
      REAL(KIND=r8) :: xi,yi,zi,teng,peng
      REAL(KIND=r8) :: rr,rr2,rr3,rr5,rr7,frr2,tfrr2,srr2
      REAL(KIND=r8) :: xy,xz,yz,xyz,fx,fy,fz,fix,fiy,fiz
      REAL(KIND=r8) :: imp1,imp2,imp3,imp4
      REAL(KIND=r8) :: jmp1,jmp2,jmp3,jmp4,jmp5,jmp6,jmp7,jmp8,jmp9,jmp10
      REAL(KIND=r8) :: b0,b1,b2,b3,b4,b5,rsq
      REAL(KIND=r8) :: dpx,dpy,dpz,dpxx,dpyy,dpzz,dpxy,dpxz,     &
                       dpyz,dpxxx,dpyyy,dpzzz,dpxxy,dpxxz,       &
                       dpxyy,dpxzz,dpxyz,dpyyz,dpyzz,            &
                       dpxzzz,dpxxxx,dpyyyy,dpzzzz,              &
                       dpxxxy,dpxxxz,dpxxyy,dpxxzz,dpxxyz,       &
                       dpxyyy,dpyzzz,dpyyzz,dpxyyz,dpxyzz,       &
                       dpyyyz,dpxxxxx,dpyyyyy,dpzzzzz,dpxxxxy,   &
                       dpxxxxz,dpyyyyz,dpxzzzz,dpyzzzz,dpxyyyy,  &
                       dpxxxyy,dpxxxzz,dpyyyzz,dpyyzzz,dpxxyyy,  &
                       dpxxzzz,dpxyyzz,dpxxyyz,dpxxyzz,dpxyyyz,  &
                       dpxyzzz,dpxxxyz,                          &
                       thrb2,thrb3,thrb4,sixb4,tenb4,fifb3,      &
                       edq,eqq,ijmp,bijmp,xx,yy,zz,xx2,yy2,      &
                       zz2,tm1,tm2,tm3,tm4,tm5,tm6,tm7,tm8,      &
                       tm9,tm10,tm11,tm12,tm13,tm14,tm15

      i=id
      xi=x(i); yi=y(i); zi=z(i)
      imp1=mpole(i,1); imp2=mpole(i,2); imp3=mpole(i,3); imp4=mpole(i,4)
      teng=0.0_r8; fix=0.0_r8; fiy=0.0_r8; fiz=0.0_r8

         DO j=i+1,numpars 
            mord=ord(j)
            peng=0.0_r8; fx=0.0_r8; fy=0.0_r8; fz=0.0_r8

            xx=xi-x(j)
            yy=yi-y(j)
            zz=zi-z(j)

            xx2=xx*xx; yy2=yy*yy; zz2=zz*zz; xy=xx*yy; xz=xx*zz; yz=yy*zz; xyz=xx*yy*zz
            rsq=1.0_r8/(xx2+yy2+zz2); rr=SQRT(rsq); ! frr2=5.0_r8*rsq; srr2=7.0_r8*rsq; tfrr2=35.0_r8*rsq
            rr3=rsq*rr; rr5=3.0_r8*rsq*rr3;! rr7=5.0_r8*rr5*rsq

! compute recurrence terms
 
            b0 = rr 
            b1 = b0*rsq
            b2 = 3.0_r8*b1 * rsq
            b3 = 5.0_r8*b2 * rsq
            b4 = 7.0_r8*b3 * rsq
            b5 = 9.0_r8*b4 * rsq

! charge-charge interaction
     
            jmp1 = mpole(j,1)

   ! potential on i due to j (L_j operator)

            teng = teng + jmp1*rr

   ! potential on j due to i (L_i operator)

            denergy(j) = denergy(j)+imp1*rr

   ! forces

            ijmp= imp1*jmp1

            bijmp = ijmp*b1
            fx  = bijmp*xx
            fy  = bijmp*yy
            fz  = bijmp*zz

! charge-dipole interactions and dipole-dipole interactions

   ! charge(j)-dipole(i) interaction
   
            tm1 = imp2*jmp1
            tm2 = imp3*jmp1
            tm3 = imp4*jmp1

   ! potential on j due to i (L_i operator - dipole part)

            denergy(j) = denergy(j) - rr3*(imp2*xx+imp3*yy+imp4*zz)

            IF (mord .GE. 1) THEN

               jmp2=mpole(j,2); jmp3=mpole(j,3); jmp4=mpole(j,4)

   ! potential on i due to j (L_j operator - dipole part)

               teng = teng + rr3*(jmp2*xx+jmp3*yy+jmp4*zz)


   ! charge(i)-dipole(j) interaction

               tm1 = tm1-imp1*jmp2
               tm2 = tm2-imp1*jmp3
               tm3 = tm3-imp1*jmp4
           
            END IF

               dpxx= xx2*b2-b1
               dpyy= yy2*b2-b1
               dpzz= zz2*b2-b1
               dpxy= xx*yy*b2
               dpxz= xx*zz*b2
               dpyz= yy*zz*b2
 

               fx = fx - (dpxx*tm1 + dpxy*tm2 + dpxz*tm3)
               fy = fy - (dpxy*tm1 + dpyy*tm2 + dpyz*tm3)
               fz = fz - (dpxz*tm1 + dpyz*tm2 + dpzz*tm3)

   ! forces (dipole-dipole)

            IF (mord .GE. 1) THEN

               dpxxx = xx*(3.0_r8*b2-xx2*b3)
               dpyyy = yy*(3.0_r8*b2-yy2*b3)
               dpzzz = zz*(3.0_r8*b2-zz2*b3)
               dpxxy = yy*(b2-xx2*b3)
               dpxxz = zz*(b2-xx2*b3)
               dpxyy = xx*(b2-yy2*b3)
               dpyyz = zz*(b2-yy2*b3)
               dpxzz = xx*(b2-zz2*b3)
               dpyzz = yy*(b2-zz2*b3)
               dpxyz = -xx*yy*zz*b3

               tm1 = -imp2*jmp2
               tm2 = -imp2*jmp3-imp3*jmp2
               tm3 = -imp2*jmp4-imp4*jmp2
               tm4 = -imp3*jmp3
               tm5 = -imp3*jmp4-imp4*jmp3
               tm6 = -imp4*jmp4

              fx = fx - (tm1*dpxxx+tm2*dpxxy+tm3*dpxxz+tm4*dpxyy+tm5*dpxyz+tm6*dpxzz)
              fy = fy - (tm1*dpxxy+tm2*dpxyy+tm3*dpxyz+tm4*dpyyy+tm5*dpyyz+tm6*dpyzz)
              fz = fz - (tm1*dpxxz+tm2*dpxyz+tm3*dpxzz+tm4*dpyyz+tm5*dpyzz+tm6*dpzzz)

            END IF

! charge-quadrupole and dipole-quadrupole interactions

            IF (mord .EQ. 2) THEN

              jmp5=mpole(j,5); jmp6=mpole(j,6); jmp7=mpole(j,7)
              jmp8=mpole(j,8); jmp9=mpole(j,9); jmp10=mpole(j,10)

   ! potential on i due to j (L_j operator)

              teng = teng - rr3*(jmp5+jmp8+jmp10)+&
              rr5*(xx*(jmp5*xx+jmp6*yy+jmp7*zz)+yy*(jmp8*yy+jmp9*zz)+zz*jmp10*zz)

   ! forces (charge-quadrupole)

              dpxxx = xx*(3.0_r8*b2-xx2*b3)
              dpyyy = yy*(3.0_r8*b2-yy2*b3)
              dpzzz = zz*(3.0_r8*b2-zz2*b3)
              dpxxy = yy*(b2-xx2*b3)
              dpxxz = zz*(b2-xx2*b3)
              dpxyy = xx*(b2-yy2*b3)
              dpyyz = zz*(b2-yy2*b3)
              dpxzz = xx*(b2-zz2*b3)
              dpyzz = yy*(b2-zz2*b3)
              dpxyz = -xx*yy*zz*b3

              tm1 = imp1*jmp5; tm2=imp1*jmp6
              tm3 = imp1*jmp7; tm4=imp1*jmp8
              tm5 = imp1*jmp9; tm6=imp1*jmp10


              fx = fx - (tm1*dpxxx+tm2*dpxxy+tm3*dpxxz+tm4*dpxyy+tm5*dpxyz+tm6*dpxzz)
              fy = fy - (tm1*dpxxy+tm2*dpxyy+tm3*dpxyz+tm4*dpyyy+tm5*dpyyz+tm6*dpyzz)
              fz = fz - (tm1*dpxxz+tm2*dpxyz+tm3*dpxzz+tm4*dpyyz+tm5*dpyzz+tm6*dpzzz)

! dipole-quadrupole interactions

              thrb2  = 3.0_r8*b2; thrb3 = 3.0_r8*b3

              dpxxxx = thrb2 - xx2*(6.0_r8*b3 - xx2*b4)
              dpyyyy = thrb2 - yy2*(6.0_r8*b3 - yy2*b4)
              dpzzzz = thrb2 - zz2*(6.0_r8*b3 - zz2*b4)
              dpxxxy = yy*xx*(xx2*b4-thrb3)
              dpxxxz = zz*xx*(xx2*b4-thrb3)
              dpxyyy = xx*yy*(yy2*b4-thrb3)
              dpyyyz = yy*zz*(yy2*b4-thrb3)
              dpyzzz = yy*zz*(zz2*b4-thrb3)
              dpxzzz = xx*zz*(zz2*b4-thrb3)
              dpxxyy = b2-(xx2+yy2)*b3+xx2*yy2*b4
              dpxxzz = b2-(xx2+zz2)*b3+xx2*zz2*b4
              dpyyzz = b2-(yy2+zz2)*b3+yy2*zz2*b4
              dpxxyz = yy*zz*(xx2*b4-b3)
              dpxyyz = xx*zz*(yy2*b4-b3)
              dpxyzz = xx*yy*(zz2*b4-b3)

              tm1 = imp2*jmp5
              tm2 = imp2*jmp6+imp3*jmp5
              tm3 = imp2*jmp7+imp4*jmp5
              tm4 = imp2*jmp8+imp3*jmp6
              tm5 = imp2*jmp9+imp3*jmp7+imp4*jmp6
              tm6 = imp2*jmp10+imp4*jmp7
              tm7 = imp3*jmp8
              tm8 = imp3*jmp9+imp4*jmp8
              tm9 = imp3*jmp10+imp4*jmp9
              tm10= imp4*jmp10

              fx = fx - (tm1*dpxxxx+tm2*dpxxxy+tm3*dpxxxz+tm4*dpxxyy+tm5*dpxxyz + &
                         tm6*dpxxzz+tm7*dpxyyy+tm8*dpxyyz+tm9*dpxyzz+tm10*dpxzzz)

              fy = fy - (tm1*dpxxxy+tm2*dpxxyy+tm3*dpxxyz+tm4*dpxyyy+tm5*dpxyyz + &
                         tm6*dpxyzz+tm7*dpyyyy+tm8*dpyyyz+tm9*dpyyzz+tm10*dpyzzz)

              fz = fz - (tm1*dpxxxz+tm2*dpxxyz+tm3*dpxxzz+tm4*dpxyyz+tm5*dpxyzz + &
                         tm6*dpxzzz+tm7*dpyyyz+tm8*dpyyzz+tm9*dpyzzz+tm10*dpzzzz)

            END IF

            fix = fix + fx
            fiy = fiy + fy
            fiz = fiz + fz
 
            fxx(j) = fxx(j) - fx
            fyy(j) = fyy(j) - fy
            fzz(j) = fzz(j) - fz

         END DO
 
         denergy(i)=denergy(i)+teng
         fxx(i) = fxx(i)+fix
         fyy(i) = fyy(i)+fiy
         fzz(i) = fzz(i)+fiz

         END SUBROUTINE COMPUTE_ORD_ONE
!!!!!!!!!!!!!!!!!!
      SUBROUTINE COMPUTE_ORD_TWO(x,y,z,ord,mpole,id,numpars,mpdim,denergy,&
                            fxx,fyy,fzz,dpeng) 

      IMPLICIT NONE

      INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(12)
      INTEGER,INTENT(IN) :: numpars,mpdim,id
      REAL(KIND=r8),DIMENSION(numpars),INTENT(IN) :: x,y,z
      REAL(KIND=r8),DIMENSION(numpars,1:mpdim),INTENT(IN) :: mpole
      REAL(KIND=r8),DIMENSION(numpars),INTENT(INOUT) :: denergy,fxx,fyy,fzz
      REAL(KIND=r8),INTENT(INOUT) :: dpeng
      INTEGER,DIMENSION(numpars),INTENT(IN) :: ord

! local variables 
   
      INTEGER :: i,j,mord
      REAL(KIND=r8) :: xi,yi,zi,teng,peng
      REAL(KIND=r8) :: rr,rr2,rr3,rr5,rr7,frr2,tfrr2,srr2
      REAL(KIND=r8) :: xy,xz,yz,xyz,fx,fy,fz,fix,fiy,fiz
      REAL(KIND=r8) :: imp1,imp2,imp3,imp4,imp5,imp6,imp7,imp8,imp9,imp10
      REAL(KIND=r8) :: jmp1,jmp2,jmp3,jmp4,jmp5,jmp6,jmp7,jmp8,jmp9,jmp10
      REAL(KIND=r8) :: b0,b1,b2,b3,b4,b5,rsq
      REAL(KIND=r8) :: dpx,dpy,dpz,dpxx,dpyy,dpzz,dpxy,dpxz,     &
                       dpyz,dpxxx,dpyyy,dpzzz,dpxxy,dpxxz,       &
                       dpxyy,dpxzz,dpxyz,dpyyz,dpyzz,            &
                       dpxzzz,dpxxxx,dpyyyy,dpzzzz,              &
                       dpxxxy,dpxxxz,dpxxyy,dpxxzz,dpxxyz,       &
                       dpxyyy,dpyzzz,dpyyzz,dpxyyz,dpxyzz,       &
                       dpyyyz,dpxxxxx,dpyyyyy,dpzzzzz,dpxxxxy,   &
                       dpxxxxz,dpyyyyz,dpxzzzz,dpyzzzz,dpxyyyy,  &
                       dpxxxyy,dpxxxzz,dpyyyzz,dpyyzzz,dpxxyyy,  &
                       dpxxzzz,dpxyyzz,dpxxyyz,dpxxyzz,dpxyyyz,  &
                       dpxyzzz,dpxxxyz,                          &
                       thrb2,thrb3,thrb4,sixb4,tenb4,fifb3,      &
                       edq,eqq,ijmp,bijmp,xx,yy,zz,xx2,yy2,      &
                       zz2,tm1,tm2,tm3,tm4,tm5,tm6,tm7,tm8,      &
                       tm9,tm10,tm11,tm12,tm13,tm14,tm15

      i=id
      xi=x(i); yi=y(i); zi=z(i)
      imp1=mpole(i,1); imp2=mpole(i,2); imp3=mpole(i,3); imp4=mpole(i,4) 
      imp5=mpole(i,5); imp6=mpole(i,6); imp7=mpole(i,7); imp8=mpole(i,8)
      imp9=mpole(i,9); imp10=mpole(i,10)

      teng=0.0_r8; fix=0.0_r8; fiy=0.0_r8; fiz=0.0_r8
         DO j=i+1,numpars 
            mord=ord(j)
            peng=0.0_r8; fx=0.0_r8; fy=0.0_r8; fz=0.0_r8

            xx=xi-x(j)
            yy=yi-y(j)
            zz=zi-z(j)

            xx2=xx*xx; yy2=yy*yy; zz2=zz*zz; xy=xx*yy; xz=xx*zz; yz=yy*zz; xyz=xx*yy*zz
            rsq=1.0_r8/(xx2+yy2+zz2); rr=SQRT(rsq); ! frr2=5.0_r8*rsq; srr2=7.0_r8*rsq; tfrr2=35.0_r8*rsq
            rr3=rsq*rr; rr5=3.0_r8*rsq*rr3;! rr7=5.0_r8*rr5*rsq

! compute recurrence terms
 
            b0 = rr 
            b1 = b0*rsq
            b2 = 3.0_r8*b1 * rsq
            b3 = 5.0_r8*b2 * rsq
            b4 = 7.0_r8*b3 * rsq
            b5 = 9.0_r8*b4 * rsq

! THIS SUBROUTINE IS A LITTLE DIFFERENT FROM THE OTHER TWO (*_ZERO and *_ONE)
! HERE THE PERSPECTIVE IS ON THE SOURCE, ATOM J, NOT ON ATOM I

! charge-charge interaction
            
            jmp1 = mpole(j,1)

   ! potential on i due to j (L_j operator)

            teng = teng + jmp1*rr

   ! potential on j due to i (L_operator)

            denergy(j) = denergy(j)+imp1*rr-rr3*(imp2*xx+imp3*yy+imp4*zz)-rr3*(imp5+imp8+imp10)+&
                 rr5*(xx*(imp5*xx+imp6*yy+imp7*zz)+yy*(imp8*yy+imp9*zz)+zz*imp10*zz) 

   ! forces

            ijmp= imp1*jmp1

            bijmp = ijmp*b1
            fx  = bijmp*xx
            fy  = bijmp*yy
            fz  = bijmp*zz

! charge - dipole interactions

            tm1 = imp2*jmp1
            tm2 = imp3*jmp1
            tm3 = imp4*jmp1

            IF (mord .GE. 1) THEN

               jmp2=mpole(j,2); jmp3=mpole(j,3); jmp4=mpole(j,4)

               teng = teng +rr3*(jmp2*xx+jmp3*yy+jmp4*zz)

               tm1 = tm1-imp1*jmp2
               tm2 = tm2-imp1*jmp3
               tm3 = tm3-imp1*jmp4

            END IF

            dpx = -xx*b1; dpy = -yy*b1; dpz = -zz*b1
            dpxx= xx2*b2-b1
            dpyy= yy2*b2-b1
            dpzz= zz2*b2-b1
            dpxy= xx*yy*b2
            dpxz= xx*zz*b2
            dpyz= yy*zz*b2


            fx = fx - (dpxx*tm1 + dpxy*tm2 + dpxz*tm3)
            fy = fy - (dpxy*tm1 + dpyy*tm2 + dpyz*tm3)
            fz = fz - (dpxz*tm1 + dpyz*tm2 + dpzz*tm3)

! charge-quadrupole and dipole-dipole interactions

              dpxxx = xx*(3.0_r8*b2-xx2*b3)
              dpyyy = yy*(3.0_r8*b2-yy2*b3)
              dpzzz = zz*(3.0_r8*b2-zz2*b3)
              dpxxy = yy*(b2-xx2*b3)
              dpxxz = zz*(b2-xx2*b3)
              dpxyy = xx*(b2-yy2*b3)
              dpyyz = zz*(b2-yy2*b3)
              dpxzz = xx*(b2-zz2*b3)
              dpyzz = yy*(b2-zz2*b3)
              dpxyz = -xx*yy*zz*b3


              ! charge(j)-quadrupole(i)
              tm1 = imp5*jmp1; tm2 = imp6*jmp1; tm3 = imp7*jmp1
              tm4 = imp8*jmp1; tm5 = imp9*jmp1; tm6 = imp10*jmp1

              ! add dipole-dipole if present
              IF (mord .GE. 1) THEN

                 tm1 = tm1-imp2*jmp2
                 tm2 = tm2-imp2*jmp3-imp3*jmp2
                 tm3 = tm3-imp2*jmp4-imp4*jmp2
                 tm4 = tm4-imp3*jmp3
                 tm5 = tm5-imp3*jmp4-imp4*jmp3
                 tm6 = tm6-imp4*jmp4

              END IF

              ! now add in charge(i)-quadrupole(j) if present
              IF (mord .EQ. 2) THEN

                jmp5=mpole(j,5); jmp6=mpole(j,6); jmp7=mpole(j,7)
                jmp8=mpole(j,8); jmp9=mpole(j,9); jmp10=mpole(j,10)

                tm1 = tm1+imp1*jmp5; tm2=tm2+imp1*jmp6
                tm3 = tm3+imp1*jmp7; tm4=tm4+imp1*jmp8
                tm5 = tm5+imp1*jmp9; tm6=tm6+imp1*jmp10
        
                teng = teng -rr3*(jmp5+jmp8+jmp10)+&
                 rr5*(xx*(jmp5*xx+jmp6*yy+jmp7*zz)+yy*(jmp8*yy+jmp9*zz)+zz*jmp10*zz)

              END IF

              fx = fx - (tm1*dpxxx+tm2*dpxxy+tm3*dpxxz+tm4*dpxyy+tm5*dpxyz+tm6*dpxzz)
              fy = fy - (tm1*dpxxy+tm2*dpxyy+tm3*dpxyz+tm4*dpyyy+tm5*dpyyz+tm6*dpyzz)
              fz = fz - (tm1*dpxxz+tm2*dpxyz+tm3*dpxzz+tm4*dpyyz+tm5*dpyzz+tm6*dpzzz)

! dipole-quadrupole interactions

              IF (mord .GE. 1) THEN

                 thrb2  = 3.0_r8*b2; thrb3 = 3.0_r8*b3

                 dpxxxx = thrb2 - xx2*(6.0_r8*b3 - xx2*b4)
                 dpyyyy = thrb2 - yy2*(6.0_r8*b3 - yy2*b4)
                 dpzzzz = thrb2 - zz2*(6.0_r8*b3 - zz2*b4)
                 dpxxxy = yy*xx*(xx2*b4-thrb3)
                 dpxxxz = zz*xx*(xx2*b4-thrb3)
                 dpxyyy = xx*yy*(yy2*b4-thrb3)
                 dpyyyz = yy*zz*(yy2*b4-thrb3)
                 dpyzzz = yy*zz*(zz2*b4-thrb3)
                 dpxzzz = xx*zz*(zz2*b4-thrb3)
                 dpxxyy = b2-(xx2+yy2)*b3+xx2*yy2*b4
                 dpxxzz = b2-(xx2+zz2)*b3+xx2*zz2*b4
                 dpyyzz = b2-(yy2+zz2)*b3+yy2*zz2*b4
                 dpxxyz = yy*zz*(xx2*b4-b3)
                 dpxyyz = xx*zz*(yy2*b4-b3)
                 dpxyzz = xx*yy*(zz2*b4-b3)

                 tm1 = -imp5*jmp2
                 tm2 = -imp6*jmp2-imp5*jmp3
                 tm3 = -imp7*jmp2-imp5*jmp4
                 tm4 = -imp8*jmp2-imp6*jmp3
                 tm5 = -imp9*jmp2-imp7*jmp3-imp6*jmp4
                 tm6 = -imp10*jmp2-imp7*jmp4
                 tm7 = -imp8*jmp3
                 tm8 = -imp9*jmp3-imp8*jmp4
                 tm9 = -imp10*jmp3-imp9*jmp4
                 tm10= -imp10*jmp4

                 IF (mord .EQ. 2) THEN

                    tm1 = tm1+imp2*jmp5
                    tm2 = tm2+imp2*jmp6+imp3*jmp5
                    tm3 = tm3+imp2*jmp7+imp4*jmp5
                    tm4 = tm4+imp2*jmp8+imp3*jmp6
                    tm5 = tm5+imp2*jmp9+imp3*jmp7+imp4*jmp6
                    tm6 = tm6+imp2*jmp10+imp4*jmp7
                    tm7 = tm7+imp3*jmp8
                    tm8 = tm8+imp3*jmp9+imp4*jmp8
                    tm9 = tm9+imp3*jmp10+imp4*jmp9
                    tm10= tm10+imp4*jmp10

                 END IF
                 fx = fx - (tm1*dpxxxx+tm2*dpxxxy+tm3*dpxxxz+tm4*dpxxyy+tm5*dpxxyz + &
                            tm6*dpxxzz+tm7*dpxyyy+tm8*dpxyyz+tm9*dpxyzz+tm10*dpxzzz)

                 fy = fy - (tm1*dpxxxy+tm2*dpxxyy+tm3*dpxxyz+tm4*dpxyyy+tm5*dpxyyz + &
                         tm6*dpxyzz+tm7*dpyyyy+tm8*dpyyyz+tm9*dpyyzz+tm10*dpyzzz)

                 fz = fz - (tm1*dpxxxz+tm2*dpxxyz+tm3*dpxxzz+tm4*dpxyyz+tm5*dpxyzz + &
                            tm6*dpxzzz+tm7*dpyyyz+tm8*dpyyzz+tm9*dpyzzz+tm10*dpzzzz)

              END IF   
! quadrupole-quadrupole interactions
 
              IF (mord .EQ. 2) THEN

                 fifb3 = 15.0_r8*b3; thrb4 = 3.0_r8*b4; sixb4 = 6.0_r8*b4; tenb4 = 10.0_r8*b4
                 xyz   = xx*yy*zz

                 dpxxxxx = xx*(xx2*(tenb4-xx2*b5)-fifb3)
                 dpyyyyy = yy*(yy2*(tenb4-yy2*b5)-fifb3)
                 dpzzzzz = zz*(zz2*(tenb4-zz2*b5)-fifb3)
                 dpxxxxy = yy*(xx2*(sixb4-xx2*b5)-thrb3)
                 dpxxxxz = zz*(xx2*(sixb4-xx2*b5)-thrb3)
                 dpyyyyz = zz*(yy2*(sixb4-yy2*b5)-thrb3)
                 dpyzzzz = yy*(zz2*(sixb4-zz2*b5)-thrb3)
                 dpxzzzz = xx*(zz2*(sixb4-zz2*b5)-thrb3)
                 dpxyyyy = xx*(yy2*(sixb4-yy2*b5)-thrb3)
                 dpxxxyy = xx*((xx2+3.0_r8*yy2)*b4-xx2*yy2*b5-thrb3)
                 dpxxxzz = xx*((xx2+3.0_r8*zz2)*b4-xx2*zz2*b5-thrb3)
                 dpyyyzz = yy*((yy2+3.0_r8*zz2)*b4-yy2*zz2*b5-thrb3)
                 dpyyzzz = zz*((3.0_r8*yy2+zz2)*b4-yy2*zz2*b5-thrb3)
                 dpxxyyy = yy*((3.0_r8*xx2+yy2)*b4-xx2*yy2*b5-thrb3)
                 dpxxzzz = zz*((3.0_r8*xx2+zz2)*b4-xx2*zz2*b5-thrb3)
                 dpxyyzz = xx*((yy2+zz2)*b4-yy2*zz2*b5-b3)
                 dpxxyyz = zz*((xx2+yy2)*b4-xx2*yy2*b5-b3)
                 dpxxyzz = yy*((xx2+zz2)*b4-xx2*zz2*b5-b3)
                 dpxyyyz = xyz*(thrb4-yy2*b5)
                 dpxyzzz = xyz*(thrb4-zz2*b5)
                 dpxxxyz = xyz*(thrb4-xx2*b5)

                 tm1 = imp5*jmp5
                 tm2 = imp5*jmp6+imp6*jmp5
                 tm3 = imp5*jmp7+imp7*jmp5
                 tm4 = imp5*jmp8+imp8*jmp5+imp6*jmp6
                 tm5 = imp5*jmp9+imp9*jmp5+imp6*jmp7+imp7*jmp6
                 tm6 = imp5*jmp10+imp10*jmp5+imp7*jmp7
                 tm7 = imp6*jmp8+imp8*jmp6
                 tm8 = imp6*jmp9+imp9*jmp6+imp7*jmp8+imp8*jmp7
                 tm9 = imp6*jmp10+imp10*jmp6+imp7*jmp9+imp9*jmp7
                 tm10= imp7*jmp10+imp10*jmp7
                 tm11= imp8*jmp8
                 tm12= imp8*jmp9+imp9*jmp8
                 tm13= imp8*jmp10+imp10*jmp8+imp9*jmp9
                 tm14= imp9*jmp10+imp10*jmp9
                 tm15= imp10*jmp10

                 fx = fx - (tm1*dpxxxxx+tm2*dpxxxxy+tm3*dpxxxxz+tm4*dpxxxyy+tm5*dpxxxyz+tm6*dpxxxzz +    &
                            tm7*dpxxyyy+tm8*dpxxyyz+tm9*dpxxyzz+tm10*dpxxzzz+tm11*dpxyyyy+tm12*dpxyyyz + &
                            tm13*dpxyyzz+tm14*dpxyzzz+tm15*dpxzzzz)

                 fy = fy - (tm1*dpxxxxy+tm2*dpxxxyy+tm3*dpxxxyz+tm4*dpxxyyy+tm5*dpxxyyz+tm6*dpxxyzz +    &
                            tm7*dpxyyyy+tm8*dpxyyyz+tm9*dpxyyzz+tm10*dpxyzzz+tm11*dpyyyyy+tm12*dpyyyyz + &
                            tm13*dpyyyzz+tm14*dpyyzzz+tm15*dpyzzzz)

                 fz = fz - (tm1*dpxxxxz+tm2*dpxxxyz+tm3*dpxxxzz+tm4*dpxxyyz+tm5*dpxxyzz+tm6*dpxxzzz +    &
                            tm7*dpxyyyz+tm8*dpxyyzz+tm9*dpxyzzz+tm10*dpxzzzz+tm11*dpyyyyz+tm12*dpyyyzz + &
                            tm13*dpyyzzz+tm14*dpyzzzz+tm15*dpzzzzz)

              END IF

              fix = fix + fx
              fiy = fiy + fy
              fiz = fiz + fz
 
              fxx(j) = fxx(j) - fx
              fyy(j) = fyy(j) - fy
              fzz(j) = fzz(j) - fz

         END DO
 
        denergy(i)=denergy(i)+teng
        fxx(i) = fxx(i)+fix
        fyy(i) = fyy(i)+fiy
        fzz(i) = fzz(i)+fiz
 
        END SUBROUTINE COMPUTE_ORD_TWO

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
   
