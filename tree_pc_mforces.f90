!
!   Author:  Henry A. Boateng  (boateng@umich.edu)
!   Department of Mathematics
!   University of Michigan, Ann Arbor
!
!   Copyright (c) 2013. The Regents of the University of Michigan.
!   All Rights Reserved.
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 3 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, <http://www.gnu.org/licenses/>.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      MODULE treecode_procedures
      IMPLICIT NONE

! r8 is 8-byte (double precision) real

      INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(12)

! global variables for taylor expansions

      INTEGER :: torder,torderlim
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: cf,cf1,cf2,rkf
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:,:,:) :: b1

! global variables for multipoles of target particle

      REAL(KIND=r8) :: imp1,imp2,imp3,imp4,imp5, &
                       imp6,imp7,imp8,imp9,imp10

! global variables to track tree levels

      INTEGER :: minlevel,maxlevel

! global variables used when computing potential

      INTEGER :: orderoffset
      REAL(KIND=r8),DIMENSION(3) :: tarpos
      REAL(KIND=r8) :: thetasq

! global variables for postition and charge storage

      INTEGER,ALLOCATABLE,DIMENSION(:)  :: orderarr

! node pointer and node type declarations

      TYPE tnode_pointer
           TYPE(tnode), POINTER :: p_to_tnode
      END TYPE tnode_pointer
      TYPE tnode
           INTEGER          :: numpar,ibeg,iend
           REAL(KIND=r8)    :: x_min,y_min,z_min
           REAL(KIND=r8)    :: x_max,y_max,z_max
           REAL(KIND=r8)    :: x_mid,y_mid,z_mid
           REAL(KIND=r8)    :: radius,sqradius,aspect
           INTEGER          :: level,num_children,exist_ms
           REAL(KIND=r8),DIMENSION(:,:,:,:),POINTER :: ms
           TYPE(tnode_pointer), DIMENSION(8) :: child
      END TYPE tnode

      CONTAINS
!!!!!!!!!!!!!!!
      SUBROUTINE SETUP(x,y,z,numpars,mporder,mpdim,order,theta,xyzminmax)
      IMPLICIT NONE
!
! SETUP allocates and initializes arrays needed for the Taylor expansion.
! Also, global variables are set and the Cartesian coordinates of
! the smallest box containing the particles is determined.
!
      INTEGER,INTENT(IN) :: numpars,mporder,mpdim,order
      REAL(KIND=r8),DIMENSION(numpars),INTENT(IN) :: x,y,z
      REAL(KIND=r8),INTENT(INOUT),DIMENSION(6) :: xyzminmax
      REAL(KIND=r8),INTENT(IN) :: theta

! local variables

      INTEGER :: err,i
      REAL(KIND=r8) :: t1

! global integers and reals:  TORDER, TORDERLIM and THETASQ

      torder=order
      orderoffset=1
      torderlim=2*mporder+torder+orderoffset
      thetasq=theta*theta

      WRITE(*,*)"torder, torderlim", torder,torderlim 
! allocate global Taylor expansion variables

      ALLOCATE(cf1(torderlim),cf2(torderlim),rkf(0:torder+1),&
              b1(-2:torderlim,-2:torderlim,-2:torderlim), STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error allocating Taylor variables! '
         STOP
      END IF

! initialize arrays for Taylor sums and coeffs

      DO i=1,torderlim
         t1=1.0_r8/REAL(i,KIND=r8)
         cf1(i)=2.0_r8-t1
         cf2(i)=t1-1.0_r8
      END DO

      DO i=0,torder+1
         rkf(i)=1.0_r8/exp(FACTORIAL(i))
      END DO

! find bounds of Cartesian box enclosing the particles

      xyzminmax(1)=MINVAL(x(1:numpars))
      xyzminmax(2)=MAXVAL(x(1:numpars))
      xyzminmax(3)=MINVAL(y(1:numpars))
      xyzminmax(4)=MAXVAL(y(1:numpars))
      xyzminmax(5)=MINVAL(z(1:numpars))
      xyzminmax(6)=MAXVAL(z(1:numpars))

      ALLOCATE(orderarr(numpars),STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error allocating orderarr'
         STOP
      END IF

      DO i=1,numpars
         orderarr(i)=i
      END DO

      RETURN
      END SUBROUTINE SETUP
!!!!!!!!!!!!!!!!

      FUNCTION Factorial(N)

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: N
        Real(KIND=r8) :: Factorial
        INTEGER :: i

        Factorial = 0.0_r8

        do i=2,N
          Factorial = Factorial+Log(Real(i,Kind=r8))
        end do

      END FUNCTION Factorial

!!!!!!!!!!!!!
      RECURSIVE SUBROUTINE CREATE_TREE_N0(p,ibeg,iend,x,y,z,mpole,shrink, &
                                      maxparnode,xyzmm,level,arrdim,mpdim)
      IMPLICIT NONE
!
! CREATE_TREE_N0 recursively creates the tree structure. Node P is
! input, which contains particles indexed from IBEG to IEND. After
! the node parameters are set subdivision occurs if IEND-IBEG+1 > MAXPARNODE.
! Real array XYZMM contains the min and max values of the coordinates
! of the particle in P, thus defining the box. The division of a cluster terminates
! when the number of particles in a cluster are is less or equal to maxparnode
!
      TYPE(tnode),POINTER :: p
      INTEGER,INTENT(IN) :: ibeg,iend,shrink,mpdim,level,maxparnode,arrdim
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(INOUT) :: x,y,z
      REAL(KIND=r8),DIMENSION(arrdim,1:mpdim),INTENT(INOUT) :: mpole
      REAL(KIND=r8),DIMENSION(6),INTENT(IN) :: xyzmm

! local variables

      REAL(KIND=r8) :: x_mid,y_mid,z_mid,xl,yl,zl,lmax,t1,t2,t3
      INTEGER, DIMENSION(8,2) :: ind
      REAL(KIND=r8), DIMENSION(6,8) :: xyzmms
      INTEGER :: i,j,limin,limax,err,loclev,numposchild
      REAL(KIND=r8), DIMENSION(6) ::  lxyzmm

! allocate pointer

      ALLOCATE(p,STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error allocating pointer! '
         STOP
      END IF

! set node fields: number of particles, exist_ms
! and xyz bounds

      p%numpar=iend-ibeg+1
      p%exist_ms=0

      IF (shrink .EQ. 1) THEN
         p%x_min=MINVAL(x(ibeg:iend))
         p%x_max=MAXVAL(x(ibeg:iend))
         p%y_min=MINVAL(y(ibeg:iend))
         p%y_max=MAXVAL(y(ibeg:iend))
         p%z_min=MINVAL(z(ibeg:iend))
         p%z_max=MAXVAL(z(ibeg:iend))
      ELSE
         p%x_min=xyzmm(1)
         p%x_max=xyzmm(2)
         p%y_min=xyzmm(3)
         p%y_max=xyzmm(4)
         p%z_min=xyzmm(5)
         p%z_max=xyzmm(6)
      END IF

! compute aspect ratio

      xl=p%x_max-p%x_min
      yl=p%y_max-p%y_min
      zl=p%z_max-p%z_min

      lmax=MAX(xl,yl,zl)
      t1=lmax
      t2=MIN(xl,yl,zl)
      IF (t2 .NE. 0.0_r8) THEN
         p%aspect=t1/t2
      ELSE
         p%aspect=0.0_r8
      END IF

! midpoint coordinates , RADIUS and SQRADIUS

      p%x_mid=(p%x_max+p%x_min)/2.0_r8
      p%y_mid=(p%y_max+p%y_min)/2.0_r8
      p%z_mid=(p%z_max+p%z_min)/2.0_r8
      t1=p%x_max-p%x_mid
      t2=p%y_max-p%y_mid
      t3=p%z_max-p%z_mid
      p%sqradius=t1*t1+t2*t2+t3*t3
      p%radius=SQRT(p%sqradius)

! set particle limits, tree level of node, and nullify children pointers

      p%ibeg=ibeg
      p%iend=iend
      p%level=level
      IF (maxlevel .LT. level) THEN
         maxlevel=level
      END IF
      p%num_children=0
      DO i=1,8
         NULLIFY(p%child(i)%p_to_tnode)
      END DO

      IF (p%numpar .GT. maxparnode) THEN
!
! set IND array to 0 and then call PARTITION routine.  IND array holds indices
! of the eight new subregions. Also, setup XYZMMS array in case SHRINK=1
!
         xyzmms(1,1)=p%x_min
         xyzmms(2,1)=p%x_max
         xyzmms(3,1)=p%y_min
         xyzmms(4,1)=p%y_max
         xyzmms(5,1)=p%z_min
         xyzmms(6,1)=p%z_max
         ind(1,1)=ibeg
         ind(1,2)=iend
         x_mid=p%x_mid
         y_mid=p%y_mid
         z_mid=p%z_mid

         CALL PARTITION_8(x,y,z,mpole,xyzmms,xl,yl,zl,lmax,numposchild, &
                         x_mid,y_mid,z_mid,ind,arrdim,mpdim)
!
! create children if indicated and store info in parent
!
         loclev=level+1
         DO i=1,numposchild
            IF (ind(i,1) .LE. ind(i,2)) THEN
               p%num_children=p%num_children+1
               lxyzmm=xyzmms(:,i)
               CALL CREATE_TREE_N0(p%child(p%num_children)%p_to_tnode, &
                               ind(i,1),ind(i,2),x,y,z,mpole,shrink, &
                               maxparnode,lxyzmm,loclev,arrdim,mpdim)
            END IF

         END DO
      ELSE
         IF (level .LT. minlevel) THEN
            minlevel=level
         END IF
      END IF

      END SUBROUTINE CREATE_TREE_N0
!!!!!!!!!!!!!!!
      RECURSIVE SUBROUTINE CREATE_TREE_LV(p,ibeg,iend,x,y,z,mpole,shrink, &
                                          treelevel,xyzmm,level,arrdim,mpdim)
      IMPLICIT NONE
!
! CREATE_TREE_LV recursively creates the tree structure. Node P is
! input, which contains particles indexed from IBEG to IEND. After
! the node parameters are set subdivision occurs if IEND-IBEG+1 > MAXPARNODE.
! Real array XYZMM contains the min and max values of the coordinates
! of the particle in P, thus defining the box. The subdivision terminates
! when the number of levels equals treelevel
!

      TYPE(tnode),POINTER :: p
      INTEGER,INTENT(IN) :: ibeg,iend,shrink,level,treelevel,arrdim,mpdim
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(INOUT) :: x,y,z
      REAL(KIND=r8),DIMENSION(arrdim,1:mpdim),INTENT(INOUT) :: mpole
      REAL(KIND=r8),DIMENSION(6),INTENT(IN) :: xyzmm

! local variables

      REAL(KIND=r8) :: x_mid,y_mid,z_mid,xl,yl,zl,lmax,t1,t2,t3
      INTEGER, DIMENSION(8,2) :: ind
      REAL(KIND=r8), DIMENSION(6,8) :: xyzmms
      INTEGER :: i,j,limin,limax,err,loclev,numposchild
      REAL(KIND=r8), DIMENSION(6) ::  lxyzmm

! allocate pointer

      ALLOCATE(p,STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error allocating pointer! '
         STOP
      END IF
! set node fields: number of particles, exist_ms
! and xyz bounds

      p%numpar=iend-ibeg+1
      p%exist_ms=0

      IF (shrink .EQ. 1) THEN
         p%x_min=MINVAL(x(ibeg:iend))
         p%x_max=MAXVAL(x(ibeg:iend))
         p%y_min=MINVAL(y(ibeg:iend))
         p%y_max=MAXVAL(y(ibeg:iend))
         p%z_min=MINVAL(z(ibeg:iend))
         p%z_max=MAXVAL(z(ibeg:iend))
      ELSE
         p%x_min=xyzmm(1)
         p%x_max=xyzmm(2)
         p%y_min=xyzmm(3)
         p%y_max=xyzmm(4)
         p%z_min=xyzmm(5)
         p%z_max=xyzmm(6)
      END IF

! compute aspect ratio

      xl=p%x_max-p%x_min
      yl=p%y_max-p%y_min
      zl=p%z_max-p%z_min

      lmax=MAX(xl,yl,zl)
      t1=lmax
      t2=MIN(xl,yl,zl)
      IF (t2 .NE. 0.0_r8) THEN
         p%aspect=t1/t2
      ELSE
         p%aspect=0.0_r8
      END IF

! midpoint coordinates , RADIUS and SQRADIUS

      p%x_mid=(p%x_max+p%x_min)/2.0_r8
      p%y_mid=(p%y_max+p%y_min)/2.0_r8
      p%z_mid=(p%z_max+p%z_min)/2.0_r8
      t1=p%x_max-p%x_mid
      t2=p%y_max-p%y_mid
      t3=p%z_max-p%z_mid
      p%sqradius=t1*t1+t2*t2+t3*t3
      p%radius=SQRT(p%sqradius)

! set particle limits, tree level of node, and nullify children pointers

      p%ibeg=ibeg
      p%iend=iend
      p%level=level
      p%num_children=0
      DO i=1,8
         NULLIFY(p%child(i)%p_to_tnode)
      END DO

      IF(level .LT. treelevel) THEN
!
! set IND array to 0 and then call PARTITION routine.  IND array holds indices
! of the eight new subregions. Also, setup XYZMMS array in case SHRINK=1
!
         xyzmms(1,1)=p%x_min
         xyzmms(2,1)=p%x_max
         xyzmms(3,1)=p%y_min
         xyzmms(4,1)=p%y_max
         xyzmms(5,1)=p%z_min
         xyzmms(6,1)=p%z_max
         ind(1,1)=ibeg
         ind(1,2)=iend
         x_mid=p%x_mid
         y_mid=p%y_mid
         z_mid=p%z_mid

         CALL PARTITION_8(x,y,z,mpole,xyzmms,xl,yl,zl,lmax,numposchild, &
                         x_mid,y_mid,z_mid,ind,arrdim,mpdim)
!
! create children if indicated and store info in parent
!
         loclev=level+1
         DO i=1,numposchild
            IF (ind(i,1) .LE. ind(i,2)) THEN
               p%num_children=p%num_children+1
               lxyzmm=xyzmms(:,i)
               CALL CREATE_TREE_LV(p%child(p%num_children)%p_to_tnode, &
                               ind(i,1),ind(i,2),x,y,z,mpole,shrink, &
                               treelevel,lxyzmm,loclev,arrdim,mpdim)
            END IF

         END DO
      ELSE
         IF (level .LT. minlevel) THEN
            minlevel=level
         END IF
      END IF

      END SUBROUTINE CREATE_TREE_LV


!!!!!!!!!!!!!!!
      SUBROUTINE PARTITION_8(x,y,z,mpole,xyzmms,xl,yl,zl,lmax,numposchild, &
                            x_mid,y_mid,z_mid,ind,arrdim,mpdim)
      IMPLICIT NONE
!
! PARTITION_8 determines the particle indices of the eight sub boxes
! containing the particles after the box defined by particles I_BEG
! to I_END is divided by its midpoints in each coordinate direction.
! The determination of the indices is accomplished by the subroutine
! PARTITION. A box is divided in a coordinate direction as long as the
! resulting aspect ratio is not too large. This avoids the creation of
! "narrow" boxes in which Talyor expansions may become inefficient.
! On exit the INTEGER array IND (dimension 8 x 2) contains
! the indice limits of each new box (node) and NUMPOSCHILD the number
! of possible children.  If IND(J,1) > IND(J,2) for a given J this indicates
! that box J is empty.
!
      INTEGER, INTENT(IN) :: arrdim,mpdim
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(INOUT) :: x,y,z
      REAL(KIND=r8),DIMENSION(arrdim,1:mpdim),INTENT(INOUT) :: mpole
      INTEGER, DIMENSION(8,2),INTENT(INOUT) :: ind
      REAL(KIND=r8),DIMENSION(6,8),INTENT(INOUT) :: xyzmms
      REAL(KIND=r8), INTENT(IN) :: x_mid,y_mid,z_mid,xl,yl,zl,lmax
      INTEGER,INTENT(INOUT) :: numposchild

! local variables

      INTEGER :: temp_ind,i
      REAL(KIND=r8) :: critlen

      numposchild=1
      critlen=lmax/sqrt(2.0_r8)

      IF (xl .GE. critlen) THEN
         CALL PARTITION(x,y,z,mpole,orderarr,ind(1,1),ind(1,2), &
                       x_mid,temp_ind,arrdim,mpdim)

         ind(2,1)=temp_ind+1
         ind(2,2)=ind(1,2)
         ind(1,2)=temp_ind
         xyzmms(:,2)=xyzmms(:,1)
         xyzmms(2,1)=x_mid
         xyzmms(1,2)=x_mid
         numposchild=2*numposchild
      END IF

      IF (yl .GE. critlen) THEN
         DO i=1,numposchild
            CALL PARTITION(y,x,z,mpole,orderarr,ind(i,1),ind(i,2), &
                          y_mid,temp_ind,arrdim,mpdim)
            ind(numposchild+i,1)=temp_ind+1
            ind(numposchild+i,2)=ind(i,2)
            ind(i,2)=temp_ind
            xyzmms(:,numposchild+i)=xyzmms(:,i)
            xyzmms(4,i)=y_mid
            xyzmms(3,numposchild+i)=y_mid
         END DO
         numposchild=2*numposchild
      END IF

      IF (zl .GE. critlen) THEN
         DO i=1,numposchild
            CALL PARTITION(z,x,y,mpole,orderarr,ind(i,1),ind(i,2), &
                          z_mid,temp_ind,arrdim,mpdim)
            ind(numposchild+i,1)=temp_ind+1
            ind(numposchild+i,2)=ind(i,2)
            ind(i,2)=temp_ind
            xyzmms(:,numposchild+i)=xyzmms(:,i)
            xyzmms(6,i)=z_mid
            xyzmms(5,numposchild+i)=z_mid
         END DO
         numposchild=2*numposchild
      END IF

      RETURN
      END SUBROUTINE PARTITION_8
!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE PC_TREECODE(p,x,y,z,mpole,tpeng,EnP,tfx,tfy,tfz,numpars,mpdim)
      IMPLICIT NONE
!
! PC_TREECODE is the driver routine which calls COMPUTE_PC for each
! target particle, setting the global variable TARPOS before the call.
! P is the root node of the source tree, xT,yT,zT are the coordinates of
! the target particles. X,Y,Z,Q are the coordinates and charges of the
! source particles. The energy at target 'i', is stored in EnP(i).
!

      INTEGER,INTENT(IN) :: numpars,mpdim
      TYPE(tnode),POINTER :: p
      REAL(KIND=r8),DIMENSION(numpars),INTENT(INOUT) :: x
      REAL(KIND=r8),DIMENSION(numpars),INTENT(IN) :: y,z
      REAL(KIND=r8),DIMENSION(numpars,1:mpdim),INTENT(INOUT) :: mpole
      REAL(KIND=r8),DIMENSION(numpars),INTENT(INOUT) :: EnP,tfx,tfy,tfz
      REAL(KIND=r8),INTENT(INOUT) :: tpeng

! local variables

      INTEGER :: i,j,nn
      REAL(KIND=r8) :: penglocal,peng,xtemp,mptemp(1:mpdim)
      REAL(KIND=r8) :: fx,fy,fz,flcx,flcy,flcz


      tpeng=0.0_r8
      DO i=1,numpars
         peng=0.0_r8; fx=0.0_r8; fy=0.0_r8; fz=0.0_r8
         tarpos(1)=x(i)
         tarpos(2)=y(i)
         tarpos(3)=z(i)
         xtemp=x(i)
         mptemp=mpole(i,:)
         imp1=mptemp(1); imp2=mptemp(2); imp3=mptemp(3); imp4=mptemp(4)
         imp5=mptemp(5); imp6=mptemp(6); imp7=mptemp(7); imp8=mptemp(8)
         imp9=mptemp(9); imp10=mptemp(10)

         x(i)=x(i)+100.0_r8
         mpole(i,:)=0.0_r8
         DO j=1,p%num_children
            CALL COMPUTE_PC(p%child(j)%p_to_tnode,penglocal, &
                           flcx,flcy,flcz,x,y,z,mpole,numpars,mpdim)
            peng=peng+penglocal
            fx=fx+flcx
            fy=fy+flcy
            fz=fz+flcz
         END DO
         x(i)=xtemp
         mpole(i,:)=mptemp
         nn=orderarr(i)
         EnP(nn)=peng
         tfx(nn)=fx
         tfy(nn)=fy
         tfz(nn)=fz
      END DO
      tpeng=SUM(EnP)

      RETURN
      END SUBROUTINE PC_TREECODE
!!!!!!!!!!!!!!
      RECURSIVE SUBROUTINE COMPUTE_PC(p,peng,fx,fy,fz,x,y,z,mpole,arrdim,mpdim)

      IMPLICIT NONE

! COMPUTE_PC is the recursive routine for computing the interaction
! between a target particle and a source cluster. If the MAC is
! satisfied the interaction is approximated using a multipole
! approximation. If the MAC is not satisfied then the algorithm
! descends to the children of the current cluster, unless the
! current cluster is a leaf then the interaction is done exactly
! via a call to the routine COMP_DIRECT
!

      INTEGER,INTENT(IN) :: arrdim,mpdim
      TYPE(tnode),POINTER :: p
      REAL(KIND=r8),INTENT(INOUT) :: peng,fx,fy,fz
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(IN) :: x,y,z
      REAL(KIND=r8),DIMENSION(arrdim,1:mpdim),INTENT(IN) :: mpole

! local variables

      REAL(KIND=r8) :: flcx,flcy,flcz
      REAL(KIND=r8) :: tx,ty,tz,distsq,t1,t2,penglocal,k2f,k3f,fac
      REAL(KIND=r8) :: m1,m2,m3,m4,m5,m6,m7,m8,m9,m10
      INTEGER :: i,j,k,i1,i2,i3,i4,i5,j1,j2,j3,j4,j5,k1,k2,k3,k4,k5,err

! determine DISTSQ for MAC test

      tx=tarpos(1)-p%x_mid
      ty=tarpos(2)-p%y_mid
      tz=tarpos(3)-p%z_mid
      distsq=tx*tx+ty*ty+tz*tz

! intialize potential energy and force

      peng=0.0_r8
      fx=0.0_r8; fy=0.0_r8; fz=0.0_r8

! If MAC is accepted and there is more than 1 particle in the
! box use the expansion for the approximation.

      IF ((p%sqradius .LT. distsq*thetasq) .AND. &
         (p%sqradius .NE. 0.0_r8)) THEN


         IF (p%exist_ms .EQ. 0) THEN
             ALLOCATE(p%ms(1:mpdim,0:torder,0:torder,0:torder),STAT=err)
             IF (err .NE. 0) THEN
                WRITE(6,*) 'Error allocating node moments! '
                STOP
             END IF
             CALL COMP_MS(p,x,y,z,mpole,arrdim,mpdim)
             p%exist_ms=1
         END IF
         penglocal=0.0_r8
         flcx=0.0_r8;flcy=0.0_r8;flcz=0.0_r8
         b1=0.0_r8
         CALL COMP_TCOEFF(p)
         DO k=0,torder
            k1=k+1; k2=k+2; k3=k+3; k4=k+4; k5=k+5
            k3f=rkf(k)
            DO j=0,torder-k
               j1=j+1; j2=j+2; j3=j+3; j4=j+4; j5=j+5
               k2f=k3f*rkf(j)
               DO i=0,torder-k-j
                  i1=i+1; i2=i+2; i3=i+3; i4=i+4; i5=i+5
                  fac=k2f*rkf(i)

                  m1=p%ms(1,i,j,k);m2=p%ms(2,i,j,k);m3=p%ms(3,i,j,k);m4=p%ms(4,i,j,k)
                  m5=p%ms(5,i,j,k);m6=p%ms(6,i,j,k);m7=p%ms(7,i,j,k);m8=p%ms(8,i,j,k)
                  m9=p%ms(9,i,j,k);m10=p%ms(10,i,j,k)

                  peng=peng+fac*( b1(i,j,k)*m1+b1(i1,j,k)*m2+b1(i,j1,k)*m3+b1(i,j,k1)*m4+&
                                  b1(i2,j,k)*m5+b1(i1,j1,k)*m6+b1(i1,j,k1)*m7+b1(i,j2,k)*m8+&
                                  b1(i,j1,k1)*m9+b1(i,j,k2)*m10 )

! charge on target interacting with source
! (0,0,0)

                  flcx=flcx+imp1*fac*( b1(i1,j,k)*m1+b1(i2,j,k)*m2+b1(i1,j1,k)*m3+b1(i1,j,k1)*m4+&
                                          b1(i3,j,k)*m5+b1(i2,j1,k)*m6+b1(i2,j,k1)*m7+b1(i1,j2,k)*m8+&
                                          b1(i1,j1,k1)*m9+b1(i1,j,k2)*m10 )
   
                  flcy=flcy+imp1*fac*( b1(i,j1,k)*m1+b1(i1,j1,k)*m2+b1(i,j2,k)*m3+b1(i,j1,k1)*m4+&
                                  b1(i2,j1,k)*m5+b1(i1,j2,k)*m6+b1(i1,j1,k1)*m7+b1(i,j3,k)*m8+&
                                  b1(i,j2,k1)*m9+b1(i,j1,k2)*m10 )

                  flcz=flcz+imp1*fac*( b1(i,j,k1)*m1+b1(i1,j,k1)*m2+b1(i,j1,k1)*m3+b1(i,j,k2)*m4+&
                                  b1(i2,j,k1)*m5+b1(i1,j1,k1)*m6+b1(i1,j,k2)*m7+b1(i,j2,k1)*m8+&
                                  b1(i,j1,k2)*m9+b1(i,j,k3)*m10 )


! dipole on target interacting with source
! (1,0,0)

                  flcx=flcx-imp2*fac*( b1(i2,j,k)*m1+b1(i3,j,k)*m2+b1(i2,j1,k)*m3+b1(i2,j,k1)*m4+&
                                  b1(i4,j,k)*m5+b1(i3,j1,k)*m6+b1(i3,j,k1)*m7+b1(i2,j2,k)*m8+&
                                  b1(i2,j1,k1)*m9+b1(i2,j,k2)*m10 )

                  flcy=flcy-imp2*fac*( b1(i1,j1,k)*m1+b1(i2,j1,k)*m2+b1(i1,j2,k)*m3+b1(i1,j1,k1)*m4+&
                                  b1(i3,j1,k)*m5+b1(i2,j2,k)*m6+b1(i2,j1,k1)*m7+b1(i1,j3,k)*m8+&
                                  b1(i1,j2,k1)*m9+b1(i1,j1,k2)*m10 )

                  flcz=flcz-imp2*fac*( b1(i1,j,k1)*m1+b1(i2,j,k1)*m2+b1(i1,j1,k1)*m3+b1(i1,j,k2)*m4+&
                                  b1(i3,j,k1)*m5+b1(i2,j1,k1)*m6+b1(i2,j,k2)*m7+b1(i1,j2,k1)*m8+&
                                  b1(i1,j1,k2)*m9+b1(i1,j,k3)*m10 )

!(0,1,0)
                  flcx=flcx-imp3*fac*( b1(i1,j1,k)*m1+b1(i2,j1,k)*m2+b1(i1,j2,k)*m3+b1(i1,j1,k1)*m4+&
                                  b1(i3,j1,k)*m5+b1(i2,j2,k)*m6+b1(i2,j1,k1)*m7+b1(i1,j3,k)*m8+&
                                  b1(i1,j2,k1)*m9+b1(i1,j1,k2)*m10 )

                  flcy=flcy-imp3*fac*( b1(i,j2,k)*m1+b1(i1,j2,k)*m2+b1(i,j3,k)*m3+b1(i,j2,k1)*m4+&
                                  b1(i2,j2,k)*m5+b1(i1,j3,k)*m6+b1(i1,j2,k1)*m7+b1(i,j4,k)*m8+&
                                  b1(i,j3,k1)*m9+b1(i,j2,k2)*m10 )

                  flcz=flcz-imp3*fac*( b1(i,j1,k1)*m1+b1(i1,j1,k1)*m2+b1(i,j2,k1)*m3+b1(i,j1,k2)*m4+&
                                  b1(i2,j1,k1)*m5+b1(i1,j2,k1)*m6+b1(i1,j1,k2)*m7+b1(i,j3,k1)*m8+&
                                  b1(i,j2,k2)*m9+b1(i,j1,k3)*m10 )

!(0,0,1)

                  flcx=flcx-imp4*fac*( b1(i1,j,k1)*m1+b1(i2,j,k1)*m2+b1(i1,j1,k1)*m3+b1(i1,j,k2)*m4+&
                                  b1(i3,j,k1)*m5+b1(i2,j1,k1)*m6+b1(i2,j,k2)*m7+b1(i1,j2,k1)*m8+&
                                  b1(i1,j1,k2)*m9+b1(i1,j,k3)*m10 )

                  flcy=flcy-imp4*fac*( b1(i,j1,k1)*m1+b1(i1,j1,k1)*m2+b1(i,j2,k1)*m3+b1(i,j1,k2)*m4+&
                                  b1(i2,j1,k1)*m5+b1(i1,j2,k1)*m6+b1(i1,j1,k2)*m7+b1(i,j3,k1)*m8+&
                                  b1(i,j2,k2)*m9+b1(i,j1,k3)*m10 )

                  flcz=flcz-imp4*fac*( b1(i,j,k2)*m1+b1(i1,j,k2)*m2+b1(i,j1,k2)*m3+b1(i,j,k3)*m4+&
                                  b1(i2,j,k2)*m5+b1(i1,j1,k2)*m6+b1(i1,j,k3)*m7+b1(i,j2,k2)*m8+&
                                  b1(i,j1,k3)*m9+b1(i,j,k4)*m10 )

! quadrupole on target interacting with source
! (2,0,0)

                  flcx=flcx+imp5*fac*( b1(i3,j,k)*m1+b1(i4,j,k)*m2+b1(i3,j1,k)*m3+b1(i3,j,k1)*m4+&
                                  b1(i5,j,k)*m5+b1(i4,j1,k)*m6+b1(i4,j,k1)*m7+b1(i3,j2,k)*m8+&
                                  b1(i3,j1,k1)*m9+b1(i3,j,k2)*m10 )

                  flcy=flcy+imp5*fac*( b1(i2,j1,k)*m1+b1(i3,j1,k)*m2+b1(i2,j2,k)*m3+b1(i2,j1,k1)*m4+&
                                  b1(i4,j1,k)*m5+b1(i3,j2,k)*m6+b1(i3,j1,k1)*m7+b1(i2,j3,k)*m8+&
                                  b1(i2,j2,k1)*m9+b1(i2,j1,k2)*m10 )

                  flcz=flcz+imp5*fac*( b1(i2,j,k1)*m1+b1(i3,j,k1)*m2+b1(i2,j1,k1)*m3+b1(i2,j,k2)*m4+&
                                  b1(i4,j,k1)*m5+b1(i3,j1,k1)*m6+b1(i3,j,k2)*m7+b1(i2,j2,k1)*m8+&
                                  b1(i2,j1,k2)*m9+b1(i2,j,k3)*m10 )

! (1,1,0)

                  flcx=flcx+imp6*fac*( b1(i2,j1,k)*m1+b1(i3,j1,k)*m2+b1(i2,j2,k)*m3+b1(i2,j1,k1)*m4+&
                                  b1(i4,j1,k)*m5+b1(i3,j2,k)*m6+b1(i3,j1,k1)*m7+b1(i2,j3,k)*m8+&
                                  b1(i2,j2,k1)*m9+b1(i2,j1,k2)*m10 )

                  flcy=flcy+imp6*fac*( b1(i1,j2,k)*m1+b1(i2,j2,k)*m2+b1(i1,j3,k)*m3+b1(i1,j2,k1)*m4+&
                                  b1(i3,j2,k)*m5+b1(i2,j3,k)*m6+b1(i2,j2,k1)*m7+b1(i1,j4,k)*m8+&
                                  b1(i1,j3,k1)*m9+b1(i1,j2,k2)*m10 )

                  flcz=flcz+imp6*fac*( b1(i1,j1,k1)*m1+b1(i2,j1,k1)*m2+b1(i1,j2,k1)*m3+b1(i1,j1,k2)*m4+&
                                  b1(i3,j1,k1)*m5+b1(i2,j2,k1)*m6+b1(i2,j1,k2)*m7+b1(i1,j3,k1)*m8+&
                                  b1(i1,j2,k2)*m9+b1(i1,j1,k3)*m10 )

! (1,0,1)

                  flcx=flcx+imp7*fac*( b1(i2,j,k1)*m1+b1(i3,j,k1)*m2+b1(i2,j1,k1)*m3+b1(i2,j,k2)*m4+&
                                  b1(i4,j,k1)*m5+b1(i3,j1,k1)*m6+b1(i3,j,k2)*m7+b1(i2,j2,k1)*m8+&
                                  b1(i2,j1,k2)*m9+b1(i2,j,k3)*m10 )

                  flcy=flcy+imp7*fac*( b1(i1,j1,k1)*m1+b1(i2,j1,k1)*m2+b1(i1,j2,k1)*m3+b1(i1,j1,k2)*m4+&
                                  b1(i3,j1,k1)*m5+b1(i2,j2,k1)*m6+b1(i2,j1,k2)*m7+b1(i1,j3,k1)*m8+&
                                  b1(i1,j2,k2)*m9+b1(i1,j1,k3)*m10 )

                  flcz=flcz+imp7*fac*( b1(i1,j,k2)*m1+b1(i2,j,k2)*m2+b1(i1,j1,k2)*m3+b1(i1,j,k3)*m4+&
                                  b1(i3,j,k2)*m5+b1(i2,j1,k2)*m6+b1(i2,j,k3)*m7+b1(i1,j2,k2)*m8+&
                                  b1(i1,j1,k3)*m9+b1(i1,j,k4)*m10 )

! (0,2,0)

                  flcx=flcx+imp8*fac*( b1(i1,j2,k)*m1+b1(i2,j2,k)*m2+b1(i1,j3,k)*m3+b1(i1,j2,k1)*m4+&
                                  b1(i3,j2,k)*m5+b1(i2,j3,k)*m6+b1(i2,j2,k1)*m7+b1(i1,j4,k)*m8+&
                                  b1(i1,j3,k1)*m9+b1(i1,j2,k2)*m10 )

                  flcy=flcy+imp8*fac*( b1(i,j3,k)*m1+b1(i1,j3,k)*m2+b1(i,j4,k)*m3+b1(i,j3,k1)*m4+&
                                  b1(i2,j3,k)*m5+b1(i1,j4,k)*m6+b1(i1,j3,k1)*m7+b1(i,j5,k)*m8+&
                                  b1(i,j4,k1)*m9+b1(i,j3,k2)*m10 )

                  flcz=flcz+imp8*fac*( b1(i,j2,k1)*m1+b1(i1,j2,k1)*m2+b1(i,j3,k1)*m3+b1(i,j2,k2)*m4+&
                                  b1(i2,j2,k1)*m5+b1(i1,j3,k1)*m6+b1(i1,j2,k2)*m7+b1(i,j4,k1)*m8+&
                                  b1(i,j3,k2)*m9+b1(i,j2,k3)*m10 )

! (0,1,1)

                  flcx=flcx+imp9*fac*( b1(i1,j1,k1)*m1+b1(i2,j1,k1)*m2+b1(i1,j2,k1)*m3+b1(i1,j1,k2)*m4+&
                                  b1(i3,j1,k1)*m5+b1(i2,j2,k1)*m6+b1(i2,j1,k2)*m7+b1(i1,j3,k1)*m8+&
                                  b1(i1,j2,k2)*m9+b1(i1,j1,k3)*m10 )

                  flcy=flcy+imp9*fac*( b1(i,j2,k1)*m1+b1(i1,j2,k1)*m2+b1(i,j3,k1)*m3+b1(i,j2,k2)*m4+&
                                  b1(i2,j2,k1)*m5+b1(i1,j3,k1)*m6+b1(i1,j2,k2)*m7+b1(i,j4,k1)*m8+&
                                  b1(i,j3,k2)*m9+b1(i,j2,k3)*m10 )

                  flcz=flcz+imp9*fac*( b1(i,j1,k2)*m1+b1(i1,j1,k2)*m2+b1(i,j2,k2)*m3+b1(i,j1,k3)*m4+&
                                  b1(i2,j1,k2)*m5+b1(i1,j2,k2)*m6+b1(i1,j1,k3)*m7+b1(i,j3,k2)*m8+&
                                  b1(i,j2,k3)*m9+b1(i,j1,k4)*m10 )

! (0,0,2)

                   flcx=flcx+imp10*fac*( b1(i1,j,k2)*m1+b1(i2,j,k2)*m2+b1(i1,j1,k2)*m3+b1(i1,j,k3)*m4+&
                                  b1(i3,j,k2)*m5+b1(i2,j1,k2)*m6+b1(i2,j,k3)*m7+b1(i1,j2,k2)*m8+&
                                  b1(i1,j1,k3)*m9+b1(i1,j,k4)*m10 )

                  flcy=flcy+imp10*fac*( b1(i,j1,k2)*m1+b1(i1,j1,k2)*m2+b1(i,j2,k2)*m3+b1(i,j1,k3)*m4+&
                                  b1(i2,j1,k2)*m5+b1(i1,j2,k2)*m6+b1(i1,j1,k3)*m7+b1(i,j3,k2)*m8+&
                                  b1(i,j2,k3)*m9+b1(i,j1,k4)*m10 )

                  flcz=flcz+imp10*fac*( b1(i,j,k3)*m1+b1(i1,j,k3)*m2+b1(i,j1,k3)*m3+b1(i,j,k4)*m4+&
                                  b1(i2,j,k3)*m5+b1(i1,j1,k3)*m6+b1(i1,j,k4)*m7+b1(i,j2,k3)*m8+&
                                  b1(i,j1,k4)*m9+b1(i,j,k5)*m10 )                                                                                                                                                                         
               END DO
            END DO
         END DO
         fx = fx+flcx; fy=fy+flcy; fz=fz+flcz
        
      ELSE

! If MAC fails check to see if the are children. If not, perform direct
! calculation.  If there are children, call routine recursively for each.
!
         IF (p%num_children .EQ. 0) THEN
            CALL COMP_DIRECT(penglocal,p%ibeg,p%iend, &
                              flcx,flcy,flcz,x,y,z,mpole,arrdim,mpdim)
            peng=penglocal
            fx=flcx; fy=flcy; fz=flcz
         ELSE
            DO i=1,p%num_children
               CALL COMPUTE_PC(p%child(i)%p_to_tnode,penglocal, &
                              flcx,flcy,flcz,x,y,z,mpole,arrdim,mpdim)
               peng=peng+penglocal
               fx=fx+flcx; fy=fy+flcy; fz=fz+flcz
            END DO
         END IF
      END IF

      RETURN
      END SUBROUTINE COMPUTE_PC
!!!!!!!!
      SUBROUTINE COMP_TCOEFF(p)
      IMPLICIT NONE
!
! COMP_TCOEFF computes the Taylor coefficients of the potential
! using a recurrence formula.  The center of the expansion is the
! midpoint of the node P.  TARPOS and TORDERLIM are globally defined.
!
      TYPE(tnode),POINTER :: p

! local variables

      REAL(KIND=r8) :: tx,ty,tz,tdx,tdy,tdz,fac,sqfac,t1,rk2,rk221,rk3,rk331
      INTEGER :: i,j,k,n,tp1,mm,i1,i2,j1,j2,k1,k2,k3,k21,k22,k31,k32

! setup variables

      tx=tarpos(1)-p%x_mid
      ty=tarpos(2)-p%y_mid
      tz=tarpos(3)-p%z_mid

      fac=1.0_r8/(tx*tx+ty*ty+tz*tz)
      sqfac=SQRT(fac)

! k1=0, k2=0, k3=0

      b1(0,0,0)=sqfac

!      k2=0, k3=0

      do k1=1,torderlim

         b1(k1,0,0)=fac*(cf1(k1) * Real(k1,Kind=r8)*tx*b1(k1-1,0,0)+&
                         cf2(k1) * Real(k1*(k1-1),Kind=r8)*b1(k1-2,0,0))
      end do

!      k3=0
      do k2=1,torderlim
         k21=k2-1; k22=k2-2;rk2=Real(k2,Kind=r8)*ty;rk221=Real(k2*k21,Kind=r8)
         do k1=0,torderlim-k2
            n=k1+k2
            b1(k1,k2,0)=fac*(cf1(n)*(Real(k1,Kind=r8)*tx*b1(k1-1,k2,0)+&
                                                      rk2*b1(k1,k21,0))+&
                             cf2(n)*(Real(k1*(k1-1),Kind=r8)*b1(k1-2,k2,0)+&
                                                     rk221*b1(k1,k22,0)))
         end do
      end do

      do k3=1,torderlim
         k31=k3-1; k32=k3-2;rk3=Real(k3,Kind=r8)*tz;rk331=Real(k3*k31,Kind=r8)
         do k2=0,torderlim-k3
            k21=k2-1; k22=k2-2;rk2=Real(k2,Kind=r8)*ty;rk221=Real(k2*k21,Kind=r8)
            do k1=0,torderlim-k3-k2
               n=k1+k2+k3
               b1(k1,k2,k3)=fac*(cf1(n)*(Real(k1,Kind=r8)*tx*b1(k1-1,k2,k3)+&
                                                       rk2*b1(k1,k21,k3)+&
                                                       rk3*b1(k1,k2,k31))+&
                                 cf2(n)*(Real(k1*(k1-1),Kind=r8)*b1(k1-2,k2,k3)+&
                                                            rk221*b1(k1,k22,k3)+&
                                                            rk331*b1(k1,k2,k32) ))
            end do
         end do
       end do

      RETURN
      END SUBROUTINE COMP_TCOEFF
!!!!!!!!!!!!!!!
      SUBROUTINE COMP_MS(p,x,y,z,mpole,arrdim,mpdim)
      IMPLICIT NONE
!
! COMP_MS computes the moments for node P needed in the Taylor approximation
!
      INTEGER,INTENT(IN) :: arrdim,mpdim
      TYPE(tnode),POINTER :: p
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(IN) :: x,y,z
      REAL(KIND=r8),DIMENSION(arrdim,1:mpdim),INTENT(IN) :: mpole

! local variables

      INTEGER :: i,j,k1,k2,k3
      REAL(KIND=r8) :: dx,dy,dz,tx,ty,tz,txyz,jmp(1:mpdim)
      REAL(KIND=r8) :: m1,m2,m3,m4,m5,m6,m7,m8,m9,m10

      p%ms=0.0_r8
      DO j=p%ibeg,p%iend
         dx=x(j)-p%x_mid
         dy=y(j)-p%y_mid
         dz=z(j)-p%z_mid
         jmp=mpole(j,:)
         m1=jmp(1); m2=jmp(2); m3=jmp(3); m4=jmp(4); m5=jmp(5)
         m6=jmp(6); m7=jmp(7); m8=jmp(8); m9=jmp(9); m10=jmp(10)

         tz=1.0_r8
         DO k3=0,torder
            ty=tz
            DO k2=0,torder-k3
               tx=ty
               DO k1=0,torder-k3-k2

                  p%ms(1,k1,k2,k3)=p%ms(1,k1,k2,k3)+m1*tx
                  p%ms(2,k1,k2,k3)=p%ms(2,k1,k2,k3)+m2*tx
                  p%ms(3,k1,k2,k3)=p%ms(3,k1,k2,k3)+m3*tx
                  p%ms(4,k1,k2,k3)=p%ms(4,k1,k2,k3)+m4*tx
                  p%ms(5,k1,k2,k3)=p%ms(5,k1,k2,k3)+m5*tx
                  p%ms(6,k1,k2,k3)=p%ms(6,k1,k2,k3)+m6*tx
                  p%ms(7,k1,k2,k3)=p%ms(7,k1,k2,k3)+m7*tx
                  p%ms(8,k1,k2,k3)=p%ms(8,k1,k2,k3)+m8*tx
                  p%ms(9,k1,k2,k3)=p%ms(9,k1,k2,k3)+m9*tx
                  p%ms(10,k1,k2,k3)=p%ms(10,k1,k2,k3)+m10*tx

                  tx=tx*dx
              END DO
              ty=ty*dy
            END DO
            tz=tz*dz
         END DO

      END DO

      RETURN
      END SUBROUTINE COMP_MS
!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE COMP_DIRECT(peng,ibeg,iend,fx,fy,fz,x,y,z,mpole,arrdim,mpdim)
      IMPLICIT NONE
!
! COMP_DIRECT directly computes the potential on the current target
! particle (determined by the global variable TARPOS) due to the
! current source cluster

      INTEGER,INTENT(IN) :: ibeg,iend,arrdim,mpdim
      REAL(KIND=r8),INTENT(OUT) :: peng,fx,fy,fz
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(IN) :: x,y,z
      REAL(KIND=r8),DIMENSION(arrdim,1:mpdim),INTENT(IN) :: mpole

! local variables

      INTEGER ::j
      REAL(KIND=r8) :: tx2,ty2,tz2,txy,txz,tyz,txyz
      REAL(KIND=r8) :: etmp,rr,rsq,rr2,rr3,rr5,rr7,frr2,tfrr2,srr2
      REAL(KIND=r8) :: jmp1,jmp2,jmp3,jmp4,jmp5,jmp6,jmp7,jmp8,jmp9,jmp10
      REAL(KIND=r8) :: b0,b1,b2,b3,b4,b5,xyz
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

      peng=0.0_r8
      fx=0.0_r8; fy=0.0_r8; fz=0.0_r8
      DO j=ibeg,iend
         xx=tarpos(1)-x(j)
         yy=tarpos(2)-y(j)
         zz=tarpos(3)-z(j)

         jmp1=mpole(j,1);jmp2=mpole(j,2);jmp3=mpole(j,3);jmp4=mpole(j,4)
         jmp5=mpole(j,5);jmp6=mpole(j,6);jmp7=mpole(j,7);jmp8=mpole(j,8)
         jmp9=mpole(j,9);jmp10=mpole(j,10)

         xx2=xx*xx; yy2=yy*yy; zz2=zz*zz
         rsq=1.0_r8/(xx2+yy2+zz2); rr=SQRT(rsq)
         rr3=rsq*rr; rr5=3.0_r8*rsq*rr3

         peng=peng+jmp1*rr+rr3*(jmp2*xx+jmp3*yy+jmp4*zz)-rr3*(jmp5+jmp8+jmp10)+&
              rr5*(xx*(jmp5*xx+jmp6*yy+jmp7*zz)+yy*(jmp8*yy+jmp9*zz)+zz*jmp10*zz)

! compute recurrence terms
 
            b0 = rr 
            b1 = b0*rsq
            b2 = 3.0_r8*b1 * rsq
            b3 = 5.0_r8*b2 * rsq
            b4 = 7.0_r8*b3 * rsq
            b5 = 9.0_r8*b4 * rsq

!==== compute forces ============

! charge-charge interaction

            ijmp= imp1*jmp1

            bijmp = ijmp*b1
            fx  = fx+bijmp*xx
            fy  = fy+bijmp*yy
            fz  = fz+bijmp*zz

! charge-dipole interactions

            tm1 = imp2*jmp1-imp1*jmp2
            tm2 = imp3*jmp1-imp1*jmp3
            tm3 = imp4*jmp1-imp1*jmp4

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

              tm1 = imp2*jmp2
              tm2 = imp2*jmp3+imp3*jmp2
              tm3 = imp2*jmp4+imp4*jmp2
              tm4 = imp3*jmp3
              tm5 = imp3*jmp4+imp4*jmp3
              tm6 = imp4*jmp4

              tm1 = imp1*jmp5+imp5*jmp1-tm1; tm2=imp1*jmp6+imp6*jmp1-tm2
              tm3 = imp1*jmp7+imp7*jmp1-tm3; tm4=imp1*jmp8+imp8*jmp1-tm4
              tm5 = imp1*jmp9+imp9*jmp1-tm5; tm6=imp1*jmp10+imp10*jmp1-tm6


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

              tm1 = imp2*jmp5-imp5*jmp2
              tm2 = imp2*jmp6-imp6*jmp2+imp3*jmp5-imp5*jmp3
              tm3 = imp2*jmp7-imp7*jmp2+imp4*jmp5-imp5*jmp4
              tm4 = imp2*jmp8-imp8*jmp2+imp3*jmp6-imp6*jmp3
              tm5 = imp2*jmp9-imp9*jmp2+imp3*jmp7-imp7*jmp3+imp4*jmp6-imp6*jmp4
              tm6 = imp2*jmp10-imp10*jmp2+imp4*jmp7-imp7*jmp4
              tm7 = imp3*jmp8-imp8*jmp3
              tm8 = imp3*jmp9-imp9*jmp3+imp4*jmp8-imp8*jmp4
              tm9 = imp3*jmp10-imp10*jmp3+imp4*jmp9-imp9*jmp4
              tm10= imp4*jmp10-imp10*jmp4

              fx = fx - (tm1*dpxxxx+tm2*dpxxxy+tm3*dpxxxz+tm4*dpxxyy+tm5*dpxxyz + &
                         tm6*dpxxzz+tm7*dpxyyy+tm8*dpxyyz+tm9*dpxyzz+tm10*dpxzzz)

              fy = fy - (tm1*dpxxxy+tm2*dpxxyy+tm3*dpxxyz+tm4*dpxyyy+tm5*dpxyyz + &
                         tm6*dpxyzz+tm7*dpyyyy+tm8*dpyyyz+tm9*dpyyzz+tm10*dpyzzz)

              fz = fz - (tm1*dpxxxz+tm2*dpxxyz+tm3*dpxxzz+tm4*dpxyyz+tm5*dpxyzz + &
                         tm6*dpxzzz+tm7*dpyyyz+tm8*dpyyzz+tm9*dpyzzz+tm10*dpzzzz)

! quadrupole-quadrupole interactions

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


      END DO

      RETURN
      END SUBROUTINE COMP_DIRECT
!!!!!!!!!!!!!!!!!
      SUBROUTINE CLEANUP(p)
      IMPLICIT NONE
!
! CLEANUP deallocates allocated global variables and then
! calls recursive routine REMOVE_NODE to delete the tree.
!
      TYPE(tnode),POINTER :: p

! local variables

      INTEGER :: err

      DEALLOCATE(cf1,cf2,b1, STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error deallocating Taylor variables! '
         STOP
      END IF

      DEALLOCATE(orderarr,STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error deallocating orderarr variables! '
         STOP
      END IF

      CALL REMOVE_NODE(p)
      DEALLOCATE(p, STAT=err)
      IF (err .NE. 0) THEN
         WRITE(6,*) 'Error deallocating root node! '
         STOP
      END IF
      NULLIFY(p)

      RETURN
      END SUBROUTINE CLEANUP
!!!!!!!!!!!
      RECURSIVE SUBROUTINE REMOVE_NODE(p)
      IMPLICIT NONE
!
! REMOVE_NODE recursively removes each node from the
! tree and deallocates its memory for MS array if it
! exits.
!
      TYPE(tnode),POINTER :: p

! local variables

      INTEGER :: i,err

      IF (p%exist_ms .EQ. 1) THEN
         DEALLOCATE(p%ms,STAT=err)
         IF (err .NE. 0) THEN
            WRITE(6,*) 'Error deallocating node MS! '
            STOP
         END IF
      END IF

      IF (p%num_children .GT. 0) THEN
          DO i=1,p%num_children
            CALL REMOVE_NODE(p%child(i)%p_to_tnode)
            DEALLOCATE(p%child(i)%p_to_tnode,STAT=err)
            IF (err .NE. 0) THEN
               WRITE(6,*) 'Error deallocating node child! '
               STOP
            END IF
          END DO
      END IF

      RETURN
      END SUBROUTINE REMOVE_NODE

      END MODULE treecode_procedures
!!!!!!!!!!!!!!!
      SUBROUTINE TREECODE(x,y,z,mpole,numpars,mporder,mpdim, &
                      tEn,tpeng,tfx,tfy,tfz,order,theta,shrink,&
                      maxparnode,timetree,treelevel,&
                      iflag)

!====================================================================
!
! x,y,z,mpole    :: x,y,z coordinates and multipoles of particles
! numpar         :: number of particles
! mporder        :: Highest multipole order
! mpdim          :: number of independent multipoles
! tEn            :: array of dimension numpar for storing potential
!                   at each particle
! tfx, tfy, tfz  :: arrays of dimension numpar for storing x, y and
!                   z components of forces at each particle
! tpeng          :: total potential
! shrink         :: flag for determining whether to shrink a
!                   cluster to the smallest Cartesian box that
!                   encloses the particles in the cluster.
! maxparnode     :: maximum number of particles in a leaf
! timetree       :: The total time for the treecode computation
! treelevel      :: maximum number of levels
! iflag          :: if iflag=0, the division of the source tree
!                   terminates when the number of particles in a leaf
!                   less than or equal to maxparnode. If iflag is
!                   not equal to zero, then the divison terminates
!                   when the number of levels of the tree is equal
!                   to treelevel
!=====================================================================

      USE treecode_procedures
      IMPLICIT NONE

      INTEGER,INTENT(IN) :: numpars,mporder,mpdim,order,shrink,&
                            maxparnode,treelevel,iflag
      REAL(KIND=r8),DIMENSION(numpars),INTENT(INOUT) :: x,y,z
      REAL(KIND=r8),DIMENSION(numpars,1:mpdim),INTENT(INOUT) :: mpole
      REAL(KIND=r8),DIMENSION(numpars),INTENT(OUT) :: tEn,tfx,tfy,tfz
      REAL(KIND=r8),INTENT(IN) :: theta
      REAL(KIND=r8),INTENT(OUT) :: tpeng,timetree

! local variables

      TYPE(tnode),POINTER :: troot
      INTEGER :: i,level,err
      REAL(KIND=r8), DIMENSION(6) :: xyzminmax

! variables needed for f90 DATE_AND_TIME intrinsic

      INTEGER,DIMENSION(8) :: time1,time2
      CHARACTER (LEN=8)  :: datec
      CHARACTER (LEN=10) :: timec
      CHARACTER (LEN=5)  :: zonec
      REAL(KIND=r8)      :: totaltime


! Call SETUP to allocate arrays for Taylor expansions
! and setup global variables.

      CALL SETUP(x,y,z,numpars,mporder,mpdim,order,theta,xyzminmax)

! nullify pointer to root of tree (TROOT) and create tree

      NULLIFY(troot)

      CALL DATE_AND_TIME(datec,timec,zonec,time1)

! set global variables to track tree levels during construction

      level=0
      minlevel=50000

         WRITE(6,*) ' '
         WRITE(6,*) 'Creating tree '

      IF(iflag .EQ. 0) THEN
         maxlevel=0
         CALL CREATE_TREE_N0(troot,1,numpars,x,y,z,mpole,shrink, &
                      maxparnode,xyzminmax,level,numpars,mpdim)
      ELSE
         maxlevel=treelevel
         CALL CREATE_TREE_LV(troot,1,numpars,x,y,z,mpole,shrink, &
                      treelevel,xyzminmax,level,numpars,mpdim)
      END IF
      CALL DATE_AND_TIME(datec,timec,zonec,time2)
      CALL TTIME(time1,time2,totaltime)
      timetree=totaltime

! print tree information to stdout

         WRITE(6,*) ' '
         WRITE(6,*) 'Tree created '
         WRITE(6,*) ' '
         WRITE(6,*) 'Tree creation time (secs) : ',totaltime
         WRITE(6,*) ' '
         WRITE(6,*) 'Run synopsis : '
         WRITE(6,*) ' '
         WRITE(6,*) 'numpar :',troot%numpar
         WRITE(6,*) 'x_mid  :',troot%x_mid
         WRITE(6,*) 'y_mid  :',troot%y_mid
         WRITE(6,*) 'z_mid  :',troot%z_mid
         WRITE(6,*) 'radius :',troot%radius
         WRITE(6,*) 'torder          :',torder
         WRITE(6,*) 'theta           :',theta
         WRITE(6,*) 'shrink          :',shrink
         WRITE(6,*) 'maxparnode      :',maxparnode
         WRITE(6,*) 'iflag           :',iflag
         WRITE(6,*) 'tree maxlevel   :',treelevel
         WRITE(6,*) ' '


      CALL DATE_AND_TIME(datec,timec,zonec,time1)

!Call driver routine for particle cluster

      CALL PC_TREECODE(troot,x,y,z,mpole,tpeng,tEn,tfx,tfy,tfz,numpars,mpdim)

      CALL DATE_AND_TIME(datec,timec,zonec,time2)
      CALL TTIME(time1,time2,totaltime)
      timetree=timetree+totaltime

         WRITE(6,*) ' '
         WRITE(6,*) 'Treecode time (secs) : ',timetree

! Call CLEANUP to deallocate global variables and tree structure.

         WRITE(6,*) ' '
         WRITE(6,*) 'Deallocating tree structure!'
         WRITE(6,*) ' '

    !  CALL CLEANUP(troot)

      END SUBROUTINE TREECODE
!!!!!!!!!!!!!
      SUBROUTINE PARTITION(a,b,c,mpole,indarr,ibeg,iend,val,midind,arrdim,mpdim)
      IMPLICIT NONE
!
! PARTITION determines the index MIDIND, after partitioning
! in place the  arrays A,B,C and Q,  such that
! A(IBEG:MIDIND) <= VAL and  A(MIDIND+1:IEND) > VAL.
! If on entry IBEG > IEND or  A(IBEG:IEND) > VAL then MIDIND
! is returned as IBEG-1.
!
      INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(12)
      INTEGER, INTENT(IN) :: arrdim,ibeg,iend,mpdim
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(INOUT) :: a,b,c
      REAL(KIND=r8),DIMENSION(arrdim,1:mpdim),INTENT(INOUT) :: mpole
      INTEGER,DIMENSION(arrdim),INTENT(INOUT) :: indarr
      INTEGER, INTENT(INOUT) :: midind
      REAL(KIND=r8) val

! local variables

      REAL(KIND=r8) ta,tb,tc,tmpole(1:mpdim)
      INTEGER lower,upper,tind

      IF (ibeg .LT. iend) THEN

! temporarily store IBEG entries and set A(IBEG)=VAL for
! the partitoning algorithm.

         ta=a(ibeg)
         tb=b(ibeg)
         tc=c(ibeg)
         tmpole=mpole(ibeg,:)
         tind=indarr(ibeg)
         a(ibeg)=val
         upper=ibeg
         lower=iend

         DO WHILE (upper .NE. lower)
            DO WHILE ((upper .LT. lower) .AND. (val .LT. a(lower)))
                  lower=lower-1
            END DO
            IF (upper .NE. lower) THEN
               a(upper)=a(lower)
               b(upper)=b(lower)
               c(upper)=c(lower)
               mpole(upper,:)=mpole(lower,:)
               indarr(upper)=indarr(lower)
            END IF
            DO WHILE ((upper .LT. lower) .AND. (val .GE. a(upper)))
                  upper=upper+1
            END DO
            IF (upper .NE. lower) THEN
               a(lower)=a(upper)
               b(lower)=b(upper)
               c(lower)=c(upper)
               mpole(lower,:)=mpole(upper,:)
               indarr(lower)=indarr(upper)
            END IF
         END DO
         midind=upper

! replace TA in position UPPER and change MIDIND if TA > VAL

         IF (ta .GT. val) THEN
            midind=upper-1
         END IF
         a(upper)=ta
         b(upper)=tb
         c(upper)=tc
         mpole(upper,:)=tmpole
         indarr(upper)=tind

      ELSEIF (ibeg .EQ. iend) THEN
         IF (a(ibeg) .LE. val) THEN
            midind=ibeg
         ELSE
            midind=ibeg-1
         END IF
      ELSE
         midind=ibeg-1
      END IF

      RETURN
      END SUBROUTINE PARTITION
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
