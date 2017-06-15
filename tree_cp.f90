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
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: cf,cf1,cf2,rkf,jmp
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:,:,:) :: b1

! global variables to track tree levels 
 
      INTEGER :: minlevel,maxlevel
      
! global variables used when computing potential

      INTEGER :: orderoffset,targid
      REAL(KIND=r8),DIMENSION(3) :: tarpos
      REAL(KIND=r8) :: thetasq

      REAL(KIND=r8) :: m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14
      REAL(KIND=r8) :: m15,m16,m17,m18,m19,m20,m21,m22,m23,m24,m25
      REAL(KIND=r8) :: m26,m27,m28,m29,m30,m31,m32,m33,m34,m35

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
           REAL(KIND=r8),DIMENSION(:,:,:),POINTER :: ms
           TYPE(tnode_pointer), DIMENSION(8) :: child
      END TYPE tnode

      CONTAINS
!!!!!!!!!!!!!!!
      SUBROUTINE SETUP(x,y,z,numpars,mpdim,order,theta,xyzminmax) 
      IMPLICIT NONE
!
! SETUP allocates and initializes arrays needed for the Taylor expansion.
! Also, global variables are set and the Cartesian coordinates of
! the smallest box containing the particles is determined. 
!
      INTEGER,INTENT(IN) :: numpars,mpdim,order
      REAL(KIND=r8),DIMENSION(numpars),INTENT(IN) :: x,y,z
      REAL(KIND=r8),INTENT(INOUT),DIMENSION(6) :: xyzminmax
      REAL(KIND=r8),INTENT(IN) :: theta

! local variables

      INTEGER :: err,i
      REAL(KIND=r8) :: t1

! global integers and reals:  TORDER, TORDERLIM and THETASQ

      torder=order
      orderoffset=4
      torderlim=torder+orderoffset
      thetasq=theta*theta

! allocate global Taylor expansion variables

      ALLOCATE(cf1(torderlim),cf2(torderlim),rkf(0:torder+1),&
              b1(-2:torderlim,-2:torderlim,-2:torderlim), &
              jmp(1:mpdim),STAT=err)
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
      SUBROUTINE CP_TREECODE(p,x,y,z,mpole,tpeng,EnP,numpars,mpdim)
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
      REAL(KIND=r8),DIMENSION(numpars),INTENT(INOUT) :: EnP
      REAL(KIND=r8),INTENT(INOUT) :: tpeng
 
! local variables

      INTEGER :: i,j

      EnP=0.0_r8
      DO i=1,numpars
         targid=i
         tarpos(1)=x(i)
         tarpos(2)=y(i)
         tarpos(3)=z(i)
         jmp=mpole(i,:)
!=======================================================================
         m1=jmp(1); m2=jmp(2); m3=jmp(3); m4=jmp(4); m5=jmp(5)
         m6=jmp(6); m7=jmp(7); m8=jmp(8); m9=jmp(9); m10=jmp(10)
         m11=jmp(11); m12=jmp(12); m13=jmp(13); m14=jmp(14); m15=jmp(15)
         m16=jmp(16); m17=jmp(17); m18=jmp(18); m19=jmp(19); m20=jmp(20)
         m21=jmp(21); m22=jmp(22); m23=jmp(23); m24=jmp(24); m25=jmp(25)
         m26=jmp(26); m27=jmp(27); m28=jmp(28); m29=jmp(29); m30=jmp(30)
         m31=jmp(31); m32=jmp(32); m33=jmp(33); m34=jmp(34); m35=jmp(35)
!========================================================================
         DO j=1,p%num_children
            CALL COMPUTE_CP1(p%child(j)%p_to_tnode,EnP, &
                           x,y,z,numpars)
         END DO
         
      END DO

      CALL COMPUTE_CP2(p,x,y,z,EnP,numpars)

      tpeng=SUM(EnP)

      RETURN
      END SUBROUTINE CP_TREECODE
!!!!!!!!!!!!!!
      RECURSIVE SUBROUTINE COMPUTE_CP1(p,EnP,x,y,z,arrdim)
      IMPLICIT NONE


! COMPUTE_CP1 is the recursive routine for computing the interaction
! between a target particle and a source cluster. If the MAC is
! satisfied the power series coefficients for the current target 
! cluster are updated. If the MAC is not satisfied then the algorithm 
! descends to the children of the current cluster, unless the
! current cluster is a leaf then the interaction is done exactly
! via a call to the routine COMP_DIRECT
!

      INTEGER,INTENT(IN) :: arrdim
      TYPE(tnode),POINTER :: p      
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(INOUT) :: EnP
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(IN) :: x,y,z

! local variables

      REAL(KIND=r8) :: tx,ty,tz,distsq,t1,t2,penglocal
      INTEGER :: i,j,k,err

! determine DISTSQ for MAC test

      tx=p%x_mid-tarpos(1)
      ty=p%y_mid-tarpos(2)
      tz=p%z_mid-tarpos(3)
      distsq=tx*tx+ty*ty+tz*tz

! intialize potential energy and force 

! If MAC is accepted and there is more than 1 particle in the 
! box use the expansion for the approximation.
      IF ((p%sqradius .LT. distsq*thetasq) .AND. &
         (p%sqradius .NE. 0.0_r8)) THEN

         b1=0.0_r8
         CALL COMP_TCOEFF(tx,ty,tz)
         IF (p%exist_ms .EQ. 0) THEN
             ALLOCATE(p%ms(0:torder,0:torder,0:torder),STAT=err)
             IF (err .NE. 0) THEN
                WRITE(6,*) 'Error allocating node moments! '
                STOP
             END IF
             p%ms=0.0_r8
             p%exist_ms=1
         END IF
         CALL COMP_CMS(p)   
    
      ELSE

! If MAC fails check to see if the are children. If not, perform direct 
! calculation.  If there are children, call routine recursively for
! each.
!
         IF (p%num_children .EQ. 0) THEN
            CALL COMP_DIRECT(EnP,p%ibeg,p%iend, &
                              x,y,z,arrdim)
         ELSE
            DO i=1,p%num_children
               CALL COMPUTE_CP1(p%child(i)%p_to_tnode,EnP, &
                              x,y,z,arrdim)
            END DO  
         END IF 
      END IF

      RETURN
      END SUBROUTINE COMPUTE_CP1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      RECURSIVE SUBROUTINE COMPUTE_CP2(ap,x,y,z,EnP,arrdim)

        IMPLICIT NONE

! COMPUTE_CP2 is a recursive routine that evaluates the power series  
! approximation of the potential at the targets in a cluster via 
! a 3-D Horner's rule.  
!
        INTEGER,INTENT(IN) :: arrdim
        TYPE(tnode),POINTER :: ap
        REAL(KIND=r8),DIMENSION(arrdim),INTENT(INOUT) :: EnP 
        REAL(KIND=r8),DIMENSION(arrdim),INTENT(IN) :: x,y,z

! local variables

        REAL(KIND=r8) :: tx,ty,tz,peng
        REAL(KIND=r8) :: xm,ym,zm,dx,dy,dz
        INTEGER :: i,nn,j,k1,k2,k3,porder,porder1

        porder=torder
        porder1=porder-1 
        IF(ap%exist_ms==1)THEN
          xm=ap%x_mid; ym=ap%y_mid; zm=ap%z_mid

         DO i=ap%ibeg,ap%iend
           nn=orderarr(i)
           dx=x(i)-xm; dy=y(i)-ym; dz=z(i)-zm
           peng=ap%ms(0,0,porder)           

           DO k3=porder1,0,-1

              ty=ap%ms(0,porder-k3,k3)

              DO k2=porder1-k3,0,-1

                 tx=ap%ms(porder-k3-k2,k2,k3)

                 DO k1=porder1-k3-k2,0,-1
                 
                   tx  = dx * tx + ap%ms(k1,k2,k3) 
                   
                 END DO

                 ty  = dy * ty + tx 
              END DO

               peng = dz * peng + ty  

             END DO

             EnP(nn)=EnP(nn)+peng

          END DO
        END IF

         DO j=1,ap%num_children
            CALL COMPUTE_CP2(ap%child(j)%p_to_tnode,x,y,z,EnP,arrdim)
         END DO

       RETURN
      END SUBROUTINE COMPUTE_CP2 

!!!!!!!!

!!!!!!!!!!!!!!
      SUBROUTINE COMP_TCOEFF(tx,ty,tz)
      IMPLICIT NONE
!
! COMP_TCOEFF computes the Taylor coefficients of the potential
! using a recurrence formula.  The center of the expansion is the
! midpoint of the node P.  TARPOS and TORDERLIM are globally defined.  
!

      REAL(KIND=r8),INTENT(IN) :: tx,ty,tz

! local variables

      REAL(KIND=r8) :: tdx,tdy,tdz,fac,sqfac,t1,rk2,rk221,rk3,rk331
      INTEGER :: i,j,k,n,tp1,mm,i1,i2,j1,j2,k1,k2,k3,k21,k22,k31,k32

! setup variables

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
      SUBROUTINE COMP_CMS(p)
      IMPLICIT NONE
!
! COMP_CMS computes the coefficients for node P needed in the power series
!
      TYPE(tnode),POINTER :: p 

! local variables

      INTEGER :: i,i1,i2,i3,i4,j,j1,j2,j3,j4,k,k1,k2,k3,k4
      REAL(KIND=r8) :: fac,tx,ty,tz,k2f,k3f

      tz=1.0_r8       
      DO k=0,torder
         ty=tz
         k1=k+1; k2=k+2; k3=k+3; k4=k+4
         k3f=rkf(k)
         DO j=0,torder-k
            tx=ty
            j1=j+1; j2=j+2; j3=j+3; j4=j+4
            k2f=k3f*rkf(j)
            DO i=0,torder-k-j
               i1=i+1; i2=i+2; i3=i+3; i4=i+4
               fac=k2f*rkf(i)*tx
 
               p%ms(i,j,k)=p%ms(i,j,k)+fac*(m1*b1(i,j,k)+&
               m2*b1(i1,j,k)+m3*b1(i,j1,k)+m4*b1(i,j,k1)+&
               m5*b1(i2,j,k)+m6*b1(i1,j1,k)+m7*b1(i1,j,k1)+&
               m8*b1(i,j2,k)+m9*b1(i,j1,k1)+m10*b1(i,j,k2)+&
               m11*b1(i3,j,k)+m12*b1(i2,j1,k)+m13*b1(i2,j,k1)+&
               m14*b1(i1,j2,k)+m15*b1(i1,j1,k1)+m16*b1(i1,j,k2)+&
               m17*b1(i,j3,k)+m18*b1(i,j2,k1)+m19*b1(i,j1,k2)+&
               m20*b1(i,j,k3)+m21*b1(i4,j,k)+m22*b1(i3,j1,k)+&
               m23*b1(i3,j,k1)+m24*b1(i2,j2,k)+m25*b1(i2,j1,k1)+&
               m26*b1(i2,j,k2)+m27*b1(i1,j3,k)+m28*b1(i1,j2,k1)+&
               m29*b1(i1,j1,k2)+m30*b1(i1,j,k3)+m31*b1(i,j4,k)+&
               m32*b1(i,j3,k1)+m33*b1(i,j2,k2)+m34*b1(i,j1,k3)+&
               m35*b1(i,j,k4))
  
!               p%ms(1,i,j,k)=p%ms(1,i,j,k)+m1*fac*b1(i,j,k)
!               p%ms(2,i,j,k)=p%ms(2,i,j,k)+m2*fac*b1(i1,j,k)
!               p%ms(3,i,j,k)=p%ms(3,i,j,k)+m3*fac*b1(i,j1,k)
!               p%ms(4,i,j,k)=p%ms(4,i,j,k)+m4*fac*b1(i,j,k1)
!               p%ms(5,i,j,k)=p%ms(5,i,j,k)+m5*fac*b1(i2,j,k)
!               p%ms(6,i,j,k)=p%ms(6,i,j,k)+m6*fac*b1(i1,j1,k)
!               p%ms(7,i,j,k)=p%ms(7,i,j,k)+m7*fac*b1(i1,j,k1)
!               p%ms(8,i,j,k)=p%ms(8,i,j,k)+m8*fac*b1(i,j2,k)
!               p%ms(9,i,j,k)=p%ms(9,i,j,k)+m9*fac*b1(i,j1,k1)
!               p%ms(10,i,j,k)=p%ms(10,i,j,k)+m10*fac*b1(i,j,k2)
!               p%ms(11,i,j,k)=p%ms(11,i,j,k)+m11*fac*b1(i3,j,k)
!               p%ms(12,i,j,k)=p%ms(12,i,j,k)+m12*fac*b1(i2,j1,k)
!               p%ms(13,i,j,k)=p%ms(13,i,j,k)+m13*fac*b1(i2,j,k1)
!               p%ms(14,i,j,k)=p%ms(14,i,j,k)+m14*fac*b1(i1,j2,k)
!               p%ms(15,i,j,k)=p%ms(15,i,j,k)+m15*fac*b1(i1,j1,k1)
!               p%ms(16,i,j,k)=p%ms(16,i,j,k)+m16*fac*b1(i1,j,k2)
!               p%ms(17,i,j,k)=p%ms(17,i,j,k)+m17*fac*b1(i,j3,k)
!               p%ms(18,i,j,k)=p%ms(18,i,j,k)+m18*fac*b1(i,j2,k1)
!               p%ms(19,i,j,k)=p%ms(19,i,j,k)+m19*fac*b1(i,j1,k2)
!               p%ms(20,i,j,k)=p%ms(20,i,j,k)+m20*fac*b1(i,j,k3)
!               p%ms(21,i,j,k)=p%ms(21,i,j,k)+m21*fac*b1(i4,j,k)
!               p%ms(22,i,j,k)=p%ms(22,i,j,k)+m22*fac*b1(i3,j1,k)
!               p%ms(23,i,j,k)=p%ms(23,i,j,k)+m23*fac*b1(i3,j,k1)
!               p%ms(24,i,j,k)=p%ms(24,i,j,k)+m24*fac*b1(i2,j2,k)
!               p%ms(25,i,j,k)=p%ms(25,i,j,k)+m25*fac*b1(i2,j1,k1)
!               p%ms(26,i,j,k)=p%ms(26,i,j,k)+m26*fac*b1(i2,j,k2)
!               p%ms(27,i,j,k)=p%ms(27,i,j,k)+m27*fac*b1(i1,j3,k)
!               p%ms(28,i,j,k)=p%ms(28,i,j,k)+m28*fac*b1(i1,j2,k1)
!               p%ms(29,i,j,k)=p%ms(29,i,j,k)+m29*fac*b1(i1,j1,k2)
!               p%ms(30,i,j,k)=p%ms(30,i,j,k)+m30*fac*b1(i1,j,k3)
!               p%ms(31,i,j,k)=p%ms(31,i,j,k)+m31*fac*b1(i,j4,k)
!               p%ms(32,i,j,k)=p%ms(32,i,j,k)+m32*fac*b1(i,j3,k1)
!               p%ms(33,i,j,k)=p%ms(33,i,j,k)+m33*fac*b1(i,j2,k2)
!               p%ms(34,i,j,k)=p%ms(34,i,j,k)+m34*fac*b1(i,j1,k3)
!               p%ms(35,i,j,k)=p%ms(35,i,j,k)+m35*fac*b1(i,j,k4)

               tx=-1.0*tx
           END DO
           ty=-1.0_r8*ty
         END DO
         tz=-1.0_r8*tz
      END DO

      RETURN
      END SUBROUTINE COMP_CMS

      SUBROUTINE COMP_DIRECT(EnP,ibeg,iend,x,y,z,arrdim)
      IMPLICIT NONE
!
! COMP_DIRECT directly computes the potential on the current target
! particle (determined by the global variable TARPOS) due to the 
! current source cluster 

      INTEGER,INTENT(IN) :: ibeg,iend,arrdim
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(INOUT) :: EnP
      REAL(KIND=r8),DIMENSION(arrdim),INTENT(IN) :: x,y,z

! local variables

      INTEGER ::j,nn 
      REAL(KIND=r8) :: tx,ty,tz,tx2,ty2,tz2,txy,txz,tyz,txyz
      REAL(KIND=r8) :: etmp,rr,rr2,rr3,rr5,rr7,frr2,tfrr2,srr2
      
      IF (targid >= ibeg .AND. targid <= iend) THEN

         DO j=ibeg,targid-1
            nn=orderarr(j)
            tx=x(j)-tarpos(1)
            ty=y(j)-tarpos(2)
            tz=z(j)-tarpos(3)

            tx2=tx*tx; ty2=ty*ty; tz2=tz*tz; txy=tx*ty; txz=tx*tz; tyz=ty*tz; txyz=tx*ty*tz
            rr2=1.0_r8/(tx*tx+ty*ty+tz*tz); rr=SQRT(rr2); frr2=5.0_r8*rr2; srr2=7.0_r8*rr2; tfrr2=35.0_r8*rr2
            rr3=rr2*rr; rr5=3.0_r8*rr2*rr3; rr7=5.0_r8*rr5*rr2

            ! Differentiate with respect to y (The L_j operator)

            EnP(nn)=EnP(nn)+m1*rr+rr3*(m2*tx+m3*ty+m4*tz)-rr3*(m5+m8+m10)+&
                 rr5*(tx*(m5*tx+m6*ty+m7*tz)+ty*(m8*ty+m9*tz)+tz*m10*tz-&
                 m11*tx*(3.0_r8-tx2*frr2)-m17*ty*(3.0_r8-ty2*frr2)-m20*tz*(3.0_r8-tz2*frr2)-&
                 (m12*ty+m13*tz)*(1.0_r8-tx2*frr2)-(m14*tx+m18*tz)*(1.0_r8-ty2*frr2)-&
                 (m16*tx+m19*ty)*(1.0_r8-tz2*frr2))+m15*txyz*rr7+&
                 rr5*(m21*(3.0_r8+tx2*rr2*(tfrr2*tx2-30.0_r8))+m24*(1.0_r8-frr2*(tx2+ty2-srr2*tx2*ty2))+&
                 m26*(1.0_r8-frr2*(tx2+tz2-srr2*tx2*tz2))+m31*(3.0_r8+ty2*rr2*(tfrr2*ty2-30.0_r8))+&
                 m33*(1.0_r8-frr2*(ty2+tz2-srr2*ty2*tz2))+m35*(3.0_r8+tz2*rr2*(tfrr2*tz2-30.0_r8)))+&
                 rr7*(m22*txy*(srr2*tx2-3.0_r8)+m23*txz*(srr2*tx2-3.0_r8)+m25*tyz*(srr2*tx2-1.0_r8)+&
                 m27*txy*(srr2*ty2-3.0_r8)+m28*txz*(srr2*ty2-1.0_r8)+m29*txy*(srr2*tz2-1.0_r8)+&
                 m30*txz*(srr2*tz2-3.0_r8)+m32*tyz*(srr2*ty2-3.0_r8)+m34*tyz*(srr2*tz2-3.0_r8))

         END DO   

         DO j=targid+1,iend
            nn=orderarr(j)
            tx=x(j)-tarpos(1)
            ty=y(j)-tarpos(2)
            tz=z(j)-tarpos(3)

            tx2=tx*tx; ty2=ty*ty; tz2=tz*tz; txy=tx*ty; txz=tx*tz; tyz=ty*tz; txyz=tx*ty*tz
            rr2=1.0_r8/(tx*tx+ty*ty+tz*tz); rr=SQRT(rr2); frr2=5.0_r8*rr2; srr2=7.0_r8*rr2; tfrr2=35.0_r8*rr2
            rr3=rr2*rr; rr5=3.0_r8*rr2*rr3; rr7=5.0_r8*rr5*rr2

            ! Differentiate with respect to y (The L_j operator)

            EnP(nn)=EnP(nn)+m1*rr+rr3*(m2*tx+m3*ty+m4*tz)-rr3*(m5+m8+m10)+&
                 rr5*(tx*(m5*tx+m6*ty+m7*tz)+ty*(m8*ty+m9*tz)+tz*m10*tz-&
                 m11*tx*(3.0_r8-tx2*frr2)-m17*ty*(3.0_r8-ty2*frr2)-m20*tz*(3.0_r8-tz2*frr2)-&
                 (m12*ty+m13*tz)*(1.0_r8-tx2*frr2)-(m14*tx+m18*tz)*(1.0_r8-ty2*frr2)-&
                 (m16*tx+m19*ty)*(1.0_r8-tz2*frr2))+m15*txyz*rr7+&
                 rr5*(m21*(3.0_r8+tx2*rr2*(tfrr2*tx2-30.0_r8))+m24*(1.0_r8-frr2*(tx2+ty2-srr2*tx2*ty2))+&
                 m26*(1.0_r8-frr2*(tx2+tz2-srr2*tx2*tz2))+m31*(3.0_r8+ty2*rr2*(tfrr2*ty2-30.0_r8))+&
                 m33*(1.0_r8-frr2*(ty2+tz2-srr2*ty2*tz2))+m35*(3.0_r8+tz2*rr2*(tfrr2*tz2-30.0_r8)))+&
                 rr7*(m22*txy*(srr2*tx2-3.0_r8)+m23*txz*(srr2*tx2-3.0_r8)+m25*tyz*(srr2*tx2-1.0_r8)+&
                 m27*txy*(srr2*ty2-3.0_r8)+m28*txz*(srr2*ty2-1.0_r8)+m29*txy*(srr2*tz2-1.0_r8)+&
                 m30*txz*(srr2*tz2-3.0_r8)+m32*tyz*(srr2*ty2-3.0_r8)+m34*tyz*(srr2*tz2-3.0_r8))

         END DO


      ELSE

         DO j=ibeg,iend
            nn=orderarr(j)
            tx=x(j)-tarpos(1)
            ty=y(j)-tarpos(2)
            tz=z(j)-tarpos(3)

            tx2=tx*tx; ty2=ty*ty; tz2=tz*tz; txy=tx*ty; txz=tx*tz; tyz=ty*tz; txyz=tx*ty*tz
            rr2=1.0_r8/(tx*tx+ty*ty+tz*tz); rr=SQRT(rr2); frr2=5.0_r8*rr2; srr2=7.0_r8*rr2; tfrr2=35.0_r8*rr2
            rr3=rr2*rr; rr5=3.0_r8*rr2*rr3; rr7=5.0_r8*rr5*rr2

            ! Differentiate with respect to y (The L_j operator)
 
            EnP(nn)=EnP(nn)+m1*rr+rr3*(m2*tx+m3*ty+m4*tz)-rr3*(m5+m8+m10)+&
                 rr5*(tx*(m5*tx+m6*ty+m7*tz)+ty*(m8*ty+m9*tz)+tz*m10*tz-&
                 m11*tx*(3.0_r8-tx2*frr2)-m17*ty*(3.0_r8-ty2*frr2)-m20*tz*(3.0_r8-tz2*frr2)-&
                 (m12*ty+m13*tz)*(1.0_r8-tx2*frr2)-(m14*tx+m18*tz)*(1.0_r8-ty2*frr2)-&
                 (m16*tx+m19*ty)*(1.0_r8-tz2*frr2))+m15*txyz*rr7+&
                 rr5*(m21*(3.0_r8+tx2*rr2*(tfrr2*tx2-30.0_r8))+m24*(1.0_r8-frr2*(tx2+ty2-srr2*tx2*ty2))+&
                 m26*(1.0_r8-frr2*(tx2+tz2-srr2*tx2*tz2))+m31*(3.0_r8+ty2*rr2*(tfrr2*ty2-30.0_r8))+&
                 m33*(1.0_r8-frr2*(ty2+tz2-srr2*ty2*tz2))+m35*(3.0_r8+tz2*rr2*(tfrr2*tz2-30.0_r8)))+&
                 rr7*(m22*txy*(srr2*tx2-3.0_r8)+m23*txz*(srr2*tx2-3.0_r8)+m25*tyz*(srr2*tx2-1.0_r8)+&
                 m27*txy*(srr2*ty2-3.0_r8)+m28*txz*(srr2*ty2-1.0_r8)+m29*txy*(srr2*tz2-1.0_r8)+&
                 m30*txz*(srr2*tz2-3.0_r8)+m32*tyz*(srr2*ty2-3.0_r8)+m34*tyz*(srr2*tz2-3.0_r8))

         END DO

      END IF

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

      DEALLOCATE(cf1,cf2,b1,jmp, STAT=err)
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
      SUBROUTINE TREECODE(x,y,z,mpole,numpars,mpdim, &
                      tEn,tpeng,order,theta,shrink,&
                      maxparnode,timetree,treelevel,&
                      iflag)

!====================================================================
!                                                                   
! x,y,z,mpole    :: x,y,z coordinates and multipoles of particles        
! numpar         :: number of particles 
! mpdim          :: number of independent multipoles
! tEn            :: array of dimension numparT for storing potential 
!                   at each target
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

      INTEGER,INTENT(IN) :: numpars,mpdim,order,shrink,maxparnode,treelevel,iflag
      REAL(KIND=r8),DIMENSION(numpars),INTENT(INOUT) :: x,y,z
      REAL(KIND=r8),DIMENSION(numpars,1:mpdim),INTENT(INOUT) :: mpole
      REAL(KIND=r8),DIMENSION(numpars),INTENT(OUT) :: tEn
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

      CALL SETUP(x,y,z,numpars,mpdim,order,theta,xyzminmax)

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

      CALL CP_TREECODE(troot,x,y,z,mpole,tpeng,tEn,numpars,mpdim) 

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



