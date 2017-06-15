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

      PROGRAM TREEDRIVER
      IMPLICIT NONE

      INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(12)

! runtime parameters

      INTEGER :: numpars,order,maxparnode
      INTEGER :: shrink,treelevel,iflag
      REAL(KIND=r8) :: theta,t1,t2,rmsE,temp 

! arrays for coordinates and charges and energy of target particles

      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: x,y,z  !source particles
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:,:) :: mpole
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: denergy,tenergy !exact energy and energy via treecode

! variables for potential energy computation

      REAL(KIND=r8) :: tpeng,dpeng

! variables needed for f90 DATE_AND_TIME intrinsic

      INTEGER,DIMENSION(8) :: time1,time2 
      CHARACTER (LEN=8)  :: datec
      CHARACTER (LEN=10) :: timec
      CHARACTER (LEN=5)  :: zonec
      CHARACTER (LEN=20) :: sampin1,sampin2,sampin3,sampout
      REAL(KIND=r8)      :: totaltime,timedirect,timetree

! variables for error calculations

      REAL(KIND=r8) :: inferr,relinferr

! local variables

      INTEGER :: i,j,err,mpdim

! EXECUTABLE STATEMENTS
      WRITE(6,*) 'Enter the name of input file 1'
      READ(5,*) sampin1
      WRITE(6,*) 'Enter the name of input file 2'
      READ(5,*) sampin2
      WRITE(6,*) 'Enter the name of output file'
      READ(5,*) sampout
      WRITE(6,*) 'Enter NUMPARS : '
      READ(5,*) numpars
      WRITE(6,*) 'Enter THETA : ' ! The multipole acceptability criterion
      READ(5,*) theta         
      WRITE(6,*) 'Enter ORDER : ' ! The order of the approximation 
      READ(5,*) order
      WRITE(6,*) 'Enter MAXPARNODE : '  ! maximum number of particles in a leaf
      READ(5,*) maxparnode 
      WRITE(6,*) 'Enter TREELEVEL : '   ! maximum number of levels of tree
      READ(5,*) treelevel
      WRITE(6,*) 'Enter SHRINK (0-No 1-Yes) : ' ! Flag to decide whether to shrink 
      READ(5,*) shrink                 ! cluster to minimal box containing sources/targets
      WRITE(6,*) 'Enter IFLAG (0-divide tree till number of particles in leaf is less or equal to maxparnode)'
      WRITE(6,*) '            (1-divide tree till maxlevel attained )'
      READ(5,*) iflag   ! Flag to decide whether tree divides till number of particles in a leaf is
                              ! at most maxparnode or till number of levels in a tree is maxlevel

      mpdim=35
      ALLOCATE(x(numpars),y(numpars),z(numpars),mpole(numpars,1:mpdim),STAT=err)
      IF (err .NE. 0) THEN
          WRITE(6,*) 'Error allocating arrays for x, y, z and mpole! '
          STOP
      END IF

      ALLOCATE(tenergy(numpars),denergy(numpars),STAT=err)
      IF (err .NE. 0) THEN
          WRITE(6,*) 'Error allocating tenergy or denergy! '
          STOP
      END IF

      OPEN(unit=82,file=sampin1,status='old',action='read')
 

! Read in coordinates and charges for the sources
      DO i=1,numpars
         READ(82,*) x(i),y(i),z(i)
      END DO
     
      CLOSE(82)

      CALL RANDOM_SEED()
      CALL RANDOM_NUMBER(mpole); mpole = -1.0_r8 + 2.0_r8*mpole

      OPEN(unit=84,file=sampin2,status='old',action='read')

! Read in the values for the exact energy at each target       
      READ(84,13)timedirect,dpeng
      DO i=1,numpars
         READ(84,17)denergy(i)
      END DO

      CLOSE (84)

      OPEN(unit=85,file=sampout,status='unknown',action='write',&
           position='append')

! Calling main subroutine to approximate the energy

      CALL TREECODE(x,y,z,mpole,numpars,mpdim,tenergy,tpeng,&
                    order,theta,shrink,maxparnode,timetree,&
                    treelevel,iflag)


      WRITE(6,*) ' '
      WRITE(6,*) 'Direct time (secs) : ',timedirect
      WRITE(6,*) ' '
      WRITE(6,*) 'Tree time (secs)   : ',timetree
      WRITE(6,*) ' '

      WRITE(6,*) ' '
      WRITE(6,*) 'Direct potential energy :',dpeng
      WRITE(6,*) 'Tree   potential energy :',tpeng
      WRITE(6,*) ' '
      WRITE(6,*) 'Absolute and Relative Inf norm error for total potential : ' 
      WRITE(6,*) ' '
      WRITE(6,14) ABS(tpeng-dpeng),ABS((tpeng-dpeng)/dpeng)
      WRITE(6,*) ' '


! compute energy errors

         inferr=0.0_r8
         relinferr=0.0_r8
         
         t1=0.0_r8; t2=0.0_r8

         DO j=1,numpars
            temp=ABS( (denergy(j)-tenergy(j))/denergy(j) )
            IF (temp .GT. relinferr)THEN
              relinferr=temp
              i=j
            END IF
            WRITE(200,*)j
            WRITE(200,14)denergy(j),tenergy(j)
         END DO

         inferr=MAXVAL(ABS(denergy-tenergy))
         rmsE=SQRT(DOT_PRODUCT(denergy-tenergy,denergy-tenergy)/DOT_PRODUCT(denergy,denergy))

! output errors to standard out

         WRITE(6,*) ' '
         WRITE(6,*) 'Absolute Inf norm error in energy :'      
         WRITE(6,*) ' '
         WRITE(6,13) inferr
         WRITE(6,*) ' '
         WRITE(6,*) 'Relative Inf norm error in energy :'     
         WRITE(6,*) ' ' 
         WRITE(6,13) relinferr 
         WRITE(6,*) 'Relative 2 norm error in energy :'
         WRITE(6,*) ' '
         WRITE(6,13)rmsE
         WRITE(6,*)'max location = ',i
         
         write(85,15)numpars,iflag,maxparnode,treelevel,order,theta
         write(85,16)ABS((tpeng-dpeng)/dpeng),rmsE,timedirect,timetree

      CLOSE(unit=85)

 13   FORMAT(2X,E24.16,2X,E24.16)
 14   FORMAT(E24.16,1X,E24.16)
 15   FORMAT(I8,2X,I4,2X,I4,2X,I4,2X,I2,2X,F12.8)
 16   FORMAT(4(2X,F24.16)) 
 17   FORMAT(14X,E24.16)
      STOP
      END PROGRAM TREEDRIVER

!!!!!!!!!!!!!!

