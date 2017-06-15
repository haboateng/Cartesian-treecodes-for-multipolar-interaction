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

      include "comms.inc"

      INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(12)

! runtime parameters

      INTEGER :: numpars,order,maxparnode
      INTEGER :: shrink,treelevel,iflag
      REAL(KIND=r8) :: theta,t1,t2,rmsE,temp 

! arrays for coordinates and charges and energy of target particles

      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: x,y,z  !source particles
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:,:) :: mpole
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: denergy,tenergy !exact energy and energy via treecode
      REAL(KIND=r8),ALLOCATABLE,DIMENSION(:) :: energy

! variables for potential energy computation

      REAL(KIND=r8) :: tpeng,dpeng,peng

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

      INTEGER :: i,j,err,ierr,mpdim

! Parameters for parallel processing
      INTEGER :: idnode, mxnode

!Begin MPI Communications

      CALL INITCOMMS()
      CALL GSYNC()

!Find the id of processor

      CALL MACHINE(idnode,mxnode)

! Let master node (idnode=0) Read in input
      IF(idnode==0)THEN

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

      END IF

!Broadcast number of number of particles-numpars,
!theta, order, maxparnode, treelevel,
!shrink, iflag, to each node

      CALL MPI_BCAST(numpars,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL GSYNC()
      CALL MPI_BCAST(theta,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      CALL GSYNC()
      CALL MPI_BCAST(order,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL GSYNC()
      CALL MPI_BCAST(maxparnode,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL GSYNC()
      CALL MPI_BCAST(treelevel,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL GSYNC()
      CALL MPI_BCAST(shrink,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL GSYNC()
      CALL MPI_BCAST(iflag,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL GSYNC()


      mpdim=35
      ALLOCATE(x(numpars),y(numpars),z(numpars),mpole(numpars,1:mpdim),STAT=err)
      IF (err .NE. 0) THEN
          WRITE(6,*) 'Error allocating arrays for x, y, z and mpole! '
          STOP
      END IF

      ALLOCATE(tenergy(numpars),denergy(numpars),energy(numpars),STAT=err)
      IF (err .NE. 0) THEN
          WRITE(6,*) 'Error allocating tenergy or denergy! '
          STOP
      END IF

      CALL GSYNC()
!On the master node, open the file containing the coordinates
!and charges of the sources and read them in

      IF(idnode==0) THEN 

         OPEN(unit=82,file=sampin1,status='old',action='read')
 

! Read in coordinates 
         DO i=1,numpars
            READ(82,*) x(i),y(i),z(i)
         END DO
     
         CLOSE(82)

         OPEN(unit=84,file=sampin2,status='old',action='read')

! Read in the values for the exact energy at each target       
         READ(84,13)timedirect,dpeng
         DO i=1,numpars
            READ(84,17)denergy(i)
         END DO

         CLOSE (84)

      END IF

! Generate multipoles on all processors

      CALL RANDOM_SEED()
      CALL RANDOM_NUMBER(mpole); mpole = -1.0_r8 + 2.0_r8*mpole


!Broadcast the coordinates and multipoles to all processors

      CALL MPI_BCAST(x,numpars,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      CALL GSYNC()
      CALL MPI_BCAST(y,numpars,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      CALL GSYNC()
      CALL MPI_BCAST(z,numpars,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      CALL GSYNC()


!Call treecode and get the computation time on master node
       
      IF(idnode==0)timetree=MPI_WTIME()

! Calling main subroutine to approximate the energy

      CALL TREECODE(x,y,z,mpole,numpars,mpdim,tenergy,tpeng,&
                    order,theta,shrink,maxparnode,treelevel,&
                    iflag,idnode,mxnode)

      CALL GSYNC
     
!Sum up energies over all processors and send to master
      CALL MPI_REDUCE(tenergy,energy,numpars,MPI_DOUBLE_PRECISION,&
           MPI_SUM,0,MPI_COMM_WORLD,ierr)
!      CALL GSYNC()
      CALL MPI_REDUCE(tpeng,peng,1,MPI_DOUBLE_PRECISION,&
           MPI_SUM,0,MPI_COMM_WORLD,ierr)
!      CALL GSYNC()

      IF(idnode==0)timetree=MPI_WTIME()-timetree

!On the master node
      IF(idnode==0)THEN
         !Find the time for the treecode 
!         timetree=MPI_WTIME()-timetree

         WRITE(6,*) ' '
         WRITE(6,*) 'Tree time (secs)   : ',timetree
         WRITE(6,*) ' '

         WRITE(6,*) ' '
         WRITE(6,*) 'Direct potential energy :',dpeng
         WRITE(6,*) 'Tree   potential energy :',peng
         WRITE(6,*) ' '
         WRITE(6,*) 'Abs and Rel Inf norm error for total potential : ' 
         WRITE(6,*) ' '
         WRITE(6,14) ABS(peng-dpeng),ABS((peng-dpeng)/dpeng)
         WRITE(6,*) ' '

! compute energy errors

         inferr=0.0_r8
         relinferr=0.0_r8
         
         t1=0.0_r8; t2=0.0_r8

         DO j=1,numpars
            temp=ABS( (denergy(j)-energy(j))/energy(j) )
            IF (temp .GT. relinferr)THEN
              relinferr=temp
              i=j
            END IF
!            WRITE(200,*)j
!            WRITE(200,14)denergy(j),energy(j)
         END DO

         inferr=MAXVAL(ABS(denergy-energy))
         rmsE=SQRT(DOT_PRODUCT(denergy-energy,denergy-energy)/DOT_PRODUCT(denergy,denergy))

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

         OPEN(unit=85,file=sampout,status='unknown',action='write',&
              position='append')

         write(85,15)mxnode,numpars,maxparnode,treelevel,order,theta
         write(85,16)ABS((peng-dpeng)/dpeng),rmsE,timedirect,timetree

         CLOSE(unit=85)

      END IF

      DEALLOCATE(x,y,z,mpole,denergy,tenergy,energy)

      CALL EXITCOMMS()

 13   FORMAT(2X,E24.16,2X,E24.16)
 14   FORMAT(E24.16,1X,E24.16)
 15   FORMAT(I4,2X,I8,2X,I4,2X,I4,2X,I2,2X,F12.8)
 16   FORMAT(4(2X,F24.16)) 
 17   FORMAT(14X,E24.16)
      STOP
      END PROGRAM TREEDRIVER

!!!!!!!!!!!!!!

