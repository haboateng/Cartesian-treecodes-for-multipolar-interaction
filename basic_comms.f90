      subroutine initcomms()
      
!*********************************************************************
!     
!     communication harness initialisation
!     
!     MPI version - t.forester may 1995
!     CPP version - w.smith may 1995
!     
!     wl
!     2008/01/14 13:33:07
!     1.4
!     Exp
!     
!*********************************************************************
      
      implicit none
      
      include "comms.inc"
      
      integer :: ierr


      call MPI_init(ierr)

      return
      end

      subroutine machine(idnode,mxnode)

!*********************************************************************
!     
!     dl_poly subroutine for obtaining charcteristics of
!     the computer on which the program is being run
!     
!     copyright daresbury laboratory 1992
!     author - w.smith july 1992
!     
!     MPI version - t.forester may 1995
!
!     wl
!     1.4
!     Exp
!*********************************************************************

      implicit none

      integer :: idnode,mxnode,mynode,numnodes

      idnode=mynode()
      mxnode=numnodes()

      return
      end

      integer function mynode()

!*********************************************************************
!
!     routine to determine identity of processing node 
!
!     MPI version - t.forester may 1995
!
!     wl
!     2008/01/14 13:33:07
!     1.4
!     Exp
!
!*********************************************************************

      implicit none

      include "comms.inc"

      integer :: ierr


      call MPI_COMM_RANK(MPI_COMM_WORLD, mynode ,ierr)

      return
      end

      integer function nodedim()

!*********************************************************************
!
!     calculate dimension of hypercube
!
!     MPI version - t.forester may 1995
!
!     wl
!     2008/01/14 13:33:07
!     1.4
!     Exp
!
!*********************************************************************

      implicit none

      include "comms.inc"

      integer :: i,n,ierr,mxnode


      call MPI_COMM_SIZE(MPI_COMM_WORLD, mxnode ,ierr)
      n=1
      nodedim = -1
      do i=0,16

         if(n.eq.mxnode)nodedim=i
         n=2*n

      enddo

      return
      end

      integer function numnodes()

!*********************************************************************
!
!     calculate number of nodes
!
!     MPI version - t.forester may 1995
!
!     wl
!     2008/01/14 13:33:07
!     1.4
!     Exp
!
!*********************************************************************

      implicit none

      include "comms.inc"

      integer :: ierr


      call MPI_COMM_SIZE(MPI_COMM_WORLD, numnodes, ierr)

      return
      end

      subroutine csend(tagmsg,buf,length,pe,idum)

!*********************************************************************
!
!     Intel-like  csend (double precision)
!
!     MPI version - t.forester may 1995
!     CPP version - w.smith may 1995
!
!     wl
!     2008/01/14 13:33:07
!     1.4
!     Exp
!
!*********************************************************************

      implicit none

      include "comms.inc"

      INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(12)
      integer :: tagmsg,length,pe,idum

      integer :: ierr
      real(kind=r8) buf(*)


      call MPI_send(buf,length,MPI_DOUBLE_PRECISION,pe,tagmsg,&
          MPI_COMM_WORLD,ierr)

      return
      end

      subroutine crecv(tagmsg,buf,length)

!*********************************************************************
!
!     Intel-like  crecv (double precision)
!
!     MPI version - t.forester may 1995
!     CPP version - w.smith may 1995
!
!     wl
!     2008/01/14 13:33:07
!     1.4
!     Exp
!
!*********************************************************************

      implicit none

      include "comms.inc"

      INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(12)
      integer :: tagmsg,length

      integer :: ierr
      integer :: status(MPI_STATUS_SIZE)
      real(kind=r8) buf(*)


      call MPI_RECV(buf,length,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,&
          tagmsg,MPI_COMM_WORLD,status,ierr)

      return 
      end

      subroutine gisum(aaa,nnn,bbb)

!***********************************************************************
!     
!     dl_poly global summation subroutine for hypercube - MPI version
!     integer version
!     
!     copyright - daresbury laboratory 1992
!     author    - w. smith march 1992.
!     MPI version - t.forester may 1995
!     CPP version - w.smith may 1995
!     
!     wl
!     2008/01/14 13:33:07
!     1.4
!     Exp
!
!***********************************************************************
      
      implicit none

      integer :: nnn,i,ierror,iii,kk,k,k0,k1,k2,msg1,msg2
      integer :: aaa(nnn),bbb(nnn)

      include "comms.inc"

      integer :: status(MPI_STATUS_SIZE)


      call MPI_allreduce(aaa,bbb,nnn,MPI_INTEGER, &
       MPI_SUM,MPI_COMM_WORLD,ierror)

      do i = 1,nnn
        aaa(i) = bbb(i)
      enddo

      return
      end

      subroutine gdsum(aaa,nnn,bbb)

!***********************************************************************
!     
!     dl_poly global summation subroutine for MPI - hypercube assumed
!     double precision version
!     
!     copyright - daresbury laboratory 1995
!     author    - w. smith march 1992.
!     MPI version - t.forester may 1995
!     CPP version - w.smith may 1995
!     
!     wl
!     2008/01/14 13:33:07
!     1.4
!     Exp
!
!***********************************************************************

      implicit none

      INTEGER,PARAMETER :: r8=SELECTED_REAL_KIND(12)
      integer :: nnn,i,iii,kk,k1,k2,ierror
      real(kind=8) aaa(nnn),bbb(nnn)

      include "comms.inc"

      integer :: status(MPI_STATUS_SIZE)


      call MPI_allreduce(aaa,bbb,nnn,MPI_DOUBLE_PRECISION, &
       MPI_SUM,MPI_COMM_WORLD,ierror)

        do i = 1,nnn
          aaa(i) = bbb(i)
        enddo

      return
      end

      subroutine gstate(check)

!***********************************************************************
!     
!     dl_poly global status subroutine : gisum version
!     
!     copyright - daresbury laboratory 1992
!     author    - w. smith       march 1992
!     MPI version -  t. forester may 1995
!     
!     wl
!     1.4
!     Exp
!***********************************************************************


      implicit none

      logical :: check
      integer :: i,idum

      i = 0
      if(.not.check) i = 1

      call gisum(i,1,idum)
      
      check = (i.eq.0)

      return
      end

      subroutine gsync()

!*********************************************************************
!     
!     barrier / synchronization routine
!
!     MPI version - t.forester may 1995
!     CPP version - w.smith
!
!     wl
!     2008/01/14 13:33:07
!     1.4
!     Exp
!
!*********************************************************************

      implicit none

      integer :: ierr

      include "comms.inc"


      call  MPI_BARRIER(MPI_COMM_WORLD,ierr)

      return
      end

      subroutine exitcomms()

!*********************************************************************
!
!     exitcomms: exit from communication harness
!
!     MPI version - t.forester may 1995
!     CPP version - w.smith may 1995
!
!     wl
!     2008/01/14 13:33:07
!     1.4
!     Exp
!
!*********************************************************************

      implicit none

      include "comms.inc"

      integer :: ierr

      call MPI_FINALIZE(ierr)
      call exit(0)

      return
      end
