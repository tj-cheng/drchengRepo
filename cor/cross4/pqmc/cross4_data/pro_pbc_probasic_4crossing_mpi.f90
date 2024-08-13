!-- systems with open boundary conditions --!
!
!   two intersecting chains
!  ------*-x-*-x---
!          !\
!  ------x-* *-x----- 
!   
!-------------!
module system
   !-------------!
   use mpi
   save

   integer :: comm,ierror  !MPI
   integer :: lx,lchain
   integer :: ly
   integer :: ly0
   integer :: lly    ! length of y in the real lattice 
   integer :: nn
   integer :: nb
   integer :: nd
   integer :: mm,mm0
   integer :: dd=1    ! disorder strength

   real(8), allocatable :: jj(:)
   real(8), allocatable :: pj(:)
   real(8), allocatable :: jj0(:)

   integer, allocatable :: ab(:)   ! A-/B- subsystems
   integer, allocatable :: oper(:)
   integer, allocatable :: lbnd(:)
   integer, allocatable :: rbnd(:)
   integer, allocatable :: lspn(:)
   integer, allocatable :: rspn(:)
   integer, allocatable :: frst(:)
   integer, allocatable :: last(:)
   integer, allocatable :: vrtx(:)
   integer, allocatable :: bsite(:,:)
   integer, allocatable :: dsite(:,:)
   integer, allocatable :: pos(:,:)  ! positions in the real lattice
   integer, allocatable :: ipos(:,:)
   integer, allocatable :: xy(:,:)
   integer, allocatable :: xyi(:,:)
   real(8), allocatable :: hamp(:,:)
   integer, allocatable :: dscale(:,:) !MH distance between two sites
   integer, allocatable :: endsite(:) !save the 8 end sites

   integer, allocatable :: ss(:) ! 4 end-spins


   real(8) :: enrg=0.d0
   real(8) :: mag2=0.d0

end module system
!-----------------!

!------------!
module mdata
   !------------!
   save

   integer :: nmsr=0.d0
   integer, allocatable :: plbnd(:)
   integer, allocatable :: prbnd(:)
   integer, allocatable :: loopnmbr(:)
   integer, allocatable :: loopsize(:)
   real(8), allocatable :: scor(:) !for loop (periodic chain) in the middle part of cross-4
   real(8), allocatable :: eescor(:,:) !there are 8 end sites and this store their correlation
   real(8), allocatable :: pscor(:) !:=1,2,3,...,128/2

end module mdata
!----------------!

!================!
program probasic
   !=====================!
   ! Anders Sandvik, 2012
   !---------------------!
   use system; use mdata; implicit none

   integer :: i,j,init,bins,mstps,istps,samples
   integer :: rank, numproc   !! MPI
   integer :: isize,num_size

   call MPI_INIT(ierror)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)
   call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)

   open(10,file='read.in',status='old')
   read(10,*)lx,ly0,mm0,dd
   read(10,*)init,samples,mstps,istps
   close(10)

   lchain=2*lx+4
   num_size = 0
   mm=mm0
   call initran(1,rank) 
   call allocateall()    ! set  mm=2**mm, mm=mm*nn/2, nn=lx*ly
   call makelattice()
   call dijkstra(rank)   !hamp

   do j=1,samples
      call initsystem(init)

      if (istps>0) then
         do i=1,istps
            call montecarlosweep()
         enddo
         !call writeconf()
      endif

      do i=1,mstps
         call montecarlosweep()
         call measure() 
      enddo
      call writebindata(mstps,rank)
      !call writeconf()
   enddo ! end of sample

   call deallocateall()
   !end do
   !open(777,file='numbins.dat',status='unknown')
   !write(777,*)samples*numproc, num_size
   !close(777)

   call MPI_FINALIZE(ierror)
end program probasic
!====================!

!----------------------------!
subroutine montecarlosweep()
   !----------------------------!
   use system; implicit none

   call diagonalupdate()
   call makevertexlist()
   call loopupdate()
   call stateupdate()

end subroutine montecarlosweep
!------------------------------!

!---------------------------!
subroutine diagonalupdate()
   !---------------------------!
   use system; implicit none

   integer :: i,b
   real(8) :: ran

   external :: ran

   lspn(:)=rspn(:)
   do i=0,mm-1
      if (mod(oper(i),2)==0) then       
         !10 b=int(ran()*dble(nb))+1
         10 call generatebond(b)
         !write(*,*) b
         if (lspn(bsite(1,b))/=lspn(bsite(2,b))) then
            oper(i)=2*b
         else
            goto 10
         endif
      else
         b=oper(i)/2
         lspn(bsite(1,b))=-lspn(bsite(1,b))
         lspn(bsite(2,b))=-lspn(bsite(2,b))     
      endif
   enddo

contains

   !--------------------------!
   subroutine generatebond(b)
      !--------------------------!
      integer :: b,b1,b2
      real(8) :: ran,p
      external :: ran

      b1=1
      b2=nb
      p=ran()
      do
         b=(b1+b2)/2               !binary search!
         if (p<=pj(b)) then        !pj=cumulative probability!
            b2=b
         else
            b1=b+1
         endif
         if (b1==b2) then
            b=b1
            exit
         endif
      enddo

   end subroutine generatebond
   !---------------------------!

end subroutine diagonalupdate
!-----------------------------!

!---------------------------!
subroutine makevertexlist()
   !---------------------------!
   use system; implicit none

   integer :: i,b,p,s1,s2,v0,v1,v2,v3,v4

   last(:)=-1
   frst(:)=-1
   do p=0,mm-1
      v0=4*p
      b=oper(p)/2
      s1=bsite(1,b)
      s2=bsite(2,b)
      v1=last(s1)
      v2=last(s2)
      if (v1/=-1) then
         vrtx(v1)=v0
         vrtx(v0)=v1
      else
         frst(s1)=v0
      endif
      if (v2/=-1) then
         vrtx(v2)=v0+1
         vrtx(v0+1)=v2
      else
         frst(s2)=v0+1
      endif
      last(s1)=v0+2
      last(s2)=v0+3       
   enddo

   v0=4*mm
   do i=1,nn
      v1=v0+1; v2=v0+2; v3=v0+3
      if (frst(i)>=0) then
         vrtx(frst(i))=v0
         vrtx(v0)=frst(i)
         vrtx(last(i))=v2
         vrtx(v2)=last(i)
      else
         vrtx(v0)=v2
         vrtx(v2)=v0
      endif
      vrtx(v1)=4*mm+4*rbnd(i)-3
      vrtx(v3)=4*mm+4*lbnd(i)-1
      v0=v0+4
   enddo  


end subroutine makevertexlist
!-----------------------------!

!-----------------------!
subroutine loopupdate()
   !-----------------------!
   use system; implicit none

   integer :: b,p,v0,v1,v2
   real(8) :: ran
   external :: ran

   do v0=0,4*mm+4*nn-2,2
      if (vrtx(v0)>=0) then
         v1=v0
         if (ran()<0.5d0) then
            call visitloop()
         else
            call fliploop()
         endif
      endif
   enddo

contains

   !------------------------!
   subroutine visitloop()
      !------------------------!

      do
         vrtx(v1)=-1
         v2=ieor(v1,1)
         v1=vrtx(v2)
         vrtx(v2)=-1
         if (v1==v0) exit
      enddo

   end subroutine visitloop
   !------------------------! 

   !---------------------!
   subroutine fliploop()
      !---------------------!

      do
         p=v1/4
         if (p<mm) then
            oper(p)=ieor(oper(p),1)
         else
            p=p-mm+1
            if (mod(v1,4)<2) then
               rspn(p)=-rspn(p)
            else
               lspn(p)=-lspn(p)
            endif
         endif
         vrtx(v1)=-1
         v2=ieor(v1,1)
         v1=vrtx(v2)
         vrtx(v2)=-1
         if (v1==v0) exit
      enddo

   end subroutine fliploop
   !-----------------------!

end subroutine loopupdate
!-------------------------!

!------------------------!
subroutine stateupdate()
   !------------------------!
   use system; implicit none

   integer :: i,d,x1,x2,y1,y2,s1,s2,s3,s4,nnn
   real(8) :: ran,w1,w2

   external :: ran


   do i=1,nn 
      d=int(ran()*nd)+1
      s1=dsite(1,d)
      s2=dsite(2,d)
      s3=lbnd(s1)
      s4=lbnd(s2)
      if (lspn(s1)==lspn(s4)) cycle
      !x1=abs(pos(1,s3)-pos(1,s1)); !x1=min(x1,lx-x1)
      !y1=abs(pos(2,s3)-pos(2,s1)); !y1=min(y1,ly-y1)
      !x2=abs(pos(1,s4)-pos(1,s2)); !x2=min(x2,lx-x2)
      !y2=abs(pos(2,s4)-pos(2,s2)); !y2=min(y2,ly-y2)
      !w1=hamp(x1,y1)*hamp(x2,y2)
      w1=hamp(s1,s3)*hamp(s2,s4)
      !x1=abs(pos(1,s4)-pos(1,s1)); !x1=min(x1,lx-x1)
      !y1=abs(pos(2,s4)-pos(2,s1)); !y1=min(y1,ly-y1)
      !x2=abs(pos(1,s3)-pos(1,s2)); !x2=min(x2,lx-x2)
      !y2=abs(pos(2,s3)-pos(2,s2)); !y2=min(y2,ly-y2)
      !w2=hamp(x1,y1)*hamp(x2,y2)
      w2=hamp(s1,s4)*hamp(s2,s3)
      if (ran()<w2/w1) then
         lbnd(s1)=s4
         lbnd(s4)=s1
         lbnd(s2)=s3
         lbnd(s3)=s2
      endif
   enddo      
   do i=1,nn
      d=int(ran()*nd)+1
      s1=dsite(1,d)
      s2=dsite(2,d)
      s3=rbnd(s1)
      s4=rbnd(s2)
      if (rspn(s1)==rspn(s4)) cycle
      !x1=abs(pos(1,s3)-pos(1,s1)); !x1=min(x1,lx-x1)
      !y1=abs(pos(2,s3)-pos(2,s1)); !y1=min(y1,ly-y1)
      !x2=abs(pos(1,s4)-pos(1,s2)); !x2=min(x2,lx-x2)
      !y2=abs(pos(2,s4)-pos(2,s2)); !y2=min(y2,ly-y2)
      !w1=hamp(x1,y1)*hamp(x2,y2)
      w1=hamp(s1,s3)*hamp(s2,s4)
      !x1=abs(pos(1,s4)-pos(1,s1)); !x1=min(x1,lx-x1)
      !y1=abs(pos(2,s4)-pos(2,s1)); !y1=min(y1,ly-y1)
      !x2=abs(pos(1,s3)-pos(1,s2)); !x2=min(x2,lx-x2)
      !y2=abs(pos(2,s3)-pos(2,s2)); !y2=min(y2,ly-y2)
      !w2=hamp(x1,y1)*hamp(x2,y2)
      w2=hamp(s1,s4)*hamp(s2,s3)
      if (ran()<w2/w1) then
         rbnd(s1)=s4
         rbnd(s4)=s1
         rbnd(s2)=s3
         rbnd(s3)=s2
      endif
   enddo      

end subroutine stateupdate
!--------------------------!

!---------------------------!
subroutine propagatebonds()
   !---------------------------!
   use system; use mdata; implicit none

   integer :: i,b,s1,s2,s3,s4

   prbnd(:)=rbnd(:)
   do i=0,mm/2-1
      b=oper(i)/2
      s1=bsite(1,b)
      s2=bsite(2,b)
      s3=prbnd(s1)      
      s4=prbnd(s2)
      prbnd(s1)=s2
      prbnd(s2)=s1
      prbnd(s3)=s4
      prbnd(s4)=s3
   enddo
   plbnd(:)=lbnd(:)
   do i=mm-1,mm/2,-1
      b=oper(i)/2
      s1=bsite(1,b)
      s2=bsite(2,b)
      s3=plbnd(s1)
      s4=plbnd(s2)
      plbnd(s1)=s2
      plbnd(s2)=s1
      plbnd(s3)=s4
      plbnd(s4)=s3
   enddo

end subroutine propagatebonds
!-----------------------------!

!----------------------!
subroutine makeloops(nl)
   !----------------------!
   use system; use mdata; implicit none

   integer :: i,s,s1,nl

   nl=0
   loopnmbr(:)=0
   loopsize(:)=0
   do s1=1,nn
      if (loopnmbr(s1)==0) then    
         nl=nl+1
         s=s1
         do 
            loopsize(nl)=loopsize(nl)+2
            loopnmbr(s)=nl
            s=prbnd(s)
            loopnmbr(s)=nl
            s=plbnd(s)
            if (s==s1) exit
         enddo
      endif
   enddo   

end subroutine makeloops
!------------------------!

!--------------------!
subroutine measure()
   !--------------------!
   use system; use mdata; implicit none

   integer :: i,j,k,k1,r,x,y,lj,nl,m2
   real(8) :: e1
   integer, allocatable :: np(:)
   real(8), allocatable :: pcor(:)

   !integer :: ss(4)   ! 4 end-spins
   allocate(np(1:lchain/2)); np=0
   allocate(pcor(1:lchain/2)); pcor=0.d0

   call propagatebonds()
   call makeloops(nl)

   do x=1,lchain/2
      y=x+lchain/2
      lj=loopnmbr(y)
      if (lj == loopnmbr(x)) scor(1)=scor(1)+0.75/(lchain/2) !-3/4
   end do

   do i=1,7
      do j=i+1,8
         x=endsite(i)
         y=endsite(j)
         lj=loopnmbr(y)
         if (lj == loopnmbr(x)) eescor(i,j)=eescor(i,j)+0.75 !-3/4
      end do
   end do


   if (lx .eq. 128) then
      do x=1,lchain/2-1
         do y=x+1,lchain
            r=dscale(x,y)
            np(r)=np(r)+1
            lj=loopnmbr(y)
            !print*,r,x,y
            if (lj == loopnmbr(x)) pcor(r)=pcor(r)+0.75 !-3/4
         end do
      end do

   do r=1,lchain/2
      pscor(r)=pcor(r)/np(r)
   end do
   end if

   nmsr=nmsr+1
   deallocate(np)
   deallocate(pcor)


end subroutine measure
!----------------------!

!-------------------------!
subroutine writebindata(steps,rank)
   !-------------------------!
   use system; use mdata; implicit none

   integer :: r,steps,rank
   integer :: i,j,x,y
   character(len=1024)  :: numasstring

   write(numasstring,'(I3)') rank

   !scor=scor/(dble(nn)*dble(nmsr))
   scor=scor/dble(nmsr)
   eescor=eescor/dble(nmsr)
   pscor=pscor/dble(nmsr)

   open(10,file='loopcor'//trim(adjustl(numasstring)) //'.dat',status='unknown',position='append')
   write(10,'(i8,f16.10)') lchain,scor(1)
   close(10)

   open(10,file='eecor'//trim(adjustl(numasstring)) //'.dat',status='unknown',position='append')
   do i=1,7
      do j=i+1,8
         x=endsite(i)
         y=endsite(j)
         write(10,'(i8,i6,i6,i6,f16.10)') lchain, x, y, dscale(x,y), eescor(i,j)
      end do
   end do
   close(10)
   if (lx .eq. 128) then
      open(10,file='paircor'//trim(adjustl(numasstring)) //'.dat',status='unknown',position='append')
      do x=1,lchain/2-1
         do y=x+1,lchain
            r=dscale(x,y)
            write(10,'(i8,i6,f16.10)') lchain, r, pscor(r)
         end do
      end do
      close(10)
   end if
   nmsr=0
   scor=0.d0
   pscor=0.d0
   eescor=0.d0

end subroutine writebindata
!!---------------------------!
!
!!----------------------!
subroutine writeconf()
   !----------------------!
   use system; implicit none

   integer :: i

   open(10,file='conf',status='replace')
   do i=1,nn
      write(10,'(2i3,2i7)')lspn(i),rspn(i),lbnd(i),rbnd(i)
   enddo
   do i=0,mm-1
      write(10,*)oper(i) 
   enddo
   close(10)

end subroutine writeconf
!------------------------!

!---------------------!
subroutine readconf()
   !---------------------!
   use system; implicit none

   integer :: i

   open(10,file='conf',status='old')
   do i=1,nn
      read(10,*)lspn(i),rspn(i),lbnd(i),rbnd(i)
   enddo
   do i=0,mm-1
      read(10,*)oper(i)
   enddo
   close(10)

end subroutine readconf
!-----------------------!

!---------------------!
subroutine readbond()
   !---------------------!
   use system; implicit none

   integer :: i

   open(10,file='bond.in',status='old')
   do i=1,nb
      read(10,*)jj0(i)
   enddo
   close(10)

end subroutine readbond
!-----------------------!

!---------------------------!
subroutine initsystem(init)
   !---------------------------!
   use system; use mdata; implicit none

   integer :: i,j,a,x,y,init,b,x0,y0
   real(8) :: ran

   external :: ran

   !-- random couplings --!
   if (ly0==0) then        ! 1D chain for test
      call readbond()
      jj = jj0
   else
      do b=1,nb
         jj(b) = (ran())**dd
      enddo
   endif

   pj(1)=jj(1)
   !write(*,*) 1, jj(1), pj(1)
   do b=2,nb
      pj(b)=pj(b-1)+jj(b)
      !write(*,*) b, jj(b), pj(b)
   enddo
   pj=pj/pj(nb)

   !-- operator-string -------!

   if (init==0) then

      !if(ly0<2) then                    ! single chain
      !   do x=0,lx-1
      !      i=xyi(x,0)
      !      lspn(i)=(-1)**(xy(1,i)+xy(2,i))
      !      ab(i)=(-1)**(xy(1,i)+xy(2,i))
      !   enddo
      !elseif (ly0==2) then              ! two intersecting chains
      !   x0=lx/2-1
      !   y=1                              ! second chain
      !   do x=0,lx-1        
      !      i=xyi(x,y)
      !      lspn(i)=(-1)**(xy(1,i)+xy(2,i))
      !      ab(i)=(-1)**(xy(1,i)+xy(2,i))
      !   enddo 
      !   y=0
      !   do x=0,x0
      !      i=xyi(x,y)
      !      lspn(i)=(-1)**(xy(1,i)+xy(2,i))
      !      ab(i)=(-1)**(xy(1,i)+xy(2,i))
      !   enddo
      !   do x=x0+1,lx-1
      !      i=xyi(x,y)
      !      j=xyi(x-1,y)                     ! sign change
      !      lspn(i)=(-1)**(xy(1,j)+xy(2,j))
      !      ab(i)=(-1)**(xy(1,j)+xy(2,j)) 
      !   enddo
      !elseif (ly0==4) then              ! four intersecting chains
      !   x0=lx/2-1
      !   y=0
      !   do x=0,lx-1
      !      i=xyi(x,y)
      !      if (x .gt. lx/4-1 .and. x .lt. 3*lx/4) then
      !         j=xyi(x-1,y)
      !         lspn(i)=(-1)**(xy(1,j)+xy(2,j))
      !         ab(i)=(-1)**(xy(1,j)+xy(2,j)) 
      !      else 
      !         lspn(i)=(-1)**(xy(1,i)+xy(2,i))
      !         ab(i)=(-1)**(xy(1,i)+xy(2,i))
      !      end if 
      !   end do
      !   y=2
      !   do x=0,lx-1
      !      i=xyi(x,y)
      !      j=xyi(x,y-2)
      !      ab(i)=(-1)*ab(j)
      !      lspn(i)=(-1)*lspn(j)
      !   end do
      !   y=1
      !   do x=0,lx-1        
      !      i=xyi(x,y)
      !      lspn(i)=(-1)**(xy(1,i)+xy(2,i))
      !      ab(i)=(-1)**(xy(1,i)+xy(2,i))
      !   enddo 
      !   y=3
      !   do x=0,lx-1
      !      i=xyi(x,y)
      !      j=xyi(x,y-2)
      !      ab(i)=(-1)*ab(j)
      !      lspn(i)=(-1)*lspn(j)
      !   end do
      !endif

      !pbc new
      lbnd=0
      lspn=0
      b=lchain
      do
         b=b+1
         x=bsite(1,b)
         y=bsite(2,b)
         if (b .gt. nb) exit
         if (lbnd(x) /= 0 .or. lbnd(y) /= 0) cycle
         lbnd(x)=y
         lbnd(y)=x
         lspn(x)=1
         lspn(y)=-1
         !print*,x,y,lspn(x),lspn(y)
      end do
      b=0
      do
         b=b+1
         x=bsite(1,b)
         y=bsite(2,b)
         if (b .gt. lchain) exit
         if (lbnd(x) /= 0 .or. lbnd(y) /= 0) cycle
         lbnd(x)=y
         lbnd(y)=x
         lspn(x)=1
         lspn(y)=-1
         !print*,x,y,lspn(x),lspn(y)
      end do
      rbnd=lbnd
      rspn=lspn    
      !pbc new

      !!check
      !do i=1,nn
      !   print*,i,ab(i),lspn(i)
      !end do

      !if(ly0<2) then        ! single chain
      !   do y=0,ly-1
      !      do x=0,lx-2,2
      !         i=xyi(x,y)
      !         j=xyi(x+1,y)
      !         lbnd(i)=j
      !         lbnd(j)=i
      !      enddo
      !   enddo
      !elseif(ly0==2) then  ! two intersecting chains
      !   do y=0,ly-1
      !      do x=0,lx-1,2
      !         i=xyi(x,y)
      !         j=xyi(x+1,y)
      !         lbnd(i)=j
      !         lbnd(j)=i
      !      enddo
      !   enddo
      !elseif(ly0==4) then  ! four intersecting chains
      !   do y=0,ly-1
      !      do x=0,lx-1,2
      !         i=xyi(x,y)
      !         j=xyi(x+1,y)
      !         lbnd(i)=j
      !         lbnd(j)=i
      !      enddo
      !   enddo
      !endif 
      !rbnd=lbnd

      do i=0,mm-1
         oper(i)=2*(int(ran()*nb)+1)
      enddo
   else
      call readconf()
   endif

end subroutine initsystem
!-------------------------!

!------------------------!
subroutine makelattice()
   !----------------------------------------------------------------------------------------!
   ! Constructs the lattice, in the form of the list of sites connected by nearest-neighbor
   ! bonds 'bsite' (defining the hamiltonian) and next-nearest-neighbor sites 'dsite' (uses.
   ! in the sampling of the valence-bond trial state.
   ! There are three types of lattices and boundary conditions, depending on ly:
   ! ly0=0  : 1D chain with a given set of couplings for test
   ! ly0=1  : 1D chain (periodic only in x-direction)
   ! ly0=2  : Two chains with one crossing 
   ! Note that the system is frustrated for ly>2 and odd and the program would not give
   ! correct results in that case (unless one changes to open boundaries in the y-direction)
   !----------------------------------------------------------------------------------------!
   use system; implicit none

   integer :: s,a,x1,x2,y1,y2,x,y,b,d,x0,y0,s0,s1,s2,s3  ! (x0,y0),(s0,s1,s2) the intersection 
   integer :: k,px1 !ly0=4
   integer :: i,j !check
   integer :: xa,xb,leftan,rightan,numan,ls !ls label of s
   integer, allocatable :: oldseed(:),dseed(:),numd
   integer :: record_la,record_ra
   integer :: num_end
   !real(8) :: ran

   !external :: ran
   num_end=0 !max=8

   b=0
   do x=1,lchain
      s=x
      b=b+1
      bsite(1,b)=s
      bsite(2,b)=mod(s+1,lchain+1)+s/lchain
      !print*,bsite(1,b),bsite(2,b)
   end do

   !print*,'nb=',nb,b
   xb=(lchain-4)/4
   xa=xb/2
   leftan=lx/4-1  !the length of antenna
   rightan=lx/4

   !!(xa+1),(xa+1)+(xb+1)*1,(xa+1)+(xb+1)*2,(xa+1)+(xb+1)*3 天線位置

   s0=xa+1
   s1=s0+xb+1
   s2=s1+xb+1
   s3=s2+xb+1
   numan=1
   !print*,'b=',b

   s=lchain+1
   ls=lchain
   k=0
   record_la=0
   record_ra=0
   do
      s=s+1
      ls=ls+1
      if (s .eq. 2*lchain-1) exit
      if (s+1 .eq. s0+lchain+(xb+1)*k) then
         b=b+1
         !bsite(1,b)=s
         bsite(1,b)=ls
         bsite(2,b)=s0+(xb+1)*k
         !print*,bsite(1,b),bsite(2,b),'lv'
         b=b+1
         !bsite(1,b)=s+1
         bsite(1,b)=ls+1
         bsite(2,b)=s0+(xb+1)*k
         !print*,bsite(1,b),bsite(2,b),'rv'
         k=k+1
         numan=0
         do
            numan=numan+1
            s=s+1
            ls=ls+1
            b=b+1
            !bsite(1,b)=s
            !bsite(2,b)=s+1
            bsite(1,b)=ls
            bsite(2,b)=ls+1
            !print*,bsite(1,b),bsite(2,b),'ra'
            if (numan .eq. rightan-1) then
               num_end=num_end+1
               endsite(num_end)=ls+1
              ! print*,'r annt',ls+1
               if (k <= 3) then
                  s=s+3
                  ls=ls+1
               end if 
               exit
            end if
         end do
         numan=0
      else
         numan=numan+1
         b=b+1
         !bsite(1,b)=s
         !bsite(2,b)=s+1
         bsite(1,b)=ls
         bsite(2,b)=ls+1
         !print*,bsite(1,b),bsite(2,b),'ot'
         !if (record_la .eq. 0) then
         !   print*,'left annt',ls
         !   record_la=ls+1
         !else if (record_la .eq. ls) then
         !   record_la=ls+1
         !else 
         !   print*,'left annt',ls
         !   record_la=ls+1
         !end if
         if (record_la .eq. ls) then
            record_la=ls+1
         else 
            num_end=num_end+1
            endsite(num_end)=ls
            !print*,'left annt',ls
            record_la=ls+1
         end if 
      end if 
   end do
   !print*,'b=',b,'rightan',rightan

   if(b/=nb) then
      write(*,*) "error: nb"
   endif

   !print*,endsite
   !new dsite
   d=0
   do s=1,2*lchain
      !s=3
      numd=0
      allocate(dseed(1:numd))
      do b=1,nb
         s1=bsite(1,b)
         s2=bsite(2,b)
         if (s1 .eq. s) then
            allocate(oldseed(numd))
            oldseed=dseed
            deallocate(dseed)
            numd=numd+1
            allocate(dseed(1:numd))
            dseed(1:numd-1)=oldseed(1:numd-1)
            dseed(numd)=s2
            deallocate(oldseed)
         else if (s2 .eq. s) then
            allocate(oldseed(numd))
            oldseed=dseed
            deallocate(dseed)
            numd=numd+1
            allocate(dseed(1:numd))
            dseed(1:numd-1)=oldseed(1:numd-1)
            dseed(numd)=s1
            deallocate(oldseed)
         end if 
      end do
      !print*,'dseed',dseed,'numd',numd
      do x=1,numd
         do y=x+1,numd
            d=d+1
            dsite(1,d)=dseed(x)
            dsite(2,d)=dseed(y)
            !print*,'dsite',d,dsite(1,d),dsite(2,d)
         end do
      end do
      deallocate(dseed)
   end do

   if(d/=nd) then
      write(*,*) "error: nd"
   endif
   !do i=1,nb
   !   print*,bsite(1,i),bsite(2,i)
   !end do

end subroutine makelattice
!--------------------------!
subroutine dijkstra(rank)
   use system; implicit none

   integer :: rank
   character(len=1024)  :: numasstring
   real(8) :: rgprob

  ! Constants
   integer, parameter :: infinity = 100000

  ! Variables
  integer :: n, m, s, t, u, v, w, i, j,  next_node
  real(8), allocatable :: visited(:), distance(:), previous(:)
  integer :: starts,ends
  real(8), allocatable :: graph(:,:)
  real(8) :: min_distance

  write(numasstring,'(I3)') rank
  !Note this is the----- chain case --------
  ! Read in the graph
  allocate(graph(1:2*lchain,1:2*lchain))
  allocate(visited(1:2*lchain))
  allocate(distance(1:2*lchain))
  allocate(previous(1:2*lchain))
  graph = infinity
  do i=1,nb
     u=bsite(1,i)
     v=bsite(2,i)
     graph(u,v)=1
     graph(v,u)=graph(u,v)
     !print*,u,v,scor(i)
  end do

  
  ! Read in the start and end nodes
  ! open chain
  !starts=6
  !ends=24
  do starts=1,2*lchain
     do ends=starts,2*lchain
        ! Initialize the visited, distance, and previous arrays
        do i = 1, 2*lchain
           visited(i) = 0
           distance(i) = infinity
           previous(i) = -1
        end do
        distance(starts) = 0

        ! Find the shortest path
        do i = 1, 2*lchain
           min_distance = infinity
           do j = 1, 2*lchain 
              if (visited(j)==0 .and. distance(j)<min_distance) then
                 min_distance = distance(j)
                 next_node = j
              end if
           end do
           visited(next_node) = 1
           if (min_distance == infinity) then
              exit
           end if
           do j = 1, 2*lchain
              if (graph(next_node,j) /= infinity) then
                 if (distance(next_node) + graph(next_node,j) < distance(j)) then
                    distance(j) = distance(next_node) + graph(next_node,j)
                    previous(j) = next_node
                 end if
              end if
           end do
        end do
        hamp(starts,ends)=int(distance(ends))
        hamp(ends,starts)=hamp(starts,ends)
        !print*,starts,'->',ends,'hamp=',hamp(starts,ends),hamp(ends,starts)
     end do
  end do
  dscale=hamp
  !print*,hamp(37,43),37,43,dscale(37,43)
  !print*,hamp(37,44),37,44,dscale(37,44)
  !print*,hamp(37,50),37,50,dscale(37,50)
  !print*,hamp(37,57),37,57,dscale(37,57)
  !print*,hamp(37,64),37,64,dscale(37,64)
  !print*,hamp(43,50),43,50,dscale(43,50)
  !print*,hamp(43,44),43,44,dscale(43,44)
  hamp=1/hamp**2
  !print*,hamp(1,24)

  !do i=1,7
  !   do j=i+1,8
  !      u=endsite(i)
  !      v=endsite(j)
  !      print*,u,v,dscale(u,v)
  !   end do
  !end do
  
  ! Print the shortest path and distance
  !print*, distance(ends), ends
  !i = ends
  !do while (previous(i) /= -1)
  !  write(*,*) i
  !  i = previous(i)
  !end do
  !write(*,*) starts


  deallocate(graph)
  deallocate(visited)
  deallocate(distance)
  deallocate(previous)
end subroutine dijkstra
!--------------------------!
subroutine allocateall()
   !--------------------------!
   use system; use mdata; implicit none

   !integer :: i,j,b,s1,s2,ml,init
   integer :: s,x1,x2,y1,y2
   real(8) :: ran 
   external :: ran

   if(ly0==0) then
      ly=1
      lly=1
      nn=lx
   endif 
   if (ly0==1) then
      ly=1
      lly=1
      nn=lx
   endif 
   if (ly0==2) then    ! one crossing; lx must be a multiple of 4  
      ly=2
      lly=lx+1
      nn=ly*lx          ! the total number of the spins = 2*lx
   end if 
   if (ly0==4) then    ! four crossing; lx must be a multiple of 8  
      !ly=4
      lly=lx+2
      !nn=ly*lx          ! the total number of the spins = 4*lx
      nn=2*lchain-8
   endif   


   mm=2**mm
   !mm=mm*(2*lx) 
   mm=mm*nn
   !write(*,*) "m:", mm

   allocate(ss(0:7))
   allocate(oper(0:mm))
   allocate(vrtx(0:4*mm+4*nn))
   allocate(ab(nn))
   allocate(lspn(nn))
   allocate(rspn(nn))
   allocate(lbnd(nn))
   allocate(rbnd(nn))
   allocate(plbnd(nn))
   allocate(prbnd(nn))
   allocate(frst(nn))
   allocate(last(nn))
   allocate(loopnmbr(nn))
   allocate(loopsize(nn/2))

   allocate(pos(2,nn))
   allocate(ipos(0:lx-1,0:lly-1))
   allocate(xy(2,nn))
   allocate(xyi(0:lx-1,0:ly-1))

   !open(99,file='d.dat',status='unknown')
   !do s=1,nn
   !   x1=mod(s-1,lx)
   !   y1=(s-1)/lx
   !   xy(1,s)=x1           ! x coordinate
   !   xy(2,s)=y1           ! y coordinate
   !   xyi(x1,y1)=s         ! site indices
   !   write(99,*)x1,y1
   !enddo
   !close(99)

   !if(ly0<2) then
   !   pos=xy
   !   ipos=xyi
   !endif

   if (ly==1) then        ! 1d open chain
      nb=lx-1
      nd=lx-2
   elseif (ly0==2) then    ! 2 intersecting chains
      nb=2*(lx-1)+1
      !nd=lx+(lx-2)+(lx-3)  !!
      !nd=lx
      nd=(lx-2)+(lly-2)
   elseif (ly0==4) then    ! 4 intersecting chains
      !nb=4*lx
      nb=2*lchain-8
      nd=2*(lx-2)+2*(lly-2)+16
   endif

   allocate(bsite(2,nb))
   allocate(dsite(2,nd))

   allocate(jj0(1:nb))
   allocate(jj(1:nb))
   allocate(pj(1:nb))
   !allocate(hamp(0:lx/2,0:ly/2))
   !allocate(hamp(0:lx,0:ly))
   !allocate(hamp(0:lx,0:lly)) 
   allocate(hamp(1:2*lchain,1:2*lchain)) 
   allocate(dscale(1:2*lchain,1:2*lchain)) 

   allocate(scor(1:1)); scor=0.d0
   !allocate(scor(1:28)); scor=0.d0 
   allocate(eescor(1:8,1:8)); eescor=0.d0
   allocate(pscor(1:lchain/2)); pscor=0.d0
   allocate(endsite(1:8)); endsite=0


end subroutine allocateall
!--------------------------!




!--------------------------!
subroutine deallocateall()
   !--------------------------! 
   use system; use mdata; implicit none

   deallocate(ab)
   deallocate(oper)
   deallocate(lspn)
   deallocate(rspn)
   deallocate(lbnd)
   deallocate(rbnd)
   deallocate(plbnd)
   deallocate(prbnd)
   deallocate(loopnmbr)
   deallocate(loopsize)
   deallocate(frst)
   deallocate(last)
   deallocate(vrtx)
   deallocate(hamp)
   deallocate(dscale)
   deallocate(bsite)
   deallocate(dsite)
   deallocate(xyi)
   deallocate(xy)
   deallocate(pos)
   deallocate(ipos)
   deallocate(jj0)
   deallocate(jj)
   deallocate(pj)
   deallocate(scor)
   deallocate(eescor)
   deallocate(pscor)
   deallocate(ss)

end subroutine deallocateall
!----------------------------!

!----------------------!
real(8) function ran()
   !----------------------------------------------!
   ! 64-bit congruental generator                 !
   ! iran64=oran64*2862933555777941757+1013904243 !
   !----------------------------------------------!
   implicit none

   real(8)    :: dmu64
   integer(8) :: ran64,mul64,add64
   common/bran64/dmu64,ran64,mul64,add64

   ran64=ran64*mul64+add64
   ran=0.5d0+dmu64*dble(ran64)

end function ran
!----------------!

!---------------------!
subroutine initran(w,rank)
   !---------------------!
   implicit none

   integer(8) :: irmax
   integer(4) :: w,nb,b
   integer :: rank
   character(len=1024)  :: numasstring 

   real(8)    :: dmu64
   integer(8) :: ran64,mul64,add64
   common/bran64/dmu64,ran64,mul64,add64

   irmax=2_8**31
   irmax=2*(irmax**2-1)+1
   mul64=2862933555777941757_8
   add64=1013904243
   dmu64=0.5d0/dble(irmax)

   open(10,file='seed.in',status='old')
   read(10,*)ran64
   close(10)

   ran64=ran64+rank*5386
   if (w.ne.0) then
      open(10,file='seed.in',status='unknown')
      write(10,*)abs((ran64*mul64)/5+5265361)
      close(10)
   endif

end subroutine initran
!----------------------!
