!--------------------!
module configuration
   !----------------------------------------------!
   ! Most important parameters and data structures
   !----------------------------------------------!
   use mpi
   save

   integer :: comm,ierror  !MPI
   integer :: lx       ! system length in x direction
   integer :: ly       ! system length in x direction
   integer :: nn       ! number of sites; nn=lx*ly
   integer :: nb       ! number of bonds (depends on lx,ly and boundary conditions)
   integer :: nh       ! number of H-operators in string
   integer :: mm       ! maximum string length (self-determined)
   integer :: dd=1    ! disorder strength
   integer :: door    !when door=0 open constructing bsite along x

   real(8) :: beta     ! inverse temperature
   real(8) :: aprob    ! part of the acceptance probability for adding operator
   real(8) :: dprob    ! part of the acceptance probability for removing operator

   integer, allocatable :: spin(:)      ! spin state
   integer, allocatable :: bsites(:,:)  ! list of sites bsites(1,b),bsites(2,b) at bond b
   integer, allocatable :: vbsites(:,:)  
   integer, allocatable :: opstring(:)  ! operator string
   integer, allocatable :: flip(:)

   integer, allocatable :: frstspinop(:) ! first operation on each site in linked vertex list
   integer, allocatable :: lastspinop(:) ! last operation on each site in linked vertex list
   integer, allocatable :: vertexlist(:) ! list of vertex links
   real(8), allocatable :: jj(:)
   real(8), allocatable :: pj(:)
   real(8), allocatable :: jj0(:)
   real(8), allocatable :: sumspin(:)
   !real(8), allocatable :: sumspinii(:)


   integer, allocatable :: t(:)
   real(8), allocatable :: ab(:),szspin(:),mab(:)

end module configuration
!------------------------!

!----------------------!
module measurementdata
   !----------------------------------------------!
   ! Data were measured quantities are accumulated
   ! - See 'measure' for explanation of variables
   !----------------------------------------------!
   save

   real(8) :: enrg1=0.d0
   real(8) :: enrg2=0.d0
   real(8) :: amag2=0.d0 
   real(8) :: ususc=0.d0 
   !real(8) :: corete=0.d0  !correlation end-to-end
   real(8) :: locsus=0.d0
   real(8) :: mlocsus=0.d0

   !--------------------------!
end module measurementdata
!--------------------------!

!============================!
program basic_heisenberg_sse
   !=====================================================!
   ! Stochastic series QMC for the S=1/2 Heisenberg model 
   ! on a lattice wit Lx*Ly sites. 
   ! - See 'makelattice' for boundary conditions
   ! - See below for input data
   ! - See 'writeresults' for output data
   !------------------------------------------------------
   ! Anders Sandvik, 2012   
   !------------------------------------------------------
   use configuration; implicit none

   integer :: i,j,nbins,msteps,isteps
   integer :: b
   integer :: rank, numproc   !! MPI
   real(8) :: rgprob !the probability of removing rung

   ! Reads input data:
   ! lx,ly  : lattice size in x and y direction
   ! beta   : inverse dimensionless temperature J/T
   ! nbins  : Number of bins (averages written to file 'res.dat' after each bin
   ! msteps : Number of MC sweeps in each bin (measurements after each sweep)
   ! isteps : Number of MC sweeps for equilibration (no measurements)

   call MPI_INIT(ierror)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)
   call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)

   open(10,file='read.in',status='old')
   read(10,*)lx,ly,beta,dd,rgprob
   read(10,*)nbins,msteps,isteps
   close(10)

   ! initializes random number generator; needs file 'seed.in' with a seed integer
   !call initran(1)   
   call initran(1,rank) 

   ! constructs array with sites connected by bonds (defines the lattice)
   !call makelattice()


   ! probabilities used in diagonal update
!   aprob=0.5d0*beta*nb
!   dprob=1.d0/(0.5d0*beta*nb)
   allocate(vbsites(2,2*(lx-1)))
   do j=1,nbins
      ! Initialization of arrays and configuration
      call initconfig(rgprob)

      ! Do isteps equilibration sweeps
      do i=1,isteps
         call diagonalupdate()
         call linkvertices()
         call loopupdate()
         ! this subroutine increases the cutoff if needed (writes cutoff to 'cut.dat')
         call adjustcutoff(i)
      enddo

      ! Do nbins bins with msteps MC sweeps in each, measure after each
      do i=1,msteps
         call diagonalupdate()
         call linkvertices()
         call loopupdate()
         call measure()
      enddo
      ! write bin averages to file 'res.dat'
      call writeresults(msteps,rank,rgprob)
      call deallocateall()
   enddo
   deallocate (vbsites)

   !deallocate(sumspinii)

   call MPI_FINALIZE(ierror)

end program basic_heisenberg_sse
!================================!

!---------------------------!
subroutine diagonalupdate()
   !------------------------------------------!
   ! Carries out one sweep of diagonal updates
   !------------------------------------------!
   use configuration; implicit none

   integer :: i,b,op
   real(8) :: ran

   external :: ran

   do i=0,mm-1
      op=opstring(i)
      if (op==0) then       
         !b=int(ran()*nb)+1
         10 call generatebond(b)
         if (spin(bsites(1,b))/=spin(bsites(2,b))) then
            if (ran()*(mm-nh)<=aprob) then
               opstring(i)=2*b
               nh=nh+1 
            endif
         else
            goto 10
         endif
      elseif (mod(op,2)==0) then        
         if (ran()<=dprob*(mm-nh+1)) then
            opstring(i)=0
            nh=nh-1
         endif
      else
         b=op/2
         spin(bsites(1,b))=-spin(bsites(1,b))
         spin(bsites(2,b))=-spin(bsites(2,b))
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

!-------------------------!
subroutine linkvertices()
   !------------------------------!
   ! Makes the linkes vertex list.
   !------------------------------!
   use configuration; implicit none

   integer :: b,op,s1,s2,v0,v1,v2

   frstspinop(:)=-1
   lastspinop(:)=-1
   do v0=0,4*mm-1,4
      op=opstring(v0/4)
      if (op/=0) then
         b=op/2
         s1=bsites(1,b)
         s2=bsites(2,b)
         v1=lastspinop(s1)
         v2=lastspinop(s2)
         if (v1/=-1) then
            vertexlist(v1)=v0
            vertexlist(v0)=v1
         else
            frstspinop(s1)=v0
         endif
         if (v2/=-1) then
            vertexlist(v2)=v0+1
            vertexlist(v0+1)=v2
         else
            frstspinop(s2)=v0+1
         endif
         lastspinop(s1)=v0+2
         lastspinop(s2)=v0+3
      else
         vertexlist(v0:v0+3)=-1
      endif
   enddo
   do s1=1,nn
      v1=frstspinop(s1)
      if (v1/=-1) then
         v2=lastspinop(s1)
         vertexlist(v2)=v1
         vertexlist(v1)=v2
      endif
   enddo

end subroutine linkvertices
!---------------------------!

!-----------------------!
subroutine loopupdate()
   !----------------------------------------------------------------------------!
   ! Carries out loop updates; each loop is constructed and flipped with
   ! probability 1/2. At the end, spins not connected to any operator
   ! (corresponding to purely time-like loops) are flipped with probability 1/2
   !----------------------------------------------------------------------------!
   use configuration; implicit none

   integer :: i,v0,v1,v2
   real(8) :: ran
   external :: ran

   do v0=0,4*mm-1,2       
      if (vertexlist(v0)<0) cycle
      v1=v0
      if (ran()<0.5d0) then
         do 
            opstring(v1/4)=ieor(opstring(v1/4),1)
            vertexlist(v1)=-2
            v2=ieor(v1,1)
            v1=vertexlist(v2)
            vertexlist(v2)=-2
            if (v1==v0) exit
         enddo
      else
         do 
            vertexlist(v1)=-1
            v2=ieor(v1,1)
            v1=vertexlist(v2)
            vertexlist(v2)=-1
            if (v1==v0) exit
         enddo
      endif
   enddo

   do i=1,nn
      if (frstspinop(i)/=-1) then
         if (vertexlist(frstspinop(i))==-2) spin(i)=-spin(i)
      else
         if (ran()<0.5) spin(i)=-spin(i)
      endif
   enddo

end subroutine loopupdate
!-------------------------!

!--------------------!
subroutine measure()
   !-----------------------------------------------------!
   ! Computes expectation values and accumulates them:
   ! enrg1 : energy
   ! enrg2 : squared energy for specific heat calculation
   ! amag2 : Staggered magnetization
   ! ususc : Uniform susceptibility
   !-----------------------------------------------------!

   use configuration; use measurementdata; implicit none

   integer :: i,b,op,s1,s2,am,j,am1
   real(8) :: am2
   integer :: x,b0,xflip    !calculate loc susceptibility
   real(8) :: xsus,const
   real(8) :: myres, rightres,ax1



   !---------------------loc-----------------------------!
!   !xsus=0.d0
!   flip=0
!   sumspin=0
!   sumspinii=0
!
!   do i=0,mm-1
!      op=opstring(i)
!      if (mod(op,2)==1) then        
!         b=op/2
!         s1=bsites(1,b)
!         s2=bsites(2,b)
!         flip(s1)=flip(s1)+1
!         flip(s2)=flip(s2)+1
!         do x=1,nn
!            sumspin(x)=sumspin(x)+(-1)**flip(x)*spin(x)/2.
!            sumspinii(x)=sumspinii(x)+1./4.
!         end do
!      else 
!         do x=1,nn
!            sumspin(x)=sumspin(x)+(-1)**flip(x)*spin(x)/2.
!            sumspinii(x)=sumspinii(x)+1./4.
!         end do
!      endif
!   end do
!
!   const=beta/mm/(mm+1)
!
!   myres=0
!   do x=1,nn
!      locsus=locsus+const*(sumspin(x)**2+sumspinii(x))/nn
!      myres=myres+const*(sumspin(x)**2+sumspinii(x))/nn
!   end do


   !end--------------------------------------------------!
   !------------------------------------------------!
   !---code which I adjust comes from Prof. Lin---!
   am=0
   do i=1,nn
      am=am+spin(i)*(-1)**(mod(i-1,lx)+(i-1)/lx)
   enddo      
   am=am/2
   am2=0.d0
   sumspin=0.d0
   t=-1
   do i=0,mm-1
      op=opstring(i)
      !  if (op==0) cycle
      if (mod(op,2)==1) then        
         b=op/2
         s1=bsites(1,b)
         s2=bsites(2,b)
         sumspin(s1)=sumspin(s1)+(i-t(s1))*spin(s1)
         sumspin(s2)=sumspin(s2)+(i-t(s2))*spin(s2)
         t(s1)=i
         t(s2)=i
         spin(s1)=-spin(s1)
         spin(s2)=-spin(s2)
         am=am+2*spin(s1)*(-1)**(mod(s1-1,lx)+(s1-1)/lx)
      endif
      am2=am2+dfloat(am)**2
   enddo

   do j=1,nn
      sumspin(j)=0.5d0*(sumspin(j)+(mm-1-t(j))*spin(j))
   enddo

   ab=0
   mab=0
   szspin=0
   do i=1,nn
      ab(i)=ab(i)+sumspin(i)**2+0.25*mm
      szspin(i)=szspin(i)+sumspin(i)/dble(mm)
   enddo

   do i=1,nn
      !ab(i)=ab(i)*beta/mm/(mm+1)-beta*szspin(i)**2
      ab(i)=ab(i)*beta/mm/(mm+1)!-beta*szspin(i)**2
      mab(i)=-1*beta*szspin(i)**2
   enddo

   locsus=locsus+sum(ab)
   mlocsus=mlocsus+sum(mab)
   !asusc=asusc+ax1
   !rightres=sum(ab)/nn
   !print*,'r:' ,rightres,'my:',myres

   !---code which I adjust comes from Prof. Lin---!
   !-----------------------------------------------------!
   !-----------------------------------------------------!
!   am=0
!   do i=1,nn
!      am=am+spin(i)*(-1)**(mod(i-1,lx)+(i-1)/lx)
!   enddo      
!   am=am/2
!   am2=0.d0
!   do i=0,mm-1
!      op=opstring(i)
!      if (mod(op,2)==1) then        
!         b=op/2
!         s1=bsites(1,b)
!         s2=bsites(2,b)
!         spin(s1)=-spin(s1)
!         spin(s2)=-spin(s2)
!         am=am+2*spin(s1)*(-1)**(mod(s1-1,lx)+(s1-1)/lx)
!      endif
!      am2=am2+dfloat(am)**2
!      corete=corete+spin(1)*spin(nn)
!   enddo

   am2=am2/dble(mm)
!   corete=corete/dble(mm)

   !   enrg1=enrg1+dble(nh)
   !   enrg2=enrg2+dble(nh)**2
   amag2=amag2+am2
   ususc=ususc+dble(sum(spin)/2)**2
   !-----------------------------------------------------!

end subroutine measure
!----------------------!

!-------------------------------!
subroutine writeresults(msteps,rank,rgprob)
   !--------------------------------------------------------------!
   ! Writes bin results to the file 'res.dat'; four columns with
   ! 1 : energy per spin
   ! 2 : specific heat per spin
   ! 3 : squared magnetization per spin
   ! 4 : uniform susceptibility per spin
   !--------------------------------------------------------------!
   use configuration; use measurementdata; implicit none

   integer :: msteps,rank
   real(8) :: rgprob
   character(len=1024)  :: numasstring

   write(numasstring,'(I3)') rank
   

   !   enrg1=enrg1/msteps
   !   enrg2=enrg2/msteps
   amag2=amag2/msteps
   ususc=ususc/msteps
   !corete=corete/msteps
   locsus=locsus/msteps
   mlocsus=mlocsus/msteps

   !   enrg2=(enrg2-enrg1*(enrg1+1.d0))/nn
   !   enrg1=-(enrg1/(beta*nn)-0.25d0*dble(nb)/dble(nn))
   amag2=3.d0*amag2/dble(nn)**2
   ususc=beta*ususc/nn
   !corete=3*corete/nn
   locsus=locsus/nn
   mlocsus=mlocsus/nn

   !   open(10,file='res.dat',position='append')
   !   write(10,*)enrg1,enrg2,amag2,ususc,corete
   !   close(10)
   !open(10,file='res.dat',position='append')
   open(10,file='res'//trim(adjustl(numasstring)) //'.dat',status='unknown',position='append')
   write(10,'(i4,i8,i8,f18.10,3f18.10,3f18.10,f18.10)')lx,ly,int(beta),amag2,ususc,locsus,mlocsus,rgprob
   close(10)

   enrg1=0.d0
   enrg2=0.d0
   amag2=0.d0
   ususc=0.d0
   !corete=0.d0
   locsus=0.d0
   mlocsus=0.d0

end subroutine writeresults
!---------------------------!

!-----------------------------!
subroutine adjustcutoff(step)
   !--------------------------------------------------------------------!
   ! Increases the cut-off mm of the expansion in case the current mm is
   ! less than nh+nh/3 (nh being the number of H operators in the string).
   ! The arrays of size mm must be deallocated and allocated again.
   !--------------------------------------------------------------------!

   use configuration; implicit none

   integer, allocatable :: stringcopy(:)
   integer :: mmnew,step

   mmnew=nh+nh/3
   if (mmnew<=mm) return

   allocate(stringcopy(0:mm-1))
   stringcopy(:)=opstring(:)

   deallocate(opstring)
   allocate(opstring(0:mmnew-1))
   opstring(0:mm-1)=stringcopy(:)
   opstring(mm:mmnew-1)=0
   deallocate(stringcopy)

   mm=mmnew
   deallocate (vertexlist)
   allocate(vertexlist(0:4*mm-1))

   !   ! Writes the current step number and cutoff to a file
   !   open(unit=10,file='cut.dat',position='append')
   !   write(10,*)step,mm
   !   close(10)

end subroutine adjustcutoff
!---------------------------!

!-----------------------!
subroutine initconfig(rgprob)
   !--------------------------------------------------------------!
   ! Allocates the most important arrays an initializes the stored
   ! spin configuration and the empty operator string.
   !--------------------------------------------------------------!
   use configuration; implicit none

   integer :: i,b
   integer :: vnb,rnb,rank
   integer :: s
   integer :: checknb,x1,y1,x2,y2
   integer, allocatable :: remainrung(:), rungsurv(:)
   integer :: testnb
   real(8) :: rgprob !the probability of removing rung
   real(8) :: ran
   external :: ran
   integer :: numrung,nrung,rungs1

   nn=lx*ly
   allocate(spin(nn))
   do i=1,nn
      spin(i)=(-1)**(mod(i-1,lx)+(i-1)/lx)
   enddo

   mm=max(4,nn/4)
   allocate(opstring(0:mm-1))
   opstring(:)=0
   nh=0

   allocate(frstspinop(nn))
   allocate(lastspinop(nn))
   allocate(vertexlist(0:4*mm-1))
   allocate(sumspin(1:nn))
   allocate(flip(1:nn))
   allocate(t(1:nn))
   allocate(ab(1:nn))
   allocate(mab(1:nn))
   allocate(szspin(1:nn))




   vnb = 2*(lx-1)
   rnb = lx
   allocate(remainrung(rnb)); remainrung=0
   allocate(rungsurv(1:lx)); rungsurv=0
   numrung=int(lx*(1-rgprob)+0.5)
   nrung=0

   !save rung which remains
   nb=0
   !do i=1,rnb
   !   if (i .eq. 1 .or. i .eq. rnb) cycle !to let end sites no rung
   !   nb = nb+1
   !   if (ran() .gt. rgprob) then  !1-rgprob:the prob of remainrung
   !      remainrung(i)=nb
   !   end if 
   !end do
   if (numrung .lt. lx-2 .and. numrung .gt. 0) then
      do
         rungs1=int(lx*ran())
         if (rungs1 .eq. 0) cycle
         if (rungsurv(rungs1) .eq. 1) cycle
         if (rungs1 .eq. 1 .or. rungs1 .eq. lx) cycle
         rungsurv(rungs1)=1
         nrung=nrung+1
         if (nrung .eq. numrung) exit
      end do
   else if (numrung .gt. lx-2) then
      numrung = lx-2
      do rungs1=2,lx-1
         rungsurv(rungs1)=1
         nrung=nrung+1
      end do
      if (nrung /=  numrung) print*,'error'
   end if
   !---------------------

   checknb = 2*(lx-1)+lx
   ! NOTE this nb is not the finale nb!
   !nb = vnb+rnb
   nb = vnb+nrung
   allocate(bsites(2,nb))
   allocate(jj0(1:nb))
   allocate(jj(1:nb))
   allocate(pj(1:nb))
   ! NOTE this nb is not the finale nb!

   !save the bonds which don't contain rung
   !here, we just need to calculate once
   !and save them. In the every sample, 
   !we just need to call out the vbsite.
   !--------------------------------- door=0 vnb++
   !| | | | | | | | | | | | | | | | |
   !x x x x x x x x x x x x x x x x x door/=0 
   !| | | | | | | | | | | | | | | | |
   !--------------------------------- door=0 vnb++
   if (door == 0) then
      nb = 0
      do y1=0,1
         do x1=1,lx-1
            nb = nb+1
            vbsites(1,nb) = y1*lx+x1
            vbsites(2,nb) = y1*lx+x1+1
         end do
      end do
      bsites(1:2,1:vnb)=vbsites(1:2,1:vnb)
      if (nb /= vnb) print*,'error nb/=vnb'
      !print*,'open the door'
   else
      bsites(1:2,1:vnb)=vbsites(1:2,1:vnb)
   end if 
   !-------------------------------------

   !label the rungs which remain
   nb = vnb
   do i=1,rnb
      !x1 = remainrung(i)
      !if (x1 /= 0) then
      x1 = rungsurv(i)
      if (x1 == 1) then
         nb = nb+1
         !print*,nb,checknb,x1
         !bsites(1,nb) = x1
         !bsites(2,nb) = x1+lx
         bsites(1,nb) = i
         bsites(2,nb) = i+lx
      end if 
   end do 
   !in this moment, nb is real!
   !nb control all, thus, no matter how large size of allocating 
   !can't influence the result. only focus on nb


   aprob=0.5d0*beta*nb
   dprob=1.d0/(0.5d0*beta*nb)

   deallocate(remainrung)
   deallocate(rungsurv)
   door = 1

   do b=1,nb
      jj(b) = (ran())**dd
   enddo
   pj(1)=jj(1)
   !write(*,*) 1, jj(1), pj(1)
   do b=2,nb
      pj(b)=pj(b-1)+jj(b)
      !write(*,*) b, jj(b), pj(b)
   enddo
   pj=pj/pj(nb)


   !-----------------------------
   !!check
   !if (nb /= checknb) print*,'error'
   !print*,nb,checknb,(lx-(nb-2*(lx-1)))*1./lx,'rgprob=',rgprob

   !print*,'bsite'
   !print*,rgprob
   !do i=1,nb
   !   print*,bsites(1,i),bsites(2,i)
   !end do
   !print*,'-----'

end subroutine initconfig
!-------------------------!

!------------------------!
subroutine makelattice()
   !----------------------------------------------------------------------------------------!
   ! Constructs the lattice, in the form of the list of sites connected by bonds 'bsites'.
   ! There are three options here, depending on ly:
   ! ly=1  : 1D chain (periodic only in x-direction)
   ! ly=2  : 2-leg ladder (periodic only in x-direction)
   ! ly>2  : generic 2d system (periodic in both directions)
   ! Note that the system is frustrated for ly>2 and odd and the program would not give
   ! correct results in that case (unless one changes to open boundaries in the y-direction)
   ! change ly=2 and ly =1 to OBC by Ta-Jung
   !----------------------------------------------------------------------------------------!
   use configuration; implicit none

   integer :: s,x1,x2,y1,y2

   nn=lx*ly
   if (ly==1) then
      ! this is for a 1D chain (periodic only in x-direction)
      !nb=lx      !pbc
      nb=lx-1      !obc
      allocate(bsites(2,nb))
      do x1=1,lx-1
         s=x1+1
         bsites(1,x1)=x1
         bsites(2,x1)=s
      end do
   elseif (ly==2) then
      ! this is for a 2-leg ladder chain (periodic only in x-direction)
      !nb=nn+lx   !pbc
!      nb=2*(lx-1)+lx   !obc
!      allocate(bsites(2,nb))
      !allocate(sumspinii(1:nb))
      !do y1=0,ly-1
      !   !do x1=0,lx-1 !pbc
      !   do x1=0,lx-2  !obc
      !      s=1+x1+y1*lx
      !      x2=mod(x1+1,lx)
      !      y2=y1
      !      bsites(1,s)=s
      !      bsites(2,s)=1+x2+y2*lx
      !      if (y1==0) then
      !         x2=x1
      !         y2=mod(y1+1,ly)
      !         bsites(1,s+nn)=s
      !         bsites(2,s+nn)=1+x2+y2*lx       
      !      endif
      !   enddo
      !enddo
!      s=0
!      do y1=0,ly-1
!         do x1=1,lx-1
!            s=s+1
!            bsites(1,s)=y1*lx+x1
!            bsites(2,s)=y1*lx+x1+1
!         end do
!      end do
!      do x1=1,lx
!         s=s+1
!         bsites(1,s)=x1
!         bsites(2,s)=x1+lx
!      end do
!      if (s /= nb) print*, 'error'
   else
      ! this is for a general fully periodic 2D system
      nb=2*nn
      allocate(bsites(2,nb))
      do y1=0,ly-1
         do x1=0,lx-1
            s=1+x1+y1*lx
            x2=mod(x1+1,lx)
            y2=y1
            bsites(1,s)=s
            bsites(2,s)=1+x2+y2*lx
            x2=x1
            y2=mod(y1+1,ly)
            bsites(1,s+nn)=s
            bsites(2,s+nn)=1+x2+y2*lx       
         enddo
      enddo
   endif

   !do y1=1,nb
   !   print*,bsites(1,y1),bsites(2,y1)
   !end do


end subroutine makelattice
!--------------------------!

!--------------------------!
subroutine deallocateall()
   !--------------------------! 
   use configuration; implicit none

   deallocate (spin)
   !deallocate (bsites)
   deallocate (opstring)
   deallocate (frstspinop)
   deallocate (lastspinop)
   deallocate (vertexlist)
   deallocate(sumspin)
   deallocate(flip)
   deallocate(t)
   deallocate(ab)
   deallocate(mab)
   deallocate(szspin)
   deallocate(jj0)
   deallocate(jj)
   deallocate(pj)
   deallocate (bsites)

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
