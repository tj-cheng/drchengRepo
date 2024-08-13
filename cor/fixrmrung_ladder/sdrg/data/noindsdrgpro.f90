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
   integer :: lx
   integer :: ly
   integer :: ly0
   integer :: lly    ! length of y in the real lattice 
   integer :: nn
   integer :: nb
   integer :: nd
   integer :: dd=1    ! disorder strength


   integer, allocatable :: ab(:)   ! A-/B- subsystems
   real(8), allocatable :: lspn(:)
   integer, allocatable :: pos(:,:)  ! positions in the real lattice
   integer, allocatable :: ipos(:,:)
   integer, allocatable :: xy(:,:)
   integer, allocatable :: xyi(:,:)

   integer, allocatable :: ss(:) ! 4 end-spins

   real(8), allocatable :: jj(:),jj0(:)
   integer, allocatable :: bsite(:,:)
   real(8), allocatable :: scor(:)  !log(C(L))
   integer, allocatable :: bondsurv(:),spinsurv(:)
   integer, allocatable :: spair(:,:),bpair(:,:),indexpair(:,:)
   integer, allocatable :: sglabel(:,:)
   integer :: killbond !when npair =0 happen at open boundary spin pair
   !indexpair store bsite(index,:) so index=1 or 2 
   !index is the selected spin e.g. bsite(1,b)=s1,bsite(2,b)=s99 and we set index=1
   integer :: npair
   integer, allocatable :: initialbsite(:,:)
   real(8), allocatable :: initiallspn(:)
   integer :: initialnb
   real(8) :: maxj
   integer :: fbond,afbond
   integer :: alpha   !jj=jj**alpha

end module system
!-----------------!


!================!
program probasic
   !=====================!
   !---------------------!
   use system; implicit none

   integer :: i,j,init,samples
   integer :: rank, numproc   !! MPI
   integer :: isize,num_size
   integer :: rgstep
   real(8) :: rgprob,lrgprob

   call MPI_INIT(ierror)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,numproc,ierror)
   call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)

   open(10,file='read.in',status='old')
   read(10,*)lx,ly0,lrgprob
   read(10,*)init,samples
   close(10)

   rgprob=1-(lrgprob/lx)
   !print*,rgprob
   alpha=1
   num_size = 0
   call initran(1,rank) 
   call allocateall()
   call makelattice()

   if(ly0 == 2) then 
      do i=1,samples
         sglabel=0
         rgstep=0
         nb=initialnb
         call initsystem(init,rgprob)
         lspn=initiallspn
         scor=-999
         init=1
         do
            rgstep=rgstep+1
            call sdrgprocedure()
         !   !--- num F bonds and AF bonds ---!
         !   if (sum(spinsurv) .gt. nn*0.03) then
         !      call measure()
         !      call writebindata(rgstep,rank,rgprob)
         !   end if
            !--------------------------------!
            if (sum(bondsurv) .eq. 0) exit
         end do
         ! --- spins correlation --- !
         call dijkstra(rank,lrgprob)
         ! ------------------------- !
         deallocate(bsite)
         deallocate(jj)
         deallocate(scor)
         deallocate(bondsurv)
         !---------------------------!
      end do
   else
      do i=1,samples
         rgstep=0
         nb=initialnb
         call initsystem(init,rgprob)
         lspn=initiallspn
         bsite=initialbsite
         scor=-999
         init=1
         do
            rgstep=rgstep+1
            call sdrgprocedure()
         !   !--- num F bonds and AF bonds ---!
         !   if (sum(spinsurv) .gt. nn*0.03) then
         !      call measure()
         !      call writebindata(rgstep,rank,rgprob)
         !   end if
            !--------------------------------!
            if (sum(bondsurv) .eq. 0) exit
         end do
         ! --- spins correlation --- !
         call dijkstra(rank,rgprob)
         ! ------------------------- !
         deallocate(bsite)
         deallocate(jj)
         deallocate(scor)
         deallocate(bondsurv)
         allocate(bsite(2,initialnb))
         allocate(jj(1:initialnb))
         allocate(scor(initialnb))
         allocate(bondsurv(initialnb))
         !---------------------------!
      end do
   end if



   deallocate(sglabel)
   call deallocateall()
   call MPI_FINALIZE(ierror)
end program probasic
!====================!

!----------------------------!

!--------------------!
subroutine sdrgprocedure()
   use system; implicit none
   integer :: nbmaxj

   call findmaxbond(nbmaxj)
   !print*,jj(nbmaxj), nbmaxj
   call construct_effective_coup(nbmaxj)

end subroutine sdrgprocedure
!--------------------!
!---------------------------!
subroutine construct_effective_coup(nbmaxj)
   use system; implicit none
   integer :: nbmaxj
   integer :: s1,s2,bp1,bp2,sp1,sp2,index1,index2
   integer :: cons1,cons2
   integer :: i,j,inib
   integer :: coeffsign
   real(8) :: omega,effspin,efjj,coeff,spin1,spin2
   integer :: case1,case2

   !call find_new_bond_pair(nbmaxj)
   omega=jj(nbmaxj)
   maxj=omega
   s1=bsite(1,nbmaxj)
   s2=bsite(2,nbmaxj)
   if (abs(lspn(s1)) .lt. abs(lspn(s2))) then
      s2=bsite(1,nbmaxj)
      s1=bsite(2,nbmaxj)
   end if 

   case1=0
   case2=0
   spin1=abs(lspn(s1))
   spin2=abs(lspn(s2))
   !if (spin1 .eq. 0 .or. spin2 .eq. 0) print*,'spin1=0 or spin2=0 error',spin1,spin2,'s1,s2',s1,s2
   !---AF---!
   if (omega .gt. 0) then
      !print*,omega
      if (spin1 .gt. spin2) then
         !print*,'AF RG',sum(bondsurv)
         !------------------------------!
         !case1:   need to set case1==1 , case2=0
         !find s3  ; ==== : maxj  s1>s2
         !        s3                s3
         !       / \               /. 
         !      /   \     --->    /.       
         !     /     \           /.      
         !  s1 ====== s2      effspin 
         !S1 and S2 belong to different sublattices 
         !coupled antiferro.
         !S1 >> S2 --> effective spin effs=S1-S2
         !3 is a connection of S1 or S2
         !effj3=j13+(j13-j23)*(S2/(S+1))
         !form effective spin
         !------------------------------!
         case1=1
         case2=0
         call find_new_bond_pair(nbmaxj,case1,case2)
         bondsurv(nbmaxj)=0
         !if (ab(s1) .eq. ab(s2)) print*,'error, we only care at cross4'
         effspin=spin1-spin2
         lspn(s1)=effspin
         lspn(s2)=0.d0
         scor(nbmaxj)=0
         spinsurv(s2)=0
         call RG_effspin(s1,s2,spin1,spin2,effspin,omega,nbmaxj)
         spinsurv(s2)=0
         deallocate(spair)
         deallocate(bpair)
         deallocate(indexpair)
      else if (spin1 .eq. spin2) then
         !print*,'SINGLET',sum(bondsurv)
         !------------------------------!
         !case2:   need to set case1==0 , case2=1
         !find s3,s4   ; == : maxj
         !     s3        s4  
         !      \       /    
         !       \     /       
         !        \   /               
         !         s1 === s2                   
         !                       --->   s3 -.-.-. s4            
         !     s3     s4                         
         !      \       \                 
         !       \       \
         !        \       \
         !         s1 === s2
         !S1=S2, antiferro
         !form singlet
         !------------------------------!
         case1=0
         case2=1
         call find_new_bond_pair(nbmaxj,case1,case2)
         bondsurv(nbmaxj)=0
         spinsurv(s1)=0
         spinsurv(s2)=0
         lspn(s1)=0
         lspn(s2)=0
         !scor(nbmaxj)=log(3./4.)
         if (spin1 .eq. 0.5) then
            sglabel(s1,s2)=1
            sglabel(s2,s1)=1
         end if 
         scor(nbmaxj)=0
         if (npair .eq. 0) then
            if (killbond /= -1) then
               bondsurv(killbond)=0
               !print*,killbond,'kill'
            end if
         else 
            call RG_singlet(omega,spin1)
         end if 
         deallocate(spair)
         deallocate(bpair)
         deallocate(indexpair)
      end if 

      !---F---!
   else if (omega .lt. 0) then
      !print*,'F RG',sum(bondsurv)
      !print*,omega
      !------------------------------!
      !case3   find pair way is same as case1 -> case1=1, case2=0
      !S1 and S2 = same subllatices
      !coupled ferro.
      !effs=S1+S2  put on initial S1 !I am not sure the location.
      !sign(j13)=sign(j23)
      !effj3=j13*c(S1,S2,S)+j23*c(S2,S1,S)
      !------------------------------!
      case1=1
      case2=0
      call find_new_bond_pair(nbmaxj,case1,case2)
      bondsurv(nbmaxj)=0
      !if (ab(s1) /= ab(s2)) print*,'error,we only care at cross4'
      effspin=spin1+spin2
      lspn(s1)=effspin
      lspn(s2)=0.d0
      spinsurv(s2)=0
      scor(nbmaxj)=0
      call RG_effspin(s1,s2,spin1,spin2,effspin,omega,nbmaxj)
      spinsurv(s2)=0
      deallocate(spair)
      deallocate(bpair)
      deallocate(indexpair)
   end if 

end subroutine construct_effective_coup
!---------------------------!
subroutine adjust_bsite_jj_bondsurv()
   use system; implicit none
      integer, allocatable :: oldbsite(:,:)
      integer, allocatable :: oldbondsurv(:)
      real(8), allocatable :: oldjj(:),oldscor(:)
      allocate(oldbsite(2,nb))
      allocate(oldjj(nb))
      allocate(oldscor(nb))
      allocate(oldbondsurv(nb))
      !oldbsite(:,1:nb)=bsite(:,1:nb)
      oldbsite=bsite
      oldscor=scor
      oldjj=jj
      oldbondsurv=bondsurv
      nb=nb+1
      deallocate(bsite)
      deallocate(scor)
      deallocate(jj)
      deallocate(bondsurv)
      allocate(bsite(2,nb))
      allocate(scor(nb))
      allocate(jj(nb))
      allocate(bondsurv(nb))
      bsite(:,1:nb-1)=oldbsite(:,1:nb-1)
      scor(1:nb-1)=oldscor(1:nb-1)
      jj(1:nb-1)=oldjj(1:nb-1)
      bondsurv(1:nb-1)=oldbondsurv(1:nb-1)
      bondsurv(nb)=1
      jj(nb)=0
      deallocate(oldscor)
      deallocate(oldbsite)
      deallocate(oldjj)
      deallocate(oldbondsurv)
end subroutine adjust_bsite_jj_bondsurv
!--------------------!
subroutine findmaxbond(nbmaxj)
   use system; implicit none
   integer :: i,nbmaxj,g1,g2
   real(8) :: oldmaxj,gap,oldgap,gs1,gs2

   gap=0
   maxj=0
   oldgap=0
   gs1=0
   gs2=0
   do i=1,nb
      if (bondsurv(i) /= 1) cycle
      g1=bsite(1,i)
      g2=bsite(2,i)
      gs1=abs(lspn(g1))
      gs2=abs(lspn(g2))
      if (jj(i) .lt. 0) then
         oldgap=abs(jj(i))*(gs1+gs2)
      else if (jj(i) .gt. 0) then
         oldgap=jj(i)*(abs(gs1-gs2)+1)
      end if 
      !oldmaxj=abs(jj(i))
      !if (oldmaxj .gt. maxj) then
      if (oldgap .gt. gap) then
         !maxj=oldmaxj
         gap=oldgap
         nbmaxj=i
      end if 
   end do
   !g1=bsite(1,nbmaxj)
   !g2=bsite(2,nbmaxj)
   !gs1=abs(lspn(g1))
   !gs2=abs(lspn(g2))
   !print*,jj(nbmaxj),gap,gs1,gs2

end subroutine findmaxbond
!--------------------!
!--------------------!
!--------------------------!
subroutine find_new_bond_pair(nbmaxj,case1,case2)
   use system; implicit none
   integer :: nbmaxj
   integer, allocatable :: lists1(:),lists2(:)
   integer, allocatable :: bondlists1(:),bondlists2(:)
   integer, allocatable :: totalslist(:),totalbondlist(:)
   integer, allocatable :: index1(:),index2(:)
   integer, allocatable :: indexbondlist(:) ! 1 or 2
   integer :: ncons1,ncons2,b
   integer :: s1,s2,s3,s4
   integer :: ns1,ns2,nss
   integer :: cons1,cons2,counts3
   integer :: i,j,k
   integer :: case1,case2,b1,bindex1,b2,bindex2
   real(8) :: omega
   !----
   ! here, s1 and s2 mean the sites of maxbond
   ! lists1 contains total connected sites with s1, lists2 ...
   ! bondlists1 stores total connected bonds with s1 instead of s1, ...
   ! ns1 is the total number of elements in the lists1
   !----
   omega=jj(nbmaxj)
   ncons1=0
   ncons2=0
   ns1=0
   ns2=0

   allocate(lists1(nn))
   allocate(lists2(nn))
   allocate(bondlists1(nn))
   allocate(bondlists2(nn))
   allocate(index1(nn))
   allocate(index2(nn))

   lists1=-1
   lists2=-1
   bondlists1=-1
   bondlists2=-1
   index1=-1
   index2=-1

   s1=bsite(1,nbmaxj)
   s2=bsite(2,nbmaxj)

   !print*,'nb=',nb
   do b=1,nb
      if(bondsurv(b) /= 1) cycle
      if(spinsurv(bsite(1,b)) .eq. 0 .or. spinsurv(bsite(2,b)) .eq. 0) print*,'surving problem'
      if (jj(b) .eq. 0) print*,'error jjb',jj(nbmaxj)
      !---connection of s1---!
      if (bsite(1,b) == s1) then
         if (bsite(2,b) /= s2) then
            ncons1=ncons1+1
            lists1(ncons1)=bsite(2,b)
            bondlists1(ncons1)=b
            ns1=ns1+1
            index1(ncons1)=1
         end if
      else if (bsite(2,b) == s1) then
         ncons1=ncons1+1
         lists1(ncons1)=bsite(1,b)
         bondlists1(ncons1)=b
         ns1=ns1+1
         index1(ncons1)=2
      end if 
      !---connection of s2---!
      if (bsite(2,b) == s2) then
         if (bsite(1,b) /= s1) then
            ncons2=ncons2+1
            lists2(ncons2)=bsite(1,b)
            bondlists2(ncons2)=b
            ns2=ns2+1
            index2(ncons2)=2
         end if
      else if (bsite(1,b) == s2) then
         ncons2=ncons2+1
         lists2(ncons2)=bsite(2,b)
         bondlists2(ncons2)=b
         ns2=ns2+1
         index2(ncons2)=1
      end if 
   end do


   nss=ns1+ns2
   allocate(totalslist(nss))
   allocate(totalbondlist(nss))
   allocate(indexbondlist(nss))

   do i=1,ns1
      totalslist(i)=lists1(i)
      totalbondlist(i)=bondlists1(i)
      indexbondlist(i)=index1(i)
   end do

   do i=1,ns2
      totalslist(ns1+i)=lists2(i)
      totalbondlist(ns1+i)=bondlists2(i)
      indexbondlist(ns1+i)=index2(i)
   end do

   
   
   npair=0
   allocate(spair(2,(nss*(nss-1))/2))
   allocate(bpair(2,(nss*(nss-1))/2))
   allocate(indexpair(2,(nss*(nss-1))/2))
   spair=-1
   bpair=-1
   indexpair=-1
   !---case1---!
   if (case1 .eq. 1 .and. case2 .eq. 0) then
      if (nss .eq. 1) then
         npair=1
         bpair(1,npair)=totalbondlist(1)
         bpair(2,npair)=0
         spair(1,npair)=totalslist(1)
         spair(2,npair)=0
         indexpair(1,npair)=indexbondlist(1)
         indexpair(2,npair)=0
      else
         do i=1,nss
            b1=totalbondlist(i)
            bindex1=indexbondlist(i)
            bindex1=2-mod(bindex1+1,2)
            s3=bsite(bindex1,b1)    ! s3
            do j=i+1,nss
               b2=totalbondlist(j)
               bindex2=indexbondlist(j)
               bindex2=2-mod(bindex2+1,2)
               s4=bsite(bindex2,b2)    !s3 
               if (s3 .eq. s4) then
                  ! find s3 
                  !        s3      
                  !       / \      
                  !      /   \        
                  !     /     \        
                  !  s1 ====== s2   
                  npair=npair+1
                  spair(1,npair)=s3
                  spair(2,npair)=s4
                  bpair(1,npair)=b1
                  bpair(2,npair)=b2
                  bindex1=2-mod(bindex1+1,2)  !s1 or s2
                  bindex2=2-mod(bindex2+1,2)  !s1 or s2
                  indexpair(1,npair)=bindex1  !note here store the index of s1 or s2 same as case2
                  indexpair(2,npair)=bindex2
               end if 
            end do
         end do
         do i=1,nss
            ! find s3 
            !        s3      
            !       /       
            !      /           
            !     /             
            !  s1 ====== s2   
            b1=totalbondlist(i)
            bindex1=indexbondlist(i)
            bindex1=2-mod(bindex1+1,2)
            s3=bsite(bindex1,b1)    ! s3
            counts3=0
            do j=1,nss
               if (s3 .eq. totalslist(j)) counts3=counts3+1
            end do
            if (counts3 .eq. 1) then
               npair=npair+1
               spair(1,npair)=s3
               spair(2,npair)=0
               bpair(1,npair)=b1
               bpair(2,npair)=0
               bindex1=2-mod(bindex1+1,2)
               indexpair(1,npair)=bindex1  !note here store the index of s1 or s2 same as case2
               indexpair(2,npair)=0
            end if 
         end do
      end if
   end if 
   !---case1---!

   !---case2---!
   if (case2 .eq. 1 .and. case1 .eq. 0) then
      npair=(nss*(nss-1))/2

      k=0
      do i=1,nss
         if (i == nss) exit
         do j=i+1,nss
            k=k+1
            spair(1,k)=totalslist(i)
            spair(2,k)=totalslist(j)
            bpair(1,k)=totalbondlist(i)
            bpair(2,k)=totalbondlist(j)
            indexpair(1,k)=indexbondlist(i)  !note here store the index of s1 or s2
            indexpair(2,k)=indexbondlist(j)
         end do
      end do
      if (npair .eq. 0) then
         killbond=totalbondlist(1)
         if (killbond .gt. nb .or. killbond .lt. 1) then 
            killbond=-1
         else
            scor(killbond) = log(abs(jj(killbond)/omega))
         end if 
      end if 
   end if 
   !---case2---!


   deallocate(lists1)
   deallocate(lists2)
   deallocate(bondlists1)
   deallocate(bondlists2)
   deallocate(totalslist)
   deallocate(totalbondlist)
end subroutine find_new_bond_pair
!---------------------------!
!---------------------------!
subroutine RG_effspin(s1,s2,spin1,spin2,effspin,omega,nbmaxj)
   use system; implicit none
   integer :: i,s1,s2
   integer :: bp1,bp2,index1,index2,cons1,cons2
   integer :: nbmaxj
   real(8) :: spin1,spin2,effspin,omega
   real(8) :: sp1,sp2

   if (omega .gt. 0) then


      do i=1,npair
         bp1=bpair(1,i)
         bp2=bpair(2,i)
         index1=indexpair(1,i)
         index2=indexpair(2,i)
         sp1=spair(1,i)
         sp2=spair(2,i)
         !if (sp1 .eq. -1 .or. sp2 .eq. -1 .or. bp1 .eq. -1 .or. bp2 .eq. -1) print*,'allocate error'
         if (bp2 .eq. 0) then
            cons1=bsite(index1,bp1)
            if (cons1 .eq. s1) then
               !s1>s2
               !     s3                s3              
               !      \                 \.          
               !       \                 \.   
               !        \                 \.
               !         s1 === s2         effs=s1-s2 
               jj(bp1)=const(spin1,spin2,effspin)*jj(bp1)
               !if (jj(bp1) .eq. 0) print*,'effjj=0'
               !print*,const(spin1,spin2,effspin),(effspin*(effspin+1)+spin1*(spin1+1)-spin2*(spin2+1))/(2.*effspin*(effspin+1))
            else if (cons1 .eq. s2) then
               !s1>s2
               !            s3                  s3                  
               !              \                 /.          
               !               \    --->       /.  adjust nb -> nb+1  
               !                \             /.    
               !         s1 === s2          effs=s1-s2
               !scor(bp1)=log(abs(jj(bp1)/omega))
               bondsurv(bp1)=0
               call adjust_bsite_jj_bondsurv() !nb -> nb+1
               bsite(1,nb)=s1                   !nb=nb+1
               bsite(2,nb)=spair(1,i)
               jj(nb)=const(spin2,spin1,effspin)*jj(bp1)
               !if (jj(nb) .eq. 0) print*,'effjj=0'
               jj(bp1)=0
            else
               print*,'error AF1'
            end if 
         else
            cons1=bsite(index1,bp1)
            cons2=bsite(index2,bp2)
            if (spair(1,i) .eq. spair(2,i)) then
               if (cons1 .eq. s1) then
                  !scor(bp2)=log(abs(jj(bp2)/omega))
                  jj(bp1)=jj(bp1)+(jj(bp1)-jj(bp2))*spin2/(effspin+1)
                  !if (jj(bp1) .eq. 0) print*,'effjj=0'
                  bondsurv(bp2)=0
                  jj(bp2)=0
                  !s1>s2
                  !        s3                   s3
                  !       / \                  /. 
                  !   bp1/   \ bp2    --->    /. bp1      
                  !     /     \              /.      
                  !  s1 ====== s2      effspin=s1-s2
               else if (cons2 .eq. s1) then
                  !scor(bp1)=log(abs(jj(bp1)/omega))
                  jj(bp2)=jj(bp2)+(jj(bp2)-jj(bp1))*spin2/(effspin+1)
                  !if (jj(bp2) .eq. 0) print*,'effjj=0'
                  bondsurv(bp1)=0
                  jj(bp1)=0
                  !s1>s2
                  !        s3                   s3
                  !       / \                  /. 
                  !   bp2/   \ bp1    --->    /. bp2      
                  !     /     \              /.      
                  !  s1 ====== s2      effspin=s1-s2
               else
                  print*,'error AF2'
               end if
            end if
         end if
      end do
   else if (omega .lt. 0) then


      do i=1,npair
         bp1=bpair(1,i)
         bp2=bpair(2,i)
         index1=indexpair(1,i)
         index2=indexpair(2,i)
         sp1=spair(1,i)
         sp2=spair(2,i)
         !if (sp1 .eq. -1 .or. sp2 .eq. -1 .or. bp1 .eq. -1 .or. bp2 .eq. -1) print*,'allocate error'
         if (bp2 .eq. 0) then
            cons1=bsite(index1,bp1)
            if (cons1 .eq. s1) then
               !s1>s2
               !     s3                s3              
               !      \                 \.          
               !       \                 \.   
               !        \                 \.
               !         s1 === s2         effs=s1+s2 
               jj(bp1)=const(spin1,spin2,effspin)*jj(bp1)
               !if (jj(bp1) .eq. 0) print*,'effjj=0'
            else if (cons1 .eq. s2) then
               !s1>s2
               !            s3                  s3                  
               !              \                 /.          
               !               \    --->       /.  adjust nb -> nb+1  
               !                \             /.    
               !         s1 === s2          effs=s1+s2
               bondsurv(bp1)=0
               !scor(bp1)=log(abs(jj(bp1)/omega))
               call adjust_bsite_jj_bondsurv()
               bsite(1,nb)=s1
               bsite(2,nb)=spair(1,i)
               if (lspn(s1) .eq. 0 .or. lspn(spair(1,i)) .eq. 0) then
                  print*,'s1',s1,lspn(s1)
                  print*,'spair(1,i)',spair(1,i),lspn(spair(1,i))
                  print*,'ajfoijsdif'
               end if 
               jj(nb)=const(spin2,spin1,effspin)*jj(bp1)
               !if (jj(nb) .eq. 0) print*,'effjj=0'
            else
               print*,'error F1',cons1,cons2
            end if 

         else
            cons1=bsite(index1,bp1)
            cons2=bsite(index2,bp2)
            !print*,'bp1,bp2',bp1,bp2
            if (spair(1,i) .eq. spair(2,i)) then
               if (cons1 .eq. s1) then
                  jj(bp1)=const(spin1,spin2,effspin)*jj(bp1)+const(spin2,spin1,effspin)*jj(bp2)
                  !if (jj(bp1) .eq. 0) print*,'effjj=0'
                  !scor(bp2)=log(abs(jj(bp2)/omega))
                  bondsurv(bp2)=0
                  !s1>s2
                  !        s3                   s3
                  !       / \                  /. 
                  !   bp1/   \ bp2    --->    /. bp1      
                  !     /     \              /.      
                  !  s1 ====== s2      effspin=s1+s2
               else if (cons2 .eq. s1) then
                  !note spin1 >= spin2
                  jj(bp2)=const(spin1,spin2,effspin)*jj(bp2)+const(spin2,spin1,effspin)*jj(bp1)
                  !if (jj(bp2) .eq. 0) print*,'effjj=0'
                  !scor(bp1)=log(abs(jj(bp1)/omega))
                  bondsurv(bp1)=0
                  !s1>s2
                  !        s3                   s3
                  !       / \                  /. 
                  !   bp2/   \ bp1    --->    /. bp2      
                  !     /     \              /.      
                  !  s1 ====== s2      effspin=s1+s2
               else
                  print*,'error F2',cons1,cons2,'s1,s2',s1,s2
               end if
            end if
         end if 
      end do
   end if 
contains 
   !---------------------------!
   real(8) function const(sp1,sp2,effsp)
      real(8), intent(in) :: sp1
      real(8), intent(in) :: sp2
      real(8), intent(in) :: effsp

      const = (effsp*(effsp+1)+sp1*(sp1+1)-sp2*(sp2+1))/(2*effsp*(effsp+1))

   end function const
   !---------------------------!
end subroutine RG_effspin
!---------------------------!
subroutine RG_singlet(omega,spin1)
   use system; implicit none
   integer :: i,j,inib
   integer :: bp1,bp2,sp1,sp2,index1,index2
   real(8) :: coeffsign,coeff,omega,spin1,efjj

   do i=1,npair
      sp1=spair(1,i)  !sp1 is the spins labels connect with s1 or s2
      sp2=spair(2,i)
      bp1=bpair(1,i)
      bp2=bpair(2,i)
      !if (sp1 .eq. -1 .or. sp2 .eq. -1 .or. bp1 .eq. -1 .or. bp2 .eq. -1) print*,'allocate error'
      index1=indexpair(1,i)
      index2=indexpair(2,i)
      scor(bp1)=log(abs(jj(bp1)/omega))
      scor(bp2)=log(abs(jj(bp2)/omega))
      if (sp1 .eq. sp2) then
         bondsurv(bp1)=0
         bondsurv(bp2)=0
      else
         !find the connection bond of sp1 and sp2
         inib=0
         do j=1,nb  !need to use the newest nb
            if (bsite(1,j) .eq. sp1 .and. bsite(2,j) .eq. sp2) then
               inib=j
               exit
            else if (bsite(1,j) .eq. sp2 .and. bsite(2,j) .eq. sp1) then
               inib=j
               exit
            end if 
         end do
         if (inib .eq. 0) then
            call adjust_bsite_jj_bondsurv()  !nb --> nb+1
            inib=nb
            bondsurv(inib)=1
            bsite(1,inib)=sp1
            bsite(2,inib)=sp2
         end if 
         !---------
         !if (bsite(1,bp1) .eq. s1 .or. bsite(2,bp1) .eq. s1)
         if (bsite(index1,bp1) .eq. bsite(index2,bp2)) then
            coeffsign=-1
         else
            coeffsign=1
         end if 
         coeff=(2*spin1*(spin1+1))/(3.*omega)
         efjj=coeffsign*coeff*jj(bp1)*jj(bp2)
         jj(inib)=jj(inib)+efjj
         !if (jj(inib) .eq. 0) print*,'effjj=0',coeff,spin1,jj(bp1),jj(bp2),'spinsurv',sum(spinsurv)

         !spinsurv(bsite(1,bp1))=0
         !spinsurv(bsite(2,bp2))=0
      end if
   end do

   do i=1,npair
      bp1=bpair(1,i)
      bp2=bpair(2,i)
      jj(bp1)=0
      jj(bp2)=0
      bondsurv(bp1)=0
      bondsurv(bp2)=0
   end do
end subroutine RG_singlet
!---------------------------!

!--------------------!
!--------------------!
subroutine measure()
   !--------------------!
   use system; implicit none

   integer :: i,j,k,k1,r,x,y,lj,nl,m2,b
   real(8) :: e1
   !integer :: ss(4)   ! 4 end-spins



   fbond=0
   afbond=0
   do b=1,nb
      if (bondsurv(b) /= 1) cycle
      if (jj(b) .gt. 0) then
         afbond=afbond+1
      else
         fbond=fbond+1
      end if
   end do

   !! two-spin correlation
   !if(ly0<2) then          ! single chain
   !   y=0
   !   j=xyi(0,y)
   !   i=xyi(lx-1,y)
   !elseif(ly0==2) then     ! 2 intersecting chains
   !   ss(0)=xyi(0,0)
   !   ss(1)=xyi(lx-1,0)
   !   ss(2)=xyi(0,1)
   !   ss(3)=xyi(lx-1,1) 

   !   !------------------
   !   ! scor1 = C(0,1)
   !   ! scor2 = C(1,2)
   !   ! scor3 = C(2,3)
   !   ! scor4 = C(3,0)
   !   ! scor5 = C(0,2)
   !   ! scor6 = C(1,3)
   !   !------------------ 
   !else if (ly0==4) then
   !   ss(0)=xyi(0,0)
   !   !print*,pos(1,ss(0)),pos(2,ss(0)),0
   !   ss(1)=xyi(lx-1,0)
   !   !print*,pos(1,ss(1)),pos(2,ss(1)),1
   !   ss(2)=xyi(0,1)
   !   !print*,pos(1,ss(2)),pos(2,ss(2)),2
   !   ss(3)=xyi(lx-1,1) 
   !   !print*,pos(1,ss(3)),pos(2,ss(3)),3
   !   ss(4)=xyi(0,3)
   !   !print*,pos(1,ss(4)),pos(2,ss(4)),4
   !   ss(5)=xyi(lx-1,3)
   !   !print*,pos(1,ss(5)),pos(2,ss(5)),5
   !   ss(6)=xyi(lx/2-1,3)
   !   !print*,pos(1,ss(6)),pos(2,ss(6)),6
   !   ss(7)=xyi(lx/2,3)
   !   !print*,pos(1,ss(7)),pos(2,ss(7)),7
   !   !------------------
   !   !(0,1)(1,2)(2,3)(3,4)(4,5)(5,6)(6,7)  7
   !   !(0,2)(1,3)(2,4)(3,5)(4,6)(5,7)       6
   !   !(0,3)(1,4)(2,5)(3,6)(4,7)            5
   !   !(0,4)(1,5)(2,6)(3,7)                 4
   !   !(0,5)(1,6)(2,7)                      3
   !   !(0,6)(1,7)                           2
   !   !(0,7)                                1
   !   !------------------ 
   !endif


end subroutine measure
!----------------------!
subroutine dijkstra(rank,lrgprob)
   use system; implicit none

   integer :: rank
   character(len=1024)  :: numasstring
   real(8) :: lrgprob


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
 ! allocate(graph(1:nn,1:nn))
 ! allocate(visited(1:nn))
 ! allocate(distance(1:nn))
 ! allocate(previous(1:nn))
 ! graph = infinity
 ! do i=1,nb
 !    u=bsite(1,i)
 !    v=bsite(2,i)
 !    graph(u,v)=abs(scor(i))
 !    graph(v,u)=graph(u,v)
 !    !print*,u,v,scor(i)
 ! end do

 ! 
 ! ! Read in the start and end nodes
 ! ! open chain
 ! starts=1
 ! ends=lx
 !  if(ly0 .lt. 4) then          ! single chain
 !     starts=1
 !     ends=lx
 !  else if (ly0 .eq. 4) then
 !     ss(0)=xyi(0,0)
 !     !print*,pos(1,ss(0)),pos(2,ss(0)),0
 !     ss(1)=xyi(lx-1,0)
 !     !print*,pos(1,ss(1)),pos(2,ss(1)),1
 !     ss(2)=xyi(0,1)
 !     !print*,pos(1,ss(2)),pos(2,ss(2)),2
 !     ss(3)=xyi(lx-1,1) 
 !     !print*,pos(1,ss(3)),pos(2,ss(3)),3
 !     ss(4)=xyi(0,3)
 !     !print*,pos(1,ss(4)),pos(2,ss(4)),4
 !     ss(5)=xyi(lx-1,3)
 !     !print*,pos(1,ss(5)),pos(2,ss(5)),5
 !     ss(6)=xyi(lx/2-1,3)
 !     !print*,pos(1,ss(6)),pos(2,ss(6)),6
 !     ss(7)=xyi(lx/2,3)
 !     !print*,pos(1,ss(7)),pos(2,ss(7)),7
 !     !------------------
 !     !(0,1)(1,2)(2,3)(3,4)(4,5)(5,6)(6,7)  7
 !     !(0,2)(1,3)(2,4)(3,5)(4,6)(5,7)       6
 !     !(0,3)(1,4)(2,5)(3,6)(4,7)            5
 !     !(0,4)(1,5)(2,6)(3,7)                 4
 !     !(0,5)(1,6)(2,7)                      3
 !     !(0,6)(1,7)                           2
 !     !(0,7)                                1
 !     !------------------ 
 !     !--(2,3)-- horizental chain ---!
 !     starts=ss(2)
 !     ends=ss(3)
 !     !print*,starts,ends
 !  endif
 ! !-----------
 ! 
 ! ! Initialize the visited, distance, and previous arrays
 ! do i = 1, nn
 !   visited(i) = 0
 !   distance(i) = infinity
 !   previous(i) = -1
 ! end do
 ! distance(starts) = 0
 ! 
 ! ! Find the shortest path
 ! do i = 1, nn
 !   min_distance = infinity
 !   do j = 1, nn
 !     if (visited(j)==0 .and. distance(j)<min_distance) then
 !       min_distance = distance(j)
 !       next_node = j
 !     end if
 !   end do
 !   visited(next_node) = 1
 !   if (min_distance == infinity) then
 !     exit
 !   end if
 !   do j = 1, nn
 !     if (graph(next_node,j) /= infinity) then
 !       if (distance(next_node) + graph(next_node,j) < distance(j)) then
 !         distance(j) = distance(next_node) + graph(next_node,j)
 !         previous(j) = next_node
 !       end if
 !     end if
 !   end do
 ! end do
 ! 
 ! ! Print the shortest path and distance
 ! !print*, distance(ends), ends
 ! !i = ends
 ! !do while (previous(i) /= -1)
 ! !  write(*,*) i
 ! !  i = previous(i)
 ! !end do
 ! !write(*,*) starts
 ! 

 ! if (ly0 .eq. 1) then
 !    open(10,file='chainlogcor'//trim(adjustl(numasstring)) //'.dat',status='unknown',position='append')
 !    write(10,*) lx, ly0, -1.d0*distance(ends)+log(3./4.)  !|C(L)|
 !    close(10)
 ! else if (ly0 .eq. 2) then
 !    open(10,file='2legladderlogcor'//trim(adjustl(numasstring)) //'.dat',status='unknown',position='append')
 !    !write(10,'(i10, i2, 10f14.10)') lx, ly0, distance(ends)+log(3./4.)  !|C(L)|
 !    write(10,*) lx, ly0, -1.d0*distance(ends)+log(3./4.), rgprob  !|C(L)|
 !    close(10)
 ! else if (ly0 .eq. 4) then
 !    open(10,file='cross4logcor'//trim(adjustl(numasstring)) //'.dat',status='unknown',position='append')
 !    write(10,*) lx, ly0, -1.d0*distance(ends)+log(3./4.)  !|C(L)|
 !    close(10)
 ! end if
 if (sglabel(1,lx) .eq. 1) then
     open(10,file='2legladderlogcor'//trim(adjustl(numasstring)) //'.dat',status='unknown',position='append')
     write(10,'(i10, i10, f20.10, f20.10)') lx, ly0, 3./4., lrgprob  !|C(L)|
     !write(10,*) lx, ly0, -1.d0*distance(ends)+log(3./4.), rgprob  !|C(L)|
     close(10)
 else 
     open(10,file='2legladderlogcor'//trim(adjustl(numasstring)) //'.dat',status='unknown',position='append')
     write(10,'(i10, i10, f20.10, f20.10)') lx, ly0, 0.d0, lrgprob  !|C(L)|
     !write(10,*) lx, ly0, -1.d0*distance(ends)+log(3./4.), rgprob  !|C(L)|
     close(10)
 end if

 if (sglabel(lx+1,2*lx) .eq. 1) then
     open(10,file='2legladderlogcor'//trim(adjustl(numasstring)) //'.dat',status='unknown',position='append')
     write(10,'(i10, i10, f20.10, f20.10)') lx, ly0, 3./4., lrgprob  !|C(L)|
     !write(10,*) lx, ly0, -1.d0*distance(ends)+log(3./4.), rgprob  !|C(L)|
     close(10)
 else 
     open(10,file='2legladderlogcor'//trim(adjustl(numasstring)) //'.dat',status='unknown',position='append')
     write(10,'(i10, i10, f20.10, f20.10)') lx, ly0, 0.d0, lrgprob  !|C(L)|
     !write(10,'(i10, i10, f20.10)') lx, ly0, 0.d0  !|C(L)|
     !write(10,*) lx, ly0, -1.d0*distance(ends)+log(3./4.), rgprob  !|C(L)|
     close(10)
 end if
 

 ! deallocate(graph)
 ! deallocate(visited)
 ! deallocate(distance)
 ! deallocate(previous)
end subroutine dijkstra
!-------------------------!
subroutine writebindata(rgstep,rank,rgprob)
   !-------------------------!
   use system; implicit none

   integer :: r,rgstep,rank
   real(8) :: rgprob
   character(len=1024)  :: numasstring

   write(numasstring,'(I3)') rank


   !print*,'rgprob',rgprob
   if (ly0 .eq. 2) then
      open(10,file='2legnumfaf'//trim(adjustl(numasstring)) //'.dat',status='unknown',position='append')
      write(10,'(f14.10,i10,i10,i10,i10,1f14.10)') abs(maxj),afbond,fbond,lx,ly,rgprob
      close(10)
   else if (ly0 .eq. 4) then
      open(10,file='cross4numfaf'//trim(adjustl(numasstring)) //'.dat',status='unknown',position='append')
      write(10,'(f14.10,i10,i10,i10,i10)') abs(maxj),afbond,fbond,lx,ly
      close(10)
   end if 
   !if(ly0<2) then
   !   write(10,'(f14.10)') scor(1)
   !elseif (ly0==2) then
   !   write(10,'(6f14.10)') scor
   !else if (ly0==4) then
   !   write(10,'(i3,28f14.10)') lx,scor
   !endif

!

end subroutine writebindata
!!---------------------------!
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
subroutine initsystem(init,rgprob)
   !---------------------------!
   use system; implicit none

   integer :: i,j,a,x,y,init,b,x0,y0
   integer :: nrung,rungs1,rungs2,numrung
   real(8) :: ran
   real(8) :: rgprob
   integer, allocatable :: rungsurv(:)

   external :: ran

   !--------------------!
   if (ly0 /= 2) then
      !-- random couplings --!
      if (ly0==0) then        ! 1D chain for test
         call readbond()
         jj = jj0
      else
         do b=1,nb
            jj(b) = (ran())**dd
         enddo
      endif

      jj=jj**alpha
   end if


   !-- operator-string -------!

   if (init==0 .and. ly0 /= 2) then

      if(ly0<2) then                    ! single chain
         do x=0,lx-1
            i=xyi(x,0)
            lspn(i)=(-1)**(xy(1,i)+xy(2,i))
            ab(i)=(-1)**(xy(1,i)+xy(2,i))
         enddo
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
      else if (ly0 == 2) then !spare 2-leg ladder
         do x=0,lx-1
            i=xyi(x,0)
            lspn(i)=(-1)**(xy(1,i)+xy(2,i))
            ab(i)=(-1)**(xy(1,i)+xy(2,i))
         enddo
         do x=0,lx-1
            i=xyi(x,1)
            lspn(i)=(-1)**(xy(1,i)+xy(2,i))
            ab(i)=(-1)**(xy(1,i)+xy(2,i))
         end do
      elseif (ly0==4) then              ! four intersecting chains
         x0=lx/2-1
         y=0
         do x=0,lx-1
            i=xyi(x,y)
            if (x .gt. lx/4-1 .and. x .lt. 3*lx/4) then
               j=xyi(x-1,y)
               lspn(i)=(-1)**(xy(1,j)+xy(2,j))
               ab(i)=(-1)**(xy(1,j)+xy(2,j)) 
            else 
               lspn(i)=(-1)**(xy(1,i)+xy(2,i))
               ab(i)=(-1)**(xy(1,i)+xy(2,i))
            end if 
         end do
         y=2
         do x=0,lx-1
            i=xyi(x,y)
            j=xyi(x,y-2)
            ab(i)=(-1)*ab(j)
            lspn(i)=(-1)*lspn(j)
         end do
         y=1
         do x=0,lx-1        
            i=xyi(x,y)
            lspn(i)=(-1)**(xy(1,i)+xy(2,i))
            ab(i)=(-1)**(xy(1,i)+xy(2,i))
         enddo 
         y=3
         do x=0,lx-1
            i=xyi(x,y)
            j=xyi(x,y-2)
            ab(i)=(-1)*ab(j)
            lspn(i)=(-1)*lspn(j)
         end do
      endif
      !!check
      !do i=1,nn
      !   print*,i,ab(i),lspn(i)
      !end do
      lspn=0.5d0*lspn
      initiallspn=lspn

   else if (ly0 == 2) then
      allocate(rungsurv(1:lx))
      !print*,rgprob
      numrung=int(lx*(1-rgprob)+0.5)
      rungsurv=0
      nrung=0
      !print*,0
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
      !print*,1

      !deallocate(bsite)
      nb=nb+numrung
      allocate(bsite(2,nb))
      bsite(:,1:nb-numrung)=initialbsite(:,1:nb-numrung)
      if (numrung .gt. 0) then
         rungs1=1
         b=initialnb
         do
            rungs1=rungs1+1
            if (rungs1 .eq. lx) exit
            if (rungsurv(rungs1) .eq. 1) then
               b=b+1
               bsite(1,b)=rungs1
               bsite(2,b)=rungs1+lx
            end if
         end do
      end if

      if (init == 0) then
         do x=0,lx-1
            i=xyi(x,0)
            lspn(i)=(-1)**(xy(1,i)+xy(2,i))
            ab(i)=(-1)**(xy(1,i)+xy(2,i))
         enddo
         do x=0,lx-1
            i=xyi(x,1)
            lspn(i)=(-1)**(xy(1,i)+xy(2,i))
            ab(i)=(-1)**(xy(1,i)+xy(2,i))
         end do
         lspn=0.5d0*lspn
         initiallspn=lspn
      end if
      allocate(jj(1:nb))
      allocate(bondsurv(1:nb))
      allocate(scor(nb))
      bondsurv=1
      spinsurv=1
      scor=-999
      do b=1,nb
         jj(b) = (ran())**dd
      enddo
      jj=jj**alpha
      deallocate(rungsurv)
   endif

   !----surving label---!
   spinsurv=1
   bondsurv=1
   !do b=1,nb
   !   print*,bsite(1,b), bsite(2,b)
   !end do


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

   integer :: s,a,x1,x2,y1,y2,x,y,b,d,x0,y0,s0,s1,s2  ! (x0,y0),(s0,s1,s2) the intersection 
   integer, allocatable :: xx0(:),ss0(:),ss1(:),ss2(:),s3(:),s4(:),s5(:)   !ly0=4
   integer :: k,px1 !ly0=4
   integer :: i,j !check
   !real(8) :: ran

   !external :: ran


   b=0
   d=0
   if (ly0==0) then        ! 1D chain for test
      !call readbond()
      do x1=0,lx-2
         b=b+1
         s=xyi(x1,0)
         x2=x1+1
         bsite(1,b)=s
         bsite(2,b)=xyi(x2,0)
      enddo
   elseif (ly0==1) then        ! 1D open chain 
      !do x1=0,lx-2
      !   b=b+1
      !   s=xyi(x1,0)
      !   x2=x1+1
      !   bsite(1,b)=s
      !   bsite(2,b)=xyi(x2,0)
      !enddo
      do x1=0,lx-2
         b=b+1
         s=x1+1
         x2=s+1
         bsite(1,b)=s
         bsite(2,b)=x2
      enddo
   elseif (ly0==2) then    ! spare two-leg ladder
      b=0
      do y1=0,1
         do x1=1,lx-1
            b=b+1
            initialbsite(1,b)=y1*lx+x1
            initialbsite(2,b)=y1*lx+x1+1
         end do
      end do
      !do x1=1,lx
      !   b=b+1
      !   bsite(1,b)=x1
      !   bsite(2,b)=x1+lx
      !end do
      !--- four intersection chains ---!
   elseif (ly0==4) then    ! two intersecting chains
      allocate(xx0(0:1))
      allocate(ss0(0:1))
      allocate(ss1(0:1))
      allocate(ss2(0:1))
      allocate(s3(0:1))
      allocate(s4(0:1))
      allocate(s5(0:1))
      xx0(0)=lx/4-1
      xx0(1)=3*lx/4
      ss0(0)=xyi(xx0(0),0)    !s0
      ss0(1)=xyi(xx0(0),2)    !s0'
      ss1(0)=xyi(xx0(0)+1,0)  !s1
      ss1(1)=xyi(xx0(0)+1,2)  !s1'
      ss2(0)=xyi(xx0(0),1)    !s2
      ss2(1)=xyi(xx0(0),3)    !s2'
      s3(0)=xyi(xx0(1),0)    !s3
      s3(1)=xyi(xx0(1),2)    !s3'
      s4(0)=xyi(xx0(1)-1,0)  !s4
      s4(1)=xyi(xx0(1)-1,2)  !s4'
      s5(0)=xyi(xx0(1),1)    !s5
      s5(1)=xyi(xx0(1),3)    !s5'
      do y1=0,ly-1
         if (mod(y1,2) == 0) then
            do x1=0,lx-2
               if (x1 /= xx0(0) .and. x1 /= xx0(1)-1) then
                  s=xyi(x1,y1)
                  x2=x1+1
                  y2=y1
                  b=b+1
                  bsite(1,b)=s
                  bsite(2,b)=xyi(x2,y2)
               end if 
            end do
            do k=0,1
               px1=y1/2*(lly/2)+k*(2*lx/8-1)
               do x1=0+k*3*lx/4,(3*k+1)*lx/4-1           !Y-chain
                  s=xyi(x1,y1)
                  pos(1,s)=xx0(k)
                  pos(2,s)=px1
                  ipos(xx0(k),px1)=s
                  if (k == 0) px1=px1+1
                  if (k == 1) px1=px1-1
               end do 
               px1=(y1/2)*(lly-1)+(-1)**(y1/2)*(lx/4)
               do x1=(k+1)*lx/4,(k+2)*lx/4-1               !X-chain
                  s=xyi(x1,y1)
                  pos(1,s)=x1
                  pos(2,s)=px1
                  ipos(x1,px1)=s
               end do
            end do
         else if (mod(y1,2) == 1) then
            do x1=0,lx-2
               if (x1 /= lx/2-1) then
                  s=xyi(x1,y1)
                  x2=x1+1
                  y2=y1
                  b=b+1
                  bsite(1,b)=s
                  bsite(2,b)=xyi(x2,y2)
               end if 
            end do
            do k=0,1
               px1=k*s4(0)
               do x1=0+k*3*lx/4,(3*k+1)*lx/4-1           !X-chain
                  s=xyi(x1,y1)
                  pos(1,s)=px1
                  pos(2,s)=(y1/2)*(lly-1)+(-1)**(y1/2)*(lx/4)
                  ipos(px1,pos(2,s))=s
                  px1=px1+1
               end do
               px1=1+pos(2,ss2(y1/(ly0-1)))+k*(2*lx/8-1)
               do x1=(k+1)*lx/4,(k+2)*lx/4-1               !Y-chain
                  s=xyi(x1,y1)
                  pos(1,s)=xx0(k)
                  pos(2,s)=px1
                  ipos(xx0(k),px1)=s
                  if (k == 0) px1=px1+1
                  if (k == 1) px1=px1-1
               end do
            end do
         end if 
      end do
      !--- bond connecting up 2-crossing and down 2-crossing---!
      b=b+1
      bsite(1,b)=xyi(lx/2-1,1)
      bsite(2,b)=xyi(0,2)
      b=b+1
      bsite(1,b)=xyi(lx/2,1)
      bsite(2,b)=xyi(lx-1,2)
      !--------------------------------------------------------!
      do k=0,1
         b=b+1
         bsite(1,b)=ss0(k)
         bsite(2,b)=ss2(k)
         b=b+1
         bsite(1,b)=ss1(k)
         bsite(2,b)=ss2(k)
         b=b+1
         bsite(1,b)=s4(k)
         bsite(2,b)=s5(k)
         b=b+1
         bsite(1,b)=s3(k)
         bsite(2,b)=s5(k)
      end do

      deallocate(xx0)
      deallocate(ss0)
      deallocate(ss1)
      deallocate(ss2)
      deallocate(s3)
      deallocate(s4)
      deallocate(s5)
   end if 

   if(b/=nb) then
      write(*,*) "error: nb (note in trial we don't care)",b,nb
   endif

   if (ly0 /= 2) then
      initialbsite=bsite !store the initial bsite
   end if

   !! check

   !! check
   !open(123,file='lattice.dat',status='unknown')
   !do b=1,nb
   !  !print*, bsite(1,b), bsite(2,b)
   !  j=bsite(1,b)
   !  i=bsite(2,b)
   !  print*,'(',pos(1,j),',',pos(2,j),')' &
   !     & ,'(',pos(1,i),',',pos(2,i),')',b
   !  write(123,*)pos(1,j),pos(2,j)
   !  write(123,*)pos(1,i),pos(2,i)
   !  write(123,*)''
   !enddo 
   !close(123)

   !! check
   !open(123,file='lattice.dat',status='unknown')
   !do s=1,nn
   !   write(123,*)pos(1,s),pos(2,s)
   !   print*, s, xy(1,s), xy(2,s), pos(1,s), pos(2,s)
   !enddo
   !close(123)

end subroutine makelattice
!--------------------------!

!--------------------------!
subroutine allocateall()
   !--------------------------!
   use system; implicit none

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
   if (ly0==2) then    !   
      ly=2
      nn=ly*lx          ! the total number of the spins = 2*lx
      allocate(sglabel(1:nn,1:nn))
   end if 
   if (ly0==4) then    ! four crossing; lx must be a multiple of 8  
      ly=4
      lly=lx+2
      nn=ly*lx          ! the total number of the spins = 4*lx
   endif   



   allocate(ss(0:7))
   allocate(ab(nn))
   allocate(lspn(nn))
   allocate(initiallspn(nn))
   allocate(spinsurv(nn))

   allocate(pos(2,nn))
   allocate(ipos(0:lx-1,0:lly-1))
   allocate(xy(2,nn))
   allocate(xyi(0:lx-1,0:ly-1))

   open(99,file='d.dat',status='unknown')
   do s=1,nn
      x1=mod(s-1,lx)
      y1=(s-1)/lx
      xy(1,s)=x1           ! x coordinate
      xy(2,s)=y1           ! y coordinate
      xyi(x1,y1)=s         ! site indices
      write(99,*)x1,y1
   enddo
   close(99)

   if(ly0<2) then
      pos=xy
      ipos=xyi
   endif

   if (ly==1) then        ! 1d open chain
      nb=lx-1
      nd=lx-2
   elseif (ly0==2) then    ! spare 2-leg ladder
      !nb=2*(lx-1)+lx
      nb=2*(lx-1)
   elseif (ly0==4) then    ! 4 intersecting chains
      nb=4*lx
      !nd=lx+(lx-2)+(lx-3)  !!
      !nd=lx
      nd=2*(lx-2)+2*(lly-2)
   endif

   allocate(initialbsite(2,nb))
   if (ly0 /= 2) then
      allocate(bsite(2,nb))
      allocate(bondsurv(nb))
      allocate(scor(nb))
      allocate(jj0(1:nb))
      allocate(jj(1:nb))
   end if
   initialnb=nb



end subroutine allocateall
!--------------------------!




!--------------------------!
subroutine deallocateall()
   !--------------------------! 
   use system; implicit none

   deallocate(ab)
   deallocate(lspn)
   deallocate(spinsurv)
   !deallocate(bsite)
   !deallocate(bondsurv)
   deallocate(xyi)
   deallocate(xy)
   deallocate(pos)
   deallocate(ipos)
   !deallocate(jj0)
   !deallocate(jj)
   !deallocate(scor)
   deallocate(ss)
   deallocate(initialbsite)
   deallocate(initiallspn)
   if (ly0 /= 2) then
      deallocate(bsite)
      deallocate(jj)
      deallocate(scor)
      deallocate(bondsurv)
      deallocate(jj0)
   end if

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
