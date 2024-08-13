!--------------------!
 subroutine measure()
!--------------------!
 use configuration; use measurementdata; implicit none

 integer :: i,b,op,s1,s2,am,j,k,l
 real(8) :: am2,scorone,ax1,am1

 k=xyi(lx/4-1,ly/2)
 l=xyi(3*lx/4-2,ly/2)
 am=0
 am1=0.d0
ax1=0.d0
scorone=0.d0
 do i=1,nn
    am=am+spin(i)*(-1)**(xy(1,i)+xy(2,i))
 enddo
!------------------------------------------------!
sumspin=0.d0
t=-1
!-----------------------------------------!
 am=am/2
 am2=0.d0
 do i=0,mm-1
   op=opstring(i)
!  if (op==0) cycle
   if(op/=0)then
     if(btest(op,0)) then
       b=op/8
       s1=bsite(1,b)
       s2=bsite(2,b)
       sumspin(s1)=sumspin(s1)+(i-t(s1))*spin(s1)
       sumspin(s2)=sumspin(s2)+(i-t(s2))*spin(s2)
       t(s1)=i
       t(s2)=i
       spin(s1)=-spin(s1)
       spin(s2)=-spin(s2)
       am=am+2*spin(s1)*(-1)**(xy(1,s1)+xy(2,s1))
    endif
     if(btest(op,1)) then
       b=op/8
       s1=bsite(3,b)
       s2=bsite(4,b)
       sumspin(s1)=sumspin(s1)+(i-t(s1))*spin(s1)
       sumspin(s2)=sumspin(s2)+(i-t(s2))*spin(s2)
       t(s1)=i
       t(s2)=i
       spin(s1)=-spin(s1)
       spin(s2)=-spin(s2)
       am=am+2*spin(s1)*(-1)**(xy(1,s1)+xy(2,s1))
    endif
     ax1=ax1+dfloat(am)
     am1=am1+dfloat(abs(am))
    am2=am2+dfloat(am)**2
endif
 enddo
     do j=1,nn
     sumspin(j)=0.5d0*(sumspin(j)+(mm-1-t(j))*spin(j))
     enddo
 if (nh/=0) then
    ax1=(ax1**2+am2)/(dfloat(nh)*dfloat(nh+1))
    am1=am1/nh
    am2=am2/nh
 else
    am1=dfloat(abs(am))
    am2=dfloat(am)**2
    ax1=am2
 endif

do i=1,nn
 ab(i)=ab(i)+sumspin(i)**2+0.25*mm
 szspin(i)=szspin(i)+sumspin(i)/dble(mm)
enddo
 enrg1=enrg1+dfloat(nh)
 enrg2=enrg2+dfloat(nh)**2
 amag2=amag2+am2
 ususc=ususc+(dble(sum(spin)/2.0))**2
 asusc=asusc+ax1

if(sum(spin)/=0)then
! print*,sum(spin)
endif
 end subroutine measure
!----------------------!

!------------------------------------!
 subroutine writeresults(msteps,bins)
!------------------------------------!
 use configuration; use measurementdata; implicit none

 integer :: i,msteps,bins
 real(8) :: wdata1(5),wdata2(5)

! scor=scor/msteps
 enrg1=enrg1/msteps
 enrg2=enrg2/msteps
 amag2=amag2/msteps
 ususc=ususc/msteps
asusc=asusc/msteps
  szspin=szspin/msteps
 sumspin=sumspin/msteps
 ab=ab/msteps
 do i=1,nn
 ab(i)=ab(i)*beta/mm/(mm+1)-beta*szspin(i)**2
 enddo
open(12,file='local-sus.dat',status='unknown',position='append')
 do i=1,nn
 write(12,*) ab(i),xy(1,i),xy(2,i)
 enddo
!  write(12,*) beta,ab(1)
!close(12)


 enrg2=(enrg2-enrg1*(enrg1+1.d0))/nn
 enrg1=enrg1/(beta*nn)-0.5d0
 amag2=3.d0*amag2/dble(nn)**2
 ususc=beta*ususc/nn
 asusc=beta*asusc/nn

 data1(1)=data1(1)+enrg1
 data1(2)=data1(2)+enrg2
 data1(3)=data1(3)+amag2
 data1(4)=data1(4)+ususc
 data1(5)=data1(5)+asusc

 data2(1)=data2(1)+enrg1**2
 data2(2)=data2(2)+enrg2**2
 data2(3)=data2(3)+amag2**2
 data2(4)=data2(4)+ususc**2
 data2(5)=data2(5)+asusc**2

 print*,bins
 wdata1(:)=data1(:)/bins
 wdata2(:)=data2(:)/bins
 wdata2(:)=sqrt(abs(wdata2(:)-wdata1(:)**2)/bins)

open(11,file='susproc.dat',status='replace')
 write(11,'(i8,2f14.6,5f18.12)')lx,beta,qqq,wdata1
close(11)
open(11,file='sus.dat',status='unknown',position='append')
 write(11,'(i8,2f14.6,3f18.12)')lx,beta,qqq,amag2,ususc,asusc
close(11)
 scor=0.d0
 enrg1=0.d0
 ab=0.d0
 enrg2=0.d0
 amag2=0.d0
 ususc=0.d0
 asusc=0.d0
 szspin=0.d0
 sumspin=0.d0
 end subroutine writeresults
!---------------------------!


