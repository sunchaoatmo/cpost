MODULE CS_STAT
  IMPLICIT NONE
  public :: QsortC
  real,parameter::log2o3=log10(2.0)/log10(3.0)
  real,parameter::sqrt2=sqrt(2.0)
  CONTAINS
subroutine kendallS(x,n,pvalue)
integer, intent(in)  ::n
real, intent(in), dimension(n)::x
real, intent(out):: pvalue
integer::i,j,k,tg
logical::belongold
integer, dimension(30):: tp_num !!assume less then 30 tie group is enough for this case
integer, dimension(30):: tp_val !!assume less then 30 tie group is enough for this case
real::s,zs,sigma
s=0
tg=0
tp_val=0
tp_num=0

do i=1,n-1
  do j=i+1,n
   if ((x(j)-x(i)).gt.0) then
    s=s+1
   elseif  ((x(j)-x(i)).lt.0) then
    s=s-1
   else
    belongold=.False.
    do k=1,tg
     if (tp_val(k).eq.x(j)) then
       tp_num(k)=tp_num(k)+1
       belongold=.True.
       exit
     endif
    end do
    if (.not.belongold) then
      tg=tg+1
      tp_num(tg)=1
      tp_val(tg)=x(j)
    end if
   endif 
  enddo 
enddo 
sigma=n*(n-1)*(2.0*n+5.0)
do i=1,tg
  sigma=sigma-tp_val(i)*(tp_val(i)-1.0)*(2.0*tp_val(i)+5.0)
enddo
sigma=sqrt(sigma/18.0)
if (s.gt.0) then 
 zs=(s-1)/sigma
elseif (s.lt.0) then
 zs=(s+1)/sigma
else
 zs=0
endif
!pvalue=1-0.5*(1+erf(zs/sqrt2)) 
if (zs.gt.0) then
  pvalue=1-0.5*(1+erf(zs/sqrt2))
else
  pvalue=0.5*(1+erf(zs/sqrt2))
endif
end subroutine kendallS

!#####################################
recursive subroutine QsortC(A)
! Recursive Fortran 95 quicksort routine
! sorts real numbers into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms,
! 1997 printing
! Made F conformant by Walt Brainerd
  real, intent(in out), dimension(:) :: A
  integer :: iq

  if(size(A) > 1) then
     call Partition(A, iq)
     call QsortC(A(:iq-1))
     call QsortC(A(iq:))
  endif
end subroutine QsortC

subroutine Partition(A, marker)
  real, intent(in out), dimension(:) :: A
  integer, intent(out) :: marker
  integer :: i, j
  real :: temp
  real :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition
!#####################################
!#################GEV parameters
real function gevppf(q,k,alpha,xi)
real,intent(in)::k,alpha,xi
real,intent(in)::q
if (k.ne.0) then
  gevppf=alpha/k*( (-log(q))**(-k)-1.0 )+xi
else
  gevppf=-alpha*log( -log(q) )+xi
end if

end function
subroutine gevfit(x,n,k,alpha,xi)
real::ub0,ub1,ub2
real,intent(out)::k,alpha,xi
real,intent(inout),dimension(n)::x
integer,intent(in)::n
call QsortC(x)
ub0=b0(x,n)
ub1=b1(x,ub0,n)
ub2=b2(x,ub0,ub1,n)

call gevp(k,alpha,xi,ub0,ub1,ub2)
end subroutine

function b0(x,n)
integer:: n
real,dimension(n)::x
real::b0
b0=1.0/n*sum(x(:))
end function

function b1(x,b0,n)
integer:: n
real,dimension(n)::x
real::b1
real::b0
real:: U
integer::i,j
U=0
do j=1,n-1
 do i=j+1,n
   U=U+(x(i)-x(j))
 end do
end do
U=1.0/n/(n-1.0)*U
b1=0.5*(U+b0)
end function

function b2(x,b0,b1,n)
integer:: n
real,dimension(n)::x
real::b2
real::b1
real::b0
real:: U
integer::i,j,k
U=0
do k=1,n-2
 do j=1+k,n-1
  do i=1+j,n
   U=U+(x(i)-2*x(j)+x(k))
  end do
 end do
end do
U=2.0/n/(n-1.0)/(n-2.0)*U
b2=(U-b0+6.0*b1)/6.0
end function

subroutine gevp(k,alpha,xi,b0,b1,b2)
real,intent(in)::b0,b1,b2
real,intent(out)::k,alpha,xi
real::c
real::gamma1k
c=(2*b1-b0)/(3*b2-b0)-log2o3
k=7.8590*c+2.9554*c*c
gamma1k=gamma(1.0+k)
alpha=(2.0*b1-b0)*k/gamma1k/(1.0-2.0**(-k))
xi=b0+alpha*(gamma1k-1.0)/k

end subroutine

SUBROUTINE convetive2total(pr,raincv,ratio_bin,&
                 maskval,mask,&
                 low_bin,up_bin,&
                 nx,ny,nt,nbin &
                 )
  ! USESindexx
  !Sorts an array ra(1:n) into ascending numerical order 
  !while making the corresponding rearrangements of the arrays rb(1:n) and rc(1:n). 
  !An index table is constructed via the routine indexx.
  INTEGER,intent(in)                   ::nx,ny,nt,nbin
  REAL,dimension(nt,nx,ny),intent(in)  ::pr,raincv
  REAL,intent(in),dimension(nx,ny)     ::mask
  INTEGER,intent(in),dimension(nbin)   ::low_bin
  INTEGER,intent(in),dimension(nbin)   ::up_bin
  REAL,intent(in)                      ::maskval
  REAL,dimension(nbin),intent(out)     ::ratio_bin
  !loc
  REAL                                 ::pr_1d(nbin),raincv_1d(nbin)
  INTEGER                              ::ibin,i,j,t,t2,iday,nday,npoints
  INTEGER                              ::up_ind,low_ind

  pr_1d(:)=0.0
  raincv_1d(:)=0.0
  do i=1,nx
    do j=1,ny
      if (mask(i,j).eq.maskval) then
        do t=1,nt
          do ibin=1,nbin
            if (pr(t,i,j).le.up_bin(ibin).and.up_bin(ibin).gt.low_bin(ibin))then
               pr_1d(ibin)=pr_1d(ibin)+pr(t,i,j)
               raincv_1d(ibin)=raincv_1d(ibin)+raincv(t,i,j)
               exit
            endif
          enddo
        enddo
      endif
    enddo
  enddo

  ratio_bin(:)=raincv_1d+pr_1d
END SUBROUTINE convetive2total



!#####################################

SUBROUTINE sort3(pr,raincv,ratio_bin,&
                 maskval,mask,&
                 low_bin,up_bin,&
                 nx,ny,nt,nbin, &
                 qctl)
  ! USESindexx
  !Sorts an array ra(1:n) into ascending numerical order 
  !while making the corresponding rearrangements of the arrays rb(1:n) and rc(1:n). 
  !An index table is constructed via the routine indexx.
  INTEGER,intent(in)                   ::nx,ny,nt,nbin
  REAL,dimension(nt,nx,ny),intent(in)  ::pr,raincv
  REAL,intent(in),dimension(nx,ny)     ::mask
  INTEGER,intent(in),dimension(nbin)   ::low_bin
  INTEGER,intent(in),dimension(nbin)   ::up_bin
  REAL,intent(in)                      ::maskval
  REAL,intent(in)                      ::qctl
  REAL,dimension(nbin),intent(out)     ::ratio_bin
  !loc
  REAL,dimension(nt      )             ::pr_temp,ratio_1d
  REAL,dimension(nbin,nx,ny)           ::ratio_bin_2d
  INTEGER                              ::iwksp(nt)
  REAL                                 ::wksp(nt)
  INTEGER                              ::ibin,i,j,t,t2,iday,nday,npoints
  INTEGER                              ::up_ind,low_ind

  ratio_bin_2d(:,:,:)=0.0
  do i=1,nx
    do j=1,ny
      if (mask(i,j).eq.maskval) then
        iday=0
        do t=1,nt
          if (pr(t,i,j)>qctl)then 
            iday=iday+1
            pr_temp(iday)=pr(t,i,j)
            ratio_1d(iday)=raincv(t,i,j)/pr_temp(iday)
          endif
        enddo
        nday=iday
        call indexx(iday,pr_temp(1:nday),iwksp(1:nday)) 
        do iday=1,nday
          wksp(iday)=ratio_1d(iday) 
        enddo 
        do  iday=1,nday
          ratio_1d(iday)=wksp(iwksp(iday))
        enddo 
        do ibin=1,nbin
          low_ind=floor(low_bin(ibin)/100.0*nday) 
          low_ind=max(low_ind,1)
          up_ind=ceiling(up_bin(ibin)/100.0*nday)
          up_ind=min(nday,up_ind)
          do t2=low_ind,up_ind
            ratio_bin_2d(ibin,i,j)=ratio_bin_2d(ibin,i,j)+ratio_1d(t2)
          enddo
          if (low_ind<up_ind) then 
            ratio_bin_2d(ibin,i,j)=ratio_bin_2d(ibin,i,j)/(up_ind-low_ind)
          else
            ratio_bin_2d(ibin,i,j)=0.0
          endif
        enddo
      endif
    enddo
  enddo
  ratio_bin(:)=0.0
  npoints=0
  do i=1,nx
    do j=1,ny
      if (mask(i,j).eq.maskval) then
        ratio_bin(:)=ratio_bin(:)+ratio_bin_2d(:,i,j)
        npoints=npoints+1
      endif
    enddo
  enddo
  ratio_bin(:)=ratio_bin(:)/npoints
END SUBROUTINE sort3

  subroutine tendif3dgev(pr , qctl,season_periods_loc ,   &
  mask                   , maskval            , pr_ctl,&
  season                 , &
  trend_i                , trend_f            , &
  nyear                  , nday               , nx         , ny)
   real, intent(in),dimension(4,nyear*2) ::season_periods_loc
   real, intent(in),dimension(nday,nx,ny)::pr
   real,intent(in),dimension(nx,ny)::mask
   real,intent(in)                 ::maskval
   real,intent(in)                 ::qctl
   integer,intent(in)::nday,nyear,nx,ny
   real, intent(in)                       ::pr_ctl
   integer, intent(in)              ::season
!   integer,parameter                ::nout=4
!   real,intent(out),dimension(:,:,:)::trend_i,trend_f
   real,intent(out),dimension(nx,ny)::trend_i,trend_f
   !real,intent(out),dimension(nout,nx,ny)::trend_i,trend_f
   !local
   real::prthresh
   real, dimension (4,nyear)   ::intense
   real, dimension (nyear)   ::pr_max
   real, dimension (4,nyear)  :: freq
   real, dimension (4,nyear) ::tend_f,tend_i
   real, dimension (nyear) ::total
   integer ::i,j,k,s,nt
   integer ::bs, es
   real,dimension(nyear)::x,sig
   real::a,b1,b2,siga,sigb,chi2,q
   integer :: mwt
   real :: pvalue,meanintense,meantotal
   real::gevk,alpha,xi,ppf
   do k=1,nyear
    x(k)=k
   end do
   mwt=0
   sig=0
   total=-999
   s=season+1 ! from python which is starting from 0
!         bs=season_periods_loc(s,2*1-1)
!         es=season_periods_loc(s,2*1)
!   print*,"season=",s,bs,es
   do j=1,ny
     do i=1,nx
     if (mask(i,j).ne.maskval) then
       bs=1
       do k=1,nyear
         bs=season_periods_loc(s,2*k-1)
         es=season_periods_loc(s,2*k)
         pr_max(k)=maxval(pr(bs:es,i,j))
       end do
       call  gevfit(pr_max,nyear,gevk,alpha,xi)
       prthresh=gevppf(0.2,gevk,alpha,xi)
!       bs=1
       do k=1,nyear
         bs=season_periods_loc(s,2*k-1)
         es=season_periods_loc(s,2*k)
         nt=es-bs+1
         call freqintensegev(pr(bs:es,i,j),prthresh,intense(s,k),freq(s,k),total(k),nt)
       end do
       call fit(x, freq(s,:),nyear,sig,mwt,a,b1,siga,sigb,chi2,q)
       call kendallS(freq(s,:),nyear,pvalue)
!       call lsqtest(x,freq(s,:),a,b1,nyear,pvalue)
!      if (b1<0) then
!       print*,"x=",freq(s,:)
!       print*,"pvalue=",pvalue
!       print*,"coef=",b1
!       pause
!      end if
       meanintense=sum(intense(s,:))/nyear  !/nintervals
       meantotal=sum(total(:))/nyear  !/nintervals
       if (pvalue<qctl) then
         !trend_f(i,j) =maxval(total)/nyear!  1e4=100(percent)*100(centrury)
         trend_f(i,j) =100*b1*meanintense/meantotal!  1e4=100(percent)*100(centrury)
       else
         trend_f(i,j)=-99999
       end if
       call fit(x, total(:),nyear,sig,mwt,a,b2,siga,sigb,chi2,q)

       if (pvalue<qctl ) then
        !trend_i(i,j)=sum(freq(s,:))/nyear
        trend_i(i,j)=(b2/meantotal*1e2-trend_f(i,j)) !sum(total(interval,s,:)) !b-trend_f(i,j)
       else
        trend_i(i,j)=-99999
       end if
     else
       trend_i(i,j)=-99999
       trend_f(i,j)=-99999
     end if

     end do
   end do

  end subroutine tendif3dgev





  subroutine tendif3d(pr , season_periods_loc , nintervals , selectintval , nout,&
  mask                   , maskval            , pr_ctl,&
  season                 , &
  trend_i                , trend_f            , &
  nyear                  , nday               , nx         , ny)
   integer,intent(in):: nintervals        !=10
   integer,dimension(nintervals),intent(in):: selectintval  !=(/0,0,0, 0,0,0, 0,0,0, 0/)
   real, intent(in),dimension(4,nyear*2) ::season_periods_loc
   real, intent(in),dimension(nday,nx,ny)::pr
   real,intent(in),dimension(nx,ny)::mask
   real,intent(in)                 ::maskval
   integer,intent(in)::nday,nyear,nx,ny
   real, intent(in)                       ::pr_ctl
   integer, intent(in)              ::season
   integer, intent(in)              ::nout
!   integer,parameter                ::nout=4
!   real,intent(out),dimension(:,:,:)::trend_i,trend_f
   real,intent(out),dimension(nx,ny)::trend_i,trend_f
   !real,intent(out),dimension(nout,nx,ny)::trend_i,trend_f
   !local
   real, dimension (nintervals,4,nyear)   ::intense
   real, dimension (4,nyear)  :: freq
   real, dimension (nintervals,4,nyear) ::tend_f,tend_i
   real, dimension (nintervals,nyear) ::total
   integer ::i,j,k,s,nt
   integer ::bs, es
   real,dimension(nyear)::x,sig
   real::a,b1,b2,siga,sigb,chi2,q
   integer :: mwt,interval,iout
   real :: pvalue,meanintense,meantotal
   do k=1,nyear
    x(k)=k
   end do
   mwt=0
   sig=0
   iout=0
   total=-999
   do interval=1,nintervals
    if (selectintval(interval)==1) then
     iout=iout+1
     do j=1,ny
       do i=1,nx
       if (mask(i,j).ne.maskval) then
         bs=1
         do k=1,nyear
           !do s=1,4
           s=season+1 ! from python which is starting from 0
           bs=season_periods_loc(s,2*k-1)
           es=season_periods_loc(s,2*k)
           nt=es-bs+1
           !print*,bs,es
!           print*,"bef",total(:,k)
           call freqintense(pr(bs:es,i,j),intense(:,s,k),freq(s,k),total(:,k),&
                               nt,nintervals, pr_ctl)
!           print*,"aft",total(:,k)



           !  bs=es+1
           !end do
         end do
!         do s=1,4
!         trend_i(i,j)=sum(tend_i(interval,s,:))/(nyear-1)
!         trend_f(i,j)=sum(tend_f(interval,s,:))/(nyear-1)
!         print*,tend_i(interval,s,:)
!         print*,trend_f(i,j)
!         end do
         call fit(x, freq(s,:)/nintervals,nyear,sig,mwt,a,b1,siga,sigb,chi2,q)
         call lsqtest(x,freq(s,:)/nintervals,a,b1,nyear,pvalue)
         meanintense=sum(intense(interval,s,:))/nyear  !/nintervals
         meantotal=sum(total(interval,:))/nyear  !/nintervals
!         call fit(x, tend_f(interval,s,:),nyear,sig,mwt,a,b,siga,sigb,chi2,q)
         if (pvalue<3) then
         !if (q>0.1) then
           trend_f(i,j) =100*b1*meanintense/meantotal!  1e4=100(percent)*100(centrury)
         else
           trend_f(i,j)=-99999
         end if
         call fit(x, total(interval,:),nyear,sig,mwt,a,b2,siga,sigb,chi2,q)

         if (pvalue<3 ) then
          trend_i(i,j)=(b2/meantotal*1e2-trend_f(i,j)) !sum(total(interval,s,:)) !b-trend_f(i,j)
         else
          trend_i(i,j)=-99999
         end if
       else
         trend_i(i,j)=-99999
         trend_f(i,j)=-99999
       end if

       end do
     end do
    endif
   end do

  end subroutine tendif3d

  subroutine tendif_1d(freq,intense,nyear,nintervals,s,interval,b,trend_f)
   real, intent(in) ,dimension (nintervals,4,nyear) ::intense
   real, intent(in) ,dimension (4,nyear)  :: freq
   real, intent(out)::trend_f
!   real, intent(out),dimension (nintervals,4,nyear) ::tend_f,tend_i
   integer,intent(in)::nintervals,nyear,interval,s
   real, intent(in)   :: b
!   real, intent(in) ,            dimension (nintervals,4,nyear) ::total
   real, dimension (nintervals,4) ::meanintense
   integer :: k,j!,s
   !tend_i=0
   !do k=1,nyear
   !  do s=1,4
   !    total(:,s,k)=freq(s,k)*intense(:,s,k) !/nintervals
   !  end do
   !end do
!   do j=1,nintervals
!     do s=1,4
       meanintense(interval,s)=sum(intense(interval,s,:))/nyear  !/nintervals
!     end do
!   end do
   trend_f=b*meanintense(interval,s)
!   do k=1,nyear
!     do s=1,4
!       tend_f(:,s,k)=b*meanintense(:,s)/nintervals
!       tend_f(:,s,k)=freq(s,k)*meanintense(:,s)/nintervals
!     end do
!   end do
   
  end subroutine tendif_1d
  subroutine freqintensegev(pr,prthresh,intense,freq,total,nt)
   integer,intent(in)::nt
   real, intent(in),dimension(nt)         ::pr
   real, intent(in)                       ::prthresh
   real, intent(out)::intense
   !real, intent(out),dimension (nintervals) ::total
   real, intent(out)::total
   real, intent(out):: freq
   !local
   integer ::s,j,k,npoints
   real, dimension(nt)         ::temp
   real ::x
   integer ::bi,ei,rainydays
   rainydays=0
   freq=0
   total=0
   do k=1,nt
    if (pr(k).gt.prthresh) then
     freq=freq+1
     total=total+pr(k)
    end if
   enddo 
   if(freq>0) then
     intense=total/freq
   else
     intense =0
   end if
  end subroutine

  subroutine freqintense(pr,intense,freq,total,nt,nintervals,pr_ctl)
   integer,intent(in)::nt
   integer,intent(in)::nintervals
   real, intent(in),dimension(nt)         ::pr
   real, intent(in)                       ::pr_ctl
   real, intent(out),dimension (nintervals) ::intense
   !real, intent(out),dimension (nintervals) ::total
   real, intent(out) ,            dimension (nintervals) ::total
   real, intent(out):: freq
   !local
   integer ::s,j,k,npoints
   real, dimension(nt)         ::temp
   real ::x
   integer ::bi,ei,rainydays
   rainydays=0
   do k=1,nt
    if (pr(k).gt.pr_ctl) then
     rainydays=rainydays+1
     temp(rainydays)=pr(k)
    end if
   enddo 
   call QsortC(temp(1:rainydays))
   freq=1.0*rainydays!/nt
   bi=1
   do j=1,nintervals
    ei=j*rainydays/nintervals!begin of interval
    npoints=ei-bi+1
    if (npoints>0) then
      x=sum(temp(bi:ei))
      total(j)=x !/nt
      intense(j)=x/npoints
    else
      intense(j)=0
      total(j)=0 !/nt
    end if
    bi=ei+1
   end do
  end subroutine
  subroutine fit_3d(y_3d,mask,maskval,trend,nt,nx,ny)
    integer,intent(in)               ::nt,nx,ny
    real,intent(in),dimension(nt,nx,ny)::y_3d
    real,intent(in),dimension(nx,ny)::mask
    real,intent(in)                 ::maskval
    real,intent(out),dimension(nx,ny)::trend
    
    !local 
    real,dimension(nt)::x,sig
    real::a,b,siga,sigb,chi2,q
    integer :: mwt,i,j,k
    do k=1,nt
     x(k)=k
    end do
    mwt=0
    sig=0
    do j=1,ny
      do i=1,nx
        if (mask(i,j).ne.maskval) then
          call fit(x, y_3d(:,i,j),nt,sig,mwt,a,b,siga,sigb,chi2,q)
          call lsqtest(x,y_3d(:,i,j),a,b,nt,q)
          trend(i,j)=b/sum(y_3d(:,i,j))*nt*100
          !if (q<0.5) then
          ! trend(i,j)=b/sum(y_3d(:,i,j))*nt*100
          !else
          ! trend(i,j)=0 !-99999
          !end if
        else
          trend(i,j)=-99999
        end if
             
      end do
    end do
  end subroutine fit_3d
  subroutine corrcoef_mask(obs,sim,mask,maskval,nx,ny,cor)
    integer,intent(in)               ::nx,ny
    real,intent(in),dimension(nx,ny)::obs,sim
    real,intent(in),dimension(nx,ny)::mask
    real,intent(in)                 ::maskval
    real,intent(out)                ::cor
    !local 
    real,dimension(nx*ny)           ::x_1d,y_1d
    integer :: i,j,k
    k=0
    do j=1,ny
      do i=1,nx
        if (mask(i,j).ne.maskval) then
          k=k+1
          x_1d(k)=obs(i,j)
          y_1d(k)=sim(i,j)
        end if
      end do
    end do
    cor=corrcoef_1d(x_1d(1:k),y_1d(1:k),k)
  end subroutine corrcoef_mask

  subroutine mean_2d_mask(idata_3d,mask,nt,nx,ny,odata_1d)
    integer,intent(in)               ::nt,nx,ny
    real,intent(in),dimension(nt,nx,ny)::idata_3d
    real,intent(in),dimension(nx,ny)::mask
    real,intent(out),dimension(nt)::odata_1d
    !local 
    integer :: i,j,k,landpoints
    landpoints =0
    do i=1,nx
      do j=1,ny
        if (mask(i,j).lt.0.5) then
          do k=1,nt
            if (idata_3d(k,i,j).gt.1000) then
            !  print*,idata_3d(k,i,j),mask(i,j)
            endif
          end do
          odata_1d(:)=odata_1d(:)+idata_3d(:,i,j)
          landpoints =landpoints+1
        endif
      end do
    end do
    print*,landpoints
    odata_1d(:)=odata_1d(:)/landpoints
  end subroutine mean_2d_mask

  subroutine mean_2d(idata_3d,nt,nx,ny,odata_1d)
    integer,intent(in)               ::nt,nx,ny
    real,intent(in),dimension(nt,nx,ny)::idata_3d
    real,intent(out),dimension(nt)::odata_1d
    !local 
    integer :: i,j,k
    do k=1,nt
      odata_1d(k)=sum(idata_3d(k,:,:))/nx/ny
    end do
  end subroutine mean_2d



  subroutine corrcoef_2d(obs,sim,nt,nx,ny,cor)
    integer,intent(in)               ::nt,nx,ny
    real,intent(in),dimension(nt,nx,ny)::obs,sim
    real,intent(out),dimension(nx,ny)::cor
    !local 
    integer :: i,j,k
    do j=1,ny
      do i=1,nx
        cor(i,j)=corrcoef_1d(obs(:,i,j),sim(:,i,j),nt)
      end do
    end do
  end subroutine corrcoef_2d

  subroutine rmse_2d(obs,sim,nt,nx,ny,rmse)
    integer,intent(in)               ::nt,nx,ny
    real,intent(in),dimension(nt,nx,ny)::obs,sim
    real,intent(out),dimension(nx,ny)::rmse
    !local 
    integer :: i,j,k
    do j=1,ny
      do i=1,nx
        rmse(i,j)=rmse_1d(obs(:,i,j),sim(:,i,j),nt)
      end do
    end do
  end subroutine rmse_2d

  real function rmse_1d(obs,sim,nt)
    implicit none
    integer,intent(in)               ::nt
    real,intent(in),dimension(nt)::obs,sim
    real,dimension(nt):: diff,diffsqre
    diff=(obs-sim)
    diffsqre=diff*diff
    rmse_1d=sqrt(sum(diffsqre)/nt)
  end function rmse_1d


  real function corrcoef_1d(x,y,n)
    implicit none
    integer,intent(in) :: n
    real,intent(in),dimension(n) :: x,y
    real                         :: x_sum,x_mean
    real                         :: y_sum,y_mean
    real                         :: xp,yp
    real                         :: sx,sy,sxsy,sxy
    integer                      :: i
    sx=0
    sy=0
    sxy=0
    x_sum=0
    y_sum=0
    do i=1,n
     x_sum=x_sum+x(i)
     y_sum=y_sum+y(i)
    enddo
    x_mean=x_sum/n
    y_mean=y_sum/n
    do i=1,n
     xp=(x(i)-x_mean)
     sx=sx+xp*xp
     yp=(y(i)-y_mean)
     sy=sy+yp*yp
     sxy=sxy+xp*yp
    enddo
     sxsy=sx*sy
     if (sxsy.le.1.0d-10) then
       corrcoef_1d = 0.0
     else
       corrcoef_1d =sxy/(sqrt(sxsy))
     endif
  end function corrcoef_1d

      FUNCTION betacf(a,b,x)
      INTEGER MAXIT
      REAL betacf,a,b,x,EPS,FPMIN
      PARAMETER (MAXIT=1000,EPS=3.e-7,FPMIN=1.e-30)
      INTEGER m,m2
      REAL aa,c,d,del,h,qab,qam,qap
      qab=a+b
      qap=a+1.
      qam=a-1.
      c=1.
      d=1.-qab*x/qap
      if(abs(d).lt.FPMIN)d=FPMIN
      d=1./d
      h=d
      do 11 m=1,MAXIT
        m2=2*m
        aa=m*(b-m)*x/((qam+m2)*(a+m2))
        d=1.+aa*d
        if(abs(d).lt.FPMIN)d=FPMIN
        c=1.+aa/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        h=h*d*c
        aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
        d=1.+aa*d
        if(abs(d).lt.FPMIN)d=FPMIN
        c=1.+aa/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        del=d*c
        h=h*del
        if(abs(del-1.).lt.EPS)goto 1
11    continue
!      pause 'a or b too big, or MAXIT too small in betacf'
1     betacf=h
      return
      END FUNCTION


      FUNCTION betai(a,b,x)
      REAL betai,a,b,x
!C    USES betacf,gammln
      REAL bt
      if(x.lt.0..or.x.gt.1.)pause 'bad argument x in betai'
      if(x.eq.0..or.x.eq.1.)then
        bt=0.
      else
        bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.-x))
      endif
      if(x.lt.(a+1.)/(a+b+2.))then
        betai=bt*betacf(a,b,x)/a
        return
      else
        betai=1.-bt*betacf(b,a,1.-x)/b
        return
      endif
      END FUNCTION
  function tcdf(t,mu)
  real::tcdf, t
  integer:: mu
  tcdf=1.0-betai(mu/2.0,0.5,mu/(mu+t*t))
  end function

  subroutine lsqtest(x,y,A,B,n,pvalue)
! Inputs     : A:      The y-intercept of the fitted straight line.
!
!              B:      The slope of the fitted straight line.
   REAL ,intent(in)::A, B
   INTEGER,intent(in):: n
   REAL ,intent(in),dimension(n)::x,y
   REAL ,dimension(n)::err
   REAL ,intent(out)::pvalue
   real  :: xmean,xstd,se,t

   err(:)=y(:)-(x(:)*B+A)
   xmean=sum(x)/n
   xstd=sum ( (x(:)-xmean)*(x(:)-xmean) )
   se=sum(err(:)*err(:))/(n-2.0)/xstd
   se=sqrt(se)
   t=B/se
   pvalue=(1-tcdf(t,n-2))  ! one tail
  end subroutine
!
! Name        : FIT
!
! Purpose     : Numerical Recipes linear least squares fitting subroutine.
!
! Explanation : Given a set of NDATA points X(I), Y(I), with standard 
!               deviations SIG(I), fit them to a straight line y = a + bx by 
!               minimizing chi**2.  Returned are A, B and their respective 
!               probable uncertainties SIGA, SIGB, the chi-squared (CHI2), 
!               and the goodness-of-fit probability Q (that the fit would 
!               have chi**2 this large or larger).  If MWT=0 on input, then 
!               the standard deviations are assumed to be unavailable: 
!               Q is returned as 1.0 and the normalization of CHI2 is to unit 
!               standard deviation on all points.
!
! Use         : CALL FIT(X, Y, NDATA, SIG, MWT, A, B, SIGA, SIGB, CHI2, Q)
!
! Inputs      : X:      One-dimensional array of single precision reals of 
!                       size NDATA containing the independent variable data.
!
!       Y:  One-dimensional array of single precision reals of 
!                       size NDATA containing the dependent variable data.  
!                       There must be a one-to-one correspondance between the 
!           array elements in X and Y.
!
!               NDATA:  Integer containing the number of elements in the X 
!                       and Y arrays.
!
!               SIG:    The uncertainty in the Y data.
!
!               MWT:    A switch to let the routine know if standard deviations
!                       (i.e., uncertainties) are known on input for the Y data.
!
!                       MWT = 0:  Standard deviations are unavailable.
!                       MWT > 0:  Standard deviations are available and stored 
!                                 in SIG.
!
! Outputs     : A:      The y-intercept of the fitted straight line.
!
!               B:      The slope of the fitted straight line.
!
!               SIGA:   The uncertainty in the value of A:  A+/-SIGA.
!
!               SIGB:   The uncertainty in the value of B:  B+/-SIGB.
!
!               CHI2:   The chi-square of the goodness-of-fit.
!
!               Q:      The probability that the goodness-of-fit of the data to  
!                       a straight line is believable.  Note that if 
!
!                       Q > 0.1:    The goodness-of-fit is believable.
!                       0.001 < Q < 0.1:  The fit may be acceptable if the errors 
!                                   are nonnormal or have been moderately 
!                                   underestimated.
!                       Q < 0.001:  The goodness-of-fit is questionable.
!
! Calls       : Function GAMMQ.
!
! Common      : None.
!
! Restrictions: FORTRAN 77 coding.
!
! Side effects: None.
!
! Category    : Data fitting.
!
! Prev. Hist. : Based on the FIT subroutine of Numerical Recipes for FORTRAN 77.
!
! Written     : Donald G. Luttermoser, ETSU/Physics, 2 October 2013.
!
! Modified    : Version 1, Donald G. Luttermoser, ETSU/Physics, 2 Oct 2013
!           Initial program.
!
! Version     : Version 1,  2 October 2013.
!
!-
      SUBROUTINE FIT(X, Y, NDATA, SIG, MWT, A, B, SIGA, SIGB, CHI2, Q)
!
      INTEGER MWT, NDATA
      REAL A, B, CHI2, Q, SIGA, SIGB, SIG(NDATA), X(NDATA), Y(NDATA)
!
      INTEGER I
      REAL SIGDAT, SS, ST2, SX, SXOSS, SY, T, WT
!
! Initialize sums to zero.
!
      SX = 0.
      SY = 0.
      ST2 = 0.
      B = 0.
!
! Accumulate sums
!
      IF (MWT .NE. 0) THEN
          SS = 0.
!
!  with weights:
!
          DO 11 I = 1, NDATA
              WT = 1. / (SIG(I)**2)
              SS = SS + WT
              SX = SX + X(I)*WT
              SY = SY + Y(I)*WT
 11       CONTINUE
      ELSE
!
!  or without weights:
!
          DO 12 I = 1, NDATA
              SX = SX + X(I)
              SY = SY + Y(I)
 12       CONTINUE
          SS = FLOAT(NDATA)
      ENDIF
!
      SXOSS = SX / SS
!
      IF (MWT .NE. 0) THEN
          DO 13 I = 1, NDATA
              T = (X(I) - SXOSS) / SIG(I)
              ST2 = ST2 + T*T
              B = B + T*Y(I) / SIG(i)
13        CONTINUE
      ELSE
          DO 14 I = 1, NDATA
              T = X(I) - SXOSS
              ST2 = ST2 + T*T
              B = B + T*Y(I)
 14       CONTINUE
      ENDIF
!
! Solve for A, B, SIGA, and SIGB.
!
      B = B / ST2
      A = (SY - SX*B) / SS
      SIGA = SQRT((1. + SX*SX / (SS*ST2)) / SS)
      SIGB = SQRT(1. /ST2)
!
! Calculate chi-square.
!
      CHI2 = 0.
      Q = 1.
      IF (MWT .EQ. 0) THEN
          DO 15 I = 1, NDATA
              CHI2 = CHI2 + (Y(I) - A - B*X(I))**2
 15       CONTINUE
          SIGDAT = SQRT(CHI2 / (FLOAT(NDATA-2)))
          SIGA = SIGA * SIGDAT
          SIGB = SIGB * SIGDAT
      ELSE
          DO 16 I = 1, NDATA
              CHI2 = CHI2 + ((Y(I) - A - B*X(I)) / SIG(I))**2
 16       CONTINUE
          IF (NDATA .GT. 2) Q = gammq(0.5*(FLOAT(NDATA-2)), 0.5*CHI2)
      ENDIF
!
      RETURN
      END SUBROUTINE FIT

      FUNCTION gammq(a,x)
      REAL a,gammq,x
!U    USES gcf,gser
      REAL gammcf,gamser,gln
      if(x.lt.0..or.a.le.0.)pause 'bad arguments in gammq'
      if(x.lt.a+1.)then
        call gser(gamser,a,x,gln)
        gammq=1.-gamser
      else
        call gcf(gammcf,a,x,gln)
        gammq=gammcf
      endif
      return
      END FUNCTION
      SUBROUTINE gcf(gammcf,a,x,gln)
      INTEGER ITMAX
      REAL a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)
!U    USES gammln
      INTEGER i
      REAL an,b,c,d,del,h
      gln=gammln(a)
      b=x+1.-a
      c=1./FPMIN
      d=1./b
      h=d
      do 11 i=1,ITMAX
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b+an/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        del=d*c
        h=h*del
        if(abs(del-1.).lt.EPS)goto 1
11    continue
      pause 'a too large, ITMAX too small in gcf'
1     gammcf=exp(-x+a*log(x)-gln)*h
      return
      END SUBROUTINE
      SUBROUTINE gser(gamser,a,x,gln)
      INTEGER ITMAX
      REAL a,gamser,gln,x,EPS
      PARAMETER (ITMAX=100,EPS=3.e-7)
!CU    USES gammln
      INTEGER n
      REAL ap,del,sum
      gln=gammln(a)
      if(x.le.0.)then
        if(x.lt.0.)pause 'x < 0 in gser'
        gamser=0.
        return
      endif
      ap=a
      sum=1./a
      del=sum
      do 11 n=1,ITMAX
        ap=ap+1.
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*EPS)goto 1
11    continue
      pause 'a too large, ITMAX too small in gser'
1     gamser=sum*exp(-x+a*log(x)-gln)
      return
      END SUBROUTINE

      FUNCTION gammln(xx)
      REAL gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,        &
     &24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,&
     &-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return

      END FUNCTION

 subroutine  consective_dry(fields,cdd,dry_lim,ntimes,nx,ny)
   IMPLICIT NONE
   INTEGER, intent(in)::ntimes,nx,ny
   REAL, intent(in)   ::dry_lim
   REAL,DIMENSION(ntimes,nx,ny),intent(in)::fields 
   REAL,DIMENSION(nx,ny),intent(out)::cdd
   !local
   INTEGER::i,k,j
   INTEGER::numcdd,maxnumcdd
   do i=1,nx
     do j=1,ny
       maxnumcdd=0
       numcdd=0
       do k=1,ntimes
         IF(fields(k,i,j)<dry_lim) THEN
           numcdd=numcdd+1
         ELSE
           IF(numcdd>maxnumcdd)THEN
             maxnumcdd=numcdd
           ENDIF 
           numcdd=0
         ENDIF
       end do
       cdd(i,j)=maxnumcdd
     end do
   end do


 end subroutine  consective_dry

 SUBROUTINE  precp_extrem(fields,r95t_hist,rd,r10,r5d,sdii,r95t,dry_lim,ntimes,nx,ny)
   IMPLICIT NONE
   REAL,PARAMETER     :: R_10=10
   INTEGER, intent(in)::ntimes,nx,ny
   REAL, intent(in)   ::dry_lim
   REAL,DIMENSION(ntimes,nx,ny),intent(in)::fields 
   REAL,DIMENSION(nx,ny),intent(in ),optional::R95T_hist
   REAL,DIMENSION(nx,ny),intent(out)::RD
   REAL,DIMENSION(nx,ny),intent(out)::R10 !No. of days with precipitation larger 10 mm/day
   REAL,DIMENSION(nx,ny),intent(out)::R5D !Maximum 5 d precipitation total
   REAL,DIMENSION(nx,ny),intent(out)::SDII!Simple daily intensity index: annual total/number of Rday
   REAL,DIMENSION(nx,ny),intent(out)::R95T!Fraction of annual total precipitation due to events exceeding the 1961-1990 95th percentile
   !local
   REAL   ::temp
   INTEGER::i,k,j
   INTEGER::rd_pt,rd_10
   SDII=0.0
   R95T=0.0
   do i=1,nx
     do j=1,ny
       rd_pt=0
       rd_10=0
       temp=sum(fields(1:5,i,j))
       r5d(i,j)=temp
       do k=1,ntimes

         IF(fields(k,i,j)>dry_lim) THEN
           rd_pt=rd_pt+1
           SDII(i,j)=SDII(i,j)+fields(k,i,j)
         ENDIF

         IF(fields(k,i,j)>R_10) THEN
           rd_10=rd_10+1
         ENDIF

         IF(k.gt.1.and.k.le.ntimes-4) THEN
           temp=temp-fields(k-1,i,j)+fields(k+4,i,j)
           IF(temp>r5d(i,j)) THEN
             r5d(i,j)= temp
           ENDIF
         ENDIF

       end do
       rd(i,j) =rd_pt
       r10(i,j)=rd_10

       DO k=1,ntimes
         IF (present(R95T_hist)) THEN
           IF(fields(k,i,j)>R95T_hist(i,j).and.fields(k,i,j)>dry_lim) THEN
  !          print*,fields(k,i,j),R95T_hist(i,j)
             R95T(i,j)=R95T(i,j)+fields(k,i,j)
           ENDIF
         ENDIF
       END DO

       IF(rd_pt>0) THEN
         IF (present(R95T_hist)) THEN
           R95T(i,j)=R95T(i,j)/SDII(i,j)  ! convert to the percentage in whole year
         ENDIF
         SDII(i,j)=SDII(i,j)/rd(i,j)  ! calculate the average over precipiation day
       ENDIF
       r5d(i,j)= r5d(i,j)/5.0
     end do
   end do


 END SUBROUTINE  precp_extrem

 SUBROUTINE  quantile_cal(pre_quantile,quantile,dry_lim,qvalue,ntimes,nx,ny)
!f2py intent(in) qvalue,dry_lim,pre_quantile
!f2py intent(out) quantile
   IMPLICIT NONE 
   INTEGER, intent(in)::ntimes,nx,ny
   REAL, intent(in)   ::qvalue
   REAL, intent(in)   ::dry_lim
   REAL,DIMENSION(ntimes,nx,ny),intent(in)::pre_quantile 
   REAL,DIMENSION(nx,ny),intent(out)::quantile 
   !local
   INTEGER::i,k,j,rd
   REAL,DIMENSION(ntimes)::sorted

   quantile(:,:)=0.0
   DO j=1,ny
     DO i=1,nx
       rd=0
       sorted=0.0
       DO k=1,ntimes
        if (pre_quantile(k,i,j)>dry_lim) then
          rd=rd+1
          sorted(rd)=pre_quantile(k,i,j)
        endif
       END DO
       IF (rd>2) THEN
         call PIKSRT(rd,sorted(1:rd))
         quantile(i,j)=sorted(min(max(1,floor(rd*qvalue)),rd))
       ELSEIF(rd>1) then
         quantile(i,j)=sorted(1)*0.05+sorted(2)*0.95
       ELSEIF(rd>0) then
         quantile(i,j)=sorted(1)
       ENDIF
     END DO
   END DO

 END SUBROUTINE  quantile_cal
 !*****************************************************
 !* Sorts an array ARR of length N in ascending order *
 !* by straight insertion.                            *
 !* ------------------------------------------------- *
 !* INPUTS:                                           *
 !*      N   size of table ARR                  *
 !*          ARR   table to be sorted                 *
 !* OUTPUT:                                           *
 !*      ARR   table sorted in ascending order    *
 !*                                                   *
 !* NOTE: Straight insertion is a N² routine and      * 
 !*       should only be used for relatively small    *
 !*       arrays (N<100).                             *
 !*****************************************************         
 SUBROUTINE PIKSRT(N,ARR)
  IMPLICIT NONE 
  INTEGER::N
  REAL ARR(N)
  REAL:: a
  INTEGER:: i,j
  DO j=2, N
    a=ARR(j)
    DO i=j-1,1,-1
      IF (ARR(i)<=a) goto 10
      ARR(i+1)=ARR(i)
    end DO
    i=0  
10  ARR(i+1)=a
  end DO
 return
 END SUBROUTINE PIKSRT

 SUBROUTINE indexx (n,arrin,indx)

!     Use the Heapsort algorithm to index an array arrin of length n.
!     Output the array indx such that arrin(indx(j)) is in ascending
!     order for j = 1,2,...,n.  The input quantities n and arrin are
!     not changed.

!     This is a Numerical Recipes routine, but modified by one
!     line to work if n equals 1.

  integer i, n, indx(n), indxt, ir, j, l
  real arrin(n), q

  do 11 j=1,n
    indx(j)=j
11    continue
  if (n .eq. 1) return
  l=n/2+1
  ir=n
10    continue
    if(l.gt.1)then
      l=l-1
      indxt=indx(l)
      q=arrin(indxt)
    else
      indxt=indx(ir)
      q=arrin(indxt)
      indx(ir)=indx(1)
      ir=ir-1
      if(ir.eq.1)then
        indx(1)=indxt
        return
      endif
    endif
    i=l
    j=l+l
20      if(j.le.ir)then
      if(j.lt.ir)then
        if(arrin(indx(j)).lt.arrin(indx(j+1)))j=j+1
      endif
      if(q.lt.arrin(indx(j)))then
        indx(i)=indx(j)
        i=j
        j=j+j
      else
        j=ir+1
      endif
    go to 20
    endif
    indx(i)=indxt
  go to 10
  END SUBROUTINE indexx
END MODULE CS_STAT

