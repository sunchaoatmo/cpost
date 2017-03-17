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
   REAL,DIMENSION(nx,ny),intent(in )::R95T_hist
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
         IF(fields(k,i,j)>R95T_hist(i,j).and.fields(k,i,j)>dry_lim) THEN
!          print*,fields(k,i,j),R95T_hist(i,j)
           R95T(i,j)=R95T(i,j)+fields(k,i,j)
         ENDIF
       END DO

       IF(rd_pt>0) THEN
         R95T(i,j)=R95T(i,j)/SDII(i,j)  ! convert to the percentage in whole year
         SDII(i,j)=SDII(i,j)/rd(i,j)  ! calculate the average over precipiation day
       ENDIF
       r5d(i,j)= r5d(i,j)/5.0
     end do
   end do


 END SUBROUTINE  precp_extrem

 SUBROUTINE  quantile_cal(pre_quantile,quantile,dry_lim,qvalue,ntimes,nx,ny)
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



