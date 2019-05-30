!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!! norrispollock:unconditional NPMLE by Norris and Pollock 1998
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

SUBROUTINE norrispollock(n,t,MLE,p,pi,noZeroP)

  !       ========================================================================================================
  !       purpose: unconditional NPMLE by Norris and Pollock 1998 and Bonhing and Scon 2005.
  !       Input :  n:       frequency of frequency data, defined of length 50.
  !                         if the actual length of n is shorter than 50, n is padded with 0 in the end.
  !                         if n is longer than 50, it's truncated at 50 from the main program of R codes.
  !                t:       cutoff defining the relative rare species (e.g, x<=t) to be used to estimate n0
  !                mle:     MLE of N
  !                p:       vector for Poisson mean parameters sorted in ascending order 
  !                pi:      the weights of the corresponding Poisson mixture components
  !                p,pi:    defined of length 10 since given n of length<=50, 
  !                         the true number of components never exceeds 10
  !               noZeroP:  number of components in p that is non-zero.
  !
  !       diagram:
  !
  !       norrispollock: main subroutine wrapper
  !       {      call unpmle: the EM algorithm 
  !              {
  !                   call unpmle_one : the M step, a VDM/EM alogrithm for NPMLE of Poisson mixture
  !                   {      call untrunEMNP_theta : the EM step for VDM/EM
  !                          call untrunwBisectionNP_theta: bisection method for VDM part
  !                   }
  !              }  
  !       }               
  !       =======================================================================================================

  IMPLICIT NONE
  INTEGER t,noZeroP
  DOUBLE PRECISION n(50),MLE
  DOUBLE PRECISION p(10),pi(10)

  MLE=0.D0;p=0.D0*p;pi=0.D0*pi;noZeroP=0
  CALL unpmle(n,t,MLE,p,pi,noZeroP)
  MLE=INT(MLE)  
END SUBROUTINE norrispollock

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!subroutine for the outside EM loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE unpmle(n,t,MLE,p,pi,noZeroP)
  !       ========================================================================================================
  !       purpose:    subroutine for the outside EM loop ofr norrispollock
  !       called by:  norrispollock
  !       Input :     same as norrispolloc
  !
  !       diagram of routines called:

  !       unpmle: the EM algorithm 
  !       {
  !              call unpmle_one : the M step, a VDM/EM alogrithm for NPMLE of Poisson mixture
  !              {      call untrunEMNP_theta : the EM step for VDM/EM
  !                     call untrunwBisectionNP_theta: bisection method for VDM part
  !              }
  !       }  
  !                     
  !       =======================================================================================================

  IMPLICIT NONE
  INTEGER noZeroP,t,bigIte,i
  DOUBLE PRECISION  likeEnd,likeStart,p(10),pi(10),n(50),n1(50),MLE
  DOUBLE PRECISION  like,tempprob,lgamma,n00, untrunPmix

  !! n contains the observed data
  !! n1 contains estimated n0 (first element) and n

  n1(2:50)=n(1:49)
  n1(1)=n1(2)**2/n1(3)/2   !!n0 started with Chao's lower bound estimate.
  n1(1)=n(1)/2             !! make it conservative so not to overshoot

  likeEnd=-10000000.D0
  likeStart=likeEnd-10000000
  bigIte=1
  MLE=0.D0
  noZeroP=0

  !! uncondtional NPMLE can be unbounded.
  !! the algorithm is terminated if the likelihood convergence criterion is met, or
  !! iteration number reaches 50000, or 
  !! the estimate of N is greater than 20 times of the observed total.

  DO WHILE( likeEnd-likeStart>1.D-10 .AND. bigIte<50000 .AND. SUM(n1)<20*SUM(n))
     likeStart=likeEnd
     bigIte=bigIte+1

     !! M-step
     CALL unpmle_one(n1,t,p,pi,noZeroP)   !! VDM/EM to maximize unconditional likelihood

     !! E-step
     like=0.D0
     DO i=1,t+1
        tempprob=0.D0
        tempprob=untrunPmix(i-1,p,pi,noZeroP)
        IF(i==1) THEN
           n00=1.D0*INT(SUM(n1(2:(t+1)))/(1-tempprob)-SUM(n1(2:(t+1))))
           n00=INT(n00)+1
           like=like+n00*LOG(tempprob)
        ELSE
           like=like+n1(i)*LOG(tempprob)
        END IF
     END DO

     like=like+lgamma(SUM(n1(2:(t+1)))+n00+1)-lgamma(n00+1)-lgamma(SUM(n1(2:(t+1)))+1)
     likeEnd=like
     !print*,'likeEnd=',likeEnd,'MLE=',sum(n1)
     IF(like>likeStart) THEN
        n1(1)=n00
     END IF
  END DO

  MLE=SUM(n1(1:50))
END SUBROUTINE unpmle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!unpmle_one: VDM/EM algorithm to maximize the unconditiona llikelihood
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE unpmle_one(n,t,p,pi,noZeroP)

  !       ========================================================================================================
  !       purpose:    M-step of the EM in norrispollock, VDM/EM algorithm unconditional MLE given n0
  !       called by:  unpmle
  !       Input :     same as norrispolloc
  !
  !       diagram of routines called:

  !        unpmle_one 
  !        {      call untrunEMNP_theta : the EM step for VDM/EM
  !               call untrunwBisectionNP_theta: bisection method for VDM part
  !        }
  !       =======================================================================================================


  IMPLICIT NONE

  INTEGER t,i,j,noZeroP,loc,ite
  DOUBLE PRECISION p(10),pi(10),pvec(1000+t*100),drv(1000+t*100),Lq(50)
  DOUBLE PRECISION Etol, Gtol,sum_a,sum_b,p_add,gap
  DOUBLE PRECISION n(50),w,test,untrunPmix,untrunPden,currentL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  gap=0.02         !! gap defines that if two components differ by <=0.02, they are merged.
  Etol=1.D-10      !! Etol: EM convergence criterion 
  Gtol=0.005       !! Gtol: VDM-EM converge criterion for the maximum gradient. 


  IF(noZeroP==0) THEN
     p=p*0.D0
     pi=pi*0.D0
     sum_a=0.0D0; sum_b=0.0D0
     DO i=1,t+1
        sum_a=sum_a + n(i)*(i-1); sum_b = sum_b + n(i)
     END DO

     p_add=sum_a/sum_b; p(1)=p_add; pi(1)=1
     p(1)=p_add; pi(1)=1.0D0

     noZeroP=1
  END IF

  CALL untrunEMNP_theta(n,p,pi,noZeroP,Etol,t)   !! EM step for untruncted Poisson

  DO i=1,1000   
     pvec(i)=1D-3*i
  END DO

  DO i=1,t*100
     pvec(i+1000)=1.+1D-2*i
  END DO

  DO i=1,t+1
     Lq(i)=untrunPmix(i-1,p,pi,noZeroP)
  END DO

  drv=drv*0.D0
  DO i=1,SIZE(pvec)
     DO j=1,t+1
        drv(i)=drv(i)+n(j)*(untrunPden(j-1,pvec(i))/Lq(j)-1)
     END DO
  END DO

  loc=MAXLOC(drv,1)
  p_add=pvec(loc)
  test=MAXVAL(drv)

  w=1D0
  ite=1

  DO WHILE(test>Gtol .AND. ite<100)

     CALL  untrunwBisectionNP_theta (n,p,pi,noZeroP,p_add,Lq,w,t)  !! bisection for VDM step for un-truncated Poisson
     p(noZeroP+1)=p_add
     pi(1:noZeroP)=(1-w)*pi(1:noZeroP)
     pi(noZeroP+1)=w
     noZeroP=noZeroP+1

     currentL=0.d0
     DO i=1,t+1
        currentL=currentL+n(i)*LOG(untrunPmix(i-1,p,pi,noZeroP))
     END DO

     CALL untrunEMNP_theta(n,p,pi,noZeroP,Etol,t)
     currentL=0.D0

     DO i=1,t+1
        currentL=currentL+n(i)*LOG(untrunPmix(i-1,p,pi,noZeroP))
     END DO

     CALL checkGap(p,pi,gap,noZeroP) 
     DO i=1,t+1
        Lq(i)=untrunPmix(i-1,p,pi,noZeroP)
     END DO

     drv=drv*0.D0

     DO i=1,SIZE(pvec)
        DO j=1,t+1
           drv(i)=drv(i)+n(j)*(untrunPden(j-1,pvec(i))/Lq(j)-1)
        END DO
     END DO

     loc=MAXLOC(drv,1)
     p_add=pvec(loc)
     test=MAXVAL(drv)
     ite=ite+1
  END DO

  CALL sortP(p,pi,noZeroP)

END SUBROUTINE unpmle_one


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! bisection part for embedded VDM part.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE untrunwBisectionNP_theta (n,p,pi,noZeroP,p_add,Lq,w,t)

  !       ========================================================================================================
  !       purpose:     bisection method to find optimal weight for newly added component in the VDM step
  !       called by:   unpmle_one
  !     
  !       Input :  same as norrispolloc
  !                p_add  ---  newly added component that achieves maximum gradient
  !                Lq     ---  the poisson mixture density vector of length (t+1) defined by p, pi
  !                w      ---  weigth for the new component
  !
  !       =======================================================================================================

  IMPLICIT NONE
  INTEGER noZeroP,ite,t,i
  DOUBLE PRECISION p(10),pi(10),temp_p(10),temp_pi(10),n(*),p_add,Lq(t+1),w
  DOUBLE PRECISION D,L(t+1),lower, upper,c,untrunPden,untrunPmix

  lower=0.D0
  upper=1.0D0
  w=0.5D0
  D=1.0D0
  ite=1
  !!print*,'here we are'

  DO WHILE (ABS(D)>0.00000000001 .AND. ite<20000)
     temp_p=p;temp_pi=pi
     temp_p(noZeroP+1)=p_add
     temp_pi(1:noZeroP)=(1-w)*pi(1:noZeroP)
     temp_pi(noZeroP+1)=w
     ite=ite+1

     D=0.0D0
     DO i=1,t+1
        L(i)=untrunPmix(i-1,temp_p,temp_pi,noZeroP+1)
     END DO

     !print*, L
     DO i=1,t+1
        c=n(i)*(-Lq(i)+untrunPden(i-1,p_add))/L(i)
        D=D+c
     END DO


     IF (D>0)  lower=w
     IF (D<0)  upper=w
     w=(lower+upper)/2
     !!print*,'D=',D
     !print*,'w=',w
     IF (ABS(w)>500 .OR. ite>=60) THEN
        w=0.01D0
     END IF
  END DO
END SUBROUTINE untrunwBisectionNP_theta


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! EM part for the embedded VDM-EM part.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE untrunEMNP_theta(n,p,pi,noZeroP,Etol,t)

  !       ========================================================================================================
  !       purpose:     EM step for updating p and pi optimally in the VDM/EM step
  !       called by:   unpmle_one
  !     
  !       Input : same as norrispolloc
  !               Etol  ---  threshold for likelihood convergence in EM algorithm
  !       =======================================================================================================

  IMPLICIT NONE
  INTEGER noZeroP,ite,i,j,t,g
  DOUBLE PRECISION p(10),pi(10)
  DOUBLE PRECISION n(50),Etol,z(t+1,noZeroP)
  DOUBLE PRECISION currentL, last,untrunPmix,untrunPden


  currentL=0.0D0
  ite=1
  g=noZeroP

  DO i=1,t+1
     currentL=currentL+n(i)*LOG(untrunPmix(i-1,p,pi,noZeroP))
  END DO

  last=2*currentL-1000
  DO WHILE(ABS(currentL-last)>Etol .AND. ite<50000)

     last=currentL
     DO i=1,t+1
        DO j=1,noZeroP
           z(i,j)=pi(j)*untrunPden(i-1,p(j))/untrunPmix(i-1,p,pi,noZeroP)
        END DO
     END DO

     DO i=1,g
        pi(i)=0.D0
        DO j=1,t+1
           pi(i)=pi(i)+ z(j,i)*n(j)
        END DO
        pi(i)=pi(i)/SUM(n(1:(t+1)))
     END DO

     currentL=0.D0
     DO i=1,t+1
        currentL=currentL+n(i)*LOG(untrunPmix(i-1,p,pi,noZeroP))
     END DO

     DO i=1,noZeroP
        p(i)=0.D0
        DO j=1,(t+1)
           p(i)=p(i)+(j-1)*n(j)*z(j,i)
        END DO
        p(i)=p(i)/SUM(n(1:(t+1))*z(1:(t+1),i))
     END DO

     currentL=0.D0
     DO i=1,t+1
        currentL=currentL+n(i)*LOG(untrunPmix(i-1,p,pi,noZeroP))
     END DO
     currentL=0.D0
     DO i=1,t+1
        currentL=currentL+n(i)*LOG(untrunPmix(i-1,p,pi,noZeroP))
     END DO
     ite=ite+1
  END DO
END SUBROUTINE  untrunEMNP_theta
