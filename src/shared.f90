
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!! 
!! NPMLEpen: a subroutines shared by WLnpmle,pnpmle, and PPCG
!!           it does penalized conditional NPMLE with a linear penalty for the odds parameter
!!           call EMNP_theta: 
!!           call wBisectionNP_theta:
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

SUBROUTINE NPMLEpen (n,p,pi,noZeroP,t,test,Etol,Gtol,gap,gamma2)

  !       ========================================================================================================
  !       purpose: the penalized NPMLE algorithm with adaptive penalty by linearlization
  !       Input :  n:       frequency of frequency data, defined of length 50.
  !                         if the actual length of n is shorter than 50, n is padded with 0 in the end.
  !                         if n is longer than 50, it's truncated at 50 from the main program of R codes.
  !                t:       cutoff defining the relative rare species (e.g, x<=t) to be used to estimate n0
  !                p:       vector for Poisson mean parameters sorted in ascending order 
  !                pi:      the weights of the corresponding Poisson mixture components
  !                p,pi:    defined of length 10 since given n of length<=50, 
  !                         the true number of components never exceeds 10
  !               noZeroP:  number of components in p that is non-zero.
  !                test:    maximum gradient achieved at convengence of algorithm
  !                Etol:    threshold value for convergence in terms of log likelihood in the EM part
  !                Gtol:    threshold value for convergence of entire algorithm in terms of gradient
  !                gap:     threshold value to control which components to be merged if two components are close.
  !                gamma2:  the penalty factor in the penalized likelihood with linearied penalty.
  !
  !       Called by: WLnpmle, and pnpmle, and PPCG
  !
  !       diagram:

  !       NPMLEpen :penalized VDM/EM alogrithm given a linear and fixed penalty
  !       {                         
  !             call EMNP_theta : the EM step for VDM/EM
  !             call wBisectionNP_theta: bisection method for VDM part
  !       }

  !       ========================================================================================================
 
  IMPLICIT NONE
  INTEGER t,i,j,noZeroP,loc,ite
  DOUBLE PRECISION p(10),pi(10),pvec(1000+(t-1)*100),drv(1000+(t-1)*100),Lq(50)
  DOUBLE PRECISION Etol, Gtol,p_add,gap,sum_a,sum_b
  DOUBLE PRECISION n(50),w,test,Pmix,Pden,gamma2,currentL

  test=1.D0

  DO i=1,1000   
     pvec(i)=1D-3*i
  END DO

  DO i=1,(t-1)*100
     pvec(i+1000)=1.+1D-2*i
  END DO

  !IF(noZeroP==0) THEN
  p=p*0.D0
  pi=pi*0.D0
  sum_a=0.0D0; sum_b=0.0D0
  DO i=1,t+1
     sum_a=sum_a + n(i)*(i-1); sum_b = sum_b + n(i)
  END DO
  p_add=sum_a/sum_b; p(1)=p_add; pi(1)=1
  p(1)=p_add; pi(1)=1.0D0
  noZeroP=1
  !END IF

  CALL EMNP_theta(n,p,pi,noZeroP,Etol,t,gamma2)
  CALL checkGap(p,pi,gap,noZeroP) 

  DO i=1,t
     Lq(i)=Pmix(i,p,pi,noZeroP)
  END DO

  DO i=1,SIZE(pvec)
     DO j=1,t
        drv(i)=drv(i)+n(j)*(Pden(j,pvec(i))/Lq(j)-1)
     END DO
     drv(i)=drv(i)-gamma2*(1/(EXP(pvec(i))-1)-SUM(pi(1:noZeroP)*(1/(EXP(p(1:noZeroP))-1))))
  END DO

  loc=MAXLOC(drv,1)
  p_add=pvec(loc)

  test=MAXVAL(drv)
  w=1.D0
  ite=1

  DO WHILE(test>Gtol .AND. ite<50)
     CALL  wBisectionNP_theta (n,p,pi,noZeroP,p_add,Lq,w,t,gamma2)
     p(noZeroP+1)=p_add
     pi(1:noZeroP)=(1-w)*pi(1:noZeroP)
     pi(noZeroP+1)=w
     noZeroP=noZeroP+1
     CALL EMNP_theta(n,p,pi,noZeroP,Etol,t,gamma2)
     currentL=-gamma2*SUM(pi(1:noZeroP)/(EXP(p(1:noZeroP))-1))

     DO i=1,t
        currentL=currentL+n(i)*LOG(Pmix(i,p,pi,noZeroP))
     END DO
     CALL checkGap(p,pi,gap,noZeroP) 
     DO i=1,t
        Lq(i)=Pmix(i,p,pi,noZeroP)
     END DO

     drv=drv*0.D0
     DO i=1,SIZE(pvec)
        DO j=1,t
           drv(i)=drv(i)+n(j)*(Pden(j,pvec(i))/Lq(j)-1)
        END DO
        drv(i)=drv(i)-gamma2*(1/(EXP(pvec(i))-1)-SUM(pi(1:noZeroP)/(EXP(p(1:noZeroP))-1)))
     END DO

     loc=MAXLOC(drv,1)
     p_add=pvec(loc)
     test=MAXVAL(drv)
     ite=ite+1
  END DO
  CALL sortP(p,pi,noZeroP)
END SUBROUTINE NPMLEpen



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!EMNP_theta: EM algorithm in the VDM/ECM algorithm for penalized likelihood where the 
!!            penalized functional is the odds parameter with penalizing factor gamma2
!!            This subourtine is called by WLnpmle, and pnpmle, and PPCG
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EMNP_theta(n,p,pi,noZeroP,Etol,t,gamma2)

  !       ========================================================================================================
  !       purpose:   EM algorithm updating p and pi given a linear penalty with fixed penalty factor gamma2
  !       Input :    same as  NPMLEpen
  !
  !
  !       Called by:  NPMLEpen
  !
  !       ========================================================================================================


  IMPLICIT NONE
  INTEGER noZeroP,ite,ite2,i,j,t,g
  DOUBLE PRECISION p(10),pi(10),delta,delta0,count_data(noZeroP)
  DOUBLE PRECISION n(50),sum_a,sum_b,Etol,z(t,noZeroP),d1,d2
  DOUBLE PRECISION currentL, last,lambda0,lambda1,Pmix,Pden,gamma2

  currentL=0.0D0
  ite=1
  g=noZeroP
  currentL=-gamma2*SUM(pi(1:noZeroP)/(EXP(p(1:noZeroP))-1))

  DO i=1,t
     currentL=currentL+n(i)*LOG(Pmix(i,p,pi,noZeroP))
  END DO

  last=2*currentL
  DO WHILE(ABS(currentL-last)>Etol .AND. ite<5000)
     last=currentL
     DO i=1,t
        DO j=1,noZeroP
           z(i,j)=pi(j)*Pden(i,p(j))/Pmix(i,p,pi,noZeroP)
        END DO
     END DO

     delta=-SUM(n(1:t))
     delta0=delta+1

     DO i=1,g
        count_data(i)=SUM(z(:,i)*n(1:t))
     END DO

     DO WHILE (ABS(delta-delta0)>0.0001)
        delta0=delta
        d1=-SUM(count_data/(delta-gamma2/(EXP(p(1:noZeroP))-1)))-1
        d2=SUM(count_data/(delta-gamma2/(EXP(p(1:noZeroP))-1))**2)
        delta=delta-d1/d2
     END DO

     pi=-count_data/(delta-gamma2/(EXP(p(1:noZeroP))-1))
     currentL=-gamma2*SUM(pi(1:noZeroP)/(EXP(p(1:noZeroP))-1))
     DO i=1,t
        currentL=currentL+n(i)*LOG(Pmix(i,p,pi,noZeroP))
     END DO

     DO i=1,noZeroP
        IF(p(i) .GT. 1.D-5) THEN
           lambda0=p(i)+1.D-3
           lambda1=p(i)
           ite2=1
           DO WHILE(ABS(lambda1-lambda0)>0.00000000001 .AND. ite2<2000)
              sum_a=0.0D0
              sum_b=0.0D0
              lambda0=lambda1
              DO j=1,t
                 sum_a=sum_a-n(j)*z(j,i)/(1-EXP(-lambda1))+n(j)*z(j,i)*j/lambda1
                 sum_b=sum_b+n(j)*z(j,i)/(1-EXP(-lambda1))**2*EXP(-lambda1)-n(j)*z(j,i)*j/lambda1**2
              END DO
              sum_a=sum_a+gamma2*pi(i)*EXP(lambda1)*(EXP(lambda1)-1)**(-2)
              sum_b=sum_b+gamma2*pi(i)*EXP(lambda1)*(EXP(lambda1)-1)**(-2)-2*gamma2*pi(i)*EXP(2*lambda1)*(EXP(lambda1)-1)**(-3)
              lambda1=lambda0-sum_a/sum_b
              ite2=ite2+1
              IF (lambda1<0 .OR. lambda1>t) ite2=2000
           END DO

           IF (lambda1==0 .OR. lambda1<0 .OR. lambda1>100 .OR. sum_b .EQ. 0.D0) THEN
              lambda1=p(i)
              lambda0=lambda1
           END IF

           p(i)=lambda1
        END IF
     END DO

     currentL=-gamma2*SUM(pi(1:noZeroP)/(EXP(p(1:noZeroP))-1))
     DO i=1,t
        currentL=currentL+n(i)*LOG(Pmix(i,p,pi,noZeroP))
     END DO
     ite=ite+1
  END DO
END SUBROUTINE EMNP_theta



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! wBisectionNP_theta: VDM part in the VDM/ECM algorithm for penalized likelihood where the 
!!            penalized functional is the odds parameter with penalizing factor gamma2.
!!            This subroutine find optimal weight for the new component added to the mixture.
!!            This subourtine is called by WLnpmle, and pnpmle, and PPCG
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE wBisectionNP_theta(n,p,pi,noZeroP,p_add,Lq,w,t,gamma2)

  !       ========================================================================================================
  !       purpose:   optimal weight for the newly added component in VDM part of the VDM/ECM algorithm
  !       Input :    same as  PCGpen
  !                  p_add --- newly added component
  !                  Lq    --- PCG density given the current p and pi
  !                  w     --- weight of the new component p_add in the updated mixture   
  !
  !       Called by: PCGpen
  !
  !       ========================================================================================================

  IMPLICIT NONE
  INTEGER noZeroP,ite,t,i
  DOUBLE PRECISION p(10),pi(10),temp_p(10),temp_pi(10),n(50),p_add,Lq(t),w
  DOUBLE PRECISION D,L(t),lower, upper,c,Pden,Pmix,gamma2

  lower=0.D0
  upper=1.0D0
  w=0.5D0
  D=1.0D0
  ite=1

  DO WHILE (ABS(D)>0.0000000001 .AND. ite<2000)
     temp_p=p;temp_pi=pi
     temp_p(noZeroP+1)=p_add
     temp_pi(1:noZeroP)=(1-w)*pi(1:noZeroP)
     temp_pi(noZeroP+1)=w
     ite=ite+1

     D=0.0D0
     DO i=1,t
        L(i)=Pmix(i,temp_p,temp_pi,noZeroP+1)
     END DO

     !print*, L
     DO i=1,t
        c=n(i)*(-Lq(i)+Pden(i,p_add))/L(i)
        D=D+c
     END DO

     D=D+gamma2*(SUM(pi(1:noZeroP)/(EXP(p(1:noZeroP))-1))-1/(EXP(p_add)-1))

     IF (D>0)  lower=w
     IF (D<0)  upper=w
     w=(lower+upper)/2
     !!print*,'D=',D
     !print*,'w=',w
     IF (ABS(w)>500 .OR. ite>=60) THEN
        w=0.01D0
     END IF
  END DO
END SUBROUTINE wBisectionNP_theta




!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!! some common functions or subroutines
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

SUBROUTINE checkGap(p,pi,gap,noZeroP)

  !       ==========================================================================
  !       purpose: merge the mixture components if two are very close
  !       Input :  p   --- mean parameter vector for Poisson components
  !                pi  --- weight vector for mixture
  !                gap --- threshold for how close two components are.
  !                noZeroP --- number of non-zero components in p/pi vector,integer
  !       ==========================================================================


  IMPLICIT NONE
  DOUBLE PRECISION p(10), pi(10), gap
  INTEGER noZeroP,i,j,k

  DO i=1,noZeroP-1
     j=i+1
     DO WHILE (j<=noZeroP)
        IF (ABS(p(j)-p(i))<gap) THEN
           noZeroP=noZeroP-1
           pi(i)=pi(i)+pi(j)
           DO k=j+1,noZeroP+1
              pi(k-1)=pi(k)
              p(k-1)=p(k)
           END DO
           p(noZeroP+1)=0.D0
           pi(noZeroP+1)=0.D0
        END IF
        IF (ABS(p(j)-p(i))>=gap) j=j+1
     END DO
  END DO

  DO i=1,noZeroP
     IF (pi(i)<0.0001) THEN
        pi(i)=0.D0
        p(i)=0.D0
        DO j=i+1,noZeroP
           p(j-1)=p(j)
           pi(j-1)=pi(j)
        END DO
        noZeroP=noZeroP-1
     END IF
     pi=pi/SUM(pi)
  END DO

END SUBROUTINE checkGap

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE sortP(p,pi,noZeroP)

  !       ==========================================================================
  !       purpose: sort the mixture components in ascending order
  !       Input :  p  --- mean parameter vector for Poisson components
  !                pi --- weight vector for mixture
  !                noZeroP --- number of non-zero components in p/pi vector,integer
  !       ==========================================================================


  IMPLICIT NONE
  DOUBLE PRECISION p(10), pi(10), p_temp, pi_temp
  INTEGER noZeroP,i,j

  !!first sort p, pi, and then eliminate boundary component

  DO j=1,noZeroP-1
     DO i=j+1,noZeroP
        IF (p(i)<p(j)) THEN
           p_temp=p(j)
           pi_temp=pi(j)
           p(j)=p(i)
           pi(j)=pi(i)
           p(i)=p_temp
           pi(i)=pi_temp
        END IF
     END DO
  END DO

END SUBROUTINE sortP


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DOUBLE PRECISION FUNCTION untrunPmix(k,p,pi,noZeroP)

  !       ==========================================================================
  !       purpose: Compute the Poisson mixture denstiy 
  !       Input :  k  --- argument  ( integer, >=0 )
  !                p  --- mean parameter vector for Poisson components
  !                pi --- weight vector for mixture
  !                noZeroP --- number of non-zero components in p/pi vector,integer
  !    
  !       Routine called: untrunPden
  !       ==========================================================================

  IMPLICIT NONE

  DOUBLE PRECISION p(10),pi(10),untrunPden
  INTEGER noZeroP,i,k
  untrunPmix=0.D0

  DO i=1, noZeroP
     untrunPmix=untrunPmix+pi(i)*untrunPden(k,p(i))
  END DO
END FUNCTION untrunPmix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DOUBLE PRECISION FUNCTION untrunPden(k,mu)

  !       ==========================================
  !       purpose: Compute the Poisson denstiy 
  !       Input :  k  --- argument  ( integer, >=0 )
  !                mu --- mean Parameter  ( q > 0 )
  !    
  !       Routine called: lgamma for log gamma function
  !       ==========================================

  IMPLICIT NONE
  DOUBLE PRECISION mu
  DOUBLE PRECISION lgamma
  INTEGER k
  untrunPden=EXP(-mu+k*LOG(mu)-lgamma((k+1)*1.D0))
  RETURN
END FUNCTION untrunPden

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DOUBLE PRECISION FUNCTION PmixScon(k,p,pi,a,noZeroP)

  !       ==========================================================================
  !       purpose: Compute the zero-truncted Poisson compound Gamma denstiy 
  !       Input :  k   --- argument  ( integer, >=0 )
  !                p  --- mean parameter vector for the compound Gamma components
  !                pi --- weight vector for mixture
  !                noZeroP --- number of non-zero components in p/pi vector,integer
  !                a   --- shape parameter of the Gamma
  !       Routine called: PmixScon 
  !       ==========================================================================

  IMPLICIT NONE

  DOUBLE PRECISION a,p(10),pi(10),PdenScon
  INTEGER noZeroP,i,k
  PmixScon=0.D0
  DO i=1, noZeroP
     PmixScon=PmixScon+pi(i)*PdenScon(k,p(i),a)
  END DO
END FUNCTION PmixScon

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DOUBLE PRECISION FUNCTION PdenScon(k,mu,a)
  !       ==========================================================================
  !       purpose: Compute the zero-truncted Poisson Gamma denstiy 
  !       Input :  k   --- argument  ( integer, >=0 )
  !                mu  --- mean parameter of the Gamma
  !                a   --- shape parameter of the Gamma
  !       Routine called: beta
  !       ==========================================================================

  IMPLICIT NONE
  DOUBLE PRECISION a,mu,beta1
  DOUBLE PRECISION BETA
  INTEGER k
  beta1= BETA(k*1D0,a)
  !!print*,'beta=', beta1
  PdenScon=1/beta1/k*EXP(a*LOG(a/mu)-(k+a)*LOG(1+a/mu))/(1-(a/(a+mu))**a)
  RETURN
END FUNCTION PdenScon

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


DOUBLE PRECISION FUNCTION Pmix(k,p,pi,noZeroP)
  !       ==========================================================================
  !       purpose: Compute the zero-truncted Poisson mixture denstiy 
  !       Input :  k  --- argument  ( integer, >=0 )
  !                p  --- mean parameter vector for Poisson components
  !                pi --- weight vector for mixture
  !                noZeroP --- number of non-zero components in p/pi vector,integer
  !    
  !       Routine called: Pden
  !       ==========================================================================

  IMPLICIT NONE

  DOUBLE PRECISION p(10),pi(10),Pden
  INTEGER noZeroP,i,k
  Pmix=0.D0

  DO i=1, noZeroP
     Pmix=Pmix+pi(i)*Pden(k,p(i))
  END DO
END FUNCTION Pmix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DOUBLE PRECISION FUNCTION Pden(k,mu)
  !       ==========================================
  !       purpose: Compute the zero-truncted Poisson denstiy 
  !       Input :  k  --- argument  ( integer, >=0 )
  !                mu --- mean Parameter  ( q > 0 )
  !    
  !       Routine called: lgamma for log gamma function
  !       ==========================================

  DOUBLE PRECISION mu
  DOUBLE PRECISION lgamma
  INTEGER k

  Pden=EXP(-mu+k*LOG(mu)-lgamma((k+1)*1D0))/(1D0-EXP(-mu))
  RETURN
END FUNCTION Pden

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DOUBLE PRECISION  FUNCTION BETA(p,q)

  !       ==========================================
  !       purpose: Compute the beta function B(p,q)
  !       Input :  p  --- Parameter  ( p > 0 )
  !                q  --- Parameter  ( q > 0 )
  !    
  !       function called: lgamma for computing ‚(x)
  !       ==========================================
  IMPLICIT NONE

  DOUBLE PRECISION p,q,lgamma
  beta=EXP(lgamma(p)+lgamma(q)-lgamma(p+q))

END FUNCTION BETA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

 
DOUBLE PRECISION FUNCTION lgamma(X)

  !----------------------------------------------------------------------
  !     ==================================================    
  !     Adapted From: http://www.netlib.org/specfun/
  !     References:
  !     1) W. J. Cody and K. E. Hillstrom, 'Chebyshev Approximations for
  !     the Natural Logarithm of the Gamma Function,' Math. Comp. 21,
  !     1967, pp. 198-203.
  !     2) K. E. Hillstrom, ANL/AMD Program ANLC366S, DGAMMA/LGAMMA, May,
  !     1969.
  !     3) Hart, Et. Al., Computer Approximations, Wiley and sons, New
  !     York, 1968.
  !     
  !     Authors: W. J. Cody and L. Stoltz
  !     Argonne National Laboratory
  !----------------------------------------------------------------------

  INTEGER I

  DOUBLE PRECISION C,&
       & CORR,D1,D2,D4,EPS,FRTBIG,FOUR,HALF,ONE,PNT68,P1,P2,P4,&
       & Q1,Q2,Q4,RES,SQRTPI,THRHAL,TWELVE,TWO,X,XBIG,XDEN,XINF,&
       & XM1,XM2,XM4,XNUM,Y,YSQ,ZERO
  DIMENSION C(7),P1(8),P2(8),P4(8),Q1(8),Q2(8),Q4(8)

  DATA ONE,HALF,TWELVE,ZERO/1.0D0,0.5D0,12.0D0,0.0D0/, &
       &     FOUR,THRHAL,TWO,PNT68/4.0D0,1.5D0,2.0D0,0.6796875D0/,&
       &     SQRTPI/0.9189385332046727417803297D0/

  DATA XBIG,XINF,EPS,FRTBIG/2.55D305,1.79D308,2.22D-16,2.25D76/
  
  DATA D1/-5.772156649015328605195174D-1/

  DATA P1/4.945235359296727046734888D0,2.018112620856775083915565D2,&
       &        2.290838373831346393026739D3,1.131967205903380828685045D4,&
       &        2.855724635671635335736389D4,3.848496228443793359990269D4,&
       &        2.637748787624195437963534D4,7.225813979700288197698961D3/

  DATA Q1/6.748212550303777196073036D1,1.113332393857199323513008D3,&
       &        7.738757056935398733233834D3,2.763987074403340708898585D4,&
       &        5.499310206226157329794414D4,6.161122180066002127833352D4,&
       &        3.635127591501940507276287D4,8.785536302431013170870835D3/

  DATA D2/4.227843350984671393993777D-1/
  DATA P2/4.974607845568932035012064D0,5.424138599891070494101986D2,&
       &        1.550693864978364947665077D4,1.847932904445632425417223D5,&
       &        1.088204769468828767498470D6,3.338152967987029735917223D6,&
       &        5.106661678927352456275255D6,3.074109054850539556250927D6/
  
  DATA Q2/1.830328399370592604055942D2,7.765049321445005871323047D3,&
       &        1.331903827966074194402448D5,1.136705821321969608938755D6,&
       &        5.267964117437946917577538D6,1.346701454311101692290052D7,&
       &        1.782736530353274213975932D7,9.533095591844353613395747D6/

  DATA D4/1.791759469228055000094023D0/
  DATA P4/1.474502166059939948905062D4,2.426813369486704502836312D6,&
       &        1.214755574045093227939592D8,2.663432449630976949898078D9,&
       &      2.940378956634553899906876D10,1.702665737765398868392998D11,&
       &      4.926125793377430887588120D11,5.606251856223951465078242D11/

  DATA Q4/2.690530175870899333379843D3,6.393885654300092398984238D5,&
       &        4.135599930241388052042842D7,1.120872109616147941376570D9,&
       &      1.488613728678813811542398D10,1.016803586272438228077304D11,&
       &      3.417476345507377132798597D11,4.463158187419713286462081D11/
  
  DATA C/-1.910444077728D-03,8.4171387781295D-04,&
       &     -5.952379913043012D-04,7.93650793500350248D-04,&
       &     -2.777777777777681622553D-03,8.333333333333333331554247D-02,&
       &      5.7083835261D-03/


  Y = X
  IF ((Y .GT. ZERO) .AND. (Y .LE. XBIG)) THEN
     IF (Y .LE. EPS) THEN
        RES = -LOG(Y)
     ELSE IF (Y .LE. THRHAL) THEN
        IF (Y .LT. PNT68) THEN
           CORR = -LOG(Y)
           XM1 = Y
        ELSE
           CORR = ZERO
           XM1 = (Y - HALF) - HALF
        END IF
        IF ((Y .LE. HALF) .OR. (Y .GE. PNT68)) THEN
           XDEN = ONE
           XNUM = ZERO
           DO 140 I = 1, 8
              XNUM = XNUM*XM1 + P1(I)
              XDEN = XDEN*XM1 + Q1(I)
140           CONTINUE
              RES = CORR + (XM1 * (D1 + XM1*(XNUM/XDEN)))
        ELSE
           XM2 = (Y - HALF) - HALF
           XDEN = ONE
           XNUM = ZERO
           DO 220 I = 1, 8
              XNUM = XNUM*XM2 + P2(I)
              XDEN = XDEN*XM2 + Q2(I)
220           CONTINUE
              RES = CORR + XM2 * (D2 + XM2*(XNUM/XDEN))
         END IF
     ELSE IF (Y .LE. FOUR) THEN
        XM2 = Y - TWO
        XDEN = ONE
        XNUM = ZERO
        DO 240 I = 1, 8
           XNUM = XNUM*XM2 + P2(I)
           XDEN = XDEN*XM2 + Q2(I)
240        CONTINUE
           RES = XM2 * (D2 + XM2*(XNUM/XDEN))
     ELSE IF (Y .LE. TWELVE) THEN
        XM4 = Y - FOUR
        XDEN = -ONE
        XNUM = ZERO
        DO 340 I = 1, 8
           XNUM = XNUM*XM4 + P4(I)
           XDEN = XDEN*XM4 + Q4(I)
340        CONTINUE
           RES = D4 + XM4*(XNUM/XDEN)
     ELSE 
        RES = ZERO
        IF (Y .LE. FRTBIG) THEN
           RES = C(7)
           YSQ = Y * Y
           DO 450 I = 1, 6
              RES = RES / YSQ + C(I)
450           CONTINUE
        END IF
        RES = RES/Y
        CORR = LOG(Y)
        RES = RES + SQRTPI - HALF*CORR
        RES = RES + Y*(CORR-ONE)
     END IF
  ELSE
     RES = XINF
  END IF

  LGAMMA = RES
  RETURN

END FUNCTION lgamma 
