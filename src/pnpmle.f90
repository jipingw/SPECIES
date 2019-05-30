
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!! Nwl: penalized conditional NPMLE approach from Wang and Lindsay 2005
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

SUBROUTINE Nwl(n,t, MLE,p,pi,noZeroP)
! SUBROUTINE Nwl(n,t, MLE)


  !       ========================================================================================================
  !       purpose:  penalized NPMLE by Wang and Lindsay 2005
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
  !       diagram
  !
  !       Nwl: main subroutine wrapper
  !       {      
  !              call penNPMLE: the penalized NPMLE algorithm with adaptive penalty by linearlization
  !              {
  !                     call NPMLEpen :penalized VDM/EM alogrithm given a linear and fixed penalty
  !                     {                         
  !                          call EMNP_theta : the EM step for VDM/EM
  !                         call wBisectionNP_theta: bisection method for VDM part
  !                     }
  !              }  
  !       }         
  !       ========================================================================================================


  IMPLICIT NONE
  INTEGER t,noZeroP
  DOUBLE PRECISION n(50),MLE
  DOUBLE PRECISION p(10),pi(10)
  MLE=0.D0;p=0.D0*p;pi=0.D0*pi;noZeroP=0
  CALL penNPMLE(n,t,MLE,p,pi,noZeroP)
  MLE=INT(MLE)
END SUBROUTINE Nwl


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! penNPMLE: the penalized likelihood has a quadratic penalty.
!!           Maximization is done by a iterative linearization procedure 
!!           from Wang 2008, such that in each iteration the likelihood with
!!           linear penalty is maximized.
!!           For each linear penalty, subroutine NPMLEpen is use.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE penNPMLE(n,t,MLE,p,pi,noZeroP)

  !       ========================================================================================================
  !       purpose:   the penalized NPMLE algorithm with adaptive penalty by linearlization
  !       Input :    same as Nwl
  !       Called by: NwL
  !
  !       diagram:
  !       penNPMLE: linearlization subroutine 
  !       {
  !               call NPMLEpen :penalized VDM/EM alogrithm given a linear and fixed penalty
  !               {                         
  !                     call EMNP_theta : the EM step for VDM/EM
  !                     call wBisectionNP_theta: bisection method for VDM part
  !               }
  !       }         
  !       ========================================================================================================


  IMPLICIT NONE
  INTEGER t,ite,noZeroP,STOP
  DOUBLE PRECISION p(10),pi(10),n(50)
  DOUBLE PRECISION Etol,Gtol,gap,MLE,gamma2
  DOUBLE PRECISION test,chao,thetaOld,theta_0

  noZeroP=0;Gtol=5.D-3;gap=2.D-2;test=1.D0;Etol=1.D-10

  chao=n(1)**2/2/n(2)/SUM(n(1:t))  !!lower bound estimator of chao
  theta_0=chao                     !!theta_0 is the odds parameter
  STOP=0
  gamma2=(0.5/theta_0-0.5/(1+theta_0))   !!gamma2 is the linearized penalty factor
  ite=1

  DO WHILE (STOP==0)
     thetaOld=theta_0
     theta_0=0
     CALL NPMLEpen(n,p,pi,noZeroP,t,test,Etol,Gtol,gap,gamma2)
     theta_0=SUM(1/(EXP(p(1:noZeroP))-1)*pi(1:nozeroP))    
     theta_0=0.5*theta_0+0.5*thetaOld
     gamma2=(theta_0-chao)/chao
     IF(gamma2<0) THEN
        theta_0=(chao+thetaOld)/2
        gamma2= (theta_0-chao)/chao
     END IF
     ite=ite+1
     IF (ABS(thetaOld-theta_0)<.1/SUM(n) .OR. theta_0>100 .OR. ite>500) STOP=1
  END DO
  MLE=SUM(n(1:t))*(1+theta_0)+SUM(n((t+1):50))
END SUBROUTINE penNPMLE

