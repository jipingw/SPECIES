
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!! WLunpmle: using penalized method for approximate unconditional NPMLE by Wang and Lindsay 2005
!!           in this method, the marginal likelihood is approximated by a penalty 
!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

SUBROUTINE WLunpmle(n,t,MLE,p,pi,noZeroP)

  !       ========================================================================================================
  !       purpose: unconditional NPMLE using approximate method by Wang and Lindsay 2005.
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
  !       WLunpmle: main subroutine wraper
  !       {      
  !              linearization of the penalty;
  !
  !              call NPMLEpen: the EM algorithm with linear penalty
  !              {             
  !                  call untrunEMNP_theta : the EM step for VDM/EM
  !                  call untrunwBisectionNP_theta: bisection method for VDM part
  !              }
  !       }  
  !       
  !       =======================================================================================================

  IMPLICIT NONE
  INTEGER t,ite,noZeroP,STOP
  DOUBLE PRECISION p(10),pi(10), n(50)
  DOUBLE PRECISION Etol,Gtol,gap,MLE,gamma2
  DOUBLE PRECISION test,chao,thetaOld,theta_0

  noZeroP=0;Gtol=5.D-3;gap=2.D-2;test=1.D0;Etol=1.D-10
  MLE=0.D0;p=0.D0*p;pi=0.D0*pi;noZeroP=0

  !!lower estimator chao

  chao=n(1)**2/2/n(2)/SUM(n(1:t))
  theta_0=chao
  STOP=0
  gamma2=(0.5/theta_0-0.5/(1+theta_0))
  ite=1
  noZeroP=0
  DO WHILE (STOP==0)
     thetaOld=theta_0
     theta_0=0
     p=0.D0*p
     pi=0.D0*pi
     CALL NPMLEpen(n,p,pi,noZeroP,t,test,Etol,Gtol,gap,gamma2)
     !print*,'p=',p
     !print*,'pi=',pi
     theta_0=SUM(1/(EXP(p(1:noZeroP))-1)*pi(1:nozeroP))    
     !print*,'theta_0'
     gamma2=(0.5/theta_0-0.5/(1+theta_0))

     IF(gamma2<0) THEN
        theta_0=(chao+thetaOld)/2
        gamma2=(0.5/theta_0-0.5/(1+theta_0))
     END IF

     ite=ite+1
     !print*,thetaOld,theta_0

     IF (ABS(thetaOld-theta_0)<.01/SUM(n) .OR. (ite>5 .AND. theta_0>20*SUM(n)) .OR. ite>500) STOP=1
  END DO
  MLE=INT(SUM(n(1:t))*(1+theta_0))
END SUBROUTINE WLunpmle
