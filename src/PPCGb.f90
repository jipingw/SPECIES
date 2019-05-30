
SUBROUTINE PPCGb(n,t, MLE,amodel,p,pi,noZeroP)


  !       ========================================================================================================
  !       purpose: Poisson-compound Gamma for species richness estimation by Wang 2010
  !       Input :  n:       frequency of frequency data, defined of length 50.
  !                         if the actual length of n is shorter than 50, n is padded with 0 in the end.
  !                         if n is longer than 50, it's truncated at 50 from the main program of R codes.
  !                 t  :    cutoff defining the relative rare species (e.g, x<=t) to be used to estimate n0
  !               mle  :    MLE of N
  !                 p  :    vector for Poisson mean parameters sorted in ascending order 
  !                pi  :    the weights of the corresponding mixture components
  !              p,pi  :    defined of length 10 since given n of length<=50, 
  !                         the true number of components never exceeds 10
  !           noZeroP  :    number of components in p that is non-zero.
  !            amodel  :    amodel selected from cross-validation in the original data
  !
  !       diagram:
  !
  !       ppcgb: main subroutine wrapper
  !       {
  !                        call PCGpenE: NPMLE for PCG model with exponential prior,linearlization
  !                        {
  !                                 call PCGpen: NPMLE for PCG model with linear penalty 
  !                                 {
  !                                   call  wBisectionCon_theta : optimal weight for VDM part
  !                                   call  EMScon_theta:      ECM part for VDM/ECM algorithm
  !                                 }
  !                        }
  !                        call ENPMLE: NPMLE with exponential partial prior, linearization
  !                        {
  !                                 call NPMLEpen: NPMLE for fixed linear penalty
  !                        }
  !            
  !      }
  !       =======================================================================================================


  IMPLICIT NONE
  INTEGER t,noZeroP
  DOUBLE PRECISION n(50),MLE,amodel,a,gamma2,Etol,Gtol,gap,test
  DOUBLE PRECISION p(10),pistat(10),pi(10),sum,theta_threshold
 

  noZeroP=0;MLE=0.D0;p=p*0.D0;pi=0.D0*pi;noZeroP=0.D0
  theta_threshold=1./SUM(n(1:t))   !!threshold value to control the linearization  iterations
  gamma2=0.1
  Etol=1.D-10;Gtol=0.005;gap=0.01
  
  !CALL PCGone(n,t,MLE,a_model,alpha,alphaK,p,pi,theta_threshold,noZeroP)

  a=amodel
  if (a<500) then  
     CALL PCGpenE (n,p,pi,noZeroP,a,t,test,Etol,Gtol,gap,gamma2,theta_threshold)  
     pistat(1:noZeroP)=1/(1-(a/(a+p(1:noZeroP)))**a)*pi(1:noZeroP)/&
          &SUM(1/(1-(a/(a+p(1:noZeroP)))**a)*pi(1:noZeroP))
     mle= SUM(n((t+1):50))+SUM(n(1:t))*(1+SUM(pistat(1:noZeroP)*EXP(a*LOG(a/(a+p(1:noZeroP)))))&
          &/(1-SUM(pistat(1:noZeroP)*EXP(a*LOG(a/(a+p(1:noZeroP)))))))
     MLE=INT(MLE)

  else
     CALL ENPMLE(n,p,pi,noZeroP,t,test,Etol,Gtol,gap,gamma2,theta_threshold) 
     pistat(1:noZeroP)=pi(1:noZeroP)/(1-EXP(-p(1:noZeroP)))/SUM(pi(1:noZeroP)/(1-EXP(-p(1:noZeroP))))
     mle= SUM(n((t+1):50))+SUM(n(1:t))*(1+SUM(EXP(-p(1:noZeroP))*pistat(1:noZeroP)/&
          &(1-SUM(EXP(-p(1:noZeroP))*pistat(1:noZeroP)))))
     MLE=INT(MLE)
 endif   

END SUBROUTINE PPCGb

