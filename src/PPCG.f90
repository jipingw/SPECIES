
SUBROUTINE PPCG(n,t, MLE,alpha,alphaK,amodel,p,pi,noZeroP)


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
  !             alpha  :    alpha is a grid vector for tuning parameter alpha
  !            alphaK  :    length of alphaK
  !            amodel  :    amodel selected from cross-validation
  !
  !       diagram:
  !
  !       ppcg: main subroutine wrapper
  !       {
  !               call PCGone: main body of codes, cross-validation for model selection
  !               {     
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
  !                }
  !        }
  !       =======================================================================================================


  IMPLICIT NONE
  INTEGER t,alphaK,noZeroP
  DOUBLE PRECISION n(50),MLE,amodel
  DOUBLE PRECISION alpha(alphaK),p(10),pi(10),a_model,theta_threshold
  noZeroP=0;MLE=0.D0;a_model=0.D0;p=p*0.D0;pi=0.D0*pi;noZeroP=0.D0
  theta_threshold=1./SUM(n(1:t))   !!threshold value to control the linearization  iterations
  !print*,alpha, alphak
  CALL PCGone(n,t,MLE,a_model,alpha,alphaK,p,pi,theta_threshold,noZeroP)
  MLE=INT(MLE)
  !print*,a_model
  amodel=a_model

END SUBROUTINE PPCG

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! PCGone: main body of codes realizing the PCG model, cross-validation etc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PCGone(n,t,MLE,a_model,alpha,alphaK,p,pi,theta_threshold,noZeroP)

  !       ========================================================================================================
  !       purpose:   PCG model with exponential paritial prior using cross-validation for model selection
  !       Input :    same as PPCG
  !                  theta_threshold -- the threhold value for linearization procedure
  !       Called by: PPCG
  !
  !
  !       diagram:
  !
  !       PCGone: main body of codes, cross-validation for model selection
  !       {
  !                  call PCGpenE: NPMLE for PCG model with exponential prior,linearlization
  !                  {
  !                             call PCGpen: NPMLE for PCG model with linear penalty 
  !                             {
  !                               call  wBisectionCon_theta : optimal weight for VDM part
  !                               call EMScon_theta:      ECM part for VDM/ECM algorithm
  !                             }
  !                   }
  !                   call ENPMLE: NPMLE with exponential partial prior, linearization
  !                   {
  !                             call NPMLEpen: NPMLE for fixed linear penalty
  !                   }
  !        }
  !       =======================================================================================================


  IMPLICIT NONE
  INTEGER t,alphaK,cumuless,i,k,noZeroP
  DOUBLE PRECISION p(10),pi(10),pistat(10),a_model, predict(alphaK+1)
  DOUBLE PRECISION boot(50),n(50),boot1(50),a
  DOUBLE PRECISION Etol, Gtol,gap,mle,theta_threshold
  DOUBLE PRECISION test,Pmix,PmixScon,alpha(alphaK),max_predict, gamma2

  Etol=1.D-10;Gtol=0.005;gap=0.02 
  boot=n

  !print*,alpha, alphak
  !print*, n

  IF(t<50) boot((t+1):50)=0.D0
  IF(t>=50) THEN 
     t=49
     boot(50)=0.D0
  END IF
  a_model=alpha(1)

  predict=predict*0.D0
  k=1
  test=1
  DO WHILE (k<=alphaK)
     a=alpha(k)
     noZeroP=0
     CALL PCGpenE (boot,p,pi,noZeroP,a,t,test,Etol,Gtol,gap,gamma2,theta_threshold)  
     !print*,p,pi
     i=1
     !print*, boot

     DO WHILE (i<=t)
        predict(k)=predict(k)+PmixScon(i,p,pi,a,noZeroP)**2
        i=i+1
        !print*, 'predict(j,k)=   ',  predict(k)
     END DO
    
     i=1
     DO WHILE (i<=t)
        boot1=boot
        IF(boot(i)>0) THEN
           boot1(i)=boot(i)-1
           CALL PCGpenE (boot1,p,pi,noZeroP,a,t,test,Etol,Gtol,gap,gamma2,theta_threshold)  
           predict(k)=predict(k)-2*boot(i)*PmixScon(i,p,pi,a,noZeroP)/SUM(boot)
        END IF
        i=i+1
     END DO
     !print*,'         k=', k,'a=',a, 'predict=   ', predict(k)
 
    
     IF(k==1) THEN 
        max_predict=predict(k)
     END IF

     IF (k>1 .AND. predict(k)<max_predict) THEN
        max_predict=predict(k)
        a_model=alpha(k)
        a=a_model
        k=k+1
     ELSE
        k=k+1
     END IF
  END DO

  IF(k==alphaK+1) THEN
     CALL ENPMLE (boot,p,pi,noZeroP,t,test,Etol,Gtol,gap,gamma2,theta_threshold)
     i=1
     DO WHILE (i<=t)
        predict(k)=predict(k)+Pmix(i,p,pi,noZeroP)**2
        i=i+1
     END DO

     !!cross validation
     i=1
     DO WHILE (i<=t)
        !print*,'                        i=', i
        boot1=boot
        IF(boot(i)>0) THEN
           boot1(i)=boot(i)-1
           !print*, boot1
           CALL ENPMLE(boot1,p,pi,noZeroP,t,test,Etol,Gtol,gap,gamma2,theta_threshold) 
           predict(k)=predict(k)-2*boot(i)*Pmix(i,p,pi,noZeroP)/SUM(boot)
        END IF
        i=i+1
     END DO


     IF (predict(k)<max_predict) THEN
        max_predict=predict(k)
        a_model=500.
        k=k+1
        CALL ENPMLE(boot,p,pi,noZeroP,t,test,Etol,Gtol,gap,gamma2,theta_threshold) 
        pistat(1:noZeroP)=pi(1:noZeroP)/(1-EXP(-p(1:noZeroP)))/SUM(pi(1:noZeroP)/(1-EXP(-p(1:noZeroP))))
        mle= SUM(boot((t+1):50))+SUM(boot(1:t))*(1+SUM(EXP(-p(1:noZeroP))*pistat(1:noZeroP)/&
             &(1-SUM(EXP(-p(1:noZeroP))*pistat(1:noZeroP)))))
     ELSE
        a=a_model
        CALL PCGpenE (boot,p,pi,noZeroP,a,t,test,Etol,Gtol,gap,gamma2,theta_threshold)  
        pistat(1:noZeroP)=1/(1-(a/(a+p(1:noZeroP)))**a)*pi(1:noZeroP)/&
             &SUM(1/(1-(a/(a+p(1:noZeroP)))**a)*pi(1:noZeroP))
        mle= SUM(boot((t+1):50))+SUM(boot(1:t))*(1+SUM(pistat(1:noZeroP)*EXP(a*LOG(a/(a+p(1:noZeroP)))))&
             &/(1-SUM(pistat(1:noZeroP)*EXP(a*LOG(a/(a+p(1:noZeroP)))))))
     END IF
  END IF
  !print*,'a=',a_model,'mle=',mle
END SUBROUTINE PCGone

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! PCGpenE: subroutine to fit PCG model for a given alpha value
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE PCGpenE(n,p,pi,noZeroP,a,t,test,Etol,Gtol,gap,gamma2,theta_threshold)

  !       ========================================================================================================
  !       purpose:   PCG model with exponential paritial prior using cross-validation for model selection
  !       Input :    same as PCGone
  !                Etol:    threshold value for convergence in terms of log likelihood in the EM part
  !                Gtol:    threshold value for convergence of entire algorithm in terms of gradient
  !                gap:     threold value to control which components to be merged if two components are close.
  !
  !               gamma2:   the initial value for penalty factor for the linearized likelihood.
  !                         This value will be optimized as algorithm proceeds. 
  !                         As this value is about the same in different runs of leave-one out cross-validation,
  !                         Carrying this value over to the next run often shortens the run time.
  !                  a:     a given alpha value for the compound Gamma model
  !
  !       Called by: PCGone
  !
  !                  PCGpenE: NPMLE for PCG model with exponential prior,linearlization
  !                  {
  !                             call PCGpen: NPMLE for PCG model with linear penalty 
  !                             {
  !                               call  wBisectionCon_theta : optimal weight for VDM part
  !                               call EMScon_theta:      ECM part for VDM/ECM algorithm
  !                             }
  !                   }
  !       ========================================================================================================


  IMPLICIT NONE
  INTEGER t,ite,noZeroP,STOP
  DOUBLE PRECISION p(10),pi(10),n(50)
  DOUBLE PRECISION Etol, Gtol,gap,boot(50),gamma2,theta_threshold
  DOUBLE PRECISION test,a,chao,thetaOld,theta_0

  !!lower bound estimator chao

  boot=n
  chao=boot(1)**2/2/boot(2)/SUM(boot)
  theta_0=chao
  STOP=0


  !! gamma2 is the penalty factor in the likelihood with linearized penalty.
  !! if a gamma2 is too small, e.g. <0.02, we use the gamma2 from the lower bound estimator as initial value
  !! otherwise, use the gamma2 carried over from last run of the subroutine as initial value.

  IF(gamma2<0.2) THEN
     gamma2=1/theta_0
  ELSE
     theta_0=1/gamma2
  END IF

  ite=1
  noZeroP=0
  !print*,'theta_0',theta_0

  DO WHILE (STOP==0)
     thetaOld=theta_0
     theta_0=0
     CALL PCGpen(boot,p,pi,noZeroP,a,t,test,Etol,Gtol,gap,gamma2)
     theta_0=SUM(pi(1:noZeroP)/((1+p(1:noZeroP)/a)**a-1))
     theta_0=.5*theta_0+.5*thetaOld
     gamma2=1/theta_0
     IF(gamma2<0) THEN
        theta_0=(chao+thetaOld)/2
        gamma2= 1/theta_0
     END IF
     ite=ite+1

     IF (ABS(thetaOld-theta_0)<theta_threshold .OR. theta_0>100 .OR. ite>10) STOP=1
  END DO
END SUBROUTINE PCGpenE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! PCGpen: linearized penalty for fit PCG model for a given alpha value
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PCGpen (n,p,pi,noZeroP,a,t,test,Etol,Gtol,gap,gamma2)

  !       ========================================================================================================
  !       purpose:   PCG model with exponential paritial prior using cross-validation for model selection
  !       Input :    same as  PCGpenE
  !                  a ----  a given alpha value 
  !
  !       Called by: PCGpenE
  !
  !                  call PCGpen: NPMLE for PCG model with linear penalty 
  !                  {
  !                         call  wBisectionCon_theta : optimal weight for VDM part
  !                         call EMScon_theta:      ECM part for VDM/ECM algorithm
  !                  }
  !       ========================================================================================================

  IMPLICIT NONE
  INTEGER t,i,j,noZeroP,loc,ite
  DOUBLE PRECISION p(10),pi(10),pvec(1000+(t-1)*100),drv(1000+(t-1)*100),Lq(50)
  DOUBLE PRECISION Etol, Gtol,sum_a,sum_b,p_add,gap
  DOUBLE PRECISION n(50),a,w,test,like,PmixScon,PdenScon,gamma2

  like=0.D0
  test=1.D0
  DO i=1,1000   
     pvec(i)=1D-3*i
  END DO

  DO i=1,(t-1)*100
     pvec(i+1000)=1.+1D-2*i
  END DO

  IF(noZeroP==0) THEN
     p=p*0.D0
     pi=pi*0.D0
     sum_a=0.0D0; sum_b=0.0D0
     DO i=1,t
        sum_a=sum_a + n(i)*(i); sum_b = sum_b + n(i)
     END DO
     p_add=sum_a/sum_b; p(1)=p_add; pi(1)=1
     p(1)=p_add; pi(1)=1.0D0
     noZeroP=1
  END IF

  CALL EMScon_theta(n,p,pi,noZeroP,Etol,a,t,gamma2)
  CALL checkGap(p,pi,gap,noZeroP) 
  !print*, p(1:noZeroP), pi(1:noZeroP),a,noZeroP

  DO i=1,t
     Lq(i)=PmixScon(i,p,pi,a,noZeroP)
  END DO

  !print*,Lq(1:t)

  DO i=1,SIZE(pvec)
     DO j=1,t
        drv(i)=drv(i)+n(j)*(PdenScon(j,pvec(i),a)/Lq(j)-1)
     END DO
     drv(i)=drv(i)-gamma2*(1/((1+pvec(i)/a)**a-1)-SUM(pi(1:noZeroP)/((1+p(1:noZeroP)/a)**a-1)))
  END DO

  loc=MAXLOC(drv,1)
  p_add=pvec(loc)
  test=MAXVAL(drv)

  w=1.D0
  ite=1
  test=1
  DO WHILE(test>Gtol .AND. ite<10)
     CALL  wBisectionCon_theta (n,p,pi,noZeroP,p_add,Lq,w,a,t,gamma2)
     p(noZeroP+1)=p_add
     pi(1:noZeroP)=(1-w)*pi(1:noZeroP)
     pi(noZeroP+1)=w
     noZeroP=noZeroP+1

     !print*,'p=',p
     !print*,'pi=',pi
     !print*,'noZeroP=',noZeroP,'a=',a,'t=',t,'Etol=',Etol
         
     CALL EMScon_theta(n,p,pi,noZeroP,Etol,a,t,gamma2)
     CALL checkGap(p,pi,gap,noZeroP) 

     DO i=1,t
        Lq(i)=PmixScon(i,p,pi,a,noZeroP)
     END DO

     drv=drv*0D0
     DO i=1,SIZE(pvec)
        DO j=1,t
           drv(i)=drv(i)+n(j)*(PdenScon(j,pvec(i),a)/Lq(j)-1)
        END DO
        drv(i)=drv(i)-gamma2*(1/((1+pvec(i)/a)**a-1)-SUM(pi(1:noZeroP)/((1+p(1:noZeroP)/a)**a-1)))
     END DO
     loc=MAXLOC(drv,1)
     p_add=pvec(loc)
     test=MAXVAL(drv)
     ite=ite+1
  END DO
  CALL sortP(p,pi,noZeroP)
END SUBROUTINE PCGpen



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ENPMLE: subroutine for PCG model with alpha=infity, in which a descrete model
!!         with exponential partial prior is fit.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ENPMLE(n,p,pi,noZeroP,t,test,Etol,Gtol,gap,gamma2,theta_threshold)

  !       ========================================================================================================
  !       purpose:   discrete NPMLE with exponential partial prior, corresponding to PCG model alpha=infinity
  !       Input :    same as PCGone
  !                  gamma2 --- the initial value for penalty factor for the linearized likelihood.
  !                             This value will be optimized as algorithm proceeds. 
  !                             As this value is about the same in different runs of leave-one out cross-validation,
  !                             Carrying this value over to the next run often shortens the run time.
  !
  !
  !       Called by: PCGone
  !
  !                   ENPMLE: discrete NPMLE with exponential partial prior, linearization
  !                   {
  !                             call NPMLEpen: NPMLE for fixed linear penalty
  !                   }
  !      ========================================================================================================



  IMPLICIT NONE

  INTEGER t,ite,noZeroP,STOP
  DOUBLE PRECISION p(10),pi(10),n(50)
  DOUBLE PRECISION Etol, Gtol,gap,boot(50),gamma2,theta_threshold
  DOUBLE PRECISION test,chao,thetaOld,theta_0

  !!lower estimator chao
  boot=n
  chao=boot(1)**2/2/boot(2)/SUM(boot)
  theta_0=chao
  STOP=0

  !! gamma2 is the penalty factor in the likelihood with linearized penalty.
  !! if a gamma2 is too small, e.g. <0.02, we use the gamma2 from the lower bound estimator as initial value
  !! otherwise, use the gamma2 carried over from last run of the subroutine as initial value.

  IF(gamma2<0.02) THEN
     gamma2=1/theta_0
  ELSE
     theta_0=1/gamma2
  END IF

  ite=1
  noZeroP=0

  DO WHILE (STOP==0)
     thetaOld=theta_0
     theta_0=0
     CALL NPMLEpen(boot,p,pi,noZeroP,t,test,Etol,Gtol,gap,gamma2)
     theta_0=SUM(1/(EXP(p(1:noZeroP))-1)*pi(1:nozeroP))    
     theta_0=0.5*theta_0+0.5*thetaOld	
     gamma2=1/(theta_0)

     IF(gamma2<0) THEN
        theta_0=(chao+thetaOld)/2
        gamma2=1/(theta_0)
     END IF
     ite=ite+1
     IF (ABS(thetaOld-theta_0)<theta_threshold .OR. theta_0>100 .OR. ite>10) STOP=1
  END DO
END SUBROUTINE ENPMLE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EMScon_theta(n,p,pi,noZeroP,Etol,a,t,gamma2)

  !       ========================================================================================================
  !       purpose:   EM algorithm updating p and pi given a linear penalty with fixed penalty factor gamma2
  !       Input :    same as  PCGpen
  !
  !
  !       Called by: PCGpen
  !
  !       ========================================================================================================

  IMPLICIT NONE
  INTEGER noZeroP,ite,ite2,i,j,t,g
  DOUBLE PRECISION p(10),pi(10),a,delta,delta0,deltaold,count_data(noZeroP)
  DOUBLE PRECISION n(50),sum_a,sum_b,Etol,z(t,noZeroP),d1,d2
  DOUBLE PRECISION currentL, last,lambda0,lambda1,PmixScon,PdenScon,gamma2

  currentL=0.0D0
  ite=1
  g=noZeroP

  currentL=-gamma2*SUM(pi(1:noZeroP)/((1+p(1:noZeroP)/a)**a-1))

  DO i=1,t
     currentL=currentL+n(i)*LOG(PmixScon(i,p,pi,a,noZeroP))
  END DO

  !print*, currentL
  last=2*currentL

  DO WHILE(ABS(currentL-last)>Etol .AND. ite<2000)
     last=currentL
     DO i=1,t
        DO j=1,noZeroP
           z(i,j)=pi(j)*PdenScon(i,p(j),a)/PmixScon(i,p,pi,a,noZeroP)
        END DO
     END DO

     delta=-SUM(n(1:t))
     delta0=delta+1
     deltaold=delta0

     DO i=1,g
        count_data(i)=SUM(z(:,i)*n(1:t))
     END DO

     DO WHILE (ABS(delta-delta0)>0.0001)
        delta0=delta
        d1=-SUM(count_data/(delta-gamma2/((1+p(1:noZeroP)/a)**a-1)))-1
        d2=SUM(count_data/(delta-gamma2/((1+p(1:noZeroP)/a)**a-1))**2)
        delta=delta-d1/d2
     END DO

     pi=-count_data/(delta-gamma2/((1+p(1:noZeroP)/a)**a-1))

     currentL=-gamma2*SUM(pi(1:noZeroP)/((1+p(1:noZeroP)/a)**a-1))
     DO i=1,t
        currentL=currentL+n(i)*LOG(PmixScon(i,p,pi,a,noZeroP))
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
                 sum_a=sum_a+n(j)*z(j,i)*(1D0*j/lambda1-(j+a)/(lambda1+a)-1.D0/((1.D0+lambda1/a)&
                      &**(a+1.D0)-(1D0+lambda1/a)))
                 sum_b=sum_b+n(j)*z(j,i)*(-1D0*j/lambda1**2+(j+a)/(lambda1+a)**2+1D0/((1+lambda1/a)&
                      &**(a+1)-(a+lambda1)/a)*(a+1-(a/(a+lambda1))**a)/((a+lambda1)-(a+lambda1)*(a/(a+lambda1))**a))
              END DO
              sum_a=sum_a+gamma2*pi(i)*(1+lambda1/a)**(a-1)*((1+lambda1/a)**a-1)**(-2)
              sum_b=sum_b-2*gamma2*((1+lambda1/a)**a-1)**(-3)*(1+lambda1/a)**(2*a-2)+&
                   &gamma2*(a-1)/a*((1+lambda1/a)**a-1)**(-2)*(1+lambda1/a)**(a-2)
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

     currentL=-gamma2*SUM(pi(1:noZeroP)/((1+p(1:noZeroP)/a)**a-1))
     DO i=1,t
        currentL=currentL+n(i)*LOG(PmixScon(i,p,pi,a,noZeroP))
     END DO
     ite=ite+1
  END DO
  !print*,currentL

END SUBROUTINE EMScon_theta


SUBROUTINE wBisectionCon_theta(n,p,pi,noZeroP,p_add,Lq,w,a,t,gamma2)

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
  DOUBLE PRECISION p(10),pi(10),temp_p(10),temp_pi(10),n(*),p_add,Lq(t),w
  DOUBLE PRECISION D,L(t),lower, upper,c,PdenScon,PmixScon,a,gamma2

  lower=0.D0
  upper=1.0D0
  w=0.5D0
  D=1.0D0
  ite=1
  !!print*,'here we are'

  DO WHILE (ABS(D)>0.0000000001 .AND. ite<2000)
     temp_p=p;temp_pi=pi
     temp_p(noZeroP+1)=p_add
     temp_pi(1:noZeroP)=(1-w)*pi(1:noZeroP)
     temp_pi(noZeroP+1)=w
     ite=ite+1
     D=0.0D0
     DO i=1,t
        L(i)=PmixScon(i,temp_p,temp_pi,a,noZeroP+1)
     END DO

     !print*, L
     DO i=1,t
        c=n(i)*(-Lq(i)+PdenScon(i,p_add,a))/L(i)
        D=D+c
     END DO
     D=D+gamma2*(SUM(pi(1:noZeroP)/((1+p(1:noZeroP)/a)**a-1))-1/((1+p_add/a)**a-1))

     IF (D>0)  lower=w
     IF (D<0)  upper=w
     w=(lower+upper)/2
     IF (ABS(w)>500 .OR. ite>=60) THEN
        w=0.01D0
     END IF
  END DO
END SUBROUTINE wBisectionCon_theta

