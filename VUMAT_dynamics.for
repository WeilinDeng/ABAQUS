cccccccccccccccccccccccccccccccccccccccccccccc
c
c    Fan Zhang
c    June.03.09
c    Based on Hopperstad, Borvik, Berstad and Benallal, J Phys IV France 134(2006), 435
c
c    Props
c	E=props(1)
c	nu=props(2)
c	Q0=props(3)
c	SRS0=props(4)
c	H=props(5)
c	alpha=props(6)
c	td=props(7)
c	HQ1=props(8)
c	HQ2=props(9)
c	HC1=props(10)
c	HC2=props(11)
c	pR0=props(12)
c	Omega=props(13)
c	state variables
c	p0=stateOld(kb,1)
c	ta0=stateOld(kb,2)
c
cccccccccccccccccccccccccccccccccccccccccccccc
      subroutine vumat(
C Read only (unmodifiable)variables -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only (modifiable) variables -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
C
      include 'vaba_param.inc'
c	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), tempOld(nblock),
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     8  defgradNew(nblock,ndir+nshr+nshr),
     9  fieldNew(nblock,nfieldv),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     2  enerInternNew(nblock), enerInelasNew(nblock)
C
      character*80 cmname
C
c  Local variables	
	INTEGER nds,kb
c
c
c
c	WRITE(*,*) 'VUMAT',TotalTime !debug
	nds=ndir+nshr
	IF (TotalTime>0) THEN 
	  DO kb=1,nblock
c	    WRITE(*,*) 'Normal,kb=',kb !debug
		CALL TimeInc(nblock,nstatev,ndir,nshr,nds,kb,nprops,
	1     props,stateOld,stateNew,stressOld,stressNew,strainInc,
     2     dt,StepTime,TotalTime,EIOld,EInEOld,EINew,EInENew,density)
	  END DO
	ELSE ! initialize
	  DO kb=1,nblock
c	    WRITE(*,*) 'Initial,kb=',kb !debug
		CALL Initial(nblock,nstatev,ndir,nshr,nds,kb,nprops,
	1     props,stateOld,stateNew,stressOld,stressNew,strainInc,
     2     dt,StepTime,TotalTime,EIOld,EInEOld,EINew,EInENew,density)
	  END DO
	END IF
      RETURN
      END
c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  MAIN SUBROUTINES
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!***********************************************************************
!
!  Time increment updates
!
      SUBROUTINE TimeInc(nblock,nstatev,ndir,nshr,nds,kb,nprops,
	1     props,stateOld,stateNew,stressOld,stressNew,strainInc,
     2     dt,StepTime,TotalTime,EIOld,EInEOld,EINew,EInENew,dens)
c	IMPLICIT NONE
      include 'vaba_param.inc'
	INTEGER,PARAMETER :: itermax=100
	DOUBLE PRECISION,PARAMETER :: err=0.0001
	INTEGER nblock,nstatev,ndir,nshr,nds,kb,nprops
c	DOUBLE PRECISION stressOld,stressNew,stateOld,stateNew,strainInc,
c	1                 props
c	DOUBLE PRECISION EIOld,EInEOld,EINew,EInENew,dens,dt,
c	1                 StepTime,TotalTime
	DIMENSION stressOld(nblock,nds),stressNew(nblock,nds),
	1          stateOld(nblock,nstatev),stateNew(nblock,nstatev),
     2          strainInc(nblock,nds),props(nprops) 
c  Local variables
	INTEGER i
	DOUBLE PRECISION Str(nds),Snew(nds)
	DOUBLE PRECISION E,nu,Q0,SRS0,H,alpha,td,HQ1,HQ2,HC1,HC2,pR0,Omega
	DOUBLE PRECISION p0,ta0,p,ta,dp,pR,dpR
	DOUBLE PRECISION QY,R,Qv,DQY,DR,DQv,RHS,DRHS,dta,tt
	DOUBLE PRECISION G,twoG,lamda,means,Qtr,Qnew,ratio
c
c	Parameters	
	E=props(1)
	nu=props(2)
	Q0=props(3)
	SRS0=props(4)
	H=props(5)
	alpha=props(6)
	td=props(7)
	HQ1=props(8)
	HQ2=props(9)
	HC1=props(10)
	HC2=props(11)
	pR0=props(12)
	Omega=props(13)
c	state variables
	p0=stateOld(kb,1)
	ta0=stateOld(kb,2)
	p=p0
	ta=ta0
c
	G=E/(2.D0*(1.D0+nu))
	twoG=2.D0*G
	lamda=twoG*nu/(1.D0-2.D0*nu)
	means=strainInc(kb,1)+strainInc(kb,2)+strainInc(kb,3)
	Str(1:nds)=stressOld(kb,1:nds)+twoG*strainInc(kb,1:ndir)
	Str(1:ndir)=Str(1:ndir)+lamda*means
	means=(Str(1)+Str(2)+Str(3))/3.D0
	Str(1:ndir)=Str(1:ndir)-means
	Qtr=DSQRT(1.5D0*(SUM(Str(1:ndir)*Str(1:ndir))
	1         +3.D0*SUM(Str(ndir+1:nds)*Str(ndir+1:nds))))
c
	QY=Q0+SRS0*H*(1.D0-DEXP(-(ta/td)**alpha))
	R=HQ1*(1.D0-DEXP(-HC1*p))+HQ2*(1.D0-DEXP(-HC2*p))
	pR=0.D0
	Qv=0.D0
c	WRITE(*,*) 'QY,R,Qtr' !debug
c	WRITE(*,*) QY,R,Qtr   !debug
	IF (Qtr>QY+R) THEN !Plastic
!	solve pR
	  DO i=1,itermax
		p=p0+pR*dt
		ta=(ta0+dt)/(1.D0+dt*pR/Omega)
c		dta=-(ta0+dt)*dt/((1.D0+dt*pR/Omega)**2*Omega)
		tt=(ta/td)**alpha
		QY=Q0+SRS0*H*(1.D0-DEXP(-tt))
		R=HQ1*(1.D0-DEXP(-HC1*p))+HQ2*(1.D0-DEXP(-HC2*p))
		Qv=SRS0*DLOG(1.D0+pR/pR0)
		DQY=-SRS0*H*DEXP(-tt)*tt*alpha*dt/(Omega+pR*dt)
		DR1=HQ1*HC1*DEXP(-HC1*p)
		DR2=HQ2*HC2*DEXP(-HC2*p)
		DQv=SRS0/(pR0+pR)
	    RHS=Qtr-3.D0*G*pR*dt-QY-R-Qv
		DRHS=dt*(-DQY+DR1+DR2+3.D0*G)+DQv
		dpR=RHS/DRHS
		pR=pR+dpR
		IF (DABS(dpR)<DABS(err*pR)) EXIT
	    IF (i>itermax/2) THEN
	      WRITE(*,*) 'RHS,DRHS,pR'
	      WRITE(*,*) RHS,DRHS,pR
	      WRITE(*,*) 'DQY,DR1,DR2,DQv,dt'
		  WRITE(*,*) DQY,DR1,DR2,DQv,dt
	    END IF
	  END DO
	  IF (i>itermax) THEN
	    WRITE(*,*) '***ERRORZF: NOT CONV'
	    WRITE(*,*) 'Qtr,QY,R,Qv,p,pR'
	    WRITE(*,*) Qtr,QY,R,Qv,p,pR
		WRITE(*,*) 'DQY,DR,DQv'
	    WRITE(*,*) DQY,DR,DQv
c		CALL XIT
	  END IF
	END IF
c
	dp=pR*dt
	p=p0+dp
	Qnew=Qtr-3.D0*G*dp
	ta=(ta0+dt)/(1.D0+dt*pR/Omega)
	IF (Qnew>1.D-20) THEN
	  ratio=1.D0/(1.D0+3.D0*G*dp/Qnew)
	ELSE
	  ratio=1.D0
	END IF
	Snew=Str*ratio
c	Update stressNew and stateNew
	stressNew(kb,1:ndir)=Snew(1:ndir)+means
	stressNew(kb,ndir+1:nds)=Snew(ndir+1:nds)
	stateNew(kb,1)=p
	stateNew(kb,2)=ta
	stateNew(kb,3)=R ! output debug
	stateNew(kb,4)=Qv ! output debug
	stateNew(kb,5)=QY-QR ! output debug
	EINew=EIOld
	1      +(0.5*SUM((stressOld(kb,1:ndir)+stressNew(kb,1:ndir))
     2        *strainInc(kb,1:ndir))+
     3        SUM((stressOld(kb,ndir+1:nds)+stressNew(kb,ndir+1:nds))
     4        *strainInc(kb,ndir+1:nds)))/dens
	EInENew=EInEOld+Qnew*dp/dens
c	WRITE(*,*) '***Final QY,R,Qv,Qnew,pR' !debug
c	WRITE(*,'(5F12.4)') QY,R,Qv,Qnew,pR   !debug
c	WRITE(*,*) 'i,RHS,pR,dpR'   !debug
c	WRITE(*,'(I4,3F12.4)') i,RHS,pR,dpR   !debug
	RETURN
	END
!***********************************************************************
!
!  Time increment updates
!
      SUBROUTINE Initial(nblock,nstatev,ndir,nshr,nds,kb,nprops,
	1     props,stateOld,stateNew,stressOld,stressNew,strainInc,
     2     dt,StepTime,TotalTime,EIOld,EInEOld,EINew,EInENew,dens)
c	IMPLICIT NONE
      include 'vaba_param.inc'
	INTEGER nblock,nstatev,ndir,nshr,nds,kb,nprops
c	DOUBLE PRECISION stressOld,stressNew,stateOld,stateNew,strainInc,
c	1                 props
c	DOUBLE PRECISION EIOld,EInEOld,EINew,EInENew,dens,dt,
c	1                 StepTime,TotalTime
	DIMENSION stressOld(nblock,nds),stressNew(nblock,nds),
	1          stateOld(nblock,nstatev),stateNew(nblock,nstatev),
     2          strainInc(nblock,nds),props(nprops) 
c  Local variables
	INTEGER i
	DOUBLE PRECISION Str(nds),Snew(nds)
	DOUBLE PRECISION E,nu,Q0,SRS0,H,alpha,td,HQ1,HQ2,HC1,HC2,pR0,Omega
	DOUBLE PRECISION p0,ta0,p,ta,dp,pR,dpR
	DOUBLE PRECISION QY,R,Qv,DQY,DR,DQv,RHS,DRHS,dta,tt
	DOUBLE PRECISION G,twoG,lamda,means,Qtr,Qnew,ratio
c
c	Parameters	
	E=props(1)
	nu=props(2)
c	state variables
	G=E/(2.D0*(1.D0+nu))
	twoG=2.D0*G
	lamda=twoG*nu/(1.D0-2.D0*nu)
	means=strainInc(kb,1)+strainInc(kb,2)+strainInc(kb,3)
	Str(1:nds)=stressOld(kb,1:nds)+twoG*strainInc(kb,1:ndir)
	Str(1:ndir)=Str(1:ndir)+lamda*means
	means=(Str(1)+Str(2)+Str(3))/3.D0
	Str(1:ndir)=Str(1:ndir)-means
	Qtr=DSQRT(1.5D0*(SUM(Str(1:ndir)*Str(1:ndir))
	1         +3.D0*SUM(Str(ndir+1:nds)*Str(ndir+1:nds))))
	p0=0.D0
	ta0=0.D0
c
	pR=0.D0
	dp=pR*dt
	p=p0+dp
	ta=0.D0
	Qnew=Qtr-3.D0*G*dp
	IF (Qnew>1.D-20) THEN
	  ratio=1.D0/(1.D0+3.D0*G*dp/Qnew)
	ELSE
	  ratio=1.D0
	END IF
	Snew=Str*ratio
c	Update stressNew and stateNew
	stressNew(kb,1:ndir)=Snew(1:ndir)+means
	stressNew(kb,ndir+1:nds)=Snew(ndir+1:nds)
	stateNew(kb,1)=p
	stateNew(kb,2)=ta
	EINew=EIOld
	1      +(0.5*SUM((stressOld(kb,1:ndir)+stressNew(kb,1:ndir))
     2        *strainInc(kb,1:ndir))+
     3        SUM((stressOld(kb,ndir+1:nds)+stressNew(kb,ndir+1:nds))
     4        *strainInc(kb,ndir+1:nds)))/dens
	EInENew=EInEOld+Qnew*dp/dens
	RETURN
	END
