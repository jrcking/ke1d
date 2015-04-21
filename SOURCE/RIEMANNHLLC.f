	SUBROUTINE RIEMANNHLLC (rhoL,momL,ENL,rhoR,momR,ENR,
     +					GAMR,PCRITR,
     +					fl1,fl2,fl3)

	IMPLICIT NONE
C	---------------------------------------------------------------
C	SUBROUTINE TO SOLVE A RIEMANN PROBLEM - HLLC SOLVER 
C	---------------------------------------------------------------
C	DIRECTLY RETURNS FLUXES. TAKEN FROM 
c	http://dx.doi.org/10.1016/j.jcp.2009.06.002 (sort of)
C	---------------------------------------------------------------
C	AUTHOR: J. R. C. KING
C	---------------------
C	CHANGE RECORD:
C		19-03-2013:	CREATED
C		28-05-2013:	MODIFIED TO RETURN FLUXES
C		29-07-2013:	MODIFIED TO DEAL WITH ONLY 1 FLUID
C		07-08-2013:	MODIFIED TO TAKE CONS VARS FOR INPUT
C	---------------------------------------------------------------
	DOUBLE PRECISION	rhoL,momL,ENL,rhoR,momR,ENR,
     +				GAMR,PCRITR
C
	DOUBLE PRECISION	uL,PL,uR,PR
	DOUBLE PRECISION	pstar,ustar,rhoLstar,rhoRstar
	DOUBLE PRECISION	ENLstar,ENRstar,eL,eR,HL,HR
	DOUBLE PRECISION	cL,cR
	DOUBLE PRECISION	bL,bR,bM
	DOUBLE PRECISION	ubar,cbar,rhobar,Povrhobar,Pbar
	DOUBLE PRECISION	hbar
	DOUBLE PRECISION	alf1,alf2,alf3
C	FLUXES!!
	DOUBLE PRECISION	fl1,fl2,fl3
	LOGICAL			subsonicflag
C	---------------------------------------------------------------
C	CALCULATE PRIMITIVE VARIABLES FROM CONSERVATIVE VARIABLES
	uL = momL/rhoL
	uR = momR/rhoR
	PL = (GAMR-1.0)*(ENL-0.5*rhoL*uL**2.0) - GAMR*PCRITR
	PR = (GAMR-1.0)*(ENR-0.5*rhoR*uR**2.0) - GAMR*PCRITR
c	CALCULATE LEFT AND RIGHT SPEEDS OF SOUND ======================
	cL = sqrt(GAMR*(PL+PCRITR)/rhoL)
	cR = sqrt(GAMR*(PR+PCRITR)/rhoR)
C	LITTLE ENERGY AND ENTHALPY
	eL = (PL+GAMR*PCRITR)/((GAMR-1.0)*rhoL)
	eR = (PR+GAMR*PCRITR)/((GAMR-1.0)*rhoR)
	HL = eL + PL/rhoL + 0.5*uL**2.0
	HR = eR + PR/rhoR + 0.5*uR**2.0
C	CALCULATE ROE AVERAGES ((.)bar ) FOR DENSITY AND PRESSURE =====
	rhobar = sqrt(rhoL*rhoR)
	call ROEAVG (rhoL,rhoR,uL,uR,ubar)
	call ROEAVG (rhoL,rhoR,HL,HR,hbar)
	Povrhobar = (PL/sqrt(rhoL) + PR/sqrt(rhoR))/
     +			(sqrt(rhoL)+sqrt(rhoR)) + 
     +			0.5*((uR-uL)/(sqrt(rhoL)+sqrt(rhoR)))**2.0
	Pbar = rhobar*Povrhobar
C	JUST ONE FLUID ------------------------------------------------
C	cbar = sqrt(GAMR*(Pbar + PCRITR)/rhobar)   !seem to have a 
C	cbar = sqrt(GAMR*Povrhobar + PCRITR/rhobar) ! choice here..
	cbar = sqrt((GAMR-1.0)*(hbar-0.5*ubar**2))
C	CALCULATE STAR STATES AND WAVE SPEEDS =========================
	bL = min(uL-cL,ubar-cbar)
	bR = max(ubar+cbar,uR+cR)
	ustar = (rhoR*uR*(bR-uR) - rhoL*uL*(bL-uL) + PL - PR)/
     +		(rhoR*(bR-uR) - rhoL*(bL-uL))
	bM = ustar
	rhoLstar = rhoL*(bL-uL)/(bL-ustar)
	rhoRstar = rhoR*(bR-uR)/(bR-ustar)
	pstar = PL + rhoL*(uL-bL)*(uL-ustar)
	ENLstar = (pstar + GAMR*PCRITR)/(GAMR-1.0) + 
     +		0.5*rhoLstar*ustar**2.0
 	pstar = PR + rhoR*(bR-uR)*(ustar-uR)
	ENRstar = (pstar + GAMR*PCRITR)/(GAMR-1.0) + 
     +		0.5*rhoRstar*ustar**2.0
C	===============================================================
C	---------------------------------------------------------------
C	CALCULATE FLUXES DIRECTLY!!
	subsonicflag = .TRUE.
	IF (bL .GT. 0.0) THEN
		fl1 = rhoL*uL
		fl2 = fl1*uL + PL
		fl3 = uL*(PL + ENL)
		subsonicflag = .FALSE.
	END IF
	IF (bR .LT. 0.0) THEN
		fl1 = rhoR*uR
		fl2 = fl1*uR + PR
		fl3 = uR*(PR + ENR)
		subsonicflag = .FALSE.
	END IF
	IF (subsonicflag) THEN
		IF (bM .GT. 0.0) THEN
			fl1 = rhoL*uL + bL*(rhoLstar - rhoL)
			fl2 = rhoL*uL*uL +
     +				bL*(rhoLstar*ustar - rhoL*uL) + PL
			fl3 = uL*(PL+ENL) + bL*(ENLstar - ENL)
		ELSE
			fl1 = rhoR*uR + bR*(rhoRstar - rhoR)
			fl2 = rhoR*uR*uR + 
     +				bR*(rhoRstar*ustar - rhoR*uR) + PR
			fl3 = uR*(PR+ENR) + bR*(ENRstar - ENR)	
		END IF
	END IF
C	===============================================================
C	===============================================================
	RETURN
	END
C	---------------------------------------------------------------
C	###############################################################
C	<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
C	###############################################################
	SUBROUTINE ROEAVG (roL,roR,ftnL,ftnR,RAVG)

	IMPLICIT NONE
C	---------------------------------------------------------------
	DOUBLE PRECISION roL,roR,ftnL,ftnR,RAVG
C	---------------------------------------------------------------
	RAVG = (sqrt(roL)*ftnL + sqrt(roR)*ftnR)/
     +			(sqrt(roL)+sqrt(roR))
C	---------------------------------------------------------------
	RETURN
	END 
C	---------------------------------------------------------------
