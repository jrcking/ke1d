	SUBROUTINE AUSM(rhoL,momL,ENL,rhoR,momR,ENR,
     +					GAMR,PCRITR,
     +					fl1,fl2,fl3)

	IMPLICIT NONE
C	---------------------------------------------------------------
C	SUBROUTINE TO CALCULATE AUSM FLUXES (Liou & Steffen 1993)
C	---------------------------------------------------------------
C	Returns the fluxes directly
C	---------------------------------------------------------------
C	AUTHOR: J. R. C. KING
C	---------------------
C	CHANGE RECORD:
C		10<11-2014:	CREATED
C	---------------------------------------------------------------
	DOUBLE PRECISION	rhoL,momL,ENL,rhoR,momR,ENR,
     +				GAMR,PCRITR,rL,rR
C
	DOUBLE PRECISION	uL,PL,uR,PR,utL,utR
	DOUBLE PRECISION	eL,eR,HL,HR,cL,cR
	DOUBLE PRECISION 	ML,MR,MLP,MRM,Mhalf
	DOUBLE PRECISION 	pLp,pRm,phalf
C	FLUXES!!
	DOUBLE PRECISION	fl1,fl2,fl3,fl4
C	---------------------------------------------------------------
C	CALCULATE THE OTHER VARIABLES FROM THE INPUTS...
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
C	The left and right Mach numbers and middle Mach number
	ML = uL/cL
	MR = uR/cR
	IF (ML.LE.1.0)THEN
	MLP = 0.25*(ML+1.0)**2.0
	MRM = -0.25*(MR-1.0)**2.0
	ELSE
	MLP = 0.5*(ML+abs(ML))
	MRM = 0.5*(MR-abs(MR))
	END IF
	Mhalf = MLP+MRM
C	---------------------------------------------------------------
c	middle pressure...
	IF (ML.LE.1.0)THEN
!	pLp = 0.25*pL*(2.0-ML)*(ML+1.0)**2.0 !alternative of L&S 93
!	pRm = 0.25*pR*(2.0+MR)*(MR-1.0)**2.0
	pLp = 0.5*pL*(1+ML)
	pRm = 0.5*pR*(1-MR)
	ELSE
	pLp = 0.5*pL*(ML+abs(ML))/ML
	pRm = 0.5*pR*(MR-abs(MR))/MR
	END IF
	phalf = pLp + pRm
C	---------------------------------------------------------------
	IF (Mhalf.GT.0)THEN
		fl1 = Mhalf*rhoL*cL
		fl2 = Mhalf*rhoL*cL*uL + phalf
		fl3 = Mhalf*rhoL*cL*HL
	ELSE
		fl1 = Mhalf*rhoR*cR
		fl2 = Mhalf*rhoR*cR*uR + phalf
		fl3 = Mhalf*rhoR*cR*HR
	END IF
C	===============================================================

	RETURN
	END
C	---------------------------------------------------------------
