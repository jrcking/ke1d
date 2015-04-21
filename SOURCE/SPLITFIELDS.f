	SUBROUTINE SPLITFIELDS

	INCLUDE "commonblock"
C	---------------------------------------------------------------
C	SUBROUTINE TO SPLIT THE FIELDS FOR THE GHOST FLUID METHOD
C	---------------------------------------------------------------
C	---------------------------------------------------------------
C	AUTHOR: J. R. C. KING
C	---------------------
C	CHANGE RECORD:
C		20-07-2013:	CREATED
C		06-08-2013:	VARIOUS CHANGES TO INTERFACE INP STATES
C		07-08-2013:	MODIFIED TO WORK ON CONSERVATIVE VARS
C		14-08-2013:	ONLY USE NECESSARY GHOST CELLS
C		06-11-2013:	OPTION FOR W-INTERP GFM INPUTS
C	---------------------------------------------------------------
C	---------------------------------------------------------------
C	---------------------------------------------------------------
	DOUBLE PRECISION	rhointL,rhointR,uintL,uintR,pintL,pintR
	DOUBLE PRECISION	spr1
C
	DOUBLE PRECISION	drobydrL,dpbydrL,dubydrL
	DOUBLE PRECISION	drobydrR,dpbydrR,dubydrR
C
	DOUBLE PRECISION	RTL,RTR
C	
C	---------------------------------------------------------------
C	===============================================================
C	CHOOSE THE INPUT STATES TO THE INTERFACE RIEMANN PROBLEM AND 
C	THE LIMITS FOR THE STATES WHICH ARE MODIFIED.
C		intL and intR are the cells either side of the 
C			interface;
C		intLin and intRin are the cells for the input states 
C			of the Riemann problem;
C		intLabs and intRabs are the indices of the last real
C			cells not modified with the GFM.
C
	rLOW(1) = 1
	rHIGH(1) = q+1+gfmbuffer
	rLOW(2) = q-gfmbuffer
	rHIGH(2) = nr+inttype-1
C
C	IN THE MAIN PROGRAM, PRESSURE IS ADJUSTED TO BE NON-NEGATIVE 
C	EVERY TIME-STEP. ENERGY IS NOT ADJUSTED. IF SPLITFIELDS
C	TAKES P FOR INPUTS RESULTS WILL DIFFER FROM IF SPLITFIELDS
C	TAKES E FOR INPUTS. 
C
	rhointL = rho(q)
	uintL = u(q)
	rhointR = rho(q+1)
	uintR = u(q+1)
	pintL = (gamma0-1.0)*(E(q)-0.5*rhointL*uintL*uintL)-
     +			gamma0*pc0
	pintR = (gamma1-1.0)*(E(q+1)-0.5*rhointR*uintR*uintR)-
     +			gamma1*pc1
C	---------------------------------------------------------------
	IF (.FALSE.) THEN
C	WEIGHTED INTERPOLATION!!! -------------------------------------
	RTL = Rint - 0.5*dr !or maybe Rint - 0.5*dr
	RTR = Rint + 0.5*dr !or maybe Rint + 0.5*dr
	call GFMWINTERP (rho(q-3),rho(q-2),rho(q-1),rho(q),
     +		r(q-3),r(q-2),r(q-1),r(q),Rint,RTL,dr,rhointL)
	call GFMWINTERP (rho(q+4),rho(q+3),rho(q+2),rho(q+1),
     +		r(q+4),r(q+3),r(q+2),r(q+1),Rint,RTR,-1.0*dr,rhointR)
	call GFMWINTERP (u(q-3),u(q-2),u(q-1),u(q),
     +		r(q-3),r(q-2),r(q-1),r(q),Rint,RTL,dr,uintL)
	call GFMWINTERP (u(q+4),u(q+3),u(q+2),u(q+1),
     +		r(q+4),r(q+3),r(q+2),r(q+1),Rint,RTR,-1.0*dr,uintR)
	call GFMWINTERP (p(q-3),p(q-2),p(q-1),p(q),
     +		r(q-3),r(q-2),r(q-1),r(q),Rint,RTL,dr,pintL)
	call GFMWINTERP (p(q+4),p(q+3),p(q+2),p(q+1),
     +		r(q+4),r(q+3),r(q+2),r(q+1),Rint,RTR,-1.0*dr,pintR)
	END IF
C	---------------------------------------------------------------
C	SOLVE THE RIEMANN PROBLEM AT THE INTERFACE
	call RIEMANNINTERFACE (rhointL,pintL,uintL,gamma0,Pc0,
     +				rhointR,pintR,uintR,gamma1,Pc1,
     +				rhoLI,rhoRI,UI,PI,SL,SR)
C	===============================================================
c	LAST GHOST CELL APPLIES A REFLECTIVE BOUNDARY IN THE GHOST 
C	BAND. THIS ALLOWS CHEAP CALCULATIONS OF CONSERVATION PROPERTIES
C	WITHIN GHOST BANDS, AND HENCE IN EACH FLUID.
C	FLUID 0 -------------------------------------------------------
C	FLUID 0: REAL CELLS -------------------------------------------
	DO i=1,q
		m0(1,1,i) = rho(i)
		m0(1,2,i) = rho(i)*u(i)
		m0(1,3,i) = E(i)
		m0(1,4,i) = p(i)
		m0(1,5,i) = u(i)
	END DO
C	FLUID 0: ORDINARY GHOST CELLS ---------------------------------
	DO i=q+1,q+gfmbuffer
		m0(1,1,i) = rhoLI
		m0(1,4,i) = pI
		m0(1,5,i) = uI
		m0(1,2,i) = m0(1,1,i)*m0(1,5,i)
		m0(1,3,i) = (m0(1,4,i)+gamma0*Pc0)/
     +			(gamma0-1.0)+0.5*(m0(1,2,i)**2)/m0(1,1,i)
	END DO
C	FLUID 0: LAST GHOST CELL - APPLICATION OF REFLECTIVE BC -------
	m0(1,1,q+1+gfmbuffer) = m0(1,1,q+gfmbuffer)
	m0(1,2,q+1+gfmbuffer) = -1.0*m0(1,2,q+gfmbuffer)
	m0(1,3,q+1+gfmbuffer) = m0(1,3,q+gfmbuffer)
	m0(1,4,q+1+gfmbuffer) = m0(1,4,q+gfmbuffer)
	m0(1,5,q+1+gfmbuffer) = -1.0*m0(1,5,q+gfmbuffer)
C	FLUID 1 -------------------------------------------------------
C	FLUID 1: REAL CELLS -------------------------------------------
	DO i=q+1,nr+inttype-1
		m0(2,1,i) = rho(i)
		m0(2,2,i) = rho(i)*u(i)
		m0(2,3,i) = E(i)
		m0(2,4,i) = p(i)
		m0(2,5,i) = u(i)
	END DO
C	FLUID 1: ORDINARY GHOST CELLS ---------------------------------
	DO i=q+1-gfmbuffer,q
		m0(2,1,i) = rhoRI
		m0(2,4,i) = pI
		m0(2,5,i) = uI
		m0(2,2,i) = m0(2,1,i)*m0(2,5,i)
		m0(2,3,i) = (m0(2,4,i)+gamma1*Pc1)/
     +			(gamma1-1.0)+0.5*(m0(2,2,i)**2)/m0(2,1,i)
	END DO
C	FLUID 1: LAST GHOST CELL - APPLICATION OF REFLECTIVE BC -------
	m0(2,1,q-gfmbuffer) = m0(2,1,q+1-gfmbuffer)
	m0(2,2,q-gfmbuffer) = -1.0*m0(2,2,q+1-gfmbuffer)
	m0(2,3,q-gfmbuffer) = m0(2,3,q+1-gfmbuffer)
	m0(2,4,q-gfmbuffer) = m0(2,4,q+1-gfmbuffer)
	m0(2,5,q-gfmbuffer) = -1.0*m0(2,5,q+1-gfmbuffer)
C	===============================================================
	RETURN
	END
C	===============================================================
C	000000000000000000000000000000000000000000000000000000000000000
C	===============================================================
	SUBROUTINE RIEMANNINTERFACE (rhoL,pL,uL,GAML,PCRITL,rhoR,pR,uR,
     +					GAMR,PCRITR,
     +					rhoLstar,rhoRstar,ustar,pstar,
     +					bL,bR)
	IMPLICIT NONE
C	---------------------------------------------------------------
C	SUBROUTINE TO SOLVE A RIEMANN PROBLEM - ROE TYPE
C	---------------------------------------------------------------
C	SIMPLY RETURNS THE STAR STATES
c	http://dx.doi.org/10.1016/j.jcp.2009.06.002 (sort of)
C	---------------------------------------------------------------
C	AUTHOR: J. R. C. KING
C	---------------------
C	CHANGE RECORD:
C		15-07-2013:	CREATED
C		25-07-2013:	EMBEDDED IN SPLITFIELDS.f
C	---------------------------------------------------------------
	DOUBLE PRECISION	rhoL,PL,uL,GAML,PCRITL,rhoR,PR,uR,
     +				GAMR,PCRITR
C
	DOUBLE PRECISION	eL,eR,HL,HR,GL,GR,PHIL,PHIR
	DOUBLE PRECISION	pstar,ustar,rhoLstar,rhoRstar
	DOUBLE PRECISION	ENLstar,ENRstar
	DOUBLE PRECISION	cL,cR
	DOUBLE PRECISION	bL,bR,bM
	DOUBLE PRECISION	ubar,cbar,rhobar,Povrhobar,Pbar,hbar
	DOUBLE PRECISION	PHIbar,Gbar
	DOUBLE PRECISION	Gdebar,PHIdrhobar
C	---------------------------------------------------------------
c	CALCULATE LEFT AND RIGHT SPEEDS OF SOUND ======================
	cL = sqrt(GAML*(PL+PCRITL)/rhoL)
	cR = sqrt(GAMR*(PR+PCRITR)/rhoR)
C	CALCULATE LEFT AND RIGHT ENERGY, ENTHALPY, etc ================
	eL = (PL+GAML*PCRITL)/((GAML-1.0)*rhoL)
	eR = (PR+GAMR*PCRITR)/((GAMR-1.0)*rhoR)
	HL = eL + PL/rhoL + 0.5*uL**2.0
	HR = eR + PR/rhoR + 0.5*uR**2.0
	GL = GAML-1.0
	GR = GAMR-1.0
	PHIL = (GAML-1.0)*eL
	PHIR = (GAMR-1.0)*eR
C	CALCULATE ROE AVERAGES ((.)bar ) FOR PROPERTIES ===============
	rhobar = sqrt(rhoL*rhoR)
	call ROEAVG (rhoL,rhoR,uL,uR,ubar)
	call ROEAVG (rhoL,rhoR,HL,HR,hbar)
	Povrhobar = (PL/sqrt(rhoL) + PR/sqrt(rhoR))/
     +			(sqrt(rhoL)+sqrt(rhoR)) + 
     +			0.5*((uR-uL)/(sqrt(rhoL)+sqrt(rhoR)))**2
	Pbar = rhobar*Povrhobar
C	WE ARE LOOKING AT 2 FLUIDS! I KNOW THIS BECAUSE WE ALWAYS ARE -
		call ROEAVG (rhoL,rhoR,GL,GR,Gbar)
		call ROEAVG (rhoL,rhoR,PHIL,PHIR,PHIbar)
		cbar = sqrt(PHIbar + Gbar*Povrhobar)
C	---------------------------------------------------------------
C	CALCULATE STAR STATES AND WAVE SPEEDS =========================
	bL = min(uL-cL,ubar-cbar)
	bR = max(ubar+cbar,uR+cR)
	ustar = (rhoR*uR*(bR-uR) + rhoL*uL*(uL-bL) + PL - PR)/
     +		(rhoR*(bR-uR) + rhoL*(uL-bL))
	bM = ustar
	pstar = PL + rhoL*(uL-bL)*(uL-ustar)
C 	pstar = PR + rhoR*(bR-uR)*(ustar-uR)
	rhoLstar = rhoL*(bL-uL)/(bL-ustar)
	rhoRstar = rhoR*(bR-uR)/(bR-ustar)
	ENLstar = (pstar + GAML*PCRITL)/(GAML-1.0) + 
     +		0.5*rhoLstar*ustar**2
	ENRstar = (pstar + GAMR*PCRITR)/(GAMR-1.0) + 
     +		0.5*rhoRstar*ustar**2
C	===============================================================
C	===============================================================
C	APPROXIMATE SOLVER - LINEAR
C	cbar = 0.5*(cL+cR)
C	rhobar = 0.5*(rhoL+rhoR)
C	pstar = 0.5*(PL+PR) - 0.5*(UR-UL)*rhobar*cbar
C	ustar = 0.5*(uL+uR) - 0.5*(PR-PL)/(rhobar*cbar)
C	rhoLstar = rhoL + (uL-ustar)*rhobar/cbar
C	rhoRstar = rhoR + (ustar-uR)*rhobar/cbar
C	===============================================================
C	---------------------------------------------------------------
	RETURN
	END
C	---------------------------------------------------------------
C	---------------------------------------------------------------
C	###############################################################
C	===============================================================
	SUBROUTINE GFMWINTERP (prop1,prop2,prop3,prop4,
     +		r1,r2,r3,r4,Ri,Rt,delr,phiT)

	IMPLICIT NONE
C	---------------------------------------------------------------
C	---------------------------------------------------------------
C	SUBROUTINE TO DO WEIGHTED INTERPOLATION FOR GFM INPUT STATES
C	THE PURPOSE OF THIS SUBROUTINE IS TO CREATE GFM INPUT STATES 
C	WHICH VARY SMOOTHLY EVEN AS THE INTERFACE (AND HENCE INPUTS TO 
C 	THIS SUBROUTINE) PASSES A CELL CENTRE. IT DOESN'T IMPROVE THE
C	RESULTS...
C	---------------------------------------------------------------
C	---------------------------------------------------------------
C	AUTHOR: J. R. C. KING
C	---------------------
C	CHANGE RECORD:
C		06-11-2013:	CREATED
C	---------------------------------------------------------------
C	INPUTS
	DOUBLE PRECISION prop1,prop2,prop3,prop4
	DOUBLE PRECISION r1,r2,r3,r4,Ri,Rt,delr
C	INTERNAL
	DOUBLE PRECISION phiE(3)
	DOUBLE PRECISION x,w1,w2,w3
	INTEGER i

C	OUTPUTS
	DOUBLE PRECISION phiT
C	---------------------------------------------------------------
C	CALCULATE THE THREE ESTIMATES
	phiE(1) = prop1 + (prop2-prop1)*(rt-r1)/delr
	phiE(2) = prop2 + (prop3-prop2)*(rt-r2)/delr
	phiE(3) = prop3 + (prop4-prop3)*(rt-r3)/delr
C	CALCULATE THE WEIGHTS!
	x = (Ri-r4)/delr
	w1 = x - x**2.0
	w3 = x - x**2.0
	w2 = 1.0-(w1+w3)
C	CALCULATE THE WEIGHTED ESTIMATE
	phiT = w1*phiE(1) + w2*phiE(2) + w3*phiE(3)
	RETURN
	END 
C	---------------------------------------------------------------
C	===============================================================
