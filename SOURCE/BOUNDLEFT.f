	SUBROUTINE BOUNDLEFT (mfl)

	INCLUDE "commonblock"
C	---------------------------------------------------------------
C	AUTHOR: J. R. C. KING
C	---------------------
C	CHANGE RECORD:
C		20-07-2013:	CREATED
C		06-08-2013:	MODIFIED TO UPSTREAM WENO-Z SCHEME
C		30-09-2013:	FOR A GENERAL NODE 'ibn' AT RIGHT BOUND
C		22-10-2013:	BODGED. NB: in some cases, especially
C					when using 5th order WENO-Z and
C					3rd order TVD, variation at or-
C					-igin is such that 1-sided 3rd
C					order difference equation is 
C					unstable. Hence the bodge.
C		11-02-2014:	MODIFIED TO FIT NEW PROGRAM STRUCTURE
C	---------------------------------------------------------------
C	THIS SUBROUTINE USES THE STANDARD CHARACTERISTIC BOUNDARY COND-
C	ITION FORMULATION PRESENTED IN Thompson, JCP 89, 439-461(1990)
C	TO APPLY A REFLECTIVE BOUNDARY AT THE LEFT EDGE OF THE DOMAIN
C	---------------------------------------------------------------
C	---------------------------------------------------------------
	INTEGER mfl
	DOUBLE PRECISION drhodrB,dudrB,dpdrB
	DOUBLE PRECISION cB,L1B,L2B,L3B
	DOUBLE PRECISION dudtB,drhodtB,dpdtB
	DOUBLE PRECISION newrho,newp,newu
C	---------------------------------------------------------------
C	---------------------------------------------------------------
C	CALCULATE THE LOCAL SPEED OF SOUND ----------------------------
	cB = sqrt(gm(mfl)*(m0(mfl,4,rLOW(mfl))+pcr(mfl))/
     +		m0(mfl,1,rLOW(mfl)))
C	---------------------------------------------------------------
C	FIRST CALCULATE DERIVATIVES OF THE PRIMITIVE VARIABLES AT THE - 
C	BOUNDARY: 1 SIDED 3rd ORDER -----------------------------------
	drhodrB = (11.0*m0(mfl,1,rLOW(mfl))- 18.0*m0(mfl,1,rLOW(mfl)+1) 
     +			+ 9.0*m0(mfl,1,rLOW(mfl)+2) - 
     +			2.0*m0(mfl,1,rLOW(mfl)+3))/(-6.0*dr)
	dudrB = (11.0*m0(mfl,5,rLOW(mfl)) - 18.0*m0(mfl,5,rLOW(mfl)+1) 
     +			+ 9.0*m0(mfl,5,rLOW(mfl)+2) - 
     +			2.0*m0(mfl,5,rLOW(mfl)+3))/(-6.0*dr)
	dpdrB = (11.0*m0(mfl,4,rLOW(mfl)) - 18.0*m0(mfl,4,rLOW(mfl)+1) 
     +			+ 9.0*m0(mfl,4,rLOW(mfl)+2) - 
     +			2.0*m0(mfl,4,rLOW(mfl)+3))/(-6.0*dr)
C	---------------------------------------------------------------
C	REFLECTIVE BOUNDARY! ------------------------------------------	
	IF (m0(mfl,5,rLOW(mfl)) - cB .LT. 0.0) THEN
C		L1 OUTGOING - FROM DEFINITION
		L1B = (m0(mfl,5,rLOW(mfl))-cB)*(dpdrB - 
     +			m0(mfl,1,rLOW(mfl))*cB*dudrB)
	ELSE
C		L1 INCOMING - SET TO ZERO (THIS SHOULD NEVER HAPPEN)
		L1B = 0.0
	END IF
	IF (m0(mfl,5,rLOW(mfl)) .LT. 0.0) THEN
C		L2 OUTGOING
		L2B = m0(mfl,5,rLOW(mfl))*(cB*drhodrB - dpdrB)
	ELSE
C		L2 INCOMING
		L2B = 0.0
	END IF
	IF (m0(mfl,5,rLOW(mfl)) + cB .LT. 0.0) THEN
C		L3 OUTGOING (THIS SHOULD NEVER HAPPEN)
		L3B = (m0(mfl,5,rLOW(mfl))+cB)*(dpdrB + 
     +			m0(mfl,1,rLOW(mfl))*cB*dudrB)
	ELSE
C		L3 INCOMING - SET TO L1B
		L3B = L1B
	END IF
C	---------------------------------------------------------------
C	CALCULATE TIME DERIVATIVES OF CHARACTERISTIC EULER EQUATIIONS -
C	USE abs(r(1)) RATHER THAN r(1) FOR SOME LONG FORGOTTEN REASON -
	dudtB = (-0.5/(m0(mfl,1,rLOW(mfl))*cB))*(L3B-L1B)
	drhodtB = -1.0*((L2B + 0.5*(L3B+L1B))/cB**2
     +			+ coordsno*m0(mfl,1,rLOW(mfl))*
     +			m0(mfl,5,rLOW(mfl))/abs(r(rLOW(mfl)))
     +			)
	dpdtB = -1.0*(0.5*(L3B+L1B)
     +			+ coordsno*m0(mfl,1,rLOW(mfl))*
     +			m0(mfl,5,rLOW(mfl))*cB*cB/abs(r(rLOW(mfl)))
     +			)
C	---------------------------------------------------------------
C	1st ORDER TIME INTEGRATION ------------------------------------
	newrho = m0(mfl,1,rLOW(mfl)) + drhodtB*dt
	newp = m0(mfl,4,rLOW(mfl)) + dpdtB*dt
	newu = m0(mfl,5,rLOW(mfl)) + dudtB*dt
C	---------------------------------------------------------------
C	THIS IS A BODGE I USE SOMETIMES (all the time) ----------------
	m0(mfl,1,rLOW(mfl)) = m0(mfl,1,rLOW(mfl)+1)
	m0(mfl,4,rLOW(mfl)) = m0(mfl,4,rLOW(mfl)+1)
	m0(mfl,5,rLOW(mfl)) = -1.0*m0(mfl,5,rLOW(mfl)+1)
	m0(mfl,2,rLOW(mfl)) = m0(mfl,5,rLOW(mfl))*m0(mfl,1,rLOW(mfl))
	m0(mfl,3,rLOW(mfl)) = m0(mfl,3,rLOW(mfl)+1)
C	---------------------------------------------------------------
	RETURN
	END 
C	---------------------------------------------------------------
