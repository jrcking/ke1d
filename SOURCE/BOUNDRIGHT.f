	SUBROUTINE BOUNDRIGHT (mfl)

	INCLUDE "commonblock"
C	---------------------------------------------------------------
C	AUTHOR: J. R. C. KING
C	---------------------
C	CHANGE RECORD:
C		20-07-2013:	CREATED
C		06-08-2013:	MODIFIED TO UPSTREAM WENO-Z SCHEME
C		30-09-2013:	FOR A GENERAL NODE 'ibn' AT RIGHT BOUND
C		15-02-2014: 	THOROUGHLY RE-SHUFFLED
C	---------------------------------------------------------------
C	THIS SUBROUTINE USES THE STANDARD CHARACTERISTIC BOUNDARY COND-
C	ITION FORMULATION PRESENTED IN Thompson, JCP 89, 439-461(1990)
C	TO APPLY BOUNDARY CONDITIONS AT THE RIGHT EDGE OF THE DOMAIN
C
C	USES SOLUTION AT RIGHT BOUNDARY CELL TO APPLY NLAA BC OR OTHER
C	BC BY PROVIDING NEW VALUES FOR THIS CELL ----------------------
C	---------------------------------------------------------------
C	---------------------------------------------------------------
	INTEGER mfl
	DOUBLE PRECISION drhodrB,dudrB,dpdrB
	DOUBLE PRECISION cB,L1B,L2B,L3B
	DOUBLE PRECISION henth
	DOUBLE PRECISION aL1,NA,NB
	DOUBLE PRECISION dudtB,drhodtB,dpdtB
	DOUBLE PRECISION newrho,newp,newu
C	---------------------------------------------------------------
C	---------------------------------------------------------------
C	CALCULATE THE LOCAL SPEED OF SOUND ----------------------------
	cB = sqrt(gm(mfl)*(m0(mfl,4,rHIGH(mfl))+pcr(mfl))/
     +		m0(mfl,1,rHIGH(mfl)))
C	---------------------------------------------------------------
C	FIRST CALCULATE DERIVATIVES OF PRIMITIVE VARIABLES AT THE -----
C	BOUNDARY: 1 SIDED 3rd ORDER -----------------------------------
	drhodrB = (-11.0*m0(mfl,1,rHIGH(mfl)) + 
     +			18.0*m0(mfl,1,rHIGH(mfl)-1) - 
     +			9.0*m0(mfl,1,rHIGH(mfl)-2) + 
     +			2.0*m0(mfl,1,rHIGH(mfl)-3))/(-6.0*dr)
	dudrB = (-11.0*m0(mfl,5,rHIGH(mfl)) + 
     +			18.0*m0(mfl,5,rHIGH(mfl)-1) - 
     +			9.0*m0(mfl,5,rHIGH(mfl)-2) + 
     +			2.0*m0(mfl,5,rHIGH(mfl)-3))/(-6.0*dr)
	dpdrB = (-11.0*m0(mfl,4,rHIGH(mfl)) + 
     +			18.0*m0(mfl,4,rHIGH(mfl)-1) - 
     +			9.0*m0(mfl,4,rHIGH(mfl)-2) + 
     +			2.0*m0(mfl,4,rHIGH(mfl)-3))/(-6.0*dr)
C	---------------------------------------------------------------
C	DETERMINE THE DIRECTIONS OF THE CHARACTERISTICS AND SET VALUES
C	ACCORDINGLY ---------------------------------------------------
	IF (m0(mfl,5,rHIGH(mfl)) + cB .LE. 0.0) THEN
C		L3 INCOMING
		L3B = 0.0
	ELSE	
C		L3 OUTGOING
		L3B = (m0(mfl,5,rHIGH(mfl))+cB)*(dpdrB + 
     +			m0(mfl,1,rHIGH(mfl))*cB*dudrB)
	END IF
C	---------------------------------------------------------------
C	THE TYPE OF BOUNDARY CONDITION IS SET BY THE obtype FLAG ------
C	IF WE'RE LOOKING AT FLUID 1 WE FORCE A REFLECTING (dudt=0) 
C	BOUNDARY. 
	IF (m0(mfl,5,rHIGH(mfl))-cB .LT. 0.0) THEN
C		L1 INCOMING -------------------------------------------
		IF (obtype .EQ. 1 .AND. mfl .EQ. 2) THEN
C			NON-REFLECTING, AS IN Thompson II, 1993 -------
			L1B = 0.0
			dudtB = (-0.5/(m0(mfl,1,rHIGH(mfl))*cB))*
     +				(L3B-L1B)
		ELSE IF (obtype .EQ. 2 .OR. mfl .EQ. 1) THEN
C			CONSTANT VELOCITY!! ALSO COURTESY OF Thompson -
			L1B = L3B
			dudtB = (-0.5/(m0(mfl,1,rHIGH(mfl))*cB))*
     +				(L3B-L1B)
		ELSE IF (obtype .EQ. 3 .AND. mfl .EQ. 2) THEN
C			NLAA BCs!!!
C			CALCULATE THE RHS AND VELOCITY TIME DERIVATIVE
			henth = (m0(mfl,4,rHIGH(mfl))-Pinf)/
     +				m0(mfl,1,rHIGH(mfl))
			aL1 = (m0(mfl,1,rHIGH(mfl))*
     +				(m0(mfl,5,rHIGH(mfl))-cB)/
     +				r(rHIGH(mfl)))*((3-m0(mfl,5,rHIGH(mfl))
     +				/cB)*0.5*m0(mfl,5,rHIGH(mfl))**2 + 
     +				2*m0(mfl,5,rHIGH(mfl))*cB - 
     +				(1+m0(mfl,5,rHIGH(mfl))/cB)*henth)
			NA = 1.0 - 0.5*(m0(mfl,5,rHIGH(mfl))-cB)*
     +				m0(mfl,5,rHIGH(mfl))/cB**2
			NB = (1/(2*m0(mfl,1,rHIGH(mfl))*cB))*(L3B-aL1)
			dudtB = -1.0*NB/NA
			L1B = aL1 + m0(mfl,1,rHIGH(mfl))*
     +				(m0(mfl,5,rHIGH(mfl))-cB)*
     +				m0(mfl,5,rHIGH(mfl))*dudtB/cB
		END IF
	ELSE
C		L1 OUTGOING -------------------------------------------
		L1B = (m0(mfl,5,rHIGH(mfl))-cB)*(dpdrB - 
     +			m0(mfl,1,rHIGH(mfl))*cB*dudrB)
		dudtB = (-0.5/(m0(mfl,1,rHIGH(mfl))*cB))*(L3B-L1B)
	END IF
	IF (m0(mfl,5,rHIGH(mfl)) .LT. 0.0) THEN
C		L2 INCOMING -------------------------------------------
		L2B = (m0(mfl,1,rHIGH(mfl))*m0(mfl,5,rHIGH(mfl))*
     +			m0(mfl,5,rHIGH(mfl))/r(rHIGH(mfl)))*( 
     +			(henth/cB) - (2*m0(mfl,5,rHIGH(mfl))) + 
     +			(0.5*m0(mfl,5,rHIGH(mfl))*
     +			m0(mfl,5,rHIGH(mfl))/cB)) + 
     +			m0(mfl,1,rHIGH(mfl))*m0(mfl,5,rHIGH(mfl))*
     +			(1-m0(mfl,5,rHIGH(mfl))/cB)*dudtB
C		THIS IS THE ZERO ENTROPY FLUX CHEAT CAUSING MINIMAL  --
C		ERRORS... ---------------------------------------------
		L2B = 0.0
	ELSE
C		L2 OUTGOING -------------------------------------------
		L2B = m0(mfl,5,rHIGH(mfl))*(cB*drhodrB - dpdrB)
	END IF
C	---------------------------------------------------------------
C	DETERMINE THE TIME DERIVATIVES --------------------------------
	drhodtB = -1.0*((L2B + 0.5*(L3B+L1B))/cB**2
     +			+ coordsno*m0(mfl,1,rHIGH(mfl))*
     +			m0(mfl,5,rHIGH(mfl))/r(rHIGH(mfl))
     +			)
	dpdtB = -1.0*(0.5*(L3B+L1B)
     +			+ coordsno*m0(mfl,1,rHIGH(mfl))*
     +			m0(mfl,5,rHIGH(mfl))*cB*cB/r(rHIGH(mfl))
     +			)
C	---------------------------------------------------------------
C	1st ORDER TIME INTEGRATION ------------------------------------
	newrho = sfrac(1,kt)*m1(mfl,1,rHIGH(mfl)) + sfrac(2,kt)*
     +		m0(mfl,1,rHIGH(mfl)) + sfrac(3,kt)*drhodtB*dt
	newp = sfrac(1,kt)*m1(mfl,4,rHIGH(mfl)) + sfrac(2,kt)*
     +		m0(mfl,4,rHIGH(mfl)) + sfrac(3,kt)*dpdtB*dt
	newu = sfrac(1,kt)*m1(mfl,5,rHIGH(mfl)) + sfrac(2,kt)*
     +		m0(mfl,5,rHIGH(mfl)) + sfrac(3,kt)*dudtB*dt
C	---------------------------------------------------------------
	m0(mfl,1,rHIGH(mfl)) = newrho
	m0(mfl,4,rHIGH(mfl)) = newp
	m0(mfl,5,rHIGH(mfl)) = newu
	m0(mfl,2,rHIGH(mfl)) = newrho*newu	
	m0(mfl,3,rHIGH(mfl)) = (newp+gm(mfl)*pcr(mfl))/(gm(mfl)-1) + 
     +			0.5*newrho*newu**2
C	---------------------------------------------------------------
	RETURN
	END 
C	---------------------------------------------------------------
