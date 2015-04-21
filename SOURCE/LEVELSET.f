	SUBROUTINE LEVELSET
	
	INCLUDE "commonblock"
C	---------------------------------------------------------------
C	AUTHOR: J. R. C. KING
C	---------------------
C	CHANGE RECORD:
C		20-07-2013:	CREATED
C		06-08-2013:	MODIFIED TO UPSTREAM WENO-Z SCHEME
C		21-11-2013:	TIME INTEGRATION SAME ORDER AS EULER
C	---------------------------------------------------------------
C	THIS SUBROUTINE UPDATES A LEVEL SET ACCORDING TO:
C	da/dt + uLS*da/dr = 0
C	A 5th ORDER UPSTREAM WENO-Z SCHEME IS USED FOR SPATIAL DISCRET-
C	ISATION. TIME INTEGRATION IS OF THE ORDER OF THE TIME INTEGRAT-
C	ION OF THE EULER SOLUTION (set by inttype).
C	---------------------------------------------------------------
C	---------------------------------------------------------------
	DOUBLE PRECISION	dadr(nrmax)
	DOUBLE PRECISION	store1,spare1,aUPD(nrmax)
	DOUBLE PRECISION	uls
C	---------------------------------------------------------------
C	DECIDE THE LEVEL SET VELOCITY! --------------------------------
C	OPTION 1 - the ghost cell velocity
	uls = uI
C	---------------------------------------------------------------
C	SET THE VALUES ------------------------------------------------
	DO i=1,nr
		aUPD(i) = alpha(i)
	END DO
C	---------------------------------------------------------------
C	ORDER OF INTEGRATION LOOP -------------------------------------
	DO k=1,inttype
C		-------------------------------------------------------
C		USE WENOZ TO CALCULATE DERIVATIVE OF ALPHA ------------
		IF (uI .GT. 0.0) THEN
			DO i=4,nr-2
C				LEFT SIDED DERIVATIVES ----------------
				call WENOZ ( (aUPD(i-2)-aUPD(i-3)),
     +				(aUPD(i-1)-aUPD(i-2)),
     +				(aUPD(i)-aUPD(i-1)),(aUPD(i+1)-aUPD(i))
     +				,(aUPD(i+2)-aUPD(i+1)),spare1,store1)
				dadr(i) = store1/dr
			END DO
		ELSE
			DO i=3,nr-3
C				RIGHT SIDED DERIVATIVES ---------------
				call WENOZ ((aUPD(i+3)-aUPD(i+2)),
     +				(aUPD(i+2)-aUPD(i+1)),
     +				(aUPD(i+1)-aUPD(i)),(aUPD(i)-aUPD(i-1))
     +				,(aUPD(i-1)-aUPD(i-2)),spare1,store1)
				dadr(i) = store1/dr
			END DO
		END IF
C		-------------------------------------------------------
C		BODGE BOUNDARY CONDITIONS(ACTUALLY EXACT FOR THIS CASE)
		DO i=1,3
			dadr(i) = dadr(4)
			dadr(i) = dadr(4)
		END DO
		DO i=nr-2,nr
			dadr(i) = dadr(nr-3)
			dadr(i) = dadr(nr-3)
		END DO
C		UPDATE THE LEVEL SET ----------------------------------
		DO i=1,nr
			aUPD(i)=sfrac(1,k)*alpha(i)+sfrac(2,k)*
     +			aUPD(i)	- sfrac(3,k)*dt*uls*dadr(i)
		END DO
	END DO
C	---------------------------------------------------------------
C	RETURN NEW VALUES TO ALPHA ------------------------------------
	DO i=1,nr
		alpha(i) = aUPD(i)
	END DO
C	---------------------------------------------------------------
	RETURN
	END
C	---------------------------------------------------------------
