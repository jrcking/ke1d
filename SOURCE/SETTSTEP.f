	SUBROUTINE SETTSTEP
	
	INCLUDE "commonblock"
C	---------------------------------------------------------------
C	AUTHOR: J. R. C. KING
C	---------------------
C	CHANGE RECORD:
C		20-11-2012:	CREATED
C	---------------------------------------------------------------
C	SUBROUTINE TO CALCULATE THE TIME STEP BASED ON CFL NUMBER, GRID
C	SIZE, MAXIMUM VELOCITY AND MAXIMUM SOUND SPEED
C	---------------------------------------------------------------
C	---------------------------------------------------------------
	DOUBLE PRECISION vsmax,vsig(nrmax)
C	---------------------------------------------------------------
C	FIND THE BIGGEST SIGNAL SPEED ---------------------------------
	DO i=1,nr
		IF (alpha(i) .LE. 0 ) THEN
			vsig(i) = u(i)+sqrt(gm(1)*(p(i)+Pcr(1))/rho(i))
		ELSE
			vsig(i) = u(i)+sqrt(gm(2)*(p(i)+Pcr(2))/rho(i))		
		END IF
	END DO
	vsmax = maxval(vsig,pmask)
C	---------------------------------------------------------------
C	SET THE TIME STEP ---------------------------------------------
	dt = CFL*dr/(vsmax)

c	REMOVE THIS WHEN FINISHED DEBUGGING w.r.t. 2D code
c	dt = 1e-5	
c
C	---------------------------------------------------------------
	RETURN
	END
C	---------------------------------------------------------------
