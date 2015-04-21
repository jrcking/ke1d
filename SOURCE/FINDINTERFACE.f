	SUBROUTINE FINDINTERFACE

	INCLUDE "commonblock"
C	---------------------------------------------------------------
C	AUTHOR: J. R. C. KING
C	---------------------
C	CHANGE RECORD:
C		??-05-2013:	CREATED
C		06-08-2013:	MODIFIED TO USE WENO-Z SCHEME DERIVS
C		21-11-2013:	SET SPATIAL DERIV TO 1 FOR NOW
C	---------------------------------------------------------------
C	SUBROUTINE TO FIND THE LOCATION OF THE ZERO LEVEL-SET
C	---------------------------------------------------------------
C	---------------------------------------------------------------
	DOUBLE PRECISION	av_alpha,dalphadr
	DOUBLE PRECISION 	fp,fm
C	---------------------------------------------------------------
C	FIND THE INDICES OF THE CELLS EITHER SIDE OF THE ZERO LEVEL SET
	DO i=1,nr
		IF (alpha(i)*alpha(i+1) .LE. 0.0) THEN
			IF (i .GT. q-4 .AND. i .LE. q+4) THEN			
				q = i
			END IF
		END IF
	END DO
C	CALCULATE SPATIAL DERIVATIVES WITH A 5th ORDER WENO-Z SCHEME --
C	EXTRAPOLATE OUT NEAR THE BOUNDARIES. IT IS A BODGE, BUT AN ----
C	ACCURATE ONE AT PRESENT. --------------------------------------
	IF (q .GE. 3 .AND. q .LE. nr-2) THEN
		dalphadr = (alpha(q-2) - 8*alpha(q-1) + 
     +			8*alpha(q+1) - alpha(q+1))/(12.0*dr)
		call WENOZ(alpha(q-2),alpha(q-1),alpha(q),alpha(q+1),
     +			alpha(q+2),fm,fp)
		dalphadr = (fp-fm)/dr
	ELSE IF (q .LT. 3) THEN
		dalphadr = (alpha(q+1) - alpha(q))/dr
	ELSE
		dalphadr = (alpha(q) - alpha(q-1))/dr
	END IF
C	NB: In the 1D case the level set is advected with a uniform
C	velocity field, so the spatial derivative is constant in time.
C	It is initially set to 1, and so will remain equal to 1. This 
C	bodge has virtually no effect on results, but MIGHT (if I 
C	remember correctly) help with analysis of conservation 
C	properties of the rGFM...
	dalphadr = 1.0  !this is a big cheat, but changes nothing in 
c			in this instance
C	INTERPOLATE SOMEHOW TO FIND Rint ------------------------------
C	NB: Need to stick a more accurate interpolation here sometime
	Rint = r(q) - alpha(q)/
     +			max(dalphadr,-1.0*alpha(q)/dr)
C	---------------------------------------------------------------
C	JUST IN CASE IT LIES EXACTLY ON A NODE!! (SEEMS UNLIKELY) -----
	DO i=1,nr
		IF (alpha(i) .EQ. 0.0) THEN
			WRITE (4,*) "CRIKEY!! INTERFACE ON NODE ",i,"!"
			q = i
			Rint = r(q)
		END IF
	END DO
C	---------------------------------------------------------------
	RETURN
	END 
C	---------------------------------------------------------------
