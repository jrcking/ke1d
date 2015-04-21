	SUBROUTINE SLOPELIMITER (sltype,x,B)

	IMPLICIT NONE
C	---------------------------------------------------------------
C	AUTHOR: J. R. C. KING
C	---------------------
C	CHANGE RECORD:
C		??-05-2013:	CREATED
C		03-02-2014: 	MODIFIED WITH GODUNOV OPTION
C	---------------------------------------------------------------
C	A SLOPE LIMITER SUBROUTINE FOR USE WITH MUSCL SCHEME
C	---------------------------------------------------------------
C	---------------------------------------------------------------
	INTEGER	sltype
	DOUBLE PRECISION x,B
C	---------------------------------------------------------------
	IF (sltype .EQ. 0) THEN
		B = 0
	ELSE
C		Monotonised Central
C		B = max(0.0,min(2.0*x,0.5*(1.0+x),2.0))
C		MINMOD
		B = max(0.0,min(1.0,x))
C		Superbee
C		B = max(0.0,min(2.0*x,1.0),min(x,2.0))
	END IF
C	---------------------------------------------------------------
	RETURN
	END
C	---------------------------------------------------------------
