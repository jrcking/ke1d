	SUBROUTINE RECONSTRUCT (mfl)

	INCLUDE "commonblock"
C	---------------------------------------------------------------
C	AUTHOR: J. R. C. KING
C	---------------------
C	CHANGE RECORD:
C		11-02-2014:	CREATED
C	---------------------------------------------------------------
C	RECONSTRUCTION SUBROUTINE. GIVEN FIELDS m0(mfl,.,.), 
C	RECONSTRUCTS TO FIND Uiph(.,.) AND Uimh(.,.). CHOICE OF PIECE-
C	WISE CONSTANT, LINEAR OR WENO-Z RECONSTRUCTION, DETERMINED BY
C	VALUE OF soltype
C	---------------------------------------------------------------
C	---------------------------------------------------------------
	INTEGER mfl,sltype
	DOUBLE PRECISION delU
C	---------------------------------------------------------------
C	PIECEWISE CONSTANT RECONSTRUCTION -----------------------------
	DO i=rLOW(mfl),rHIGH(mfl)
		DO j=1,5
			Uiph(j,i) = m0(mfl,j,i)
			Uimh(j,i) = m0(mfl,j,i)
		END DO
	END DO
C	---------------------------------------------------------------
C	OR FOR MUSCLE SCHEME IF REQUESTED! ----------------------------
C	PIECEWISE LINEAR RECONSTRUCTION OF DENSITY, VELOCITY, ---------
C	RESSURE AND ENERGY. MOMENTUM CALCULATED FROM DENSITY ----------
C	AND VELOCITY --------------------------------------------------
	IF (soltype .GT. 1) THEN
		DO i=rLOW(mfl)+1,rHIGH(mfl)-1
C		SNEAKY REDUCTION TO 1ST ORDER BESIDE INTERFACE --------
C			IF(mfl.EQ.1.AND.i.GE.q-1)THEN
C				sltype = 0
C			ELSE IF(mfl.EQ.2.AND.i.LE.q+2) THEN
C				sltype =0
C			ELSE
C				sltype = 1
C			END IF
C			OR I CAN JUST USE THE MUSCL SCHEME IN THE 
C			BUBBLE AND NOT IN THE WATER!!
			IF (mfl.EQ.1.AND.i.LT.q-2)THEN
				sltype = 1
			ELSE
				sltype = 0
			END IF
			DO j=1,5
				call SLOPELIMITER (sltype,
     +				(m0(mfl,j,i)-m0(mfl,j,i-1))/
     +				(m0(mfl,j,i+1)-m0(mfl,j,i)),delU)
				Uiph(j,i) = m0(mfl,j,i) + 
     +				delU*(m0(mfl,j,i+1)-m0(mfl,j,i))/2
				Uimh(j,i) = m0(mfl,j,i) - 
     +				delU*(m0(mfl,j,i+1)-m0(mfl,j,i))/2
			END DO
			Uiph(2,i) = Uiph(1,i)*Uiph(5,i)
			Uimh(2,i) = Uimh(1,i)*Uimh(5,i)
		END DO
	END IF
C	---------------------------------------------------------------
C	---------------------------------------------------------------
C	WHY NOT TRY A WENO-Z RECONSTRUCTION OF THE PRIMITIVE ----------
C	VARIABLES? ----------------------------------------------------
	IF (soltype .EQ. 5) THEN
		DO i=rLOW(mfl)+2,rHIGH(mfl)-2
			DO j=1,5
				call WENOZ(m0(mfl,j,i-2),
     +				m0(mfl,j,i-1),m0(mfl,j,i),
     +				m0(mfl,j,i+1),m0(mfl,j,i+2),
     +				Uimh(j,i),Uiph(j,i))
			END DO
C			MOMENTUM AND ENERGY SET FROM PRIMITIVE --------
C			VARIABLES -------------------------------------
			Uiph(2,i) = Uiph(1,i)*Uiph(5,i)
			Uimh(2,i) = Uimh(1,i)*Uimh(5,i)
			Uiph(3,i) = (Uiph(4,i)+gm(mfl)*
     +			pcr(mfl))/(gm(mfl)-1.0)+0.5*Uiph(1,i)*
     +			Uiph(5,i)**2.0
			Uimh(3,i) = (Uimh(4,i)+gm(mfl)*
     +			pcr(mfl))/(gm(mfl)-1.0)+0.5*Uimh(1,i)*
     +			Uimh(5,i)**2.0
		END DO
	END IF
C	===============================================================
C	---------------------------------------------------------------
	RETURN
	END
C	---------------------------------------------------------------
