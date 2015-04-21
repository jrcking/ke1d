	SUBROUTINE FLUXES (mfl)

	INCLUDE "commonblock"
C	---------------------------------------------------------------
C	AUTHOR: J. R. C. KING
C	---------------------
C	CHANGE RECORD:
C		11-02-2014:	CREATED
C	---------------------------------------------------------------
C	CALCULATES FLUXES FROM CELL EDGE VALUES
C	---------------------------------------------------------------
C	---------------------------------------------------------------
	INTEGER mfl
C	---------------------------------------------------------------
C	SOLVE THE RIEMANN PROBLEMS TO FIND THE FLUXES -----------------	
	DO i=rLOW(mfl),rHIGH(mfl)-1
C		CALCULATE AT RIGHT "CELL WALLS": i,i+1
C		LEFT DENSITY, MOMENTUM, ENERGY,
C		RIGHT DENSITY, MOMENTUM, ENERGY,
C		FLUID PROPERTIES (gamma and Pc),
C		FLUXES OF DENSITY, MOMENTUM, ENERGY

!!!!!		Note this little switch which lets me use AUSM when
!!!!!		if feel like it! Note AUSM unstable for big shocks
		IF (.false.)THEN 
			call AUSM
     +			(Uiph(1,i),Uiph(2,i),
     +			Uiph(3,i),
     +			Uimh(1,i+1),Uimh(2,i+1),Uimh(3,i+1),
     +			gm(mfl),pcr(mfl),
     +			fiph(1,i),fiph(2,i),fiph(3,i))
		ELSE
			call RIEMANNHLLC 
     +			(Uiph(1,i),Uiph(2,i),
     +			Uiph(3,i),
     +			Uimh(1,i+1),Uimh(2,i+1),Uimh(3,i+1),
     +			gm(mfl),pcr(mfl),
     +			fiph(1,i),fiph(2,i),fiph(3,i))
		END IF
	END DO
C	---------------------------------------------------------------
	RETURN
	END
C	---------------------------------------------------------------
