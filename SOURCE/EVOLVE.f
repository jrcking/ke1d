	SUBROUTINE EVOLVE (mfl)

	INCLUDE "commonblock"
C	---------------------------------------------------------------
C	AUTHOR: J. R. C. KING
C	---------------------
C	CHANGE RECORD:
C		20-07-2013:	CREATED
C		07-08-2013:	SPHERICALITY BY SOURCE TERMS
C		07-08-2013:	MODIFIED TO WORK ON CONSERVATIVE VARS
C		14-08-2013:	STOPPED SOLVING DUPLICATE RIEMANN PROBS
C		14-08-2013:	ONLY USE NECESSARY GHOST CELLS..
C		25-09-2013:	3rd ORDER TVD TIME INTEGRATION OPTION!
C	---------------------------------------------------------------
C	EVOLVE PROPERTIES GIVEN FLUXES. SOURCE TERMS NOT TIME-SPLIT.
C	---------------------------------------------------------------
C	---------------------------------------------------------------
	INTEGER mfl
	DOUBLE PRECISION delV
	DOUBLE PRECISION aph,amh
	DOUBLE PRECISION wr,wu,wm,we,wp,S(3)
c
	DOUBLE PRECISION dpdr,pp,pm
	DOUBLE PRECISION sp1,sp2,sp3,sp4,sp5
C	---------------------------------------------------------------
C	USE THE FLUXES TO UPDATE THE PROPERTIES -----------------------
	DO i=rLOW(mfl)+1,rHIGH(mfl)-1
C		-------------------------------------------------------
		delV = dr
		aph = 1.0
		amh = 1.0
		if(.false.)then !alternative method...
		delV = (1.0/(coordsno+1.0))*((r(i)+0.5*dr)**
     +		(coordsno+1.0)-(r(i)-0.5*dr)**(coordsno+1.0))
		aph = (r(i)+0.5*dr)**coordsno
		amh = (r(i)-0.5*dr)**coordsno
		end if
C		-------------------------------------------------------
C		HAVE NOT BOTHERED TO OPERATOR SPLIT... ----------------
C		-------------------------------------------------------
		S(1)=m0(mfl,2,i)
		S(2)=m0(mfl,2,i)*m0(mfl,2,i)/m0(mfl,1,i)
		S(3)=m0(mfl,2,i)*(m0(mfl,3,i)+m0(mfl,4,i))/m0(mfl,1,i)
		if(.false.)then ! alternative method...
		S(1)=0
		S(2)=0
		S(3)=0
		end if
		DO j=1,3
			m0(mfl,j,i) = sfrac(1,kt)*m1(mfl,j,i) +
     +			sfrac(2,kt)*m0(mfl,j,i) - sfrac(3,kt)*
     +			((dt/delV)*(aph*fiph(j,i)-amh*fiph(j,i-1))
     +			+dt*coordsno*S(j)/r(i))
		END DO
		m0(mfl,4,i) = (gm(mfl)-1.0)*(m0(mfl,3,i)-0.5*
     +			m0(mfl,2,i)*m0(mfl,2,i)/m0(mfl,1,i)) - 
     +			gm(mfl)*pcr(mfl)
		m0(mfl,5,i) = m0(mfl,2,i)/m0(mfl,1,i)
	END DO
C	===============================================================
C	---------------------------------------------------------------
	RETURN
	END
C	---------------------------------------------------------------
