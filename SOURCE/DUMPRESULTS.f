	SUBROUTINE DUMPRESULTS

	INCLUDE "commonblock"
C	---------------------------------------------------------------
C	AUTHOR: J. R. C. KING
C	---------------------
C	CHANGE RECORD:
C		20-11-2012:	CREATED
C		??<08-2013:	MODIFIED FOR FAR FIELD OUTPUT IF REQ
C	---------------------------------------------------------------
C	SUBROUTINE TO WRITE TO OUTPUTS IF IT IS THE REQUIRED TIME
C	---------------------------------------------------------------
C	---------------------------------------------------------------
	INTEGER ofspc
C	---------------------------------------------------------------
C	EITHER WRITE TO OUTPUT OR ADD 1 TO THE OUTPUT COUNTER ---------
	IF (ofc1 .GE. outfreq) THEN
C		IT IS TIME TO WRITE TO OUTPUT ------------------------
C		SET THE SPATIAL OUTPUTTING FREQUENCY ------------------
		ofspc = int(float(nr)/float(nrout))
C		WRITE TO SOME FILES!!
		DO k=1,int(float(nr)/float(ofspc))
			i=int(float(ofspc)*float(k))
			WRITE (7,*) r(i),rho(i)
			WRITE (8,*) r(i),p(i)
			WRITE (9,*) r(i),u(i)
			WRITE (10,*) r(i),gamm(i)
			WRITE (11,*) r(i),Pc(i)
			WRITE (12,*) r(i),E(i)
			WRITE (13,*) r(i),alpha(i)
		END DO
C		RESET THE OUTPUT COUNTER ------------------------------
		ofc1 = 1
	ELSE
C		IT ISN'T TIME TO WRITE TO OUTPUT SO ADD 1 TO THE ------
C		OUTPUT COUNTER ----------------------------------------
		ofc1 = ofc1 + 1
	END IF
C	WRITE OTHER THINGS TO OUTPUT FILES..
	WRITE (20,*) t,dt
	WRITE (22,*) t,Rint
	WRITE (23,*) t,P(q)
	WRITE (24,*) t,p(nr)
	WRITE (25,*) t,u(nr)
	WRITE (26,*) t,rho(nr)
C	---------------------------------------------------------------
	RETURN
	END
C	---------------------------------------------------------------
