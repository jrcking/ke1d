	PROGRAM KINGEULER1D
	
	INCLUDE	"commonblock"
C	---------------------------------------------------------------
C	AUTHOR: J. R. C. KING
C	---------------------
C	---------------------------------------------------------------
	INTEGER	spareINT
	DOUBLE PRECISION store1
C	---------------------------------------------------------------
C	OPEN FILES ----------------------------------------------------
	OPEN(unit=1,file="init.params",status='old')
	OPEN(unit=7,file="rho.out")
	OPEN(unit=8,file="P.out")
	OPEN(unit=9,file="u.out")
	OPEN(unit=10,file="gamma.out")
	OPEN(unit=11,file="Pc.out")
	OPEN(unit=12,file="E.out")
	OPEN(unit=13,file="A.out")
	OPEN(unit=20,file="dt.out")
	OPEN(unit=22,file="interface.out")
	OPEN(unit=23,file="pint.out")
	OPEN(unit=24,file="pbound.out")	
	OPEN(unit=25,file="ubound.out")
	OPEN(unit=26,file="rhobound.out")
C	---------------------------------------------------------------
C	GET EVERYTHING READY FOR THE TIME LOOP ------------------------
	call SETUP
C	SET TIME INDEX TO 1
	tn=1
	call FINDINTERFACE
C	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C	THIS IS THE TIME LOOP -----------------------------------------
	DO tn = 1,nt
!	DO WHILE (t.LE.0.0044) !alternative if want sod shock tube
!		tn=tn+1        ! ditto
C		-------------------------------------------------------
C		PROGRESS COUNTER EVERY NOW AND THEN -------------------
C		THIS IS WHAT PRINTS TO TERMINAL. IF IT FLICKERS -------
C		INCREASE VALUE OF spareINT. ---------------------------
		spareINT=200
		IF (mod(tn,spareINT).EQ.0) THEN
		call PROGRESS
		END IF
C		SET THE TIME STEP AND WRITE IT TO A FILE SOMEWHERE ----
		call SETTSTEP
C		WRITE RESULTS -----------------------------------------
		call DUMPRESULTS
C		-------------------------------------------------------
C		SET UP STATES FOR EACH MEDIUM ACCORDING TO GFM --------
C		Note that the rGFM is currently used once per complete
C		time step, rather than every sub-step of the TVD time
C		integration scheme. I cannot decide whether this is the
C		correct way to do things. Perhaps not, as if I use a
C		multi-step integration scheme it all goes wrong...
		call SPLITFIELDS
C		-------------------------------------------------------
C		UPDATE THE LEVEL SET ----------------------------------
		call LEVELSET
C		-------------------------------------------------------
C		FIND THE INTERFACE ------------------------------------
		call FINDINTERFACE
C		AN EULER SOLVER FOR EACH FLUID ------------------------
		DO k=1,2
C			PROPERTIES AT START OF TIME STEP --------------
			DO i=rLOW(k),rHIGH(k)
				DO j=1,5
					m1(k,j,i) = m0(k,j,i)
				END DO
			END DO
C			MULTI-STEP TIME INTEGRATION -------------------
			DO kt =1,inttype
C				BOUNDARY CONDITIONS AT LEFT OF FLUID --
				call BOUNDLEFT(k)
C				DOUNDARY CONDITIONS AT RIGHT OF FLUID -
				call BOUNDRIGHT(k)
C				RECONSTRUCTION FOR CELL EDGE PROPERTIES
				call RECONSTRUCT(k)
C				FIND THE FLUXES -----------------------
				call FLUXES(k)
C				EVOLVE!!! -----------------------------
				call EVOLVE(k)
			END DO
		END DO
C		-------------------------------------------------------
C		RECOMBINE FIELDS FROM GHOST FLUID METHOD --------------
		call RECOMBINE
C		UPDATE TIME -------------------------------------------
		t = t + dt
C		-------------------------------------------------------
	END DO
C	THAT'S THE END OF THE TIME LOOP -------------------------------
	! this is for getting sod results...
	OPEN(unit=56,file='final_rho')
	OPEN(unit=57,file='final_u')
	OPEN(unit=58,file='final_p')
	OPEN(unit=59,file='final_E')
	DO i=1,nr
		write(56,*)r(i),rho(i)
		write(57,*)r(i),u(i)
		write(58,*)r(i),p(i)
		write(59,*)r(i),E(i)
	END DO
	do i=56,59
		close(i)
	end do
C	AND THAT'S THE END OF THE PROGRAM. STOP. END. GOODBYE! --------
C	---------------------------------------------------------------
	STOP
	END
C	---------------------------------------------------------------
