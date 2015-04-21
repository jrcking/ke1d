	SUBROUTINE SETUP

	INCLUDE "commonblock"
C	---------------------------------------------------------------
C	AUTHOR: J. R. C. KING
C	---------------------
C	CHANGE RECORD:
C		20-11-2012:	CREATED
C		05-08-2013:	A MYRIAD OF CHANGES, UNRECORDED
C	---------------------------------------------------------------
C	SUBROUTINE TO READ THE INPUT FILES FOR PROGRAM KINGEULER1D
C	AND TO SET THE INITIAL FIELDS. THIS SUBROUTINE IS A BIT MESSY.
C	---------------------------------------------------------------
C	---------------------------------------------------------------
	DOUBLE PRECISION rmin
	CHARACTER spare
C	READ THE MAIN PARAMETER FILE ----------------------------------
	READ (1,*) spare
	READ (1,*) nr
	READ (1,*) nrout
	READ (1,*) nt
	READ (1,*) outfreq
	READ (1,*) Rd
	READ (1,*) CFL
	READ (1,*) soltype
	READ (1,*) inttype
	READ (1,*) obtype
	READ (1,*) intype
	READ (1,*) coordsno
C	---------------------------------------------------------------
C	HOW MANY CELLS FOR BOUNDARIES? 5 IS BE SAFE
C	NEEDS TO BE >= 1 + ORDER OF SCHEME
	gfmbuffer = 1 + soltype
C	---------------------------------------------------------------
C	SET THE COEFFICIENTS FOR TIME INTEGRATION ---------------------
	IF (inttype .EQ. 1) THEN
C		1st ORDER EULER
		sfrac(1,1) = 1.0
		sfrac(2,1) = 0.0
		sfrac(3,1) = 1.0
	ELSE IF (inttype .EQ. 2) THEN
C		2nd ORDER TVD (thanks to Shu & Osher)
		sfrac(1,1) = 1.0
		sfrac(2,1) = 0.0
		sfrac(3,1) = 1.0
		sfrac(1,2) = 0.5
		sfrac(2,2) = 0.5
		sfrac(3,2) = 0.5
	ELSE IF (inttype .EQ. 3) THEN
C		3rd ORDER TVD (thanks to Shu & Osher)
		sfrac(1,1) = 1.0
		sfrac(2,1) = 0.0
		sfrac(3,1) = 1.0
		sfrac(1,2) = 0.75
		sfrac(2,2) = 0.25
		sfrac(3,2) = 0.25
		sfrac(1,3) = 1.0/3.0
		sfrac(2,3) = 2.0/3.0
		sfrac(3,3) = 2.0/3.0
	END IF
C	SET THE INITIAL FIELDS ----------------------------------------
C	WANT TO MODEL SOME KIND OF SHOCK PROBLEM - WHAT TYPE ----------
	OPEN(unit=51,file="SOD.dat",status='old') !Sod shock tube
	OPEN(unit=52,file="AG.dat",status='old') !Air gun bubble
	OPEN(unit=53,file="UE.dat",status='old') !underwater explosion
	OPEN(unit=54,file="SEDOV.dat",status='old') !sedoiv explosion
C	---------------------------------------------------------------
	READ (intype,*) spare
	READ (intype,*) rho0
	READ (intype,*) u0
	READ (intype,*) p0
	READ (intype,*) gamma0
	READ (intype,*) pc0
	READ (intype,*) rho1
	READ (intype,*) u1
	READ (intype,*) p1
	READ (intype,*) gamma1
	READ (intype,*) pc1
	READ (intype,*) Rint
	Pinf = p1
	gm(1) = gamma0
	gm(2) = gamma1
	pcr(1) = pc0
	pcr(2) = pc1
C	CREATE THE MESH
	call MESHGEN
C	===============================================================
C	OUTPUTTING AND GIVING AN INITIAL VALUE TO INTERFACE INDEX -----
	ofc1 = outfreq + 1
	t=0.0
	DO i=1,nr-1
		IF (r(i) .LT. Rint .AND. r(i+1) .GE. Rint) THEN
			q = i
		END IF
	END DO
C	SET UP INITIAL FIELDS -----------------------------------------
	DO i=1,q
		rho(i) = rho0
		u(i) = u0
		P(i) = p0
		gamm(i) = gamma0
		Pc(i) = pc0
		E(i) = (P(i)+Pc(i)*gamm(i))/(gamm(i)-1)
     +			 + 0.5*rho(i)*u(i)**2
		alpha(i) = r(i) - Rint
	END DO
	DO i=q+1,nr+inttype-1
		rho(i) = rho1
		u(i) = u1
		P(i) = p1
		gamm(i) = gamma1
		Pc(i) = pc1
		E(i) = (P(i)+Pc(i)*gamm(i))/(gamm(i)-1)
     +			 + 0.5*rho(i)*u(i)**2
		alpha(i) = r(i) - Rint
	END DO
C	---------------------------------------------------------------
C	===============================================================
C	CREATE AN ARRAY CALLED pmask FOR SOME HISTORIC (FORGOTTEN)
C	REASON... COULD BE FOR SETTING CFL..
	DO i=1,nr-1
		pmask(i) = .TRUE.
	END DO
	DO i=nr+1,nrmax
		pmask(i) = .FALSE.
	END DO
C	===============================================================
	RETURN
	END
C	---------------------------------------------------------------
C	000000000000000000000000000000000000000000000000000000000000000
	SUBROUTINE MESHGEN

	INCLUDE "commonblock"
C	---------------------------------------------------------------
C	CREATE THE MESH BASED ON Rd and nr
C	---------------------------------------------------------------
C	SOMETIMES PRESSURE WAVES CONVERGING ON ORIGIN ARE TOO STRONG AT
C	r=0, IN WHICH CASE USE rmin TO SHIFT BOUNDARY OUT A BIT (CHEAT)
	rmin = 0.0
	dr = (Rd)/float(nr-1)
	DO i=1,nr+inttype-1
		r(i) = rmin + dr*float(i-1) + 0.5*dr
	END DO
C	---------------------------------------------------------------
	RETURN
	END
C	000000000000000000000000000000000000000000000000000000000000000
