	SUBROUTINE PROGRESS
	
	INCLUDE "commonblock"
C	---------------------------------------------------------------
C	AUTHOR: J. R. C. KING
C	---------------------
C	CHANGE RECORD:
C		08-01-2014:	CREATED
C	---------------------------------------------------------------
C	SUBROUTINE TO PRINT SOME THINGS TO SCREEN. HOPEFULLY THIS LOOKS
C	NICE ON A 24 LINE TERMINAL...
C	---------------------------------------------------------------
	WRITE (6,'(A24,I2)') 'KE1D running test case: ',intype
	WRITE (6,'(A43)') "==========================================="
	WRITE (6,'(A14,I3,A5,I5,A6)') 'Domain size = ',int(Rd),' with',
     +			nr,' cells'
	WRITE (6,'(A35,I1)') 'Spatial discretisation is of order ',
     +			soltype
	WRITE (6,'(A29,I1)') 'Time integration is of order ',inttype
	WRITE (6,'(A26,F3.1)') 'Coordinate system is type ',coordsno
	WRITE (6,'(A43)') "==========================================="
	WRITE (6,9) 'MAX P = ',maxval(p,pmask),'; MIN P = ',
     +				minval(p,pmask)
	WRITE (6,9) 'MAX u = ',maxval(u,pmask),'; MIN u = ',
     +				minval(u,pmask)
	WRITE (6,9) 'MAX rho = ',maxval(rho,pmask),
     +			'; MIN rho = ',minval(rho,pmask)
	WRITE (6,'(A43)') "==========================================="
	WRITE (6,'(A21,F5.2)') 'Interface location = ',Rint
	WRITE (6,'(A21,E8.2)') 'Interface pressure = ',p(q)
	WRITE (6,'(A43)') "==========================================="
	WRITE (6,'(A7,E10.4,A7,E10.4)') 'Time = ',t,'; dt = ',dt
	WRITE (6,'(A10,I6,A10,F8.4,A1)') 'time step:',tn,' Progress:',
     +			100.0*float(tn)/float(nt),'%'
	WRITE (6,10) ' '
  9	FORMAT(A8,E10.2,A10,E10.2)
  10	FORMAT(/,/,/,/,/,A1)
C	---------------------------------------------------------------
	RETURN
	END
C	---------------------------------------------------------------
