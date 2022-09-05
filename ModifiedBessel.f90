!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%** PROGRAM     MODIFIED BESSEL FUNCTION OF THE FIRST KIND OF ORDER ZERO **%%
!%%** AUTHOR      ALPIXELS                                                 **%%
!%%** LICENSE     LGPL-V3                                                  **%%
!%%**                                                                      **%%
!%%** OBS         FUNCTION COMPUTES THE FUNCTION FOR ANY X VALUE AND IS    **%%
!%%**             BASED ON POLYNOMIAL COEFFCIENTS GIVEN BY ABRAMOWITZ AND  **%%
!%%**             STEGUN FOLLOWING THE IMPLEMENTATION ON FORTRAN NUMERICAL **%%
!%%**             RECIPIES BUT ITS CONSTRUCTED IN FULL DOUBLE PRECISION    **%%
!%%**                                                                      **%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PROGRAM M_BESSEL
 IMPLICIT NONE
 INTEGER, PARAMETER:: D = KIND(1.0D0)
 INTEGER:: I,NX
 REAL(D):: X, BESSI0
 REAL(D):: IX,FX,DX

 IX=0.0D0
 FX=4.0D0
 DX=0.01D0

 NX=INT((FX-IX)/DX) + 1
 X=IX

 DO I=1,NX
    PRINT*,X,BESSI0(X)
    X=X + DX
 ENDDO

 STOP
END PROGRAM M_BESSEL
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   FUNCTIONS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FUNCTION BESSI0(X)
 INTEGER, PARAMETER:: D = KIND(1.0D0)
 REAL(D):: AX,Y,X
 REAL(D):: BESSI0
 REAL(D):: P1,P2,P3,P4,P5,P6,P7
 REAL(D):: Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9
 SAVE:: P1,P2,P3,P4,P5,P6,P7
 SAVE:: Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9

 DATA P1,P2/1.0D0,3.5156229D0/
 DATA P3,P4/3.0899424D0,1.2067492D0/
 DATA P5,P6/0.2659732D0,0.360768D-1/
 DATA P7/0.45813D-2/
 DATA Q1,Q2,Q3/0.39894228D0,0.1328592D-1,0.225319D-2/
 DATA Q4,Q5,Q6/-0.157565D-2,0.916281D-2,-0.2057706D-1/
 DATA Q7,Q8,Q9/0.2635537D-1,-0.1647633D-1,0.392377D-2/

 IF(ABS(X) .LT. 3.75D0)THEN
   Y=(X/3.75D0)**2
   BESSI0=P1 + Y*(P2 + Y*(P3 + Y*(P4 + Y*(P5 + Y*(P6 + Y*P7)))))
 ELSE
   AX=DABS(X)
   Y=3.75D0/AX
   BESSI0=Q1 + Y*(Q2 + Y*(Q3 + Y*(Q4 + Y*(Q5 + Y*(Q6 + Y*(Q7 + Y*(Q8 + Y*Q9)))))))
   BESSI0=BESSI0*(DEXP(AX)/DSQRT(AX))
 ENDIF

 RETURN
END 
