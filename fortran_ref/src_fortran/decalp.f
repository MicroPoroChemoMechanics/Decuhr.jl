C ***************************************************************************
C * All the software  contained in this library  is protected by copyright. *
C * Permission  to use, copy, modify, and  distribute this software for any *
C * purpose without fee is hereby granted, provided that this entire notice *
C * is included  in all copies  of any software which is or includes a copy *
C * or modification  of this software  and in all copies  of the supporting *
C * documentation for such software.                                        *
C ***************************************************************************
C * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED *
C * WARRANTY. IN NO EVENT, NEITHER  THE AUTHORS, NOR THE PUBLISHER, NOR ANY *
C * MEMBER  OF THE EDITORIAL BOARD OF  THE JOURNAL  "NUMERICAL ALGORITHMS", *
C * NOR ITS EDITOR-IN-CHIEF, BE  LIABLE FOR ANY ERROR  IN THE SOFTWARE, ANY *
C * MISUSE  OF IT  OR ANY DAMAGE ARISING OUT OF ITS USE. THE ENTIRE RISK OF *
C * USING THE SOFTWARE LIES WITH THE PARTY DOING SO.                        *
C ***************************************************************************
C * ANY USE  OF THE SOFTWARE  CONSTITUTES  ACCEPTANCE  OF THE TERMS  OF THE *
C * ABOVE STATEMENT.                                                        *
C ***************************************************************************
C
C Reference: T.O. Espelid and A. Genz, "DECUHR: An Algorithm for Automatic
C Integration of Singular Functions over a Hyperrectangular Region",
C Numerical Algorithms 8 (1994), pp. 201-220.
C
      SUBROUTINE DECALP ( NDIM, A, B, NUMFUN, ALPHA, SINGUL, LOGF,
     +     FUNSUB, X, FUNS ) 
C***BEGIN PROLOGUE DECALP
C***PURPOSE  DECALP computes  ALPHA and LOGF.
C***LAST MODIFICATION 94-05-04
C 
C   ON ENTRY
C 
C     NDIM   Integer.
C            Number of variables. 1 < NDIM <= MAXDIM.
C     A      Real array of dimension NDIM.
C            Lower limits of integration.
C     B      Real array of dimension NDIM.
C            Upper limits of integration.
C     NUMFUN Integer.
C            Number of components of the integral.
C     ALPHA  Real.
C            Degree of homogeneous function.
C     SINGUL Integer 
C            Dimension of the singularity.
C     LOGF   Integer.
C            LOGF = 0, then no logarithmic singularity
C            LOGF = 1, then a logarithmic singularity of order 1.
C     FUNSUB Externally declared subroutine for computing
C            all components of the integrand at the given
C            evaluation point.
C            It must have parameters (NDIM,X,NUMFUN,FUNVLS)
C            Input parameters:
C              NDIM   Integer that defines the dimension of the
C                     integral.
C              X      Real array of dimension NDIM
C                     that defines the evaluation point.
C              NUMFUN Integer that defines the number of
C                     components of I.
C     X       Real work array of length at least NDIM.
C     FUNS    Real work array of length at least NUMFUN.
C 
C   ON RETURN
C 
C     ALPHA  Real.
C            Degree of homogeneous function. The value will be estimated
C            by this routine.
C      LOGF  Integer.
C            LOGF = 0, then no logarithmic singularity
C            LOGF = 1, then a logarithmic singularity of order 1.
C            The value will be set by this routine.
C
C***END PROLOGUE DECALP
C
C   Global variables.
C 
      EXTERNAL FUNSUB
      DOUBLE PRECISION A(*), B(*), X(*), FUNS(*), ALPHA
      INTEGER NDIM, NUMFUN, LOGF, SINGUL
C
C   Local variables.
C

      INTEGER I, N, J, ITW, BDG, L, NB, INDEX, SIGNAL
      DOUBLE PRECISION H, LGTWO, TERM1, TERM2, KEST
      PARAMETER ( BDG = 10, L = 10 )
      DOUBLE PRECISION T(0:2*BDG+L), SUMX, SUMY
C
C***FIRST EXECUTABLE STATEMENT DECALP
C

      LGTWO = LOG(2D0)
      H = 2
      
C 
C   Estimate ALPHA and LOGF
C 
C   Step 1: Create the ratio table
C
C
C   Compute points along a line by halving the singular directions.
C
 10   H = H/2
      DO 20 I = 1,NDIM
         X(I) = A(I) + H*(B(I)-A(I))
 20   CONTINUE
      CALL FUNSUB ( NDIM, X, NUMFUN, FUNS )
      SUMY = 0
      DO 30 I = 1,NUMFUN
         SUMY = SUMY + ABS( FUNS(I) )
 30   CONTINUE
      IF ( SUMY .LE. 0 ) GO TO 10
      DO 60 J = 0,2*BDG+L
         SUMX = SUMY
         H = H/2
         DO 40 I = 1,SINGUL
            X(I) = A(I) + H*(B(I)-A(I))
 40      CONTINUE
         CALL FUNSUB ( NDIM, X, NUMFUN, FUNS )
         SUMY = 0
         DO 50 I = 1,NUMFUN
            SUMY = SUMY + ABS( FUNS(I) )
 50      CONTINUE
C
C     Check if we can compute the next ratio 
C     
         IF ( SUMY .LE. 0 ) GO TO 10
         T(J) = LOG( SUMX/SUMY )/LGTWO
 60   CONTINUE
C
C    Step 2: Perform linear extrapolation
C
      DO 80 I = 1,2*BDG+L
         N = 0
         DO 70 J = I-1,MAX(0,I-L),-1
            N = 2*N+1
            T(J) = T(J+1) + (T(J+1) - T(J))/N
 70      CONTINUE
 80   CONTINUE
C
C     Now T(0), T(1), ..., T(2*BDG) are all approximations to ALPHA.
C     
C     Check for logarithmic term:
C     
      IF ( ABS( T(2*BDG)-T(2*BDG-1) ) .GT. 5D-6 ) THEN
         LOGF = 1
C     
C     Step 3: Eliminate the effect of the logarithmic term on the
C     alpha estimates:
C     
C     SIGNAL is the index at which we have to stop the
C     extrapolation. This is based on the paper by Bjorstad,
C     Grosse and Dahlquist, BIT, 21, 1981. Based on our experience
C     we allow one step more than the stopping criteria indicates.
C     
         INDEX = 0
         SIGNAL = -1
         NB = 2*BDG 
 90      NB = NB - 2
         ITW = 2*BDG - NB
         IF ( NB .GT. 0) THEN
            IF ( ABS( T(NB+1)+T(NB-1)-2*T(NB) ) .GT. 0 ) THEN
               TERM1 = ( T(NB+1)-T(NB) )/( T(NB+1)+T(NB-1)-2*T(NB) )
            ELSE
               TERM1 = 0 
            END IF
            IF ( ABS( T(NB+2)+T(NB)-2*T(NB+1) ) .GT. 0) THEN
               TERM2 = ( T(NB+2)-T(NB+1) )/( T(NB+2)+T(NB)-2*T(NB+1) )
            ELSE
               TERM2 = 0
            END IF
            KEST = -1 - 1/( TERM2-TERM1 )
         ELSE
            KEST = 0
         END IF
         DO  100 J = 0,NB
            IF ( ABS( T(J+2)+T(J)-2*T(J+1) ) .GT. 0 ) THEN
               T(J) = T(J+1) - ITW*( T(J+2)-T(J+1) )*( T(J+1)-T(J) )
     +              /( ( T(J+2)+T(J)-2*T(J+1) )*( ITW-1 ) )
            ELSE
               T(J) = T(J+1)
            END IF
 100     CONTINUE
         IF (ABS(KEST-ITW+1).GT.1 .AND. SIGNAL.LT.0) SIGNAL=MAX(0,NB-2)
         IF ( NB .NE. SIGNAL .AND. NB .GT. 0 ) GO TO 90
         INDEX = SIGNAL
         ALPHA = T(INDEX)
      ELSE
         LOGF = 0
         ALPHA = T(L)
      END IF
C     
C***END DECALP
C 
      END
