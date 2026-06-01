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
      SUBROUTINE DEINHR (NDIM,KEY,WTLENG,W,G,ERRCOF,RULPTS,SCALES,NORMS)
C***BEGIN PROLOGUE DEINHR
C***PURPOSE DEINHR computes abscissas and weights of the integration
C            rule and the null rules to be used in error estimation.
C            These are computed as functions of NDIM and KEY.
C***LAST MODIFICATION 88-04-08
C***DESCRIPTION DEINHR will for given values of NDIM and KEY compute or
C            select the correct values of the abscissas and
C            corresponding weights for different
C            integration rules and null rules and assign them to
C            G and W.
C            The heuristic error coefficients ERRCOF
C            will be computed as a function of KEY.
C            Scaling factors SCALES and normalization factors NORMS
C            used in the error estimation are computed.
C 
C 
C   ON ENTRY
C 
C     NDIM   Integer.
C            Number of variables.
C     KEY    Integer.
C            Key to selected local integration rule.
C     WTLENG Integer.
C            The number of weights in each of the rules.
C 
C   ON RETURN
C 
C     W      Real array of dimension (5,WTLENG).
C            The weights for the basic and null rules.
C            W(1,1), ...,W(1,WTLENG) are weights for the basic rule.
C            W(I,1), ...,W(I,WTLENG), for I > 1 are null rule weights.
C     G      Real array of dimension (NDIM,WTLENG).
C            The fully symmetric sum generators for the rules.
C            G(1,J),...,G(NDIM,J) are the generators for the points
C            associated with the the Jth weights.
C     ERRCOF Real array of dimension 6.
C            Heuristic error coefficients that are used in the
C            error estimation in BASRUL.
C            It is assumed that the error is computed using:
C             IF (N1*ERRCOF(1) < N2 and N2*ERRCOF(2) < N3)
C               THEN ERROR = ERRCOF(3)*N1
C               ELSE ERROR = ERRCOF(4)*MAX(N1, N2, N3)
C             ERROR = ERROR + EP*(ERRCOF(5)*ERROR/(ES+ERROR)+ERRCOF(6))
C            where N1-N3 are the null rules, EP is the error for
C            the parent
C            subregion and ES is the error for the sibling subregion.
C     RULPTS Real array of dimension WTLENG.
C            A work array containing the number of points produced by
C            each generator of the selected rule.
C     SCALES Real array of dimension (3,WTLENG).
C            Scaling factors used to construct new null rules,
C            N1, N2 and N3,
C            based on a linear combination of two successive null rules
C            in the sequence of null rules.
C     NORMS  Real array of dimension (3,WTLENG).
C            2**NDIM/(1-norm of the null rule constructed by each of
C            the scaling factors.)
C 
C***ROUTINES CALLED  D132RE,D113RE,D07HRE,D09HRE
C***END PROLOGUE DEINHR
C 
C   Global variables.
C 
      INTEGER NDIM,KEY,WTLENG
      DOUBLE PRECISION G(NDIM,WTLENG),W(5,WTLENG),ERRCOF(6)
      DOUBLE PRECISION RULPTS(WTLENG),SCALES(3,WTLENG)
      DOUBLE PRECISION NORMS(3,WTLENG)
C 
C   Local variables.
C 
      INTEGER I,J,K
      DOUBLE PRECISION WE(14)
C 
C***FIRST EXECUTABLE STATEMENT DEINHR
C 
C   Compute W, G and ERRCOF.
C 
      IF (KEY.EQ.1) THEN
         CALL D132RE (WTLENG,W,G,ERRCOF,RULPTS)
      ELSE IF (KEY.EQ.2) THEN
         CALL D113RE (WTLENG,W,G,ERRCOF,RULPTS)
      ELSE IF (KEY.EQ.3) THEN
         CALL D09HRE (NDIM,WTLENG,W,G,ERRCOF,RULPTS)
      ELSE IF (KEY.EQ.4) THEN
         CALL D07HRE (NDIM,WTLENG,W,G,ERRCOF,RULPTS)
      END IF
C 
C   Compute SCALES and NORMS.
C 
      DO 40 K=1,3
         DO 30 I=1,WTLENG
            IF (ABS(W(K+1,I)).GT.0) THEN
               SCALES(K,I)=-W(K+2,I)/W(K+1,I)
            ELSE
               SCALES(K,I)=100
            END IF
            DO 10 J=1,WTLENG
               WE(J)=W(K+2,J)+SCALES(K,I)*W(K+1,J)
 10         CONTINUE
            NORMS(K,I)=0
            DO 20 J=1,WTLENG
               NORMS(K,I)=NORMS(K,I)+RULPTS(J)*ABS(WE(J))
 20         CONTINUE
            NORMS(K,I)=2**NDIM/NORMS(K,I)
 30      CONTINUE
 40   CONTINUE
C 
C***END DEINHR
C 
      END
