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
      SUBROUTINE DECHHR (MAXDIM,NDIM,NUMFUN,A,B,ALPHA,SINGUL,LOGF,
     +     MINPTS,MAXPTS,EPSABS,EPSREL,KEY,NW,RESTAR,EMAX,WRKSUB,
     +     NUM,MAXSUB,KEYF,IFAIL,WTLENG)
C***BEGIN PROLOGUE DECHHR
C***PURPOSE  DECHHR checks the validity of input parameters to DECUHR.
C            
C***LAST MODIFICATION 94-08-02
C***DESCRIPTION
C            DECHHR computes NUM, MAXSUB, KEYF, WTLENG and
C            IFAIL as functions of the input parameters to DECUHR,
C            and checks the validity the input parameters to DECUHR.
C 
C   ON ENTRY
C 
C     MAXDIM Integer.
C            The maximum allowed number of dimensions.
C     NDIM   Integer.
C            Number of variables. 1 < NDIM <= MAXDIM.
C     NUMFUN Integer.
C            Number of components of the integral.
C     A      Real array of dimension NDIM.
C            Lower limits of integration.
C     B      Real array of dimension NDIM.
C            Upper limits of integration.
C     ALPHA  Real.
C            Degree of homogeneous function.
C     SINGUL Integer.
C            Dimension of the singularity.
C     LOGF   Integer.
C            LOGF = 0, then no logarithmic singularity
C            LOGF = 1, then a logarithmic singularity of order 1.
C     MINPTS Integer.
C            Minimum number of function evaluations.
C     MAXPTS Integer.
C            Maximum number of function evaluations.
C            The number of function evaluations over each subregion
C            is NUM.
C            If (KEY = 0 or KEY = 1) and (NDIM = 2) Then
C              NUM = 65
C            Elseif (KEY = 0 or KEY = 2) and (NDIM = 3) Then
C              NUM = 127
C            Elseif (KEY = 0 and NDIM > 3) or (KEY = 3) Then
C              NUM = 1 + 2*NDIM + 6*NDIM*NDIM +
C                    4*NDIM*(NDIM-1)*(NDIM-2)/3 + 2**NDIM
C            Elseif (KEY = 4) Then
C              NUM = 1 + 2*NDIM*(NDIM+2) + 2**NDIM
C            MAXPTS >= (SINGUL+2)*NUM and MAXPTS >= MINPTS
C     EPSABS Real.
C            Requested absolute error.
C     EPSREL Real.
C            Requested relative error.
C     KEY    Integer.
C            Key to selected local integration rule.
C            KEY = 0 is the default value.
C                  For NDIM = 2 the degree 13 rule is selected.
C                  For NDIM = 3 the degree 11 rule is selected.
C                  For NDIM > 3 the degree  9 rule is selected.
C            KEY = 1 gives the user the 2 dimensional degree 13
C                  integration rule that uses 65 evaluation points.
C            KEY = 2 gives the user the 3 dimensional degree 11
C                  integration rule that uses 127 evaluation points.
C            KEY = 3 gives the user the degree 9 integration rule.
C            KEY = 4 gives the user the degree 7 integration rule.
C                  This is the recommended rule for problems that
C                  require the most adaptivity.
C     NW     Integer.
C            Defines the length of the working array WORK.
C            Let WRKSUB denote the maximum allowed number of subregions
C            NW should then be greater or equal to
C                3+17*NUMFUN+WRKSUB*(NDIM+NUMFUN+1)*2+EMAX+
C                NUMFUN*(3*WRKSUB+9+EMAX)+(EMAX+1)**2 +3*NDIM
C     RESTAR Integer.
C            If RESTAR = 0, this is the first attempt to compute
C            the integral.
C            If RESTAR = 1, then we restart a previous attempt.
C     EMAX   Integer
C            The maximum number of extrapolation steps.
C     WRKSUB Integer.
C            Maximum size of MAXSUB.
C 
C   ON RETURN
C 
C     NUM    Integer.
C            The number of function evaluations over each subregion.
C     MAXSUB Integer.
C            The maximum possible number of subregions for the
C            given values of MAXPTS, KEY and NDIM. (See the code)
C            MAXSUB is not allowed to be greater than WRKSUB.
C            The fact that MAXSUB may be smaller allows a restart
C            option with MAXPTS increased.
C     KEYF   Integer.
C            Key to selected integration rule.
C     IFAIL  Integer.
C            IFAIL =  0 for normal exit.
C            IFAIL =  2 if KEY is less than 0 or KEY greater than 4.
C            IFAIL =  3 if NDIM is less than 2 or NDIM greater than
C                       MAXDIM.
C            IFAIL =  4 if KEY = 1 and NDIM not equal to 2.
C            IFAIL =  5 if KEY = 2 and NDIM not equal to 3.
C            IFAIL =  6 if NUMFUN less than 1.
C            IFAIL =  7 if A(I) >= B(I), for a value of 0 < I < NDIM+1.
C            IFAIL =  8 if MAXPTS is less than MINPTS.
C            IFAIL =  9 if EPSABS < 0 and EPSREL < 0.
C            IFAIL = 10 if RESTAR > 1 or RESTAR < 0.
C            IFAIL = 11 if SINGUL > NDIM or SINGUL < 1.
C            IFAIL = 12 if LOG < 0 or LOG > 1.
C            IFAIL = 13 if MAXPTS is less than (SINGUL+2)*NUM.
C            IFAIL = 14 if ALPHA <=  -SINGUL.(Set outside this routine).
C            IFAIL = 15 if EMAX < 1.
C            IFAIL = 16 if NW is too small.
C            IFAIL = 17 if WRKSUB is too small.
C     WTLENG Integer.
C            The number of generators of the chosen integration rule.
C 
C***ROUTINES CALLED-NONE
C***END PROLOGUE DECHHR
C 
C   Global variables.
C 
      INTEGER NDIM,NUMFUN,SINGUL,MINPTS,MAXPTS,KEY,NW,MAXSUB,LOGF
      INTEGER RESTAR,NUM,KEYF,IFAIL,MAXDIM,WTLENG,EMAX,WRKSUB
      DOUBLE PRECISION A(NDIM),B(NDIM),ALPHA,EPSABS,EPSREL
C 
C   Local variables.
C 
      INTEGER LIMIT,J,C1,CSING
C 
C***FIRST EXECUTABLE STATEMENT DECHHR
C 
      IFAIL=0
C 
C   Check KEY.
C 
      IF (KEY.LT.0.OR.KEY.GT.4) THEN
         IFAIL=2
         RETURN
      END IF
C 
C   Check NDIM.
C 
      IF (NDIM.LT.2.OR.NDIM.GT.MAXDIM) THEN
         IFAIL=3
         RETURN
      END IF
C 
C   For KEY = 1, NDIM must be equal to 2.
C 
      IF (KEY.EQ.1.AND.NDIM.NE.2) THEN
         IFAIL=4
         RETURN
      END IF
C 
C   For KEY = 2, NDIM must be equal to 3.
C 
      IF (KEY.EQ.2.AND.NDIM.NE.3) THEN
         IFAIL=5
         RETURN
      END IF
C 
C   For KEY = 0, we point at the selected integration rule.
C 
      IF (KEY.EQ.0) THEN
         IF (NDIM.EQ.2) THEN
            KEYF=1
         ELSE IF (NDIM.EQ.3) THEN
            KEYF=2
         ELSE
            KEYF=3
         END IF
      ELSE
         KEYF = KEY
      END IF
C 
C   Compute NUM and WTLENG as a function of KEYF and NDIM.
C 
      IF (KEYF.EQ.1) THEN
         NUM = 65
         WTLENG = 14
      ELSE IF (KEYF.EQ.2) THEN
         NUM = 127
         WTLENG = 13
      ELSE IF (KEYF.EQ.3) THEN
         NUM = 1 + 2*NDIM + 6*NDIM*NDIM 
     +        + 4*NDIM*(NDIM-1)*(NDIM-2)/3 + 2**NDIM
         WTLENG = 9
         IF (NDIM.EQ.2) WTLENG = 8
      ELSE IF (KEYF.EQ.4) THEN
         NUM = 1 + 2*NDIM*(NDIM+2) + 2**NDIM
         WTLENG = 6
      END IF
C 
C   Check NUMFUN.
C 
      IF (NUMFUN.LT.1) THEN
         IFAIL=6
         RETURN
      END IF
C 
C   Check upper and lower limits of integration.
C 
      DO 10 J=1,NDIM
         IF (A(J).GE.B(J)) THEN
            IFAIL=7
            RETURN
         END IF
 10   CONTINUE
C 
C   Check MAXPTS >= MINPTS.
C 
      IF (MAXPTS.LT.MINPTS) THEN
         IFAIL=8
         RETURN
      END IF
C 
C   Check accuracy requests.
C 
      IF (EPSABS.LT.0.AND.EPSREL.LT.0) THEN
         IFAIL=9
         RETURN
      END IF
C 
C    Check RESTAR.
C 
      IF (RESTAR.NE.0.AND.RESTAR.NE.1) THEN
         IFAIL=10
         RETURN
      END IF
C 
C    Check SINGUL.
C 
      IF (SINGUL.LT.1.OR.SINGUL.GT.NDIM) THEN
         IFAIL=11
         RETURN
      END IF
C 
C    Check LOGF.
C 
      IF (LOGF.LT.0.OR.LOGF.GT.1) THEN
         IFAIL=12
         RETURN
      END IF
C 
C   Check that MAXPTS allows at least one subdivision of the
C   singular region in SINGUL + 1 new pieces.
C 
      IF (MAXPTS.LT.(SINGUL+2)*NUM) THEN
         IFAIL=13
         RETURN
      END IF
C 
C    Check EMAX.
C 
      IF (EMAX.LT.1) THEN
         IFAIL=15
         RETURN
      END IF
C 
C   Compute MAXSUB. This is a worst case computation. We assume that
C   the singular region will be cut CSING times, each time in
C   SINGUL + 1 new subregions. CSING is the maximum number of times
C   that this can happen, given MAXPTS and NUM. Then we assume that only
C   bisections (of non-singular regions) will take place. C1 is the
C   maximum number of times this can happen, given MAXPTS, NUM and CSING
C 
      CSING = (MAXPTS-NUM)/((SINGUL+1)*NUM)
      C1 = (MAXPTS-NUM*(1+(SINGUL+1)*CSING))/(2*NUM)
      IF ( C1.EQ.0 ) THEN
         C1 = -1
      END IF
      MAXSUB = 2 + CSING*SINGUL + C1
C 
C   Check size of double precision workspace.
C 
      LIMIT= 3 + 17*NUMFUN + WRKSUB*(NDIM+NUMFUN+1)*2 + EMAX 
     +     + NUMFUN*(3*WRKSUB+9+EMAX) + (EMAX+1)**2 + 3*NDIM
      IF ( NW.LT.LIMIT ) THEN
         IFAIL = 16
         RETURN
      END IF
C 
C   Check if MAXSUB can be achieved.
C 
      IF (MAXSUB.GT.WRKSUB) THEN
         IFAIL=17
         RETURN
      END IF
C 
C***END DECHHR
C 
      END
