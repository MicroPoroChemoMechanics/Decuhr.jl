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
      SUBROUTINE DEADHR (NDIM,NUMFUN,A,B,MAXSUB,FUNSUB,SINGUL,ALPHA,
     +LOGF,EPSABS,EPSREL,KEY,RESTAR,NUM,LENW,WTLENG,EMAX,MINPTS,MAXPTS,
     +NSUB,RESULT,ABSERR,NEVAL,IFAIL,VALUES,ERRORS,CENTRS,HWIDTS,GREATE,
     +DIR,WORK,Q,U,E,QOLD,QNEW,UNEW,UERR,EXTERR,NE,BETA,T,DIFF,G,W,
     +RULPTS,CENTER,HWIDTH,X,SCALES,NORMS,LIST,UPOINT,ORDER)
C***BEGIN PROLOGUE DEADHR
C***KEYWORDS automatic multidimensional integrator,
C            singular integrands,
C            n-dimensional hyper-rectangles,
C            general purpose, global adaptive with extrapolation.
C***PURPOSE  The routine calculates an approximation to a given
C            vector of definite integrals I, over a hyper-rectangular
C            region hopefully satisfying for each component of I the
C            following claim for accuracy:
C            ABS(I(K)-RESULT(K)).LE.MAX(EPSABS,EPSREL*ABS(I(K)))
C***LAST MODIFICATION 92-09-16
C***DESCRIPTION Computation of integrals over hyper-rectangular
C            regions with singularities.
C            DECUHR is a driver for the integration routine
C            DEADHR, which repeatedly subdivides the region
C            of integration. DEADHR estimates the integrals and the
C            errors over the subregions and subdivides the subregion
C            with greatest estimated errors until the error request
C            is met or MAXSUB subregions are stored.
C            The subdivision is done in such a way that we at any
C            stage have only one region containing the singularity.
C            When the singular region is picked it is subdivided in
C            SINGUL + 1 new regions (cutting each of the SINGUL
C            first directions once), creating one new singular region
C            and SINGUL non-singular regions.
C            The non-singular regions are divided in two equally
C            sized parts along the direction with greatest absolute
C            fourth order difference.
C 
C   ON ENTRY
C 
C     NDIM   Integer.
C            Number of variables. 1 < NDIM <= MAXDIM.
C     NUMFUN Integer.
C            Number of components of the integral.
C     A      Real array of dimension NDIM.
C            Lower limits of integration.
C     B      Real array of dimension NDIM.
C            Upper limits of integration.
C     MAXSUB Integer.
C            The computations proceed until there are at most
C            MAXSUB subregions in the data structure.
C 
C     FUNSUB Externally declared subroutine for computing
C            all components of the integrand in the given
C            evaluation point.
C            It must have parameters (NDIM,X,NUMFUN,FUNVLS)
C            Input parameters:
C              NDIM   Integer that defines the dimension of the
C                     integral.
C              X      Real array of dimension NDIM
C                     that defines the evaluation point.
C              NUMFUN Integer that defines the number of
C                     components of I.
C            Output parameter:
C              FUNVLS Real array of dimension NUMFUN
C                     that defines NUMFUN components of the integrand.
C 
C     SINGUL Integer.
C            Dimension of the singularity
C     ALPHA  Real
C            Degree of homogeneous function. ALPHA > -SINGUL.
C     LOGF   Integer
C            Indicating power of logarithmic term in all components.
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
C                  require great adaptivity.
C     RESTAR Integer.
C            If RESTAR = 0, this is the first attempt to compute
C            the integral.
C            If RESTAR = 1, then we restart a previous attempt.
C            (In this case the output parameters from DEADHR
C            must not be changed since the last
C            exit from DEADHR.)
C     NUM    Integer.
C            The number of function evaluations over each subregion.
C     LENW   Integer.
C            Defines the length of the working array WORK.
C            LENW should be greater or equal to
C            16*NUMFUN.
C     WTLENG Integer.
C            The number of weights in the basic integration rule.
C     EMAX
C            The maximum number of extrapolation steps.
C     MINPTS Integer.
C            Minimum number of function evaluations.
C     MAXPTS
C            Maximum number of function evaluations.
C     NSUB   Integer.
C            If RESTAR = 1, then NSUB must specify the number
C            of subregions stored in the previous call to DEADHR.
C 
C   ON RETURN
C 
C     NSUB   Integer.
C            Number of stored subregions.
C     RESULT Real array of dimension NUMFUN.
C            Approximations to all components of the integral.
C     ABSERR Real array of dimension NUMFUN.
C            Estimates of absolute accuracies.
C     NEVAL  Integer.
C            Number of function evaluations used by DEADHR.
C     IFAIL  Integer.
C            IFAIL = 0 for normal exit, when ABSERR(K) <=  EPSABS or
C              ABSERR(K) <=  ABS(RESULT(K))*EPSREL with MAXSUB or less
C              subregions processed for all values of K,
C              1 <=  K <=  NUMFUN.
C            IFAIL = 1 if MAXSUB was too small for DEADHR
C              to obtain the required accuracy. In this case DEADHR
C              returns values of RESULT with estimated absolute
C              accuracies ABSERR.
C     VALUES Real array of dimension (NUMFUN,MAXSUB+1).
C            Used to store estimated values of the integrals
C            over the subregions.
C     ERRORS Real array of dimension (NUMFUN,MAXSUB+1).
C            Used to store the corresponding estimated errors.
C     CENTRS Real array of dimension (NDIM,MAXSUB+1).
C            Used to store the centers of the stored subregions.
C     HWIDTS Real array of dimension (NDIM,MAXSUB+1).
C            Used to store the half widths of the stored subregions.
C     GREATE Real array of dimension MAXSUB+1.
C            Used to store the greatest estimated errors in
C            all subregions.
C     DIR    Real array of dimension MAXSUB+1.
C            DIR is used to store the directions for
C            further subdivision.
C     WORK   Real array of dimension LENW.
C            Used  in DERLHR and DETRHR.
C     Q      Real array of dimension (NUMFUN,0:MAXSUB)
C            contains the estimates over the singular regions
C            (Tail-estimates).
C     U      Real array of dimension (NUMFUN,MAXSUB)
C            contains the terms in series.
C     E      Real array of dimension (NUMFUN,MAXSUB)
C            contains the estimated errors in each U-term.
C     QOLD   Real array of dimension NUMFUN.
C            The last tail corrections for all functions in the vector.
C     QNEW   Real array of dimension NUMFUN.
C            The new tail correction for all functions in the vector.
C     UNEW   Real array of dimension NUMFUN.
C            This gives the next terms in the series (new extrapolation
C            step) else it is the correction to the u-values (updating).
C     UERR   Real array of dimension NUMFUN.
C            The estimated errors of all U-terms in the series.
C     EXTERR Real array of dimension NUMFUN.
C            These errors are associated with the singular region and
C            they are the pure extrapolation errors.
C     NE     Real array of dimension (0:EMAX)
C            Dummy parameter
C     BETA   Real array of dimension (0:EMAX,0:EMAX)
C            Dummy parameter
C     T      Real array of dimension (NUMFUN,0:EMAX)
C            contains the last rows in the extrapolation tableaus.
C     DIFF   Real array of dimension NDIM.
C            Work array.
C     G      Real array of dimension (NDIM,WTLENG).
C            The fully symmetric sum generators for the rules.
C            G(1,J),...,G(NDIM,J) are the generators for the
C            points associated with the Jth weights.
C     W      Real array of dimension (5,WTLENG).
C            The weights for the basic and null rules.
C            W(1,1), ..., W(1,WTLENG) are weights for the basic rule.
C            W(I,1), ..., W(I,WTLENG) , for I > 1 are null rule weights.
C     RULPTS Real array of dimension WTLENG.
C            Work array used in DEINHR.
C     CENTER Real array of dimension NDIM.
C            Work array used in DETRHR.
C     HWIDTH Real array of dimension NDIM.
C            Work array used in DETRHR.
C     X      Real array of dimension NDIM.
C            Work array used in DERLHR.
C     SCALES Real array of dimension (3,WTLENG).
C            Work array used by DEINHR and DERLHR.
C     NORMS  Real array of dimension (3,WTLENG).
C            Work array used by DEINHR and DERLHR.
C     LIST   Integer array used in DETRHR of dimension MAXSUB.
C            Is a partially sorted list, where LIST(1) is the top
C            element in a heap of subregions.
C     UPOINT Integer array of dimension MAXSUB
C            Is an array of pointers to where in the U-sequence
C            a region belongs. This information is used when updating
C            the corresponding U-term after a subdivision.
C     ORDER  Integer array of dimension NDIM.
C            Work array used to give the order in which the singular
C            region will be cut.
C***REFERENCES
C 
C   T.O.Espelid and A.Genz, DECUHR: An Algorithm for Automatic 
C   Integration of Singular Functions over a Hyperrectangular Region.
C   Numerical Algorithms 8(1994), PP. 201-220.
C 
C   T.O.Espelid, On integrating Vertex Singularities using
C   Extrapolation,  BIT 34,1(1994) 62-79.
C 
C   T.O.Espelid, On integrating Singularities using non-Uniform
C   Subdivision and Extrapolation,  in Numerical Integration IV, eds.
C   Brass and Hammerlin, Birhauser, ISNM Vol. 112(1993) 77-89.
C 
C   P. van Dooren and L. de Ridder, Algorithm 6, An adaptive algorithm
C   for numerical integration over an n-dimensional cube, J.Comput.Appl.
C   Math. 2(1976)207-217.
C 
C   A.C.Genz and A.A.Malik, Algorithm 019. Remarks on algorithm 006:
C   An adaptive algorithm for numerical integration over an
C   N-dimensional rectangular region, J.Comput.Appl.Math.6(1980)295-302.
C 
C***ROUTINES CALLED DEINHR,DERLHR,DETRHR,DESBHR,DEXTHR.
C***END PROLOGUE DEADHR
C 
C   Global variables.
C 
      EXTERNAL FUNSUB
      INTEGER NDIM,NUMFUN,MAXSUB,KEY,LENW,RESTAR,MINPTS,MAXPTS
      INTEGER NUM,NEVAL,NSUB,IFAIL,WTLENG,SINGUL
      INTEGER UPOINT(MAXSUB),LIST(MAXSUB),LOGF,EMAX,ORDER(NDIM)
      DOUBLE PRECISION A(NDIM),B(NDIM),EPSABS,EPSREL,ALPHA
      DOUBLE PRECISION RESULT(NUMFUN),ABSERR(NUMFUN),DIFF(NDIM)
      DOUBLE PRECISION VALUES(NUMFUN,0:MAXSUB),ERRORS(NUMFUN,0:MAXSUB)
      DOUBLE PRECISION CENTRS(NDIM,0:MAXSUB)
      DOUBLE PRECISION HWIDTS(NDIM,0:MAXSUB),T(NUMFUN,0:EMAX)
      DOUBLE PRECISION GREATE(0:MAXSUB),DIR(0:MAXSUB)
      DOUBLE PRECISION BETA(0:EMAX,0:EMAX)
      DOUBLE PRECISION WORK(LENW),RULPTS(WTLENG),EXTERR(NUMFUN)
      DOUBLE PRECISION G(NDIM,WTLENG),W(5,WTLENG),NE(EMAX)
      DOUBLE PRECISION CENTER(NDIM),HWIDTH(NDIM),X(NDIM)
      DOUBLE PRECISION SCALES(3,WTLENG),NORMS(3,WTLENG),QOLD(NUMFUN)
      DOUBLE PRECISION U(NUMFUN,MAXSUB),E(NUMFUN,MAXSUB),QNEW(NUMFUN)
      DOUBLE PRECISION Q(NUMFUN,0:MAXSUB),UNEW(NUMFUN),UERR(NUMFUN)
C 
C   Local variables.
C 
C   SBRGNS is the number of stored subregions.
C   NDIV   The number of subregions to be divided in each main step.
C   POINTR Pointer to the position in the data structure where
C          the new subregions are to be stored.
C   TOP    is a pointer to the top element in the heap of subregions.
C   VACANT Pointer to the vacant element.
C   DIRECT Direction of subdivision.
C   ERRCOF Heuristic error coefficient defined in DEINHR and used by
C          DERLHR and DEADHR.
C   UPDATE Pointer: either to a new element or to an old one to be updat
C          Value = 0 used to indicate a new extrapolation step.
C   FIRST  Logical: .TRUE. indicates that this is the first time
C          DEXTHR is called.
C   EXSTEP The exponent sequence's step size in the error expansion.
C          Default value = 1. Can be set in a parameter statement to a
C          different integer > 1 if that is appropriate. Advantage:
C          avoid to eliminate terms that are not present for the given
C          function. Warning: If the step size is set wrong (too big)
C          then this may result in poor performance.
C 
      LOGICAL FIRST
      INTEGER I,J,K,SBRGNS,TOP,UPDATE,EXSTEP,NUMU
      INTEGER NDIV,POINTR,DIRECT,INDEX,VACANT
      DOUBLE PRECISION ERRCOF(6)
      PARAMETER (EXSTEP=1)
C 
C***FIRST EXECUTABLE STATEMENT DEADHR
C 
C   Call DEINHR to compute the weights and abscissas of
C   the function evaluation points.
C 
      CALL DEINHR (NDIM,KEY,WTLENG,W,G,ERRCOF,RULPTS,SCALES,NORMS)
C 
C   If RESTAR = 1, then this is a restart run.
C 
      IF (RESTAR.EQ.1) THEN
         SBRGNS=NSUB
         GO TO 50
      END IF
C 
C   Initialize the SBRGNS, CENTRS,  HWIDTS and ORDER
C   Note the singular region will stay in position 0 and will not be
C   a part of the heap at any time.
C 
      SBRGNS=0
      DO 10 J=1,NDIM
         CENTRS(J,0)=(A(J)+B(J))/2
         HWIDTS(J,0)=ABS(B(J)-A(J))/2
         ORDER(J)=J
 10   CONTINUE
C 
C   Apply DERLHR over the whole region.
C   DERLHR will only change DIR(0) if the input value is non-negative.
C 
      DIR(0)=-SINGUL
      CALL DERLHR (NDIM,CENTRS(1,0),HWIDTS(1,0),WTLENG,G,W,ERRCOF,
     +     NUMFUN,FUNSUB,SCALES,NORMS,X,WORK,VALUES(1,0),ERRORS(1,0),
     +     DIR(0),GREATE(0),DIFF,ORDER)
C 
C   Initialize RESULT, ABSERR, Q(*,0), NEVAL, NUMU and FIRST.
C 
      DO 20 J=1,NUMFUN
         RESULT(J)=VALUES(J,0)
         ABSERR(J)=ERRORS(J,0)
         Q(J,0)=VALUES(J,0)
 20   CONTINUE
      NEVAL=NUM
      FIRST=.TRUE.
      NUMU=0
C 
C   Initialize E and U
C 
      DO 40 J=1,NUMFUN
         DO 30 I=1,MAXSUB
            E(J,I)=0
            U(J,I)=0
 30      CONTINUE
 40   CONTINUE
C 
C   Initialize the pointer LIST(1) to point on the singular region.
C 
      LIST(1)=0
C 
C***End initialization.
C 
C***Begin loop while the error is too large,
C 
C     First we determine if we will subdivide the singular region
C     or a regular region, and then the number, NDIV, of new subregions.
C 
 50   TOP=LIST(1)
C 
C     First we determine if we will subdivide the singular region or a
C     non-singular region, and then the number, NDIV, of new subregions.
C 
      IF (GREATE(TOP).LT.(GREATE(0))) THEN
         VACANT=0
      ELSE
         VACANT=TOP
      END IF
      DIRECT=DIR(VACANT)
      NDIV=MAX(2,-DIRECT+1)
C 
C     Check if NEVAL+NDIV*NUM is less than or equal to MAXPTS:
C     MAXPTS is the maximum number of function evaluations that are
C     allowed to be computed.
C 
      IF (NEVAL+NDIV*NUM.LE.MAXPTS) THEN
C 
C     We are allowed to divide further: prepare to remove the region
C     with greatest error estimate from the heap or replace the singular
C     region with a smaller one.
C 
C     Let POINTR point to the first free position in the data structure.
C 
         POINTR=SBRGNS+1
C 
C     Adjust, if necessary, the heap. (Reduce the size by one element.
C     Note: the information about the region we will subdivide
C     remains in the vacant position.)
C 
         IF (VACANT.GT.0) THEN
            CALL DETRHR (1,SBRGNS,GREATE,LIST,K)
         END IF
C 
C 
C     Determine the new subregions.
C 
         CALL DESBHR (NDIM,CENTRS,HWIDTS,DIR,DIRECT,POINTR,VACANT,ORDER)
C 
C     Determine if this is a new extrapolation step or an update.
C     UPDATE will point to the element in the
C     U-series to be corrected or created.
C 
         IF (VACANT.EQ.0) THEN
            NUMU=NUMU+1
            UPDATE=NUMU
         ELSE
            UPDATE=UPOINT(VACANT)
         END IF
C 
C     Apply the basic rule to the new regions and let each new region
C     (except the singular region) point to its U-element.
C     First the region in the VACANT position. This region is the only
C     one that may be the singular and thus give a new Q-element.
C 
         INDEX=VACANT
         IF (VACANT.EQ.0) THEN
            DO 60 J=1,NUMFUN
               UERR(J)=0
               UNEW(J)=0
 60         CONTINUE
         ELSE
            DO 70 J=1,NUMFUN
               UERR(J)=-ERRORS(J,INDEX)
               UNEW(J)=-VALUES(J,INDEX)
 70         CONTINUE
         END IF
         CALL DERLHR (NDIM,CENTRS(1,INDEX),HWIDTS(1,INDEX),WTLENG,G,W,
     +    ERRCOF,NUMFUN,FUNSUB,SCALES,NORMS,X,WORK,VALUES(1,INDEX),
     +    ERRORS(1,INDEX),DIR(INDEX),GREATE(INDEX),DIFF,ORDER)
         IF (VACANT.EQ.0) THEN
            DO 80 J=1,NUMFUN
               Q(J,NUMU)=VALUES(J,0)
 80         CONTINUE
         ELSE
            UPOINT(INDEX)=UPDATE
            DO 90 J=1,NUMFUN
               UERR(J)=UERR(J)+ERRORS(J,INDEX)
               UNEW(J)=UNEW(J)+VALUES(J,INDEX)
 90         CONTINUE
         END IF
C 
C     Then the rest of the regions
C 
         DO 110 I=2,NDIV
            INDEX=POINTR+I-2
            CALL DERLHR (NDIM,CENTRS(1,INDEX),HWIDTS(1,INDEX),WTLENG,G,
     +       W,ERRCOF,NUMFUN,FUNSUB,SCALES,NORMS,X,WORK,VALUES(1,INDEX),
     +       ERRORS(1,INDEX),DIR(INDEX),GREATE(INDEX),DIFF,ORDER)
            UPOINT(INDEX)=UPDATE
            DO 100 J=1,NUMFUN
               UERR(J)=UERR(J)+ERRORS(J,INDEX)
               UNEW(J)=UNEW(J)+VALUES(J,INDEX)
 100        CONTINUE
 110     CONTINUE
C 
C     Compute the E and U terms (These may be new terms or terms that
C     have to be updated), QNEW and QOLD have no influence in a true
C     update.
C 
         DO 120 J=1,NUMFUN
            QNEW(J)=Q(J,NUMU)
            QOLD(J)=Q(J,NUMU-1)
            U(J,UPDATE)=U(J,UPDATE)+UNEW(J)
            E(J,UPDATE)=E(J,UPDATE)+UERR(J)
 120     CONTINUE
C 
C     Do the extrapolation and compute the global results and errors.
C     UPDATE is used to signal an extrapolation step.
C 
         IF (VACANT.EQ.0) THEN
            UPDATE=0
         END IF
         CALL DEXTHR (NUMFUN,ALPHA,LOGF,SINGUL,EXSTEP,NUMU,T,UPDATE,
     +    UNEW,QNEW,QOLD,FIRST,EMAX,E,RESULT,ABSERR,EXTERR,NE,BETA)
         FIRST=.FALSE.
         NEVAL=NEVAL+NDIV*NUM
C 
C     Change the error estimates in position 0.
C     We define the pure extrapolation error to be the new error of
C     the singular region.
C 
         GREATE(0)=0
         DO 130 J=1,NUMFUN
            GREATE(0)=MAX(GREATE(0),EXTERR(J))
            ERRORS(J,0)=EXTERR(J)
 130     CONTINUE
C 
C     Store results in heap.
C 
         IF (VACANT.GT.0) THEN
            CALL DETRHR (2,POINTR-1,GREATE,LIST,VACANT)
         END IF
         DO 140 I=2,NDIV
            INDEX=POINTR+I-2
            CALL DETRHR (2,INDEX,GREATE,LIST,INDEX)
 140     CONTINUE
         SBRGNS=POINTR+NDIV-2
C 
C     Check for termination.
C 
         IF (NEVAL.LT.MINPTS) THEN
            GO TO 50
         END IF
         DO 150 J=1,NUMFUN
            IF (ABSERR(J).GT.EPSREL*ABS(RESULT(J)).AND.ABSERR(J).GT.
     +       EPSABS) THEN
               GO TO 50
            END IF
 150     CONTINUE
         IFAIL=0
C 
C     Else we did not succeed with the
C     given value of MAXSUB.
C 
      ELSE
         IFAIL=1
      END IF
      NSUB=SBRGNS
C 
C***END DEADHR
C 
      END
