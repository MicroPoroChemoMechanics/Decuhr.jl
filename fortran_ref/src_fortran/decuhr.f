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
      SUBROUTINE DECUHR ( NDIM, NUMFUN, A, B, MINPTS, MAXPTS, FUNSUB,
     +     SINGUL, ALPHA, LOGF, EPSABS, EPSREL, KEY, WRKSUB, NW, RESTAR,
     +     EMAX, RESULT, ABSERR, NEVAL, IFAIL, WORK, IWORK)
C***BEGIN PROLOGUE DECUHR
C***DATE WRITTEN   930816   (YYMMDD)
C***REVISION DATE  940802   (YYMMDD)
C***CATEGORY NO. H2B1A1
C***AUTHOR
C            Terje O. Espelid, Department of Informatics,
C            University of Bergen,  Hoyteknologisenteret,
C            N-5020 Bergen, Norway
C            Email..  terje@ii.uib.no
C 
C            Alan Genz,
C            Department of Mathematics
C            Washington State University
C            Pullman, WA 99164-3113, USA
C            Email..  genz@gauss.math.wsu.edu
C 
C***KEYWORDS automatic multidimensional integrator,
C            singular integrands,
C            n-dimensional hyper-rectangles,
C            general purpose, global adaptive with extrapolation.
C***PURPOSE  The routine calculates an approximation to a given
C            vector of definite integrals 
C 
C      B(1) B(2)     B(NDIM)
C     I    I    ... I       (F ,F ,...,F      ) DX(NDIM)...DX(2)DX(1),
C      A(1) A(2)     A(NDIM)  1  2      NUMFUN
C 
C       where F = F (X ,X ,...,X    ), I = 1,2,...,NUMFUN,
C              I   I  1  2      NDIM
C 
C            hopefully satisfying for each component of I the following
C            claim for accuracy:
C            ABS(I(K)-RESULT(K)).LE.MAX(EPSABS,EPSREL*ABS(I(K)))
C            The vector of definite integrals are all assumed to have
C            a singularity of dimension SINGUL at the point,
C            X(1) = A(1), X(2) = A(2), ... , X(SINGUL) = A(SINGUL).
C            The singularity is assumed to be caused by a homogeneous
C            function of degree ALPHA at the point. A logarithmic
C            singularity at the same point can also be handled.
C 
C***DESCRIPTION Computation of integrals over hyper-rectangular
C            regions with singularities.
C            DECUHR is a driver for the integration routine
C            DEADHR, which repeatedly subdivides the region
C            of integration. DEADHR estimates the integrals and the
C            errors over the subregions and subdivides the subregion
C            with greatest estimated errors until the error request
C            is met or MAXPTS function evaluations have been used.
C            The subdivision is done in such a way that we at any
C            stage have only one region containing the singularity.
C            When the singular region is picked it is subdivided in
C            SINGUL + 1 new regions (cutting each of the SINGUL
C            first directions once), creating one new singular region
C            and SINGUL non-singular regions. When a non-singular
C            region is picked we divide this in two halves.
C 
C            For NDIM = 2 the default integration rule is of
C            degree 13 and uses 65 evaluation points.
C            For NDIM = 3 the default integration rule is of
C            degree 11 and uses 127 evaluation points.
C            For NDIM greater then 3 the default integration rule
C            is of degree 9 and uses NUM evaluation points where
C              NUM = 1 + 4*2*NDIM + 2*NDIM*(NDIM-1) + 4*NDIM*(NDIM-1) +
C                    4*NDIM*(NDIM-1)*(NDIM-2)/3 + 2**NDIM
C            The degree 9 rule may also be applied for NDIM = 2
C            and NDIM = 3.
C            A rule of degree 7 is available in all dimensions.
C            The number of evaluation
C            points used by the degree 7 rule is
C              NUM = 1 + 3*2*NDIM + 2*NDIM*(NDIM-1) + 2**NDIM
C 
C            When DECUHR computes estimates to a vector of
C            integrals, all components of the vector are given
C            the same treatment. That is, I(F ) and I(F ) for
C                                            J         K
C            J not equal to K, are estimated with the same
C            subdivision of the region of integration.
C            For integrals with enough similarity, we may save
C            time by applying DECUHR to all integrands in one call.
C            For integrals that varies continuously as functions of
C            some parameter, the estimates produced by DECUHR will
C            also vary continuously when the same subdivision is
C            applied to all components. This will generally not be
C            the case when the different components are given
C            separate treatment. It is essential that all components
C            have the same type of singular behavior. Thus only one
C            value for ALPHA is allowed, which thus has to be the same
C            for all components.
C 
C            On the other hand this feature should be used with
C            caution when the different components of the integrals
C            require clearly different subdivisions.
C 
C   ON ENTRY
C 
C     NDIM   Integer.
C            Number of variables. 1 < NDIM <=  15.
C     NUMFUN Integer.
C            Number of components of the integral.
C     A      Real array of dimension NDIM.
C            Lower limits of integration.
C     B      Real array of dimension NDIM.
C            Upper limits of integration.
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
C            MAXPTS >= 3*NUM and MAXPTS >= MINPTS
C            For 3 < NDIM < 13 the minimum values for MAXPTS are:
C             NDIM =    4   5   6    7    8    9    10   11    12
C            KEY = 3:  459 819 1359 2151 3315 5067 7815 12351 20235
C            KEY = 4:  195 309  483  765 1251 2133 3795  7005 13299
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
C            Output parameter:
C              FUNVLS Real array of dimension NUMFUN
C                     that defines NUMFUN components of the integrand.
C     SINGUL Integer.
C            Dimension of the singularity
C     ALPHA  Real.
C            Degree of homogeneous function.
C            If input ALPHA <= -SINGUL, then ALPHA is estimated by the
C            code and this estimate is used for extrapolation.
C            SINGUL is assumed to be given correct.
C     LOGF   Integer.
C            LOGF = 0, then no logarithmic singularity
C            LOGF = 1, then a logarithmic singularity of order 1.
C            It's value will be estimated by the code if 
C            input ALPHA <= -SINGUL.
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
C     WRKSUB Integer.
C            The maximum allowed number of subregions.
C     NW     Integer.
C            Defines the length of the working array WORK.
C            NW should be >= 
C              3 + 17*NUMFUN + WRKSUB*(NDIM+NUMFUN+1)*2 + EMAX
C              + NUMFUN*(3*WRKSUB+9+EMAX) + (EMAX+1)**2 + 3*NDIM
C 
C     RESTAR Integer.
C            If RESTAR = 0, this is the first attempt to compute
C            the integral.
C            If RESTAR = 1, then we restart a previous attempt.
C            In this case the only parameters for DECUHR that may
C            be changed (with respect to the previous call of DECUHR)
C            are MINPTS, MAXPTS, EPSABS, EPSREL and RESTAR.
C       EMAX Integer.
C            The maximum number of extrapolation steps.
C   ON RETURN
C 
C     RESULT Real array of dimension NUMFUN.
C            Approximations to all components of the integral.
C     ABSERR Real array of dimension NUMFUN.
C            Estimates of absolute errors.
C     NEVAL  Integer.
C            Number of function evaluations used by DECUHR.
C     IFAIL  Integer.
C            IFAIL = 0 for normal exit, when ABSERR(K) <=  EPSABS or
C              ABSERR(K) <=  ABS(RESULT(K))*EPSREL with MAXPTS or less
C              function evaluations for all values of K,
C              1 <= K <= NUMFUN .
C            IFAIL = 1 if MAXPTS was too small for DECUHR
C              to obtain the required accuracy. In this case DECUHR
C              returns values of RESULT with estimated absolute
C              errors ABSERR.
C            IFAIL =  2 if KEY is less than 0 or KEY greater than 4.
C            IFAIL =  3 if NDIM is less than 2 or NDIM greater than 15.
C            IFAIL =  4 if KEY = 1 and NDIM not equal to 2.
C            IFAIL =  5 if KEY = 2 and NDIM not equal to 3.
C            IFAIL =  6 if NUMFUN is less than 1.
C            IFAIL =  7 if A(J) >= B(J), for some J.
C            IFAIL =  8 if MAXPTS is less than MINPTS.
C            IFAIL =  9 if EPSABS < 0 and EPSREL < 0.
C            IFAIL = 10 if RESTAR > 1 or RESTAR < 0.
C            IFAIL = 11 if SINGUL > NDIM or SINGUL < 1.
C            IFAIL = 12 if LOG < 0 or LOG > 1.
C            IFAIL = 13 if MAXPTS is less than (SINGUL+2)*NUM.
C            IFAIL = 14 if ALPHA <=  -SINGUL. (Integral does not exist).
C            IFAIL = 15 if EMAX < 1.
C            IFAIL = 16 if NW is too small.
C            IFAIL = 17 if WRKSUB is too small.
C     WORK   Real array of dimension NW.
C            Used as working storage.
C            WORK(NW) = NSUB, the number of subregions in the data
C            structure.
C            WORK(1),...,WORK(NUMFUN*(WRKSUB+1)) contain
C              the estimated components of the integrals over the
C              subregions.
C            WORK(NUMFUN*(WRKSUB+1)+1),...,WORK(2*NUMFUN*(WRKSUB+1))
C              contain the estimated errors over the subregions.
C            WORK(2*NUMFUN*(WRKSUB+1)+1),...,WORK(2*NUMFUN*(WRKSUB+1)+
C              NDIM*(WRKSUB+1)) contain the centers of the subregions.
C            WORK((2*NUMFUN*+NDIM)*(WRKSUB+1)+1),...,
C              WORK((2*NUMFUN+2*NDIM)*(WRKSUB+1)) contain the
C              subregion half widths.
C            WORK((2*NUMFUN+2*NDIM)*(WRKSUB+1)+1),...,
C              WORK((2*NUMFUN+2*NDIM+1)*(WRKSUB+1))
C              contain the greatest errors in each subregion.
C            WORK((2*NUMFUN+2*NDIM+1)*(WRKSUB+1)+1),...,
C              WORK(2*(NUMFUN+NDIM+1)*(WRKSUB+1)) contain the
C              direction of subdivision in each subregion.
C            WORK(2*(NDIM+NUMFUN+1)*(WRKSUB+1)+1),...,
C              WORK(2*(NDIM+NUMFUN+1)*(WRKSUB+1)+LENW) is used as
C              temporary storage in DEADHR.
C            WORK(2*(NDIM+NUMFUN+1)*(WRKSUB+1)+LENW+1),...,
C             WORK((2*NDIM+3*NUMFUN+2)*(WRKSUB+1)+LENW) contain the
C             tail estimates, Q, over the singular regions.
C            WORK((2*NDIM+3*NUMFUN+2)*(WRKSUB+1)+LENW+1),...,
C             WORK((2*NDIM+3*NUMFUN+2)*(WRKSUB+1)+LENW+NUMFUN*WRKSUB)
C             contain the U-sequence.
C            WORK((2*NDIM+3*NUMFUN+2)*(WRKSUB+1)+LENW+NUMFUN*WRKSUB+1),
C             ...,
C             WORK((2*NDIM+3*NUMFUN+2)*(WRKSUB+1)+LENW+2*NUMFUN*WRKSUB),
C             contain the error estimates to the elements in the U-seq.
C            Define
C            BASE=(2*NDIM+3*NUMFUN+2)*(WRKSUB+1)+LENW+2*NUMFUN*WRKSUB
C            then we have:
C            WORK(BASE+1),...,WORK(BASE+NUMFUN), contain old Q
C             estimates.
C            WORK(BASE+NUMFUN+1),...,WORK(BASE+2*NUMFUN), contain
C             new Q estimates.
C            WORK(BASE+2*NUMFUN+1), ...,WORK(BASE+3*NUMFUN), contain
C             new U-estimates.
C            WORK(BASE+3*NUMFUN+1),...,WORK(BASE+4*NUMFUN), contain
C             error estimates of the U-sequence.
C            WORK(BASE+4*NUMFUN+1),...,WORK(BASE+5*NUMFUN), contain
C             estimates of the extrapolation error.
C            WORK(BASE+5*NUMFUN+1),...,WORK(BASE+5*NUMFUN+1+EMAX),
C             contain extrapolation denominators.
C            WORK(BASE+5*NUMFUN+1+EMAX+1),...,
C             WORK(BASE+5*NUMFUN+(1+EMAX)**2), contain a work array
C             used to compute the global array.
C            WORK(BASE+5*NUMFUN+(1+EMAX)**2+1),...,
C             WORK(BASE+5*NUMFUN+(1+EMAX)**2+NUMFUN*(EMAX+1)), contain
C             extrapolation tableau.
C            WORK(BASE+5*NUMFUN+(1+EMAX)**2+NUMFUN*(EMAX+1)+1),...,
C             WORK(BASE+5*NUMFUN+(1+EMAX)**2+NUMFUN*(EMAX+1)+NDIM),
C             contain a working array used in the rule evaluation
C           routine.
C     IWORK  Integer array of dimension 2*WRKSUB + NDIM.
C            Used as working storage.
C 
C***LONG DESCRIPTION
C 
C   The information for each subregion is contained in the
C   data structure WORK.
C   When passed on to DEADHR, WORK is split into nineteen arrays:
C   VALUES, ERRORS, CENTRS, HWIDTS, GREATE, DIR,
C   WORK, Q,U,E,QOLD,QNEW,UNEW,UERR,EXTERR,NE,BETA,T and DIFF.
C   The first 6 are connected to each subregion:
C    VALUES contains the estimated values of the integrals.
C    ERRORS contains the estimated errors.
C    CENTRS contains the centers of the subregions.
C    HWIDTS contains the half widths of the subregions.
C    GREATE contains the greatest estimated error for each subregion.
C    DIR    contains the directions for further subdivision.
C   Number 7 is a work array for the integration rules:
C    WORK is used as work array in DEADHR.
C   The next 11 are used by the extrapolation procedure:
C    Q contains the estimates over the singular regions(Tail-estimates).
C    U contains the terms in the series.
C    E contains the estimated errors in each U-term.
C    QOLD  contains the estimates for the previous singular region.
C    QNEW  contains the estimates for the current singular region.
C    UNEW  contains the estimates for the new U-term (or update
C          corrections).
C    UERR  contains the estimates for the error in the new U-term.
C    EXERR contains the estimates for the extrapolation error.
C    NE    work array used in the extrapolation process.
C    BETA  work array used to compute the global error.
C    T     contains the last row in the extrapolation tableau.
C   And finally a work array for the singular region.
C    DIFF  work array used in the rule evaluation routine.
C 
C   The integer work array is spit in three: LIST, UPOINT and ORDER.
C    LIST   contains the pointers used in maintaining the heap.
C    UPOINT contains pointers from each subregion to the U-term
C           where this region is one part. This is used in the
C           the updating process.
C    ORDER  Pointer array giving the subdivision sequence of the
C           singular region. Connected to DIFF.
C   The data structures for the subregions are in DEADHR organized
C   as a heap, and the size of GREATE(I) defines the position of
C   region I in the heap. The heap is maintained by the program
C   DETRHR. The singular region is not part of the heap due to the
C   fact that the error estimates associated with this region will
C   have to be updated even when this region has not been subdivided.
C 
C   The subroutine for estimating the integral and the error over
C   each subregion, DERLHR, uses WORK2 as a work array.
C 
C***REFERENCES
C 
C   T.O.Espelid and A.Genz, DECUHR: An Algorithm for Automatic 
C   Integration of Singular Functions over a Hyperrectangular Region.
C   Numerical Algorithms 8(1994), PP. 201-220.
C 
C   T.O.Espelid, On integrating Vertex Singularities using
C   Extrapolation,  BIT, 34:62-79, 1994.
C 
C   T.O.Espelid, On integrating Singularities using non-Uniform
C   Subdivision and Extrapolation,  In Numerical Integration IV,
C   H. Brass and G. Hammerlin (Eds), Birkhauser Verlag, Basel, 
C   Vol. 112:77-89, 1993.
C 
C   J.Berntsen, T.O.Espelid and A.Genz, An Adaptive Algorithm
C   for the Approximate Calculation of Multiple Integrals,
C   ACM Transactions on Mathematical Software,Vol.17,No.4,
C   December 1991,Pages 437-451.
C 
C   J.Berntsen, T.O.Espelid and A.Genz, DCUHRE: An Adaptive
C   Multidimensional Integration Routine for a Vector of
C   Integrals, ACM Transactions on Mathematical Software, Vol.17,No.4,
C   December 1991,Pages 452-456.
C 
C***ROUTINES CALLED DECHHR, DEADHR, DECALP
C***END PROLOGUE DECUHR
C 
C   Global variables.
C 
      EXTERNAL FUNSUB
      INTEGER NDIM,NUMFUN,MINPTS,MAXPTS,KEY,NW,RESTAR,WRKSUB
      INTEGER NEVAL,IFAIL,SINGUL,IWORK(2*WRKSUB),EMAX,LOGF
      DOUBLE PRECISION A(NDIM),B(NDIM),EPSABS,EPSREL,ALPHA
      DOUBLE PRECISION RESULT(NUMFUN),ABSERR(NUMFUN),WORK(NW)
C 
C   Local variables.
C 
C   MAXDIM Integer.
C          The maximum allowed value of NDIM.
C   MAXWT  Integer. The maximum number of weights used by the
C          integration rule.
C   WTLENG Integer.
C          The number of generators used by the selected rule.
C   WORK2  Real work space. The length
C          depends on the parameters MAXDIM and MAXWT.
C   MAXSUB Integer.
C          The maximum allowed number of subdivisions
C          for the given values of KEY, NDIM and MAXPTS.
C   NUM    Integer. The number of integrand evaluations needed
C          over each subregion.
C 
      INTEGER MAXWT, WTLENG, MAXDIM, LENW2, MAXSUB
      INTEGER NUM, NSUB, LENW, KEYF
      PARAMETER ( MAXDIM = 15 )
      PARAMETER ( MAXWT = 14 )
      PARAMETER ( LENW2 = 2*MAXDIM*(MAXWT+1) + 12*MAXWT + 2*MAXDIM )
      INTEGER I1, I2, I3, I4, I5, I6, I7, I8, I9, I10, 
     +     I11, I12, I13, I14, I15, I16, I17, I18, I19
      INTEGER K1, K2, K3, K4, K5, K6, K7, K8
      DOUBLE PRECISION WORK2(LENW2)
C 
C***FIRST EXECUTABLE STATEMENT DECUHR
C 
C   Compute NUM, WTLENG, MAXSUB,
C   and check the input parameters.
C 
      CALL DECHHR ( MAXDIM, NDIM, NUMFUN, A, B, ALPHA, SINGUL, LOGF,
     +     MINPTS, MAXPTS, EPSABS, EPSREL, KEY, NW, RESTAR, EMAX, 
     +     WRKSUB, NUM, MAXSUB, KEYF, IFAIL, WTLENG )
      IF ( IFAIL .NE. 0 ) RETURN
C
C   Check if we have to estimate ALPHA
C
      IF (ALPHA.LE.-SINGUL) THEN
	   CALL DECALP ( NDIM, A, B, NUMFUN, ALPHA, 
     +     SINGUL, LOGF, FUNSUB, WORK(NW-NDIM-NUMFUN), WORK(NW-NUMFUN) )
C 
C    Check the computed ALPHA-value.
C 
          IF (ALPHA.LE.-SINGUL) THEN
            IFAIL=14
            RETURN
          END IF
      END IF
C     
C     Compute the size of the temporary work space needed in DEADHR.
C     
      LENW = 16*NUMFUN
C     
C     Split up the work space.
C     
      I1 = 1
      I2 = I1+(WRKSUB+1)*NUMFUN
      I3 = I2+(WRKSUB+1)*NUMFUN
      I4 = I3+(WRKSUB+1)*NDIM
      I5 = I4+(WRKSUB+1)*NDIM
      I6 = I5+WRKSUB+1
      I7 = I6+WRKSUB+1
      I8 = I7+LENW
      I9 = I8+(WRKSUB+1)*NUMFUN
      I10 = I9+WRKSUB*NUMFUN
      I11 = I10+WRKSUB*NUMFUN
      I12 = I11+NUMFUN
      I13 = I12+NUMFUN
      I14 = I13+NUMFUN
      I15 = I14+NUMFUN
      I16 = I15+NUMFUN
      I17 = I16+EMAX
      I18 = I17+(EMAX+1)**2
      I19 = I18+(EMAX+1)*NUMFUN
      K1 = 1
      K2 = K1+WTLENG*NDIM
      K3 = K2+WTLENG*5
      K4 = K3+WTLENG
      K5 = K4+NDIM
      K6 = K5+NDIM
      K7 = K6+NDIM
      K8 = K7+3*WTLENG
C 
C   On restart runs the number of subregions from the
C   previous call is assigned to NSUB.
C 
      IF ( RESTAR .EQ. 1 ) NSUB = WORK(NW)
      CALL DEADHR (NDIM,NUMFUN,A,B,MAXSUB,FUNSUB,SINGUL,ALPHA,LOGF,
     +     EPSABS,EPSREL,KEYF,RESTAR,NUM,LENW,WTLENG,EMAX,MINPTS,MAXPTS,
     +     NSUB,RESULT,ABSERR,NEVAL,IFAIL,WORK(I1),WORK(I2),WORK(I3),
     +     WORK(I4),WORK(I5),WORK(I6),WORK(I7),WORK(I8),WORK(I9),
     +     WORK(I10),WORK(I11),WORK(I12),WORK(I13),WORK(I14),WORK(I15),
     +     WORK(I16),WORK(I17),WORK(I18),WORK(I19),WORK2(K1),WORK2(K2),
     +     WORK2(K3),WORK2(K4),WORK2(K5),WORK2(K6),WORK2(K7),WORK2(K8),
     +     IWORK(1),IWORK(WRKSUB+1),IWORK(2*WRKSUB+1))
      WORK(NW) = NSUB
      RETURN
C 
C***END DECUHR
C 
      END
