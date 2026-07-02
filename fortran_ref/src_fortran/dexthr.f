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
      SUBROUTINE DEXTHR (NUMFUN,ALPHA,LOGF,SINGUL,EXSTEP,N,T,UPDATE,
     +UNEW,QNEW,QOLD,FIRST,EMAX,UERR,RESULT,ABSERR,EXTERR,NE,BETA)
C***BEGIN PROLOGUE DEXTHR
C***KEYWORDS Linear extrapolation, homogeneous functions,
C     logarithmic singularities, error estimation.
C***PURPOSE To compute better estimates to a vector of approximations
C     to multidimensional integrals and to provide new and
C     updated error estimates.
C***LAST MODIFICATION 92-09-16
C***DESCRIPTION
C            The routine uses linear extrapolation to compute better
C            approximations to each component in a vector of
C            multidimensional integrals. All components are assumed to
C            be singular due to a homogeneous function of degree ALPHA
C            in a vertex of dimension SINGUL. In addition we may have
C            a logarithmic singularity in the same vertex(dim. SINGUL).
C            A series, with tail correction, approach is used, assuming
C            that the terms are given with estimates of the error in
C            each term. New error estimates are computed too. The
C            routine have two options: either a new extrapolation term
C            is provided and we take a new extrapolation step,
C            or we update one previously computed term in the series and
C            therefore have to update the extrapolation tableau.
C 
C   ON ENTRY
C 
C     NUMFUN Integer.
C            Number of components of the integral.
C     ALPHA  Real.
C            Degree of homogeneous function: ALPHA > -SINGUL.
C            This singularity is assumed the same in all integrands.
C     LOGF    Integer
C            If LOGF = 1 then there is a logarithmic singularity in all
C            integrands, else there is no logarithmic singularity.
C            We assume the logarithm to appear only in power of 1.
C     SINGUL Integer
C            The dimension of the singularity.
C     EXSTEP Integer.
C            The exponent sequence's step size in the error expansion.
C     N      Integer
C            The number of U-terms in the series.
C     T      Real array of dimension (NUMFUN,0:EMAX)
C            Contains the last row in the extrapolation tableau for
C            each function in the vector.
C     UPDATE Integer
C            = 0 then this is a new extrapolation step.
C            > 0 then this is  a step where we have to correct the
C            existing tableau. The value of UPDATE gives the index to th
C            u-elment that has been modified.
C     UNEW   Real array of dimension NUMFUN.
C            If UPDATE = 0 then this gives the next terms in the series,
C            else it is the correction to the u-values with index UPDATE
C     QNEW   Real array of dimension NUMFUN.
C            The new tail correction for all functions in the vector.
C     QOLD   Real array of dimension NUMFUN.
C            The last tail corrections for all functions in the vector.
C     FIRST  Logical.
C            Value .TRUE. indicates that this is the first time
C            this routine is called.
C     EMAX   Integer
C            The maximum allowed number of extrapolation steps.
C     UERR   Real array of dimension (NUMFUN,N)
C            The estimated errors of all U-terms in the series.
C   ON RETURN
C 
C     T      Real array of dimension (NUMFUN,0:EMAX)
C            Contains the last row in the extrapolation tableau for
C            each function in the vector after the extrapolation. In
C            case this is an updating steps, then each element may
C            have been changed.
C     RESULT Real array of dimension NUMFUN
C            Contains the new approximations for the components
C            of the integral.
C     ABSERR Real array of dimension NUMFUN.
C            Returns the global errors for all components.
C            This includes both the pure extrapolation
C            error and the effect of not knowing the U-terms exactly.
C     EXTERR Real array of dimension NUMFUN.
C            These errors are associated with the singular region and
C            they are the pure extrapolation errors.
C     NE     Real array of dimension (0:EMAX).
C            A table of denominators to be used in the extrapolation.
C            Dummy parameter
C     BETA   Real Array of dimension (EMAX +1)(EMAX +1)
C            A table of coefficients to be used in the error estimation.
C 
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
C***ROUTINES CALLED NONE
C***END PROLOGUE DEXTHR
C 
C   Global variables.
C 
      INTEGER N,SINGUL,EMAX,UPDATE,NUMFUN,LOGF,EXSTEP
      DOUBLE PRECISION ALPHA,T(NUMFUN,0:EMAX),QOLD(NUMFUN),QNEW(NUMFUN),
     +UNEW(NUMFUN),UERR(NUMFUN,N),NE(EMAX),BETA(0:EMAX,0:EMAX),
     +RESULT(NUMFUN),ABSERR(NUMFUN),EXTERR(NUMFUN)
      LOGICAL FIRST
C 
C   Local variables.
C 
      INTEGER I,J,STEPS
      DOUBLE PRECISION SAVE1,SAVE2,CONST,ES
      PARAMETER (CONST=10)
C 
C  CONST heuristic constant used in the error estimation.
C  STEP  integer; keeping track of the number of extrapolation
C        steps actually used.
C 
C***FIRST EXECUTABLE STATEMENT DEXTHR
C 
C     Check if this is the first time the routine is called
C 
      IF (FIRST) THEN
C 
C     Initialize the extrapolation tableaus
C 
         DO 10 J=1,NUMFUN
            T(J,0)=QOLD(J)
 10      CONTINUE
C 
C     Compute all denominators and take into account if there
C     is a logarithmic term present.
C 
         NE(1)=2**(SINGUL+ALPHA)-1
         IF (LOGF.EQ.1) THEN
            DO 30 I=3,EMAX,2
               ES=NE(I-2)
               DO 20 J=1,EXSTEP
                  ES=2*ES+1
 20            CONTINUE
               NE(I)=ES
 30         CONTINUE
            DO 40 I=2,EMAX,2
               NE(I)=NE(I-1)
 40         CONTINUE
         ELSE
C 
C     assuming log = 0, (note log > 1 has not been implemented)
C 
            DO 60 I=2,EMAX,1
               ES=NE(I-1)
               DO 50 J=1,EXSTEP
                  ES=2*ES+1
 50            CONTINUE
               NE(I)=ES
 60         CONTINUE
         END IF
C 
C     Initialize the beta-factors to be used in the error estimation.
C 
         DO 90 J=0,EMAX
            BETA(0,J)=1
            DO 70 I=1,J
               BETA(I,J)=BETA(I,J-1)+(BETA(I,J-1)-BETA(I-1,J-1))/NE(J)
 70         CONTINUE
            DO 80 I=J+1,EMAX
               BETA(I,J)=0
 80         CONTINUE
 90      CONTINUE
      END IF
C 
C     The number of extrapolation steps
C 
      STEPS=MIN(N,EMAX)
C 
C     Check what kind of step this is: new extrapolation or modifying
C 
      IF (UPDATE.EQ.0) THEN
C 
C     A new extrapolation step.
C 
         DO 110 J=1,NUMFUN
            SAVE1=T(J,0)+(UNEW(J)+(QNEW(J)-QOLD(J)))
            DO 100 I=1,STEPS
               SAVE2=SAVE1+(SAVE1-T(J,I-1))/NE(I)
               T(J,I-1)=SAVE1
               SAVE1=SAVE2
 100        CONTINUE
            T(J,STEPS)=SAVE1
 110     CONTINUE
C 
C     A modification step.
C 
      ELSE IF (UPDATE.LT.N-STEPS) THEN
C 
C     Simply add the correction to all elements in the tableau.
C 
         DO 130 J=1,NUMFUN
            DO 120 I=0,STEPS
               T(J,I)=T(J,I)+UNEW(J)
 120        CONTINUE
 130     CONTINUE
      ELSE
         DO 150 J=1,NUMFUN
            DO 140 I=0,STEPS
               T(J,I)=T(J,I)+UNEW(J)*(1-BETA(N-UPDATE+1,I))
 140        CONTINUE
 150     CONTINUE
      END IF
C 
C     Then compute the error estimates.
C     First the extrapolation error and then the U-effect.
C     The error is accumulated in qnew
C 
      DO 180 J=1,NUMFUN
         EXTERR(J)=CONST*ABS(T(J,STEPS)-T(J,STEPS-1))
         QNEW(J)=EXTERR(J)
C 
C     Note: The last U-errors are effected by the extrapolation-process
C 
         DO 160 I=1,STEPS
            QNEW(J)=QNEW(J)+ABS(1-BETA(I,STEPS))*UERR(J,N+1-I)
 160     CONTINUE
         DO 170 I=STEPS+1,N
            QNEW(J)=QNEW(J)+UERR(J,N+1-I)
 170     CONTINUE
 180  CONTINUE
C 
C   Define the results and the new errors. We update only those
C   components which have an improved error estimate.
C 
      DO 190 J=1,NUMFUN
	  IF ( QNEW(J).LE.ABSERR(J)) THEN
            RESULT(J) = T(J,STEPS)
       	    ABSERR(J) = QNEW(J)
          END IF
 190  CONTINUE
C 
C***END DEXTHR
C 
      END
