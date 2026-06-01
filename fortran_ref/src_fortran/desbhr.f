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
      SUBROUTINE DESBHR (NDIM,CENTRS,HWIDTS,DIR,DIRECT,POINTR,VACANT,
     +     ORDER)
C***BEGIN PROLOGUE DESBHR
C***REFER TO DECUHR
C***PURPOSE DESBHR cuts the given region one or more times.
C***LAST MODIFICATION 92-09-16
C***DESCRIPTION The number of new subregions depends on the region.
C 
C           Each cut creates two new halves of the current region.
C 
C           The singular region will be cut -DIRECT times in
C           directions 1,2,3,...,-DIRECT. The order of these
C           -DIRECT cuts is given by ORDER.
C           This gives -DIRECT + 1 new subregions.
C 
C           A non-singular region will be cut once in direction
C           DIRECT. This gives two new subregions.
C***ROUTINES CALLED-NONE
C 
C   ON ENTRY
C 
C     NDIM   Integer.
C            Number of variables.
C     CENTRS Real array of dimension (NDIM,0:*).
C            Used to store the centers of the stored subregions.
C     HWIDTS Real array of dimension (NDIM,0:*).
C            Used to store the half widths of the stored subregions.
C     DIR    Real array of dimension (0:*).
C            DIR is used to store the directions for
C            further subdivision.
C     DIRECT If positive: the direction of subdivision.
C            If negative: cut direction 1,2,...,-DIRECT.
C            This will be the case if this region is the singular one.
C     POINTR Pointer to the position in the data structure where
C            the new subregions are to be stored.
C     VACANT Pointer to the region which has to be subdivided.
C            One of the new reions will be stored at this vacant element
C     ORDER  Integer array of pointers.
C            In case this is the singular region it will be cut in
C            directions 1,2,...,-DIRECT, in the order given by this
C            pointer array.
C 
C   ON RETURN
C 
C     CENTRS Real array of dimension (NDIM,0:*).
C            Unchanged array except for position VACANT and
C            positions POINTR, POINTR+1,... (depending on the
C            number of new subregions).
C     HWIDTS Real array of dimension (NDIM,0:*).
C            Unchanged array except for position VACANT and
C            positions POINTR, POINTR+1,... (depending of the
C            number of new subregions).
C     DIR    Real array of dimension (0:*)
C            Unchanged array except for position VACANT and
C            positions POINTR, POINTR+1,... (depending of the
C            number of new subregions). Saves the last direction
C            these subregions are cut. Except for the singular region
C            where the value is left unchanged equal to DIRECT.
C 
C***END PROLOGUE DESBHR
C
C   Global variables.
C
      INTEGER NDIM,POINTR,DIRECT,VACANT,ORDER(NDIM)
      DOUBLE PRECISION CENTRS(NDIM,0:*),HWIDTS(NDIM,0:*),DIR(0:*)
C
C   Local variables.
C
      INTEGER I,J,NEXT
C
C   FIRST EXECUTABLE STATEMENT DESBHR
C
      IF (DIRECT.GT.0) THEN
C 
C     When DIRECT is positive, indicating a non-singular region,
C     we only cut the region in two halves.
C 
         DO 10 J=1,NDIM
            HWIDTS(J,POINTR)=HWIDTS(J,VACANT)
            CENTRS(J,POINTR)=CENTRS(J,VACANT)
 10      CONTINUE
         HWIDTS(DIRECT,VACANT)=HWIDTS(DIRECT,VACANT)/2
         HWIDTS(DIRECT,POINTR)=HWIDTS(DIRECT,VACANT)
         CENTRS(DIRECT,POINTR)=CENTRS(DIRECT,VACANT)+HWIDTS(DIRECT,
     +    POINTR)
         CENTRS(DIRECT,VACANT)=CENTRS(DIRECT,VACANT)-HWIDTS(DIRECT,
     +    VACANT)
         DIR(POINTR)=DIRECT
      ELSE
C 
C     DIRECT negative signals the singular region: the absolute value
C     tells us how many co-ordinates are involved (always assumed to
C     be 1,2,...,-DIRECT).
C     In this case we slice the singular region -DIRECT times
C     in directions 1,2,...,-DIRECT. These cuts are performed according
C     to the array ORDER. We get -DIRECT+1 new pieces
C     and the reduced singular region ends up in the
C     the same position as the previous one: thus VACANT = 0.
C 
         DO 30 I=1,-DIRECT
            NEXT=ORDER(I)
            HWIDTS(NEXT,VACANT)=HWIDTS(NEXT,VACANT)/2
            DO 20 J=1,NDIM
               HWIDTS(J,POINTR+I-1)=HWIDTS(J,VACANT)
               CENTRS(J,POINTR+I-1)=CENTRS(J,VACANT)
 20         CONTINUE
            DIR(POINTR+I-1)=NEXT
            CENTRS(NEXT,POINTR+I-1)=CENTRS(NEXT,VACANT)
     +           + HWIDTS(NEXT,VACANT)
            CENTRS(NEXT,VACANT)=CENTRS(NEXT,VACANT)-HWIDTS(NEXT,VACANT)
 30      CONTINUE
      END IF
      DIR(VACANT)=DIRECT
      END
