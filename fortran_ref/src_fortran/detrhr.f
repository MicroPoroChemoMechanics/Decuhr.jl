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
      SUBROUTINE DETRHR (DVFLAG,SBRGNS,GREATE,LIST,NEW)
C***BEGIN PROLOGUE DETRHR
C***REFER TO DECUHR
C***PURPOSE  DETRHR maintains a heap of subregions.
C***AUTHOR   Terje O. Espelid, Department of Informatics,
C            University of Bergen,  Hoyteknologisenteret,
C            N-5020 Bergen, Norway
C            Email..  terje@ii.uib.no
C***LAST MODIFICATION 92-09-16
C***DESCRIPTION DETRHR maintains a heap of subregions.
C            The subregions are stored in a partially sorted
C            binary tree, ordered according to the size of the
C            greatest error estimates of each subregion(GREATE).
C            The subregion with greatest error estimate is in the
C            first position of the heap.
C 
C   PARAMETERS
C 
C     DVFLAG Integer.
C            If DVFLAG = 1, we remove the subregion with
C            greatest error from the heap.
C            If DVFLAG = 2, we insert a new subregion in the heap.
C     SBRGNS Integer.
C            Number of subregions in the heap.
C     GREATE Real array of dimension SBRGNS.
C            Used to store the greatest estimated errors in
C            all subregions.
C     LIST   Integer array of dimension SBRGNS.
C            Used as a partially ordered list of pointers to the
C            different subregions. This list is a heap where the
C            element on top of the list is the subregion with the
C            greatest error estimate.
C     NEW    Integer.
C            Index to the new region to be inserted in the heap.
C***ROUTINES CALLED-NONE
C***END PROLOGUE DETRHR
C 
C   Global variables.
C 
      INTEGER DVFLAG,NEW,SBRGNS,LIST(*)
      DOUBLE PRECISION GREATE(0:*)
C 
C   Local variables.
C 
C   GREAT  is used as intermediate storage for the greatest error of a
C          subregion.
C   SUBRGN Position of child/parent subregion in the heap.
C   SUBTMP Position of parent/child subregion in the heap.
      INTEGER SUBRGN,SUBTMP
      DOUBLE PRECISION GREAT
C 
C***FIRST EXECUTABLE STATEMENT DTRTRI
C 
C     If DVFLAG = 1, we will reduce the partial ordered list by the
C     element with greatest estimated error. Thus the element in
C     in the heap with index LIST(1) is vacant and can be used later.
C     Reducing the heap by one element implies that the last element
C     should be re-positioned.
C 
      IF (DVFLAG.EQ.1) THEN
         GREAT=GREATE(LIST(SBRGNS))
         SBRGNS=SBRGNS-1
         SUBRGN=1
 10      SUBTMP=2*SUBRGN
         IF (SUBTMP.LE.SBRGNS) THEN
            IF (SUBTMP.NE.SBRGNS) THEN
C 
C     Find max. of left and right child.
C 
               IF (GREATE(LIST(SUBTMP)).LT.GREATE(LIST(SUBTMP+1))) THEN
                  SUBTMP=SUBTMP+1
               END IF
            END IF
C 
C     Compare max.child with parent.
C     If parent is max., then done.
C 
            IF (GREAT.LT.GREATE(LIST(SUBTMP))) THEN
C 
C     Move the pointer at position SUBTMP up the heap.
C 
               LIST(SUBRGN)=LIST(SUBTMP)
               SUBRGN=SUBTMP
               GO TO 10
            END IF
         END IF
C 
C     Update the pointer.
C 
         IF (SBRGNS.GT.0) THEN
            LIST(SUBRGN)=LIST(SBRGNS+1)
         END IF
      ELSE IF (DVFLAG.EQ.2) THEN
C 
C     If DVFLAG = 2, find the position for the NEW region in the heap.
C 
         GREAT=GREATE(NEW)
         SUBRGN=SBRGNS
 20      SUBTMP=SUBRGN/2
         IF (SUBTMP.GE.1) THEN
C 
C     Compare max.child with parent.
C     If parent is max, then done.
C 
            IF (GREAT.GT.GREATE(LIST(SUBTMP))) THEN
C 
C     Move the pointer at position SUBTMP down the heap.
C 
               LIST(SUBRGN)=LIST(SUBTMP)
               SUBRGN=SUBTMP
               GO TO 20
            END IF
         END IF
C 
C     Set the pointer to the new region in the heap.
C 
         LIST(SUBRGN)=NEW
      END IF
C 
C***END DETRHR
C 
      END
