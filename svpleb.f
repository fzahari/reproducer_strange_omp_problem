C  5 JUL 05 - MWS - SELECT NEW ATOM,BASIS,EFP,PCM,DAF DIMENSIONS
C 13 FEB 05 - DMC - NEW MODULE FOR LEBEDEV ROOTS/WEIGHTS ON UNIT SPHERE
C
C*DECK SVPLEB *DECK LEBPTW
       SUBROUTINE LEBPTW(XYZ,XYZTR,W,NREQ,MXDIM,DOTRANSP)
       DOUBLE PRECISION XYZ,W,XYZTR
       LOGICAL DOTRANSP
       LOGICAL SKIP
C
C      ----- LEBEDEV INTEGRATION POINTS AND WEIGHTS ON UNIT SPHERE -----
C
       COMMON/IOFILE/IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
C
C      IF DOTRANSP IS SET TO .TRUE.:
C      XYZ IS FOR TEMPORARY STORAGE; RESULTS APPEAR IN XYZTR AND W
C      IF DOTRANSP IS SET TO .FALSE.:
C      XYZTR IS NOT USED; RESULTS APPEAR IN XYZ AND W
C
C      WEIGHTS ARE NORMALIZED TO SUM UP TO UNITY
C
C      THIS PROGRAM SKIPS 74, 230, AND 266 BECAUSE THESE HAVE NEGATIVE
C      WEIGHTS FOR SOME QUADRATURE POINTS.  THE GRIDS FOR 146, 170,
C      350 (PARTICULARLY), 4802, AND 5294 CAN GIVE RESULTS THAT SEEM
C      TO DIFFER CONSIDERABLY FROM NEARBY POINTS.  SEE D.M.CHIPMAN
C      AND M.DUPUIS, THEORET.CHEM.ACC. 107, 90-102(2002) FOR MORE INFO.
C
C      LEBEDEV GRIDS EXIST ONLY FOR ODD L.  L=3 HAS 6 POINTS, L=194
C      HAS 194 POINTS, L=29 HAS 302 POINTS, L=XX HAS 1202 POINTS, ETC.
C      SEE X.G.WANG AND T.CARRINGTON, J.THEORET.COMPUT.CHEM.
C      2, 599-608(2003) FOR A PICTURE COMPARING LEBEDEV GRIDS TO A
C      POLAR COORDINATE (AKA LEGENDRE) GRID.
C
       DIMENSION XYZ(NREQ,3),W(NREQ),XYZTR(3,MXDIM)
C
       SKIP = .TRUE.
       N=-99
       IF(NREQ.EQ.   6) CALL LD0006(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
       IF(NREQ.EQ.  14) CALL LD0014(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
       IF(NREQ.EQ.  26) CALL LD0026(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
       IF(NREQ.EQ.  38) CALL LD0038(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
       IF(NREQ.EQ.  50) CALL LD0050(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
C      OMIT N=74 BECAUSE OF NEGATIVE WEIGHTS
       IF(SKIP) GO TO 76
       IF(NREQ.EQ.  74) CALL LD0074(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
   76  CONTINUE
       IF(NREQ.EQ.  86) CALL LD0086(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
       IF(NREQ.EQ. 110) CALL LD0110(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
       IF(NREQ.EQ. 146) CALL LD0146(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
       IF(NREQ.EQ. 170) CALL LD0170(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
       IF(NREQ.EQ. 194) CALL LD0194(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
       IF(SKIP) GO TO 301
C      OMIT N=230 BECAUSE OF NEGATIVE WEIGHTS
       IF(NREQ.EQ. 230) CALL LD0230(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
C      OMIT N=266 BECAUSE OF NEGATIVE WEIGHTS
       IF(NREQ.EQ. 266) CALL LD0266(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
  301  CONTINUE
       IF(NREQ.EQ. 302) CALL LD0302(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
C      OHD0350 seems to have better accuracy
C       IF(NREQ.EQ. 350) CALL LD0350(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
       IF(NREQ.EQ. 434) CALL LD0434(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
       IF(NREQ.EQ. 590) CALL LD0590(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
       IF(NREQ.EQ. 770) CALL LD0770(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
       IF(NREQ.EQ. 974) CALL LD0974(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
       IF(NREQ.EQ.1202) CALL LD1202(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
       IF(NREQ.EQ.1454) CALL LD1454(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
       IF(NREQ.EQ.1730) CALL LD1730(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
       IF(NREQ.EQ.2030) CALL LD2030(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
       IF(NREQ.EQ.2354) CALL LD2354(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
       IF(NREQ.EQ.2702) CALL LD2702(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
       IF(NREQ.EQ.3074) CALL LD3074(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
       IF(NREQ.EQ.3470) CALL LD3470(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
       IF(NREQ.EQ.3890) CALL LD3890(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
       IF(NREQ.EQ.4334) CALL LD4334(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
       IF(NREQ.EQ.4802) CALL LD4802(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
       IF(NREQ.EQ.5294) CALL LD5294(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
       IF(NREQ.EQ.5810) CALL LD5810(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
C
CMV    Class of grids developed by A.S.Popov. They can be more accurate
CMV    than Lebedev-Laikov grids with comparable number of points,
CMV    especially those with rotational octahedral symmetry.
CMV    Rotational octahedral symmetry grids (without inversion):
       IF(NREQ.EQ. 246) CALL  OD0246(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
       IF(NREQ.EQ. 264) CALL  OD0264(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
       IF(NREQ.EQ. 342) CALL  OD0342(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
       IF(NREQ.EQ. 432) CALL  OD0432(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
CMV    Full octahedral symmetry grids (same symmetry as Lebedev grids):
       IF(NREQ.EQ. 350) CALL OHD0350(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
       IF(NREQ.EQ. 398) CALL OHD0398(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
CMV    Rotational icosahedral symmetry grids:
       IF(NREQ.EQ. 192) CALL  ID0192(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
       IF(NREQ.EQ. 212) CALL  ID0212(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
       IF(NREQ.EQ. 242) CALL  ID0242(XYZ(1,1),XYZ(1,2),XYZ(1,3),W,N)
       IF(N.NE.NREQ) THEN
          WRITE(IW,9999) NREQ
          WRITE(IW,'(A,I8)') 'N=',N
!          CALL ABRT
          CALL EXIT()
       ENDIF
       IF (DOTRANSP) THEN
           CALL LEBTR(XYZ,XYZTR,N,MXDIM)
       END IF
C
 9999  FORMAT(' LEBEDEV INTEGRATION WITH ',I4,' POINTS NOT AVAILABLE.',
     1        ' STOP')
       RETURN
       END
       SUBROUTINE LEBTR(XYZ,XYZTR,N,MXDIM)
       DOUBLE PRECISION XYZ,XYZTR,TMP
       DIMENSION XYZ(N,3),XYZTR(3,MXDIM)
       DO K=1,3
          DO I=1,N
             TMP=XYZ(I,K)
             XYZ(I,K)=XYZTR(K,I)
             XYZTR(K,I)=TMP
          ENDDO
       ENDDO
       RETURN
       END
       SUBROUTINE LEB_GEN_OH(CODE, NUM, X, Y, Z, W, A, B, V)
       IMPLICIT LOGICAL(A-Z)
       DOUBLE PRECISION X(*),Y(*),Z(*),W(*)
       DOUBLE PRECISION A,B,V
       DOUBLE PRECISION C
       INTEGER CODE
       INTEGER NUM
       INTEGER IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA
       COMMON/IOFILE/IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED FROM C TO FORTRAN77 BY HAND.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] V.I. LEBEDEV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
CHVD        ORDER OF ACCURACY"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
CHVD
CHVD   [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
CHVD       "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
CHVD
CHVD   [4] V.I. LEBEDEV
CHVD       "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
CHVD       SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
CHVD
CHVD   [5] V.I. LEBEDEV
CHVD       "QUADRATURES ON A SPHERE"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
CHVD       1976, PP. 10-24.
CHVD
CHVD   [6] V.I. LEBEDEV
CHVD       "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
CHVD        ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
CHVD        OCTAHEDRON GROUP WITH INVERSION"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
CHVD       1975, PP. 44-51.
CHVD
CVW
CVW    GIVEN A POINT ON A SPHERE (SPECIFIED BY A AND B), GENERATE ALL
CVW    THE EQUIVALENT POINTS UNDER OH SYMMETRY, MAKING GRID POINTS WITH
CVW    WEIGHT V.
CVW    THE VARIABLE NUM IS INCREASED BY THE NUMBER OF DIFFERENT POINTS
CVW    GENERATED.
CVW
CVW    DEPENDING ON CODE, THERE ARE 6...48 DIFFERENT BUT EQUIVALENT
CVW    POINTS.
CVW
CVW    CODE=1:   (0,0,1) ETC                                (  6 POINTS)
CVW    CODE=2:   (0,A,A) ETC, A=1/SQRT(2)                   ( 12 POINTS)
CVW    CODE=3:   (A,A,A) ETC, A=1/SQRT(3)                   (  8 POINTS)
CVW    CODE=4:   (A,A,B) ETC, B=SQRT(1-2 A^2)               ( 24 POINTS)
CVW    CODE=5:   (A,B,0) ETC, B=SQRT(1-A^2), A INPUT        ( 24 POINTS)
CVW    CODE=6:   (A,B,C) ETC, C=SQRT(1-A^2-B^2), A/B INPUT  ( 48 POINTS)
CMV    General octahedral points:
CMV    CODE=7:   (A,B,C) ETC, C=SQRT(1-A^2-B^2), A/B INPUT  ( 24 POINTS)
CVW
       GOTO (1,2,3,4,5,6,7) CODE
       WRITE (IW,*) 'GEN_OH: INVALID CODE'
!       CALL ABRT
       CALL EXIT()
    1  CONTINUE
       A=1.0D+00
       X(1) =  A
       Y(1) =  0.0D+00
       Z(1) =  0.0D+00
       W(1) =  V
       X(2) = -A
       Y(2) =  0.0D+00
       Z(2) =  0.0D+00
       W(2) =  V
       X(3) =  0.0D+00
       Y(3) =  A
       Z(3) =  0.0D+00
       W(3) =  V
       X(4) =  0.0D+00
       Y(4) = -A
       Z(4) =  0.0D+00
       W(4) =  V
       X(5) =  0.0D+00
       Y(5) =  0.0D+00
       Z(5) =  A
       W(5) =  V
       X(6) =  0.0D+00
       Y(6) =  0.0D+00
       Z(6) = -A
       W(6) =  V
       NUM=NUM+6
       RETURN
CVW
    2  CONTINUE
       A=SQRT(0.5D+00)
       X( 1) =  0.0D+00
       Y( 1) =  A
       Z( 1) =  A
       W( 1) =  V
       X( 2) =  0.0D+00
       Y( 2) = -A
       Z( 2) =  A
       W( 2) =  V
       X( 3) =  0.0D+00
       Y( 3) =  A
       Z( 3) = -A
       W( 3) =  V
       X( 4) =  0.0D+00
       Y( 4) = -A
       Z( 4) = -A
       W( 4) =  V
       X( 5) =  A
       Y( 5) =  0.0D+00
       Z( 5) =  A
       W( 5) =  V
       X( 6) = -A
       Y( 6) =  0.0D+00
       Z( 6) =  A
       W( 6) =  V
       X( 7) =  A
       Y( 7) =  0.0D+00
       Z( 7) = -A
       W( 7) =  V
       X( 8) = -A
       Y( 8) =  0.0D+00
       Z( 8) = -A
       W( 8) =  V
       X( 9) =  A
       Y( 9) =  A
       Z( 9) =  0.0D+00
       W( 9) =  V
       X(10) = -A
       Y(10) =  A
       Z(10) =  0.0D+00
       W(10) =  V
       X(11) =  A
       Y(11) = -A
       Z(11) =  0.0D+00
       W(11) =  V
       X(12) = -A
       Y(12) = -A
       Z(12) =  0.0D+00
       W(12) =  V
       NUM=NUM+12
       RETURN
CVW
    3  CONTINUE
       A = SQRT(1.0D+00/3.0D+00)
       X(1) =  A
       Y(1) =  A
       Z(1) =  A
       W(1) =  V
       X(2) = -A
       Y(2) =  A
       Z(2) =  A
       W(2) =  V
       X(3) =  A
       Y(3) = -A
       Z(3) =  A
       W(3) =  V
       X(4) = -A
       Y(4) = -A
       Z(4) =  A
       W(4) =  V
       X(5) =  A
       Y(5) =  A
       Z(5) = -A
       W(5) =  V
       X(6) = -A
       Y(6) =  A
       Z(6) = -A
       W(6) =  V
       X(7) =  A
       Y(7) = -A
       Z(7) = -A
       W(7) =  V
       X(8) = -A
       Y(8) = -A
       Z(8) = -A
       W(8) =  V
       NUM=NUM+8
       RETURN
CVW
    4  CONTINUE
       B = SQRT(1.0D+00 - 2.0D+00*A*A)
       X( 1) =  A
       Y( 1) =  A
       Z( 1) =  B
       W( 1) =  V
       X( 2) = -A
       Y( 2) =  A
       Z( 2) =  B
       W( 2) =  V
       X( 3) =  A
       Y( 3) = -A
       Z( 3) =  B
       W( 3) =  V
       X( 4) = -A
       Y( 4) = -A
       Z( 4) =  B
       W( 4) =  V
       X( 5) =  A
       Y( 5) =  A
       Z( 5) = -B
       W( 5) =  V
       X( 6) = -A
       Y( 6) =  A
       Z( 6) = -B
       W( 6) =  V
       X( 7) =  A
       Y( 7) = -A
       Z( 7) = -B
       W( 7) =  V
       X( 8) = -A
       Y( 8) = -A
       Z( 8) = -B
       W( 8) =  V
       X( 9) =  A
       Y( 9) =  B
       Z( 9) =  A
       W( 9) =  V
       X(10) = -A
       Y(10) =  B
       Z(10) =  A
       W(10) =  V
       X(11) =  A
       Y(11) = -B
       Z(11) =  A
       W(11) =  V
       X(12) = -A
       Y(12) = -B
       Z(12) =  A
       W(12) =  V
       X(13) =  A
       Y(13) =  B
       Z(13) = -A
       W(13) =  V
       X(14) = -A
       Y(14) =  B
       Z(14) = -A
       W(14) =  V
       X(15) =  A
       Y(15) = -B
       Z(15) = -A
       W(15) =  V
       X(16) = -A
       Y(16) = -B
       Z(16) = -A
       W(16) =  V
       X(17) =  B
       Y(17) =  A
       Z(17) =  A
       W(17) =  V
       X(18) = -B
       Y(18) =  A
       Z(18) =  A
       W(18) =  V
       X(19) =  B
       Y(19) = -A
       Z(19) =  A
       W(19) =  V
       X(20) = -B
       Y(20) = -A
       Z(20) =  A
       W(20) =  V
       X(21) =  B
       Y(21) =  A
       Z(21) = -A
       W(21) =  V
       X(22) = -B
       Y(22) =  A
       Z(22) = -A
       W(22) =  V
       X(23) =  B
       Y(23) = -A
       Z(23) = -A
       W(23) =  V
       X(24) = -B
       Y(24) = -A
       Z(24) = -A
       W(24) =  V
       NUM=NUM+24
       RETURN
CVW
    5  CONTINUE
       B=SQRT(1.0D+00-A*A)
       X( 1) =  A
       Y( 1) =  B
       Z( 1) =  0.0D+00
       W( 1) =  V
       X( 2) = -A
       Y( 2) =  B
       Z( 2) =  0.0D+00
       W( 2) =  V
       X( 3) =  A
       Y( 3) = -B
       Z( 3) =  0.0D+00
       W( 3) =  V
       X( 4) = -A
       Y( 4) = -B
       Z( 4) =  0.0D+00
       W( 4) =  V
       X( 5) =  B
       Y( 5) =  A
       Z( 5) =  0.0D+00
       W( 5) =  V
       X( 6) = -B
       Y( 6) =  A
       Z( 6) =  0.0D+00
       W( 6) =  V
       X( 7) =  B
       Y( 7) = -A
       Z( 7) =  0.0D+00
       W( 7) =  V
       X( 8) = -B
       Y( 8) = -A
       Z( 8) =  0.0D+00
       W( 8) =  V
       X( 9) =  A
       Y( 9) =  0.0D+00
       Z( 9) =  B
       W( 9) =  V
       X(10) = -A
       Y(10) =  0.0D+00
       Z(10) =  B
       W(10) =  V
       X(11) =  A
       Y(11) =  0.0D+00
       Z(11) = -B
       W(11) =  V
       X(12) = -A
       Y(12) =  0.0D+00
       Z(12) = -B
       W(12) =  V
       X(13) =  B
       Y(13) =  0.0D+00
       Z(13) =  A
       W(13) =  V
       X(14) = -B
       Y(14) =  0.0D+00
       Z(14) =  A
       W(14) =  V
       X(15) =  B
       Y(15) =  0.0D+00
       Z(15) = -A
       W(15) =  V
       X(16) = -B
       Y(16) =  0.0D+00
       Z(16) = -A
       W(16) =  V
       X(17) =  0.0D+00
       Y(17) =  A
       Z(17) =  B
       W(17) =  V
       X(18) =  0.0D+00
       Y(18) = -A
       Z(18) =  B
       W(18) =  V
       X(19) =  0.0D+00
       Y(19) =  A
       Z(19) = -B
       W(19) =  V
       X(20) =  0.0D+00
       Y(20) = -A
       Z(20) = -B
       W(20) =  V
       X(21) =  0.0D+00
       Y(21) =  B
       Z(21) =  A
       W(21) =  V
       X(22) =  0.0D+00
       Y(22) = -B
       Z(22) =  A
       W(22) =  V
       X(23) =  0.0D+00
       Y(23) =  B
       Z(23) = -A
       W(23) =  V
       X(24) =  0.0D+00
       Y(24) = -B
       Z(24) = -A
       W(24) =  V
       NUM=NUM+24
       RETURN
CVW
    6  CONTINUE
       C=SQRT(1.0D+00 - A*A - B*B)
       X( 1) =  A
       Y( 1) =  B
       Z( 1) =  C
       W( 1) =  V
       X( 2) = -A
       Y( 2) =  B
       Z( 2) =  C
       W( 2) =  V
       X( 3) =  A
       Y( 3) = -B
       Z( 3) =  C
       W( 3) =  V
       X( 4) = -A
       Y( 4) = -B
       Z( 4) =  C
       W( 4) =  V
       X( 5) =  A
       Y( 5) =  B
       Z( 5) = -C
       W( 5) =  V
       X( 6) = -A
       Y( 6) =  B
       Z( 6) = -C
       W( 6) =  V
       X( 7) =  A
       Y( 7) = -B
       Z( 7) = -C
       W( 7) =  V
       X( 8) = -A
       Y( 8) = -B
       Z( 8) = -C
       W( 8) =  V
       X( 9) =  A
       Y( 9) =  C
       Z( 9) =  B
       W( 9) =  V
       X(10) = -A
       Y(10) =  C
       Z(10) =  B
       W(10) =  V
       X(11) =  A
       Y(11) = -C
       Z(11) =  B
       W(11) =  V
       X(12) = -A
       Y(12) = -C
       Z(12) =  B
       W(12) =  V
       X(13) =  A
       Y(13) =  C
       Z(13) = -B
       W(13) =  V
       X(14) = -A
       Y(14) =  C
       Z(14) = -B
       W(14) =  V
       X(15) =  A
       Y(15) = -C
       Z(15) = -B
       W(15) =  V
       X(16) = -A
       Y(16) = -C
       Z(16) = -B
       W(16) =  V
       X(17) =  B
       Y(17) =  A
       Z(17) =  C
       W(17) =  V
       X(18) = -B
       Y(18) =  A
       Z(18) =  C
       W(18) =  V
       X(19) =  B
       Y(19) = -A
       Z(19) =  C
       W(19) =  V
       X(20) = -B
       Y(20) = -A
       Z(20) =  C
       W(20) =  V
       X(21) =  B
       Y(21) =  A
       Z(21) = -C
       W(21) =  V
       X(22) = -B
       Y(22) =  A
       Z(22) = -C
       W(22) =  V
       X(23) =  B
       Y(23) = -A
       Z(23) = -C
       W(23) =  V
       X(24) = -B
       Y(24) = -A
       Z(24) = -C
       W(24) =  V
       X(25) =  B
       Y(25) =  C
       Z(25) =  A
       W(25) =  V
       X(26) = -B
       Y(26) =  C
       Z(26) =  A
       W(26) =  V
       X(27) =  B
       Y(27) = -C
       Z(27) =  A
       W(27) =  V
       X(28) = -B
       Y(28) = -C
       Z(28) =  A
       W(28) =  V
       X(29) =  B
       Y(29) =  C
       Z(29) = -A
       W(29) =  V
       X(30) = -B
       Y(30) =  C
       Z(30) = -A
       W(30) =  V
       X(31) =  B
       Y(31) = -C
       Z(31) = -A
       W(31) =  V
       X(32) = -B
       Y(32) = -C
       Z(32) = -A
       W(32) =  V
       X(33) =  C
       Y(33) =  A
       Z(33) =  B
       W(33) =  V
       X(34) = -C
       Y(34) =  A
       Z(34) =  B
       W(34) =  V
       X(35) =  C
       Y(35) = -A
       Z(35) =  B
       W(35) =  V
       X(36) = -C
       Y(36) = -A
       Z(36) =  B
       W(36) =  V
       X(37) =  C
       Y(37) =  A
       Z(37) = -B
       W(37) =  V
       X(38) = -C
       Y(38) =  A
       Z(38) = -B
       W(38) =  V
       X(39) =  C
       Y(39) = -A
       Z(39) = -B
       W(39) =  V
       X(40) = -C
       Y(40) = -A
       Z(40) = -B
       W(40) =  V
       X(41) =  C
       Y(41) =  B
       Z(41) =  A
       W(41) =  V
       X(42) = -C
       Y(42) =  B
       Z(42) =  A
       W(42) =  V
       X(43) =  C
       Y(43) = -B
       Z(43) =  A
       W(43) =  V
       X(44) = -C
       Y(44) = -B
       Z(44) =  A
       W(44) =  V
       X(45) =  C
       Y(45) =  B
       Z(45) = -A
       W(45) =  V
       X(46) = -C
       Y(46) =  B
       Z(46) = -A
       W(46) =  V
       X(47) =  C
       Y(47) = -B
       Z(47) = -A
       W(47) =  V
       X(48) = -C
       Y(48) = -B
       Z(48) = -A
       W(48) =  V
       NUM=NUM+48
       RETURN
CMV
    7  CONTINUE
       C=SQRT(1.0D+00 - A*A - B*B)
       X( 1) =  A
       X( 2) =  A
       X( 3) = -A
       X( 4) = -A
       X( 5) =  C
       X( 6) =  C
       X( 7) = -C
       X( 8) = -C
       X( 9) =  B
       X(10) =  B
       X(11) = -B
       X(12) = -B
C
       Y( 1) =  B
       Y( 2) = -B
       Y( 3) =  B
       Y( 4) = -B
       Y( 5) =  A
       Y( 6) = -A
       Y( 7) =  A
       Y( 8) = -A
       Y( 9) =  C
       Y(10) = -C
       Y(11) =  C
       Y(12) = -C
C
       Z( 1) =  C
       Z( 2) = -C
       Z( 3) = -C
       Z( 4) =  C
       Z( 5) =  B
       Z( 6) = -B
       Z( 7) = -B
       Z( 8) =  B
       Z( 9) =  A
       Z(10) = -A
       Z(11) = -A
       Z(12) =  A
C
       X(13) =  A
       X(14) =  A
       X(15) = -A
       X(16) = -A
       X(17) = -B
       X(18) = -B
       X(19) =  B
       X(20) =  B
       X(21) =  C
       X(22) =  C
       X(23) = -C
       X(24) = -C
C
       Y(13) =  C
       Y(14) = -C
       Y(15) =  C
       Y(16) = -C
       Y(17) =  A
       Y(18) = -A
       Y(19) =  A
       Y(20) = -A
       Y(21) = -B
       Y(22) =  B
       Y(23) = -B
       Y(24) =  B
C
       Z(13) = -B
       Z(14) =  B
       Z(15) =  B
       Z(16) = -B
       Z(17) =  C
       Z(18) = -C
       Z(19) = -C
       Z(20) =  C
       Z(21) =  A
       Z(22) = -A
       Z(23) = -A
       Z(24) =  A
C
       W( 1) =  V
       W( 2) =  V
       W( 3) =  V
       W( 4) =  V
       W( 5) =  V
       W( 6) =  V
       W( 7) =  V
       W( 8) =  V
       W( 9) =  V
       W(10) =  V
       W(11) =  V
       W(12) =  V
       W(13) =  V
       W(14) =  V
       W(15) =  V
       W(16) =  V
       W(17) =  V
       W(18) =  V
       W(19) =  V
       W(20) =  V
       W(21) =  V
       W(22) =  V
       W(23) =  V
       W(24) =  V
C
       NUM=NUM+24
       RETURN
       END
       SUBROUTINE GEN_IH(CODE, NUM, X, Y, Z, W, A, B, V)
       IMPLICIT LOGICAL(A-Z)
       DOUBLE PRECISION X(*),Y(*),Z(*),W(*)
       DOUBLE PRECISION A,B,V
       DOUBLE PRECISION C
       INTEGER CODE
       INTEGER NUM
       INTEGER N
       DOUBLE PRECISION G, H, T
       PARAMETER (G=(sqrt(5.0D0)+1.0D0)*0.25D0)
       PARAMETER (H=(sqrt(5.0D0)-1.0D0)*0.25D0)
       PARAMETER (T=0.5D0)
       INTEGER IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA
       COMMON/IOFILE/IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
       N = 1
       GOTO (1,2,3,4) CODE
       WRITE (IW,*) 'GEN_IH: INVALID CODE'
    1  CONTINUE
       A = sqrt((5.0D0+sqrt(5.0D0))/10.0D0)
       B = sqrt((5.0D0-sqrt(5.0D0))/10.0D0)
       C = 0.0D0
       CALL GEN_TD( 1, N, X(N), Y(N), Z(N), W(N), A, B, C, V)
       NUM = NUM + 12
       RETURN
    2  CONTINUE
       A = sqrt((3.0D0-sqrt(5.0D0))/6.0D0)
       B = sqrt((3.0D0+sqrt(5.0D0))/6.0D0)
       C = 0.0D0
       CALL GEN_TD( 1, N, X(N), Y(N), Z(N), W(N), A, B, C, V)
       C = 1.0D0/sqrt(3.0D0)
       X(13) =  C
       X(14) =  C
       X(15) =  C
       X(16) =  C
       X(17) = -C
       X(18) = -C
       X(19) = -C
       X(20) = -C

       Y(13) =  C
       Y(14) =  C
       Y(15) = -C
       Y(16) = -C
       Y(17) =  C
       Y(18) =  C
       Y(19) = -C
       Y(20) = -C

       Z(13) =  C
       Z(14) = -C
       Z(15) =  C
       Z(16) = -C
       Z(17) =  C
       Z(18) = -C
       Z(19) =  C
       Z(20) = -C

       W(13) = V
       W(14) = V
       W(15) = V
       W(16) = V
       W(17) = V
       W(18) = V
       W(19) = V
       W(20) = V
       NUM = NUM + 20
       RETURN
    3  CONTINUE
       CALL LEB_GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A = G
       B = H
       C = T
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       NUM = NUM + 30
       RETURN
    4  CONTINUE
       C=SQRT(1.0D+00 - A*A - B*B)
       CALL GEN_TD( 1, N, X(N), Y(N), Z(N), W(N), A, B, C, V)
       CALL IH_TR(A, B, C)
       CALL GEN_TD( 1, N, X(N), Y(N), Z(N), W(N), A, B, C, V)
       CALL IH_TR(A, B, C)
       CALL GEN_TD( 1, N, X(N), Y(N), Z(N), W(N), A, B, C, V)
       CALL IH_TR(A, B, C)
       CALL GEN_TD( 1, N, X(N), Y(N), Z(N), W(N), A, B, C, V)
       CALL IH_TR(A, B, C)
       CALL GEN_TD( 1, N, X(N), Y(N), Z(N), W(N), A, B, C, V)
       NUM = NUM + 60
       RETURN
       END
       SUBROUTINE IH_TR(A, B, C)
       IMPLICIT LOGICAL(A-Z)
       DOUBLE PRECISION A, B, C
       DOUBLE PRECISION A1, B1
       DOUBLE PRECISION G, H, T
       PARAMETER (G=(sqrt(5.0D0)+1.0D0)*0.25D0)
       PARAMETER (H=(sqrt(5.0D0)-1.0D0)*0.25D0)
       PARAMETER (T=0.5D0)
       A1 = G*A + H*B - T*C
       B1 = H*A + T*B + G*C
       C  = T*A - G*B + H*C
       A = A1
       B = B1
       RETURN
       END
       SUBROUTINE GEN_TD(CODE, NUM, X, Y, Z, W, A, B, C, V)
       IMPLICIT LOGICAL(A-Z)
       DOUBLE PRECISION X(*),Y(*),Z(*),W(*)
       DOUBLE PRECISION A,B,C,V
       DOUBLE PRECISION G, H, T
       INTEGER CODE
       INTEGER NUM
       INTEGER IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA
       COMMON/IOFILE/IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
C       C=SQRT(1.0D+00 - A*A - B*B)
       X( 1) =  A
       X( 2) =  A
       X( 3) = -A
       X( 4) = -A
       X( 5) =  C
       X( 6) =  C
       X( 7) = -C
       X( 8) = -C
       X( 9) =  B
       X(10) =  B
       X(11) = -B
       X(12) = -B
C
       Y( 1) =  B
       Y( 2) = -B
       Y( 3) =  B
       Y( 4) = -B
       Y( 5) =  A
       Y( 6) = -A
       Y( 7) =  A
       Y( 8) = -A
       Y( 9) =  C
       Y(10) = -C
       Y(11) =  C
       Y(12) = -C
C
       Z( 1) =  C
       Z( 2) = -C
       Z( 3) = -C
       Z( 4) =  C
       Z( 5) =  B
       Z( 6) = -B
       Z( 7) = -B
       Z( 8) =  B
       Z( 9) =  A
       Z(10) = -A
       Z(11) = -A
       Z(12) =  A
C
       W( 1) =  V
       W( 2) =  V
       W( 3) =  V
       W( 4) =  V
       W( 5) =  V
       W( 6) =  V
       W( 7) =  V
       W( 8) =  V
       W( 9) =  V
       W(10) =  V
       W(11) =  V
       W(12) =  V
       NUM=NUM+12
       RETURN
       END
       SUBROUTINE LD0006(X,Y,Z,W,N)
       DOUBLE PRECISION X(   6)
       DOUBLE PRECISION Y(   6)
       DOUBLE PRECISION Z(   6)
       DOUBLE PRECISION W(   6)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV    6-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] V.I. LEBEDEV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
CHVD        ORDER OF ACCURACY"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
CHVD
CHVD   [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
CHVD       "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
CHVD
CHVD   [4] V.I. LEBEDEV
CHVD       "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
CHVD       SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
CHVD
CHVD   [5] V.I. LEBEDEV
CHVD       "QUADRATURES ON A SPHERE"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
CHVD       1976, PP. 10-24.
CHVD
CHVD   [6] V.I. LEBEDEV
CHVD       "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
CHVD        ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
CHVD        OCTAHEDRON GROUP WITH INVERSION"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
CHVD       1975, PP. 44-51.
CHVD
       N=1
       V=0.1666666666666667D+00
       CALL LEB_GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD0014(X,Y,Z,W,N)
       DOUBLE PRECISION X(  14)
       DOUBLE PRECISION Y(  14)
       DOUBLE PRECISION Z(  14)
       DOUBLE PRECISION W(  14)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV   14-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] V.I. LEBEDEV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
CHVD        ORDER OF ACCURACY"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
CHVD
CHVD   [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
CHVD       "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
CHVD
CHVD   [4] V.I. LEBEDEV
CHVD       "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
CHVD       SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
CHVD
CHVD   [5] V.I. LEBEDEV
CHVD       "QUADRATURES ON A SPHERE"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
CHVD       1976, PP. 10-24.
CHVD
CHVD   [6] V.I. LEBEDEV
CHVD       "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
CHVD        ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
CHVD        OCTAHEDRON GROUP WITH INVERSION"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
CHVD       1975, PP. 44-51.
CHVD
       N=1
       V=0.6666666666666667D-01
       CALL LEB_GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.7500000000000000D-01
       CALL LEB_GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD0026(X,Y,Z,W,N)
       DOUBLE PRECISION X(  26)
       DOUBLE PRECISION Y(  26)
       DOUBLE PRECISION Z(  26)
       DOUBLE PRECISION W(  26)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV   26-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] V.I. LEBEDEV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
CHVD        ORDER OF ACCURACY"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
CHVD
CHVD   [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
CHVD       "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
CHVD
CHVD   [4] V.I. LEBEDEV
CHVD       "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
CHVD       SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
CHVD
CHVD   [5] V.I. LEBEDEV
CHVD       "QUADRATURES ON A SPHERE"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
CHVD       1976, PP. 10-24.
CHVD
CHVD   [6] V.I. LEBEDEV
CHVD       "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
CHVD        ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
CHVD        OCTAHEDRON GROUP WITH INVERSION"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
CHVD       1975, PP. 44-51.
CHVD
       N=1
       V=0.4761904761904762D-01
       CALL LEB_GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3809523809523810D-01
       CALL LEB_GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3214285714285714D-01
       CALL LEB_GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD0038(X,Y,Z,W,N)
       DOUBLE PRECISION X(  38)
       DOUBLE PRECISION Y(  38)
       DOUBLE PRECISION Z(  38)
       DOUBLE PRECISION W(  38)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV   38-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] V.I. LEBEDEV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
CHVD        ORDER OF ACCURACY"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
CHVD
CHVD   [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
CHVD       "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
CHVD
CHVD   [4] V.I. LEBEDEV
CHVD       "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
CHVD       SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
CHVD
CHVD   [5] V.I. LEBEDEV
CHVD       "QUADRATURES ON A SPHERE"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
CHVD       1976, PP. 10-24.
CHVD
CHVD   [6] V.I. LEBEDEV
CHVD       "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
CHVD        ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
CHVD        OCTAHEDRON GROUP WITH INVERSION"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
CHVD       1975, PP. 44-51.
CHVD
       N=1
       V=0.9523809523809524D-02
       CALL LEB_GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3214285714285714D-01
       CALL LEB_GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4597008433809831D+00
       V=0.2857142857142857D-01
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD0050(X,Y,Z,W,N)
       DOUBLE PRECISION X(  50)
       DOUBLE PRECISION Y(  50)
       DOUBLE PRECISION Z(  50)
       DOUBLE PRECISION W(  50)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV   50-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] V.I. LEBEDEV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
CHVD        ORDER OF ACCURACY"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
CHVD
CHVD   [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
CHVD       "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
CHVD
CHVD   [4] V.I. LEBEDEV
CHVD       "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
CHVD       SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
CHVD
CHVD   [5] V.I. LEBEDEV
CHVD       "QUADRATURES ON A SPHERE"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
CHVD       1976, PP. 10-24.
CHVD
CHVD   [6] V.I. LEBEDEV
CHVD       "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
CHVD        ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
CHVD        OCTAHEDRON GROUP WITH INVERSION"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
CHVD       1975, PP. 44-51.
CHVD
       N=1
       V=0.1269841269841270D-01
       CALL LEB_GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2257495590828924D-01
       CALL LEB_GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2109375000000000D-01
       CALL LEB_GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3015113445777636D+00
       V=0.2017333553791887D-01
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD0074(X,Y,Z,W,N)
       DOUBLE PRECISION X(  74)
       DOUBLE PRECISION Y(  74)
       DOUBLE PRECISION Z(  74)
       DOUBLE PRECISION W(  74)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV   74-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] V.I. LEBEDEV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
CHVD        ORDER OF ACCURACY"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
CHVD
CHVD   [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
CHVD       "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
CHVD
CHVD   [4] V.I. LEBEDEV
CHVD       "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
CHVD       SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
CHVD
CHVD   [5] V.I. LEBEDEV
CHVD       "QUADRATURES ON A SPHERE"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
CHVD       1976, PP. 10-24.
CHVD
CHVD   [6] V.I. LEBEDEV
CHVD       "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
CHVD        ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
CHVD        OCTAHEDRON GROUP WITH INVERSION"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
CHVD       1975, PP. 44-51.
CHVD
       N=1
       V=0.5130671797338464D-03
       CALL LEB_GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.1660406956574204D-01
       CALL LEB_GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=-0.2958603896103896D-01
       CALL LEB_GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4803844614152614D+00
       V=0.2657620708215946D-01
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3207726489807764D+00
       V=0.1652217099371571D-01
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD0086(X,Y,Z,W,N)
       DOUBLE PRECISION X(  86)
       DOUBLE PRECISION Y(  86)
       DOUBLE PRECISION Z(  86)
       DOUBLE PRECISION W(  86)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV   86-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] V.I. LEBEDEV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
CHVD        ORDER OF ACCURACY"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
CHVD
CHVD   [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
CHVD       "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
CHVD
CHVD   [4] V.I. LEBEDEV
CHVD       "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
CHVD       SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
CHVD
CHVD   [5] V.I. LEBEDEV
CHVD       "QUADRATURES ON A SPHERE"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
CHVD       1976, PP. 10-24.
CHVD
CHVD   [6] V.I. LEBEDEV
CHVD       "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
CHVD        ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
CHVD        OCTAHEDRON GROUP WITH INVERSION"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
CHVD       1975, PP. 44-51.
CHVD
       N=1
       V=0.1154401154401154D-01
       CALL LEB_GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.1194390908585628D-01
       CALL LEB_GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3696028464541502D+00
       V=0.1111055571060340D-01
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6943540066026664D+00
       V=0.1187650129453714D-01
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3742430390903412D+00
       V=0.1181230374690448D-01
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD0110(X,Y,Z,W,N)
       DOUBLE PRECISION X( 110)
       DOUBLE PRECISION Y( 110)
       DOUBLE PRECISION Z( 110)
       DOUBLE PRECISION W( 110)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV  110-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] V.I. LEBEDEV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
CHVD        ORDER OF ACCURACY"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
CHVD
CHVD   [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
CHVD       "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
CHVD
CHVD   [4] V.I. LEBEDEV
CHVD       "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
CHVD       SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
CHVD
CHVD   [5] V.I. LEBEDEV
CHVD       "QUADRATURES ON A SPHERE"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
CHVD       1976, PP. 10-24.
CHVD
CHVD   [6] V.I. LEBEDEV
CHVD       "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
CHVD        ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
CHVD        OCTAHEDRON GROUP WITH INVERSION"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
CHVD       1975, PP. 44-51.
CHVD
       N=1
       V=0.3828270494937162D-02
       CALL LEB_GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.9793737512487512D-02
       CALL LEB_GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1851156353447362D+00
       V=0.8211737283191111D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6904210483822922D+00
       V=0.9942814891178103D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3956894730559419D+00
       V=0.9595471336070963D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4783690288121502D+00
       V=0.9694996361663028D-02
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD0146(X,Y,Z,W,N)
       DOUBLE PRECISION X( 146)
       DOUBLE PRECISION Y( 146)
       DOUBLE PRECISION Z( 146)
       DOUBLE PRECISION W( 146)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV  146-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] V.I. LEBEDEV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
CHVD        ORDER OF ACCURACY"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
CHVD
CHVD   [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
CHVD       "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
CHVD
CHVD   [4] V.I. LEBEDEV
CHVD       "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
CHVD       SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
CHVD
CHVD   [5] V.I. LEBEDEV
CHVD       "QUADRATURES ON A SPHERE"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
CHVD       1976, PP. 10-24.
CHVD
CHVD   [6] V.I. LEBEDEV
CHVD       "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
CHVD        ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
CHVD        OCTAHEDRON GROUP WITH INVERSION"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
CHVD       1975, PP. 44-51.
CHVD
       N=1
       V=0.5996313688621381D-03
       CALL LEB_GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.7372999718620756D-02
       CALL LEB_GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.7210515360144488D-02
       CALL LEB_GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6764410400114264D+00
       V=0.7116355493117555D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4174961227965453D+00
       V=0.6753829486314477D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1574676672039082D+00
       V=0.7574394159054034D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1403553811713183D+00
       B=0.4493328323269557D+00
       V=0.6991087353303262D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD0170(X,Y,Z,W,N)
       DOUBLE PRECISION X( 170)
       DOUBLE PRECISION Y( 170)
       DOUBLE PRECISION Z( 170)
       DOUBLE PRECISION W( 170)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV  170-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] V.I. LEBEDEV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
CHVD        ORDER OF ACCURACY"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
CHVD
CHVD   [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
CHVD       "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
CHVD
CHVD   [4] V.I. LEBEDEV
CHVD       "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
CHVD       SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
CHVD
CHVD   [5] V.I. LEBEDEV
CHVD       "QUADRATURES ON A SPHERE"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
CHVD       1976, PP. 10-24.
CHVD
CHVD   [6] V.I. LEBEDEV
CHVD       "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
CHVD        ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
CHVD        OCTAHEDRON GROUP WITH INVERSION"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
CHVD       1975, PP. 44-51.
CHVD
       N=1
       V=0.5544842902037365D-02
       CALL LEB_GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.6071332770670752D-02
       CALL LEB_GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.6383674773515093D-02
       CALL LEB_GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2551252621114134D+00
       V=0.5183387587747790D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6743601460362766D+00
       V=0.6317929009813725D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4318910696719410D+00
       V=0.6201670006589077D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2613931360335988D+00
       V=0.5477143385137348D-02
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4990453161796037D+00
       B=0.1446630744325115D+00
       V=0.5968383987681156D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD0194(X,Y,Z,W,N)
       DOUBLE PRECISION X( 194)
       DOUBLE PRECISION Y( 194)
       DOUBLE PRECISION Z( 194)
       DOUBLE PRECISION W( 194)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV  194-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] V.I. LEBEDEV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
CHVD        ORDER OF ACCURACY"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
CHVD
CHVD   [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
CHVD       "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
CHVD
CHVD   [4] V.I. LEBEDEV
CHVD       "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
CHVD       SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
CHVD
CHVD   [5] V.I. LEBEDEV
CHVD       "QUADRATURES ON A SPHERE"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
CHVD       1976, PP. 10-24.
CHVD
CHVD   [6] V.I. LEBEDEV
CHVD       "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
CHVD        ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
CHVD        OCTAHEDRON GROUP WITH INVERSION"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
CHVD       1975, PP. 44-51.
CHVD
       N=1
       V=0.1782340447244611D-02
       CALL LEB_GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.5716905949977102D-02
       CALL LEB_GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.5573383178848738D-02
       CALL LEB_GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6712973442695226D+00
       V=0.5608704082587997D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2892465627575439D+00
       V=0.5158237711805383D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4446933178717437D+00
       V=0.5518771467273614D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1299335447650067D+00
       V=0.4106777028169394D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3457702197611283D+00
       V=0.5051846064614808D-02
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1590417105383530D+00
       B=0.8360360154824589D+00
       V=0.5530248916233094D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD0230(X,Y,Z,W,N)
       DOUBLE PRECISION X( 230)
       DOUBLE PRECISION Y( 230)
       DOUBLE PRECISION Z( 230)
       DOUBLE PRECISION W( 230)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV  230-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] V.I. LEBEDEV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
CHVD        ORDER OF ACCURACY"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
CHVD
CHVD   [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
CHVD       "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
CHVD
CHVD   [4] V.I. LEBEDEV
CHVD       "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
CHVD       SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
CHVD
CHVD   [5] V.I. LEBEDEV
CHVD       "QUADRATURES ON A SPHERE"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
CHVD       1976, PP. 10-24.
CHVD
CHVD   [6] V.I. LEBEDEV
CHVD       "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
CHVD        ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
CHVD        OCTAHEDRON GROUP WITH INVERSION"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
CHVD       1975, PP. 44-51.
CHVD
       N=1
       V=-0.5522639919727325D-01
       CALL LEB_GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.4450274607445226D-02
       CALL LEB_GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4492044687397611D+00
       V=0.4496841067921404D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2520419490210201D+00
       V=0.5049153450478750D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6981906658447242D+00
       V=0.3976408018051883D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6587405243460960D+00
       V=0.4401400650381014D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4038544050097660D-01
       V=0.1724544350544401D-01
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5823842309715585D+00
       V=0.4231083095357343D-02
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3545877390518688D+00
       V=0.5198069864064399D-02
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2272181808998187D+00
       B=0.4864661535886647D+00
       V=0.4695720972568883D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD0266(X,Y,Z,W,N)
       DOUBLE PRECISION X( 266)
       DOUBLE PRECISION Y( 266)
       DOUBLE PRECISION Z( 266)
       DOUBLE PRECISION W( 266)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV  266-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] V.I. LEBEDEV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
CHVD        ORDER OF ACCURACY"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
CHVD
CHVD   [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
CHVD       "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
CHVD
CHVD   [4] V.I. LEBEDEV
CHVD       "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
CHVD       SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
CHVD
CHVD   [5] V.I. LEBEDEV
CHVD       "QUADRATURES ON A SPHERE"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
CHVD       1976, PP. 10-24.
CHVD
CHVD   [6] V.I. LEBEDEV
CHVD       "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
CHVD        ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
CHVD        OCTAHEDRON GROUP WITH INVERSION"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
CHVD       1975, PP. 44-51.
CHVD
       N=1
       V=-0.1313769127326952D-02
       CALL LEB_GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=-0.2522728704859336D-02
       CALL LEB_GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.4186853881700583D-02
       CALL LEB_GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7039373391585475D+00
       V=0.5315167977810885D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1012526248572414D+00
       V=0.4047142377086219D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4647448726420539D+00
       V=0.4112482394406990D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3277420654971629D+00
       V=0.3595584899758782D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6620338663699974D+00
       V=0.4256131351428158D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8506508083520399D+00
       V=0.4229582700647240D-02
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3233484542692899D+00
       B=0.1153112011009701D+00
       V=0.4080914225780505D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2314790158712601D+00
       B=0.5244939240922365D+00
       V=0.4071467593830964D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD0302(X,Y,Z,W,N)
       DOUBLE PRECISION X( 302)
       DOUBLE PRECISION Y( 302)
       DOUBLE PRECISION Z( 302)
       DOUBLE PRECISION W( 302)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV  302-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] V.I. LEBEDEV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
CHVD        ORDER OF ACCURACY"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
CHVD
CHVD   [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
CHVD       "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
CHVD
CHVD   [4] V.I. LEBEDEV
CHVD       "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
CHVD       SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
CHVD
CHVD   [5] V.I. LEBEDEV
CHVD       "QUADRATURES ON A SPHERE"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
CHVD       1976, PP. 10-24.
CHVD
CHVD   [6] V.I. LEBEDEV
CHVD       "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
CHVD        ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
CHVD        OCTAHEDRON GROUP WITH INVERSION"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
CHVD       1975, PP. 44-51.
CHVD
       N=1
       V=0.8545911725128148D-03
       CALL LEB_GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3599119285025571D-02
       CALL LEB_GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3515640345570105D+00
       V=0.3449788424305883D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6566329410219612D+00
       V=0.3604822601419882D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4729054132581005D+00
       V=0.3576729661743367D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9618308522614784D-01
       V=0.2352101413689164D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2219645236294178D+00
       V=0.3108953122413675D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7011766416089545D+00
       V=0.3650045807677255D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2644152887060663D+00
       V=0.2982344963171804D-02
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5718955891878961D+00
       V=0.3600820932216460D-02
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2510034751770465D+00
       B=0.8000727494073952D+00
       V=0.3571540554273387D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1233548532583327D+00
       B=0.4127724083168531D+00
       V=0.3392312205006170D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD0350(X,Y,Z,W,N)
       DOUBLE PRECISION X( 350)
       DOUBLE PRECISION Y( 350)
       DOUBLE PRECISION Z( 350)
       DOUBLE PRECISION W( 350)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV  350-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] V.I. LEBEDEV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
CHVD        ORDER OF ACCURACY"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
CHVD
CHVD   [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
CHVD       "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
CHVD
CHVD   [4] V.I. LEBEDEV
CHVD       "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
CHVD       SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
CHVD
CHVD   [5] V.I. LEBEDEV
CHVD       "QUADRATURES ON A SPHERE"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
CHVD       1976, PP. 10-24.
CHVD
CHVD   [6] V.I. LEBEDEV
CHVD       "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
CHVD        ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
CHVD        OCTAHEDRON GROUP WITH INVERSION"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
CHVD       1975, PP. 44-51.
CHVD
       N=1
       V=0.3006796749453936D-02
       CALL LEB_GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3050627745650771D-02
       CALL LEB_GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7068965463912316D+00
       V=0.1621104600288991D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4794682625712025D+00
       V=0.3005701484901752D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1927533154878019D+00
       V=0.2990992529653774D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6930357961327123D+00
       V=0.2982170644107595D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3608302115520091D+00
       V=0.2721564237310992D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6498486161496169D+00
       V=0.3033513795811141D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1932945013230339D+00
       V=0.3007949555218533D-02
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3800494919899303D+00
       V=0.2881964603055307D-02
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2899558825499574D+00
       B=0.7934537856582316D+00
       V=0.2958357626535696D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9684121455103957D-01
       B=0.8280801506686862D+00
       V=0.3036020026407088D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1833434647041659D+00
       B=0.9074658265305127D+00
       V=0.2832187403926303D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD0434(X,Y,Z,W,N)
       DOUBLE PRECISION X( 434)
       DOUBLE PRECISION Y( 434)
       DOUBLE PRECISION Z( 434)
       DOUBLE PRECISION W( 434)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV  434-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] V.I. LEBEDEV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
CHVD        ORDER OF ACCURACY"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
CHVD
CHVD   [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
CHVD       "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
CHVD
CHVD   [4] V.I. LEBEDEV
CHVD       "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
CHVD       SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
CHVD
CHVD   [5] V.I. LEBEDEV
CHVD       "QUADRATURES ON A SPHERE"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
CHVD       1976, PP. 10-24.
CHVD
CHVD   [6] V.I. LEBEDEV
CHVD       "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
CHVD        ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
CHVD        OCTAHEDRON GROUP WITH INVERSION"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
CHVD       1975, PP. 44-51.
CHVD
       N=1
       V=0.5265897968224436D-03
       CALL LEB_GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2548219972002607D-02
       CALL LEB_GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2512317418927307D-02
       CALL LEB_GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6909346307509111D+00
       V=0.2530403801186355D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1774836054609158D+00
       V=0.2014279020918528D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4914342637784746D+00
       V=0.2501725168402936D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6456664707424256D+00
       V=0.2513267174597564D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2861289010307638D+00
       V=0.2302694782227416D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7568084367178018D-01
       V=0.1462495621594614D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3927259763368002D+00
       V=0.2445373437312980D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8818132877794288D+00
       V=0.2417442375638981D-02
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9776428111182649D+00
       V=0.1910951282179532D-02
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2054823696403044D+00
       B=0.8689460322872412D+00
       V=0.2416930044324775D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5905157048925271D+00
       B=0.7999278543857286D+00
       V=0.2512236854563495D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5550152361076807D+00
       B=0.7717462626915901D+00
       V=0.2496644054553086D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9371809858553722D+00
       B=0.3344363145343455D+00
       V=0.2236607760437849D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD0590(X,Y,Z,W,N)
       DOUBLE PRECISION X( 590)
       DOUBLE PRECISION Y( 590)
       DOUBLE PRECISION Z( 590)
       DOUBLE PRECISION W( 590)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV  590-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] V.I. LEBEDEV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
CHVD        ORDER OF ACCURACY"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
CHVD
CHVD   [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
CHVD       "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
CHVD
CHVD   [4] V.I. LEBEDEV
CHVD       "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
CHVD       SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
CHVD
CHVD   [5] V.I. LEBEDEV
CHVD       "QUADRATURES ON A SPHERE"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
CHVD       1976, PP. 10-24.
CHVD
CHVD   [6] V.I. LEBEDEV
CHVD       "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
CHVD        ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
CHVD        OCTAHEDRON GROUP WITH INVERSION"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
CHVD       1975, PP. 44-51.
CHVD
       N=1
       V=0.3095121295306187D-03
       CALL LEB_GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.1852379698597489D-02
       CALL LEB_GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7040954938227469D+00
       V=0.1871790639277744D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6807744066455243D+00
       V=0.1858812585438317D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6372546939258752D+00
       V=0.1852028828296213D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5044419707800358D+00
       V=0.1846715956151242D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4215761784010967D+00
       V=0.1818471778162769D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3317920736472123D+00
       V=0.1749564657281154D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2384736701421887D+00
       V=0.1617210647254411D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1459036449157763D+00
       V=0.1384737234851692D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6095034115507196D-01
       V=0.9764331165051050D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6116843442009876D+00
       V=0.1857161196774078D-02
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3964755348199858D+00
       V=0.1705153996395864D-02
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1724782009907724D+00
       V=0.1300321685886048D-02
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5610263808622060D+00
       B=0.3518280927733519D+00
       V=0.1842866472905286D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4742392842551980D+00
       B=0.2634716655937950D+00
       V=0.1802658934377451D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5984126497885380D+00
       B=0.1816640840360209D+00
       V=0.1849830560443660D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3791035407695563D+00
       B=0.1720795225656878D+00
       V=0.1713904507106709D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2778673190586244D+00
       B=0.8213021581932511D-01
       V=0.1555213603396808D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5033564271075117D+00
       B=0.8999205842074875D-01
       V=0.1802239128008525D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD0770(X,Y,Z,W,N)
       DOUBLE PRECISION X( 770)
       DOUBLE PRECISION Y( 770)
       DOUBLE PRECISION Z( 770)
       DOUBLE PRECISION W( 770)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV  770-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] V.I. LEBEDEV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
CHVD        ORDER OF ACCURACY"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
CHVD
CHVD   [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
CHVD       "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
CHVD
CHVD   [4] V.I. LEBEDEV
CHVD       "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
CHVD       SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
CHVD
CHVD   [5] V.I. LEBEDEV
CHVD       "QUADRATURES ON A SPHERE"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
CHVD       1976, PP. 10-24.
CHVD
CHVD   [6] V.I. LEBEDEV
CHVD       "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
CHVD        ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
CHVD        OCTAHEDRON GROUP WITH INVERSION"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
CHVD       1975, PP. 44-51.
CHVD
       N=1
       V=0.2192942088181184D-03
       CALL LEB_GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.1436433617319080D-02
       CALL LEB_GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.1421940344335877D-02
       CALL LEB_GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5087204410502360D-01
       V=0.6798123511050502D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1228198790178831D+00
       V=0.9913184235294912D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2026890814408786D+00
       V=0.1180207833238949D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2847745156464294D+00
       V=0.1296599602080921D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3656719078978026D+00
       V=0.1365871427428316D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4428264886713469D+00
       V=0.1402988604775325D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5140619627249735D+00
       V=0.1418645563595609D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6306401219166803D+00
       V=0.1421376741851662D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6716883332022612D+00
       V=0.1423996475490962D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6979792685336881D+00
       V=0.1431554042178567D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1446865674195309D+00
       V=0.9254401499865368D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3390263475411216D+00
       V=0.1250239995053509D-02
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5335804651263506D+00
       V=0.1394365843329230D-02
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6944024393349413D-01
       B=0.2355187894242326D+00
       V=0.1127089094671749D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2269004109529460D+00
       B=0.4102182474045730D+00
       V=0.1345753760910670D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8025574607775339D-01
       B=0.6214302417481605D+00
       V=0.1424957283316783D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1467999527896572D+00
       B=0.3245284345717394D+00
       V=0.1261523341237750D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1571507769824727D+00
       B=0.5224482189696630D+00
       V=0.1392547106052696D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2365702993157246D+00
       B=0.6017546634089558D+00
       V=0.1418761677877656D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7714815866765732D-01
       B=0.4346575516141163D+00
       V=0.1338366684479554D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3062936666210730D+00
       B=0.4908826589037616D+00
       V=0.1393700862676131D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3822477379524787D+00
       B=0.5648768149099500D+00
       V=0.1415914757466932D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD0974(X,Y,Z,W,N)
       DOUBLE PRECISION X( 974)
       DOUBLE PRECISION Y( 974)
       DOUBLE PRECISION Z( 974)
       DOUBLE PRECISION W( 974)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV  974-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] V.I. LEBEDEV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
CHVD        ORDER OF ACCURACY"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
CHVD
CHVD   [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
CHVD       "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
CHVD
CHVD   [4] V.I. LEBEDEV
CHVD       "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
CHVD       SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
CHVD
CHVD   [5] V.I. LEBEDEV
CHVD       "QUADRATURES ON A SPHERE"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
CHVD       1976, PP. 10-24.
CHVD
CHVD   [6] V.I. LEBEDEV
CHVD       "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
CHVD        ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
CHVD        OCTAHEDRON GROUP WITH INVERSION"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
CHVD       1975, PP. 44-51.
CHVD
       N=1
       V=0.1438294190527431D-03
       CALL LEB_GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.1125772288287004D-02
       CALL LEB_GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4292963545341347D-01
       V=0.4948029341949241D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1051426854086404D+00
       V=0.7357990109125470D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1750024867623087D+00
       V=0.8889132771304384D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2477653379650257D+00
       V=0.9888347838921435D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3206567123955957D+00
       V=0.1053299681709471D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3916520749849983D+00
       V=0.1092778807014578D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4590825874187624D+00
       V=0.1114389394063227D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5214563888415861D+00
       V=0.1123724788051555D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6253170244654199D+00
       V=0.1125239325243814D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6637926744523170D+00
       V=0.1126153271815905D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6910410398498301D+00
       V=0.1130286931123841D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7052907007457760D+00
       V=0.1134986534363955D-02
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1236686762657990D+00
       V=0.6823367927109931D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2940777114468387D+00
       V=0.9454158160447096D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4697753849207649D+00
       V=0.1074429975385679D-02
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6334563241139567D+00
       V=0.1129300086569132D-02
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5974048614181342D-01
       B=0.2029128752777523D+00
       V=0.8436884500901954D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1375760408473636D+00
       B=0.4602621942484054D+00
       V=0.1075255720448885D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3391016526336286D+00
       B=0.5030673999662036D+00
       V=0.1108577236864462D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1271675191439820D+00
       B=0.2817606422442134D+00
       V=0.9566475323783357D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2693120740413512D+00
       B=0.4331561291720157D+00
       V=0.1080663250717391D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1419786452601918D+00
       B=0.6256167358580814D+00
       V=0.1126797131196295D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6709284600738255D-01
       B=0.3798395216859157D+00
       V=0.1022568715358061D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7057738183256172D-01
       B=0.5517505421423520D+00
       V=0.1108960267713108D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2783888477882155D+00
       B=0.6029619156159187D+00
       V=0.1122790653435766D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1979578938917407D+00
       B=0.3589606329589096D+00
       V=0.1032401847117460D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2087307061103274D+00
       B=0.5348666438135476D+00
       V=0.1107249382283854D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4055122137872836D+00
       B=0.5674997546074373D+00
       V=0.1121780048519972D-02
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD1202(X,Y,Z,W,N)
       DOUBLE PRECISION X(1202)
       DOUBLE PRECISION Y(1202)
       DOUBLE PRECISION Z(1202)
       DOUBLE PRECISION W(1202)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV 1202-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] V.I. LEBEDEV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
CHVD        ORDER OF ACCURACY"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
CHVD
CHVD   [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
CHVD       "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
CHVD
CHVD   [4] V.I. LEBEDEV
CHVD       "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
CHVD       SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
CHVD
CHVD   [5] V.I. LEBEDEV
CHVD       "QUADRATURES ON A SPHERE"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
CHVD       1976, PP. 10-24.
CHVD
CHVD   [6] V.I. LEBEDEV
CHVD       "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
CHVD        ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
CHVD        OCTAHEDRON GROUP WITH INVERSION"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
CHVD       1975, PP. 44-51.
CHVD
       N=1
       V=0.1105189233267572D-03
       CALL LEB_GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.9205232738090741D-03
       CALL LEB_GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.9133159786443561D-03
       CALL LEB_GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3712636449657089D-01
       V=0.3690421898017899D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9140060412262223D-01
       V=0.5603990928680660D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1531077852469906D+00
       V=0.6865297629282609D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2180928891660612D+00
       V=0.7720338551145630D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2839874532200175D+00
       V=0.8301545958894795D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3491177600963764D+00
       V=0.8686692550179628D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4121431461444309D+00
       V=0.8927076285846890D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4718993627149127D+00
       V=0.9060820238568219D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5273145452842337D+00
       V=0.9119777254940867D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6209475332444019D+00
       V=0.9128720138604181D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6569722711857291D+00
       V=0.9130714935691735D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6841788309070143D+00
       V=0.9152873784554116D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7012604330123631D+00
       V=0.9187436274321654D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1072382215478166D+00
       V=0.5176977312965694D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2582068959496968D+00
       V=0.7331143682101417D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4172752955306717D+00
       V=0.8463232836379928D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5700366911792503D+00
       V=0.9031122694253992D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9827986018263947D+00
       B=0.1771774022615325D+00
       V=0.6485778453163257D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9624249230326228D+00
       B=0.2475716463426288D+00
       V=0.7435030910982369D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9402007994128811D+00
       B=0.3354616289066489D+00
       V=0.7998527891839054D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9320822040143202D+00
       B=0.3173615246611977D+00
       V=0.8101731497468018D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9043674199393299D+00
       B=0.4090268427085357D+00
       V=0.8483389574594331D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8912407560074747D+00
       B=0.3854291150669224D+00
       V=0.8556299257311812D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8676435628462708D+00
       B=0.4932221184851285D+00
       V=0.8803208679738260D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8581979986041619D+00
       B=0.4785320675922435D+00
       V=0.8811048182425720D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8396753624049856D+00
       B=0.4507422593157064D+00
       V=0.8850282341265444D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8165288564022188D+00
       B=0.5632123020762100D+00
       V=0.9021342299040653D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8015469370783529D+00
       B=0.5434303569693900D+00
       V=0.9010091677105086D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7773563069070351D+00
       B=0.5123518486419871D+00
       V=0.9022692938426915D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7661621213900394D+00
       B=0.6394279634749102D+00
       V=0.9158016174693465D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7553584143533510D+00
       B=0.6269805509024392D+00
       V=0.9131578003189435D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7344305757559503D+00
       B=0.6031161693096310D+00
       V=0.9107813579482705D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7043837184021765D+00
       B=0.5693702498468441D+00
       V=0.9105760258970126D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD1454(X,Y,Z,W,N)
       DOUBLE PRECISION X(1454)
       DOUBLE PRECISION Y(1454)
       DOUBLE PRECISION Z(1454)
       DOUBLE PRECISION W(1454)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV 1454-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] V.I. LEBEDEV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
CHVD        ORDER OF ACCURACY"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
CHVD
CHVD   [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
CHVD       "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
CHVD
CHVD   [4] V.I. LEBEDEV
CHVD       "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
CHVD       SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
CHVD
CHVD   [5] V.I. LEBEDEV
CHVD       "QUADRATURES ON A SPHERE"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
CHVD       1976, PP. 10-24.
CHVD
CHVD   [6] V.I. LEBEDEV
CHVD       "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
CHVD        ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
CHVD        OCTAHEDRON GROUP WITH INVERSION"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
CHVD       1975, PP. 44-51.
CHVD
       N=1
       V=0.7777160743261247D-04
       CALL LEB_GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.7557646413004701D-03
       CALL LEB_GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3229290663413854D-01
       V=0.2841633806090617D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8036733271462222D-01
       V=0.4374419127053555D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1354289960531653D+00
       V=0.5417174740872172D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1938963861114426D+00
       V=0.6148000891358593D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2537343715011275D+00
       V=0.6664394485800705D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3135251434752570D+00
       V=0.7025039356923220D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3721558339375338D+00
       V=0.7268511789249627D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4286809575195696D+00
       V=0.7422637534208629D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4822510128282994D+00
       V=0.7509545035841214D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5320679333566263D+00
       V=0.7548535057718401D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6172998195394274D+00
       V=0.7554088969774001D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6510679849127481D+00
       V=0.7553147174442808D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6777315251687360D+00
       V=0.7564767653292297D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6963109410648741D+00
       V=0.7587991808518730D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7058935009831749D+00
       V=0.7608261832033027D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9955546194091857D+00
       V=0.4021680447874916D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9734115901794209D+00
       V=0.5804871793945964D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9275693732388626D+00
       V=0.6792151955945159D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8568022422795103D+00
       V=0.7336741211286294D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7623495553719372D+00
       V=0.7581866300989608D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5707522908892223D+00
       B=0.4387028039889501D+00
       V=0.7538257859800743D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5196463388403083D+00
       B=0.3858908414762617D+00
       V=0.7483517247053123D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4646337531215351D+00
       B=0.3301937372343854D+00
       V=0.7371763661112059D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4063901697557691D+00
       B=0.2725423573563777D+00
       V=0.7183448895756934D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3456329466643087D+00
       B=0.2139510237495250D+00
       V=0.6895815529822191D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2831395121050332D+00
       B=0.1555922309786647D+00
       V=0.6480105801792886D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2197682022925330D+00
       B=0.9892878979686097D-01
       V=0.5897558896594636D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1564696098650355D+00
       B=0.4598642910675510D-01
       V=0.5095708849247346D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6027356673721295D+00
       B=0.3376625140173426D+00
       V=0.7536906428909755D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5496032320255096D+00
       B=0.2822301309727988D+00
       V=0.7472505965575118D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4921707755234567D+00
       B=0.2248632342592540D+00
       V=0.7343017132279698D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4309422998598483D+00
       B=0.1666224723456479D+00
       V=0.7130871582177445D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3664108182313672D+00
       B=0.1086964901822169D+00
       V=0.6817022032112776D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2990189057758436D+00
       B=0.5251989784120085D-01
       V=0.6380941145604121D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6268724013144998D+00
       B=0.2297523657550023D+00
       V=0.7550381377920310D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5707324144834607D+00
       B=0.1723080607093800D+00
       V=0.7478646640144802D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5096360901960365D+00
       B=0.1140238465390513D+00
       V=0.7335918720601220D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4438729938312456D+00
       B=0.5611522095882537D-01
       V=0.7110120527658118D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6419978471082389D+00
       B=0.1164174423140873D+00
       V=0.7571363978689501D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5817218061802611D+00
       B=0.5797589531445219D-01
       V=0.7489908329079234D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD1730(X,Y,Z,W,N)
       DOUBLE PRECISION X(1730)
       DOUBLE PRECISION Y(1730)
       DOUBLE PRECISION Z(1730)
       DOUBLE PRECISION W(1730)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV 1730-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] V.I. LEBEDEV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
CHVD        ORDER OF ACCURACY"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
CHVD
CHVD   [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
CHVD       "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
CHVD
CHVD   [4] V.I. LEBEDEV
CHVD       "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
CHVD       SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
CHVD
CHVD   [5] V.I. LEBEDEV
CHVD       "QUADRATURES ON A SPHERE"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
CHVD       1976, PP. 10-24.
CHVD
CHVD   [6] V.I. LEBEDEV
CHVD       "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
CHVD        ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
CHVD        OCTAHEDRON GROUP WITH INVERSION"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
CHVD       1975, PP. 44-51.
CHVD
       N=1
       V=0.6309049437420976D-04
       CALL LEB_GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.6398287705571748D-03
       CALL LEB_GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.6357185073530720D-03
       CALL LEB_GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2860923126194662D-01
       V=0.2221207162188168D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7142556767711522D-01
       V=0.3475784022286848D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1209199540995559D+00
       V=0.4350742443589804D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1738673106594379D+00
       V=0.4978569136522127D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2284645438467734D+00
       V=0.5435036221998053D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2834807671701512D+00
       V=0.5765913388219542D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3379680145467339D+00
       V=0.6001200359226003D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3911355454819537D+00
       V=0.6162178172717512D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4422860353001403D+00
       V=0.6265218152438485D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4907781568726057D+00
       V=0.6323987160974212D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5360006153211468D+00
       V=0.6350767851540569D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6142105973596603D+00
       V=0.6354362775297107D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6459300387977504D+00
       V=0.6352302462706235D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6718056125089225D+00
       V=0.6358117881417972D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6910888533186254D+00
       V=0.6373101590310117D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7030467416823252D+00
       V=0.6390428961368665D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8354951166354646D-01
       V=0.3186913449946576D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2050143009099486D+00
       V=0.4678028558591711D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3370208290706637D+00
       V=0.5538829697598626D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4689051484233963D+00
       V=0.6044475907190476D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5939400424557334D+00
       V=0.6313575103509012D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1394983311832261D+00
       B=0.4097581162050343D-01
       V=0.4078626431855630D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1967999180485014D+00
       B=0.8851987391293348D-01
       V=0.4759933057812725D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2546183732548967D+00
       B=0.1397680182969819D+00
       V=0.5268151186413440D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3121281074713875D+00
       B=0.1929452542226526D+00
       V=0.5643048560507316D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3685981078502492D+00
       B=0.2467898337061562D+00
       V=0.5914501076613073D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4233760321547856D+00
       B=0.3003104124785409D+00
       V=0.6104561257874195D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4758671236059246D+00
       B=0.3526684328175033D+00
       V=0.6230252860707806D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5255178579796463D+00
       B=0.4031134861145713D+00
       V=0.6305618761760796D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5718025633734589D+00
       B=0.4509426448342351D+00
       V=0.6343092767597889D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2686927772723415D+00
       B=0.4711322502423248D-01
       V=0.5176268945737826D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3306006819904809D+00
       B=0.9784487303942695D-01
       V=0.5564840313313692D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3904906850594983D+00
       B=0.1505395810025273D+00
       V=0.5856426671038980D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4479957951904390D+00
       B=0.2039728156296050D+00
       V=0.6066386925777091D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5027076848919780D+00
       B=0.2571529941121107D+00
       V=0.6208824962234458D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5542087392260217D+00
       B=0.3092191375815670D+00
       V=0.6296314297822907D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6020850887375187D+00
       B=0.3593807506130276D+00
       V=0.6340423756791859D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4019851409179594D+00
       B=0.5063389934378671D-01
       V=0.5829627677107342D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4635614567449800D+00
       B=0.1032422269160612D+00
       V=0.6048693376081110D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5215860931591575D+00
       B=0.1566322094006254D+00
       V=0.6202362317732461D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5758202499099271D+00
       B=0.2098082827491099D+00
       V=0.6299005328403779D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6259893683876795D+00
       B=0.2618824114553391D+00
       V=0.6347722390609353D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5313795124811891D+00
       B=0.5263245019338556D-01
       V=0.6203778981238834D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5893317955931995D+00
       B=0.1061059730982005D+00
       V=0.6308414671239979D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6426246321215801D+00
       B=0.1594171564034221D+00
       V=0.6362706466959498D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6511904367376113D+00
       B=0.5354789536565540D-01
       V=0.6375414170333233D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD2030(X,Y,Z,W,N)
       DOUBLE PRECISION X(2030)
       DOUBLE PRECISION Y(2030)
       DOUBLE PRECISION Z(2030)
       DOUBLE PRECISION W(2030)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV 2030-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] V.I. LEBEDEV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
CHVD        ORDER OF ACCURACY"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
CHVD
CHVD   [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
CHVD       "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
CHVD
CHVD   [4] V.I. LEBEDEV
CHVD       "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
CHVD       SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
CHVD
CHVD   [5] V.I. LEBEDEV
CHVD       "QUADRATURES ON A SPHERE"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
CHVD       1976, PP. 10-24.
CHVD
CHVD   [6] V.I. LEBEDEV
CHVD       "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
CHVD        ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
CHVD        OCTAHEDRON GROUP WITH INVERSION"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
CHVD       1975, PP. 44-51.
CHVD
       N=1
       V=0.4656031899197431D-04
       CALL LEB_GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.5421549195295507D-03
       CALL LEB_GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2540835336814348D-01
       V=0.1778522133346553D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6399322800504915D-01
       V=0.2811325405682796D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1088269469804125D+00
       V=0.3548896312631459D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1570670798818287D+00
       V=0.4090310897173364D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2071163932282514D+00
       V=0.4493286134169965D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2578914044450844D+00
       V=0.4793728447962723D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3085687558169623D+00
       V=0.5015415319164265D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3584719706267024D+00
       V=0.5175127372677937D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4070135594428709D+00
       V=0.5285522262081019D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4536618626222638D+00
       V=0.5356832703713962D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4979195686463577D+00
       V=0.5397914736175170D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5393075111126999D+00
       V=0.5416899441599930D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6115617676843916D+00
       V=0.5419308476889938D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6414308435160159D+00
       V=0.5416936902030596D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6664099412721607D+00
       V=0.5419544338703164D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6859161771214913D+00
       V=0.5428983656630975D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6993625593503890D+00
       V=0.5442286500098193D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7062393387719380D+00
       V=0.5452250345057301D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7479028168349763D-01
       V=0.2568002497728530D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1848951153969366D+00
       V=0.3827211700292145D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3059529066581305D+00
       V=0.4579491561917824D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4285556101021362D+00
       V=0.5042003969083574D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5468758653496526D+00
       V=0.5312708889976025D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6565821978343439D+00
       V=0.5438401790747117D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1253901572367117D+00
       B=0.3681917226439641D-01
       V=0.3316041873197344D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1775721510383941D+00
       B=0.7982487607213301D-01
       V=0.3899113567153771D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2305693358216114D+00
       B=0.1264640966592335D+00
       V=0.4343343327201309D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2836502845992063D+00
       B=0.1751585683418957D+00
       V=0.4679415262318919D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3361794746232590D+00
       B=0.2247995907632670D+00
       V=0.4930847981631031D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3875979172264824D+00
       B=0.2745299257422246D+00
       V=0.5115031867540091D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4374019316999074D+00
       B=0.3236373482441118D+00
       V=0.5245217148457367D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4851275843340022D+00
       B=0.3714967859436741D+00
       V=0.5332041499895321D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5303391803806868D+00
       B=0.4175353646321745D+00
       V=0.5384583126021542D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5726197380596287D+00
       B=0.4612084406355461D+00
       V=0.5411067210798852D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2431520732564863D+00
       B=0.4258040133043952D-01
       V=0.4259797391468714D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3002096800895869D+00
       B=0.8869424306722721D-01
       V=0.4604931368460021D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3558554457457432D+00
       B=0.1368811706510655D+00
       V=0.4871814878255202D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4097782537048887D+00
       B=0.1860739985015033D+00
       V=0.5072242910074885D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4616337666067458D+00
       B=0.2354235077395853D+00
       V=0.5217069845235350D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5110707008417874D+00
       B=0.2842074921347011D+00
       V=0.5315785966280310D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5577415286163795D+00
       B=0.3317784414984102D+00
       V=0.5376833708758905D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6013060431366950D+00
       B=0.3775299002040700D+00
       V=0.5408032092069521D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3661596767261781D+00
       B=0.4599367887164592D-01
       V=0.4842744917904866D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4237633153506581D+00
       B=0.9404893773654421D-01
       V=0.5048926076188130D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4786328454658452D+00
       B=0.1431377109091971D+00
       V=0.5202607980478373D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5305702076789774D+00
       B=0.1924186388843570D+00
       V=0.5309932388325743D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5793436224231788D+00
       B=0.2411590944775190D+00
       V=0.5377419770895208D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6247069017094747D+00
       B=0.2886871491583605D+00
       V=0.5411696331677717D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4874315552535204D+00
       B=0.4804978774953206D-01
       V=0.5197996293282420D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5427337322059053D+00
       B=0.9716857199366665D-01
       V=0.5311120836622945D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5943493747246700D+00
       B=0.1465205839795055D+00
       V=0.5384309319956951D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6421314033564943D+00
       B=0.1953579449803574D+00
       V=0.5421859504051886D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6020628374713980D+00
       B=0.4916375015738108D-01
       V=0.5390948355046314D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6529222529856881D+00
       B=0.9861621540127005D-01
       V=0.5433312705027845D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD2354(X,Y,Z,W,N)
       DOUBLE PRECISION X(2354)
       DOUBLE PRECISION Y(2354)
       DOUBLE PRECISION Z(2354)
       DOUBLE PRECISION W(2354)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV 2354-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] V.I. LEBEDEV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
CHVD        ORDER OF ACCURACY"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
CHVD
CHVD   [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
CHVD       "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
CHVD
CHVD   [4] V.I. LEBEDEV
CHVD       "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
CHVD       SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
CHVD
CHVD   [5] V.I. LEBEDEV
CHVD       "QUADRATURES ON A SPHERE"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
CHVD       1976, PP. 10-24.
CHVD
CHVD   [6] V.I. LEBEDEV
CHVD       "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
CHVD        ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
CHVD        OCTAHEDRON GROUP WITH INVERSION"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
CHVD       1975, PP. 44-51.
CHVD
       N=1
       V=0.3922616270665292D-04
       CALL LEB_GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.4703831750854424D-03
       CALL LEB_GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.4678202801282136D-03
       CALL LEB_GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2290024646530589D-01
       V=0.1437832228979900D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5779086652271284D-01
       V=0.2303572493577644D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9863103576375984D-01
       V=0.2933110752447454D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1428155792982185D+00
       V=0.3402905998359838D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1888978116601463D+00
       V=0.3759138466870372D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2359091682970210D+00
       V=0.4030638447899798D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2831228833706171D+00
       V=0.4236591432242211D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3299495857966693D+00
       V=0.4390522656946746D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3758840802660796D+00
       V=0.4502523466626247D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4204751831009480D+00
       V=0.4580577727783541D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4633068518751051D+00
       V=0.4631391616615899D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5039849474507313D+00
       V=0.4660928953698676D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5421265793440747D+00
       V=0.4674751807936953D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6092660230557310D+00
       V=0.4676414903932920D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6374654204984869D+00
       V=0.4674086492347870D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6615136472609892D+00
       V=0.4674928539483207D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6809487285958127D+00
       V=0.4680748979686447D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6952980021665196D+00
       V=0.4690449806389040D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7041245497695400D+00
       V=0.4699877075860818D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6744033088306065D-01
       V=0.2099942281069176D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1678684485334166D+00
       V=0.3172269150712804D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2793559049539613D+00
       V=0.3832051358546523D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3935264218057639D+00
       V=0.4252193818146985D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5052629268232558D+00
       V=0.4513807963755000D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6107905315437531D+00
       V=0.4657797469114178D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1135081039843524D+00
       B=0.3331954884662588D-01
       V=0.2733362800522836D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1612866626099378D+00
       B=0.7247167465436538D-01
       V=0.3235485368463559D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2100786550168205D+00
       B=0.1151539110849745D+00
       V=0.3624908726013453D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2592282009459942D+00
       B=0.1599491097143677D+00
       V=0.3925540070712828D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3081740561320203D+00
       B=0.2058699956028027D+00
       V=0.4156129781116235D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3564289781578164D+00
       B=0.2521624953502911D+00
       V=0.4330644984623263D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4035587288240703D+00
       B=0.2982090785797674D+00
       V=0.4459677725921312D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4491671196373903D+00
       B=0.3434762087235733D+00
       V=0.4551593004456795D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4928854782917489D+00
       B=0.3874831357203437D+00
       V=0.4613341462749918D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5343646791958988D+00
       B=0.4297814821746926D+00
       V=0.4651019618269806D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5732683216530990D+00
       B=0.4699402260943537D+00
       V=0.4670249536100625D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2214131583218986D+00
       B=0.3873602040643895D-01
       V=0.3549555576441708D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2741796504750071D+00
       B=0.8089496256902013D-01
       V=0.3856108245249010D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3259797439149485D+00
       B=0.1251732177620872D+00
       V=0.4098622845756882D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3765441148826891D+00
       B=0.1706260286403185D+00
       V=0.4286328604268950D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4255773574530558D+00
       B=0.2165115147300408D+00
       V=0.4427802198993945D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4727795117058430D+00
       B=0.2622089812225259D+00
       V=0.4530473511488561D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5178546895819012D+00
       B=0.3071721431296201D+00
       V=0.4600805475703138D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5605141192097460D+00
       B=0.3508998998801138D+00
       V=0.4644599059958017D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6004763319352512D+00
       B=0.3929160876166931D+00
       V=0.4667274455712508D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3352842634946949D+00
       B=0.4202563457288019D-01
       V=0.4069360518020356D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3891971629814670D+00
       B=0.8614309758870850D-01
       V=0.4260442819919195D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4409875565542281D+00
       B=0.1314500879380001D+00
       V=0.4408678508029063D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4904893058592484D+00
       B=0.1772189657383859D+00
       V=0.4518748115548597D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5375056138769549D+00
       B=0.2228277110050294D+00
       V=0.4595564875375116D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5818255708669969D+00
       B=0.2677179935014386D+00
       V=0.4643988774315846D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6232334858144959D+00
       B=0.3113675035544165D+00
       V=0.4668827491646946D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4489485354492058D+00
       B=0.4409162378368174D-01
       V=0.4400541823741973D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5015136875933150D+00
       B=0.8939009917748489D-01
       V=0.4514512890193797D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5511300550512623D+00
       B=0.1351806029383365D+00
       V=0.4596198627347549D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5976720409858000D+00
       B=0.1808370355053196D+00
       V=0.4648659016801781D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6409956378989354D+00
       B=0.2257852192301602D+00
       V=0.4675502017157673D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5581222330827514D+00
       B=0.4532173421637160D-01
       V=0.4598494476455523D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6074705984161695D+00
       B=0.9117488031840314D-01
       V=0.4654916955152048D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6532272537379033D+00
       B=0.1369294213140155D+00
       V=0.4684709779505137D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6594761494500487D+00
       B=0.4589901487275583D-01
       V=0.4691445539106986D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD2702(X,Y,Z,W,N)
       DOUBLE PRECISION X(2702)
       DOUBLE PRECISION Y(2702)
       DOUBLE PRECISION Z(2702)
       DOUBLE PRECISION W(2702)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV 2702-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] V.I. LEBEDEV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
CHVD        ORDER OF ACCURACY"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
CHVD
CHVD   [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
CHVD       "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
CHVD
CHVD   [4] V.I. LEBEDEV
CHVD       "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
CHVD       SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
CHVD
CHVD   [5] V.I. LEBEDEV
CHVD       "QUADRATURES ON A SPHERE"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
CHVD       1976, PP. 10-24.
CHVD
CHVD   [6] V.I. LEBEDEV
CHVD       "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
CHVD        ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
CHVD        OCTAHEDRON GROUP WITH INVERSION"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
CHVD       1975, PP. 44-51.
CHVD
       N=1
       V=0.2998675149888161D-04
       CALL LEB_GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.4077860529495355D-03
       CALL LEB_GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2065562538818703D-01
       V=0.1185349192520667D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5250918173022379D-01
       V=0.1913408643425751D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8993480082038376D-01
       V=0.2452886577209897D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1306023924436019D+00
       V=0.2862408183288702D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1732060388531418D+00
       V=0.3178032258257357D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2168727084820249D+00
       V=0.3422945667633690D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2609528309173586D+00
       V=0.3612790520235922D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3049252927938952D+00
       V=0.3758638229818521D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3483484138084404D+00
       V=0.3868711798859953D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3908321549106406D+00
       V=0.3949429933189938D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4320210071894814D+00
       V=0.4006068107541156D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4715824795890053D+00
       V=0.4043192149672723D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5091984794078453D+00
       V=0.4064947495808078D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5445580145650803D+00
       V=0.4075245619813152D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6072575796841768D+00
       V=0.4076423540893566D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6339484505755803D+00
       V=0.4074280862251555D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6570718257486958D+00
       V=0.4074163756012244D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6762557330090709D+00
       V=0.4077647795071246D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6911161696923790D+00
       V=0.4084517552782530D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7012841911659961D+00
       V=0.4092468459224052D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7064559272410020D+00
       V=0.4097872687240906D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6123554989894765D-01
       V=0.1738986811745028D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1533070348312393D+00
       V=0.2659616045280191D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2563902605244206D+00
       V=0.3240596008171533D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3629346991663361D+00
       V=0.3621195964432943D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4683949968987538D+00
       V=0.3868838330760539D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5694479240657952D+00
       V=0.4018911532693111D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6634465430993955D+00
       V=0.4089929432983252D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1033958573552305D+00
       B=0.3034544009063584D-01
       V=0.2279907527706409D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1473521412414395D+00
       B=0.6618803044247135D-01
       V=0.2715205490578897D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1924552158705967D+00
       B=0.1054431128987715D+00
       V=0.3057917896703976D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2381094362890328D+00
       B=0.1468263551238858D+00
       V=0.3326913052452555D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2838121707936760D+00
       B=0.1894486108187886D+00
       V=0.3537334711890037D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3291323133373415D+00
       B=0.2326374238761579D+00
       V=0.3700567500783129D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3736896978741460D+00
       B=0.2758485808485768D+00
       V=0.3825245372589122D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4171406040760013D+00
       B=0.3186179331996921D+00
       V=0.3918125171518296D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4591677985256915D+00
       B=0.3605329796303794D+00
       V=0.3984720419937579D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4994733831718418D+00
       B=0.4012147253586509D+00
       V=0.4029746003338211D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5377731830445096D+00
       B=0.4403050025570692D+00
       V=0.4057428632156627D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5737917830001331D+00
       B=0.4774565904277483D+00
       V=0.4071719274114857D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2027323586271389D+00
       B=0.3544122504976147D-01
       V=0.2990236950664119D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2516942375187273D+00
       B=0.7418304388646328D-01
       V=0.3262951734212878D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3000227995257181D+00
       B=0.1150502745727186D+00
       V=0.3482634608242413D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3474806691046342D+00
       B=0.1571963371209364D+00
       V=0.3656596681700892D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3938103180359209D+00
       B=0.1999631877247100D+00
       V=0.3791740467794218D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4387519590455703D+00
       B=0.2428073457846535D+00
       V=0.3894034450156905D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4820503960077787D+00
       B=0.2852575132906155D+00
       V=0.3968600245508371D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5234573778475101D+00
       B=0.3268884208674639D+00
       V=0.4019931351420050D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5627318647235282D+00
       B=0.3673033321675939D+00
       V=0.4052108801278599D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5996390607156954D+00
       B=0.4061211551830290D+00
       V=0.4068978613940934D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3084780753791947D+00
       B=0.3860125523100059D-01
       V=0.3454275351319704D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3589988275920223D+00
       B=0.7928938987104867D-01
       V=0.3629963537007920D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4078628415881973D+00
       B=0.1212614643030087D+00
       V=0.3770187233889873D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4549287258889735D+00
       B=0.1638770827382693D+00
       V=0.3878608613694378D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5000278512957279D+00
       B=0.2065965798260176D+00
       V=0.3959065270221274D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5429785044928199D+00
       B=0.2489436378852235D+00
       V=0.4015286975463570D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5835939850491711D+00
       B=0.2904811368946891D+00
       V=0.4050866785614717D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6216870353444856D+00
       B=0.3307941957666609D+00
       V=0.4069320185051913D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4151104662709091D+00
       B=0.4064829146052554D-01
       V=0.3760120964062763D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4649804275009218D+00
       B=0.8258424547294755D-01
       V=0.3870969564418064D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5124695757009662D+00
       B=0.1251841962027289D+00
       V=0.3955287790534055D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5574711100606224D+00
       B=0.1679107505976331D+00
       V=0.4015361911302668D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5998597333287227D+00
       B=0.2102805057358715D+00
       V=0.4053836986719548D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6395007148516600D+00
       B=0.2518418087774107D+00
       V=0.4073578673299117D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5188456224746252D+00
       B=0.4194321676077518D-01
       V=0.3954628379231406D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5664190707942778D+00
       B=0.8457661551921499D-01
       V=0.4017645508847530D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6110464353283153D+00
       B=0.1273652932519396D+00
       V=0.4059030348651293D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6526430302051563D+00
       B=0.1698173239076354D+00
       V=0.4080565809484880D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6167551880377548D+00
       B=0.4266398851548864D-01
       V=0.4063018753664651D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6607195418355383D+00
       B=0.8551925814238349D-01
       V=0.4087191292799671D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD3074(X,Y,Z,W,N)
       DOUBLE PRECISION X(3074)
       DOUBLE PRECISION Y(3074)
       DOUBLE PRECISION Z(3074)
       DOUBLE PRECISION W(3074)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV 3074-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] V.I. LEBEDEV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
CHVD        ORDER OF ACCURACY"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
CHVD
CHVD   [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
CHVD       "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
CHVD
CHVD   [4] V.I. LEBEDEV
CHVD       "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
CHVD       SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
CHVD
CHVD   [5] V.I. LEBEDEV
CHVD       "QUADRATURES ON A SPHERE"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
CHVD       1976, PP. 10-24.
CHVD
CHVD   [6] V.I. LEBEDEV
CHVD       "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
CHVD        ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
CHVD        OCTAHEDRON GROUP WITH INVERSION"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
CHVD       1975, PP. 44-51.
CHVD
       N=1
       V=0.2599095953754734D-04
       CALL LEB_GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3603134089687541D-03
       CALL LEB_GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3586067974412447D-03
       CALL LEB_GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1886108518723392D-01
       V=0.9831528474385880D-04
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4800217244625303D-01
       V=0.1605023107954450D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8244922058397242D-01
       V=0.2072200131464099D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1200408362484023D+00
       V=0.2431297618814187D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1595773530809965D+00
       V=0.2711819064496707D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2002635973434064D+00
       V=0.2932762038321116D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2415127590139982D+00
       V=0.3107032514197368D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2828584158458477D+00
       V=0.3243808058921213D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3239091015338138D+00
       V=0.3349899091374030D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3643225097962194D+00
       V=0.3430580688505218D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4037897083691802D+00
       V=0.3490124109290343D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4420247515194127D+00
       V=0.3532148948561955D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4787572538464938D+00
       V=0.3559862669062833D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5137265251275234D+00
       V=0.3576224317551411D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5466764056654611D+00
       V=0.3584050533086076D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6054859420813535D+00
       V=0.3584903581373224D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6308106701764562D+00
       V=0.3582991879040586D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6530369230179584D+00
       V=0.3582371187963125D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6718609524611158D+00
       V=0.3584353631122350D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6869676499894013D+00
       V=0.3589120166517785D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6980467077240748D+00
       V=0.3595445704531601D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7048241721250522D+00
       V=0.3600943557111074D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5591105222058232D-01
       V=0.1456447096742039D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1407384078513916D+00
       V=0.2252370188283782D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2364035438976309D+00
       V=0.2766135443474897D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3360602737818170D+00
       V=0.3110729491500851D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4356292630054665D+00
       V=0.3342506712303391D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5321569415256174D+00
       V=0.3491981834026860D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6232956305040554D+00
       V=0.3576003604348932D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9469870086838469D-01
       B=0.2778748387309470D-01
       V=0.1921921305788564D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1353170300568141D+00
       B=0.6076569878628364D-01
       V=0.2301458216495632D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1771679481726077D+00
       B=0.9703072762711040D-01
       V=0.2604248549522893D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2197066664231751D+00
       B=0.1354112458524762D+00
       V=0.2845275425870697D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2624783557374927D+00
       B=0.1750996479744100D+00
       V=0.3036870897974840D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3050969521214442D+00
       B=0.2154896907449802D+00
       V=0.3188414832298066D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3472252637196021D+00
       B=0.2560954625740152D+00
       V=0.3307046414722089D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3885610219026360D+00
       B=0.2965070050624096D+00
       V=0.3398330969031360D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4288273776062765D+00
       B=0.3363641488734497D+00
       V=0.3466757899705373D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4677662471302948D+00
       B=0.3753400029836788D+00
       V=0.3516095923230054D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5051333589553359D+00
       B=0.4131297522144286D+00
       V=0.3549645184048486D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5406942145810492D+00
       B=0.4494423776081795D+00
       V=0.3570415969441392D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5742204122576457D+00
       B=0.4839938958841502D+00
       V=0.3581251798496118D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1865407027225188D+00
       B=0.3259144851070796D-01
       V=0.2543491329913348D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2321186453689432D+00
       B=0.6835679505297343D-01
       V=0.2786711051330776D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2773159142523882D+00
       B=0.1062284864451989D+00
       V=0.2985552361083679D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3219200192237254D+00
       B=0.1454404409323047D+00
       V=0.3145867929154039D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3657032593944029D+00
       B=0.1854018282582510D+00
       V=0.3273290662067609D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4084376778363622D+00
       B=0.2256297412014750D+00
       V=0.3372705511943501D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4499004945751427D+00
       B=0.2657104425000896D+00
       V=0.3448274437851510D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4898758141326335D+00
       B=0.3052755487631557D+00
       V=0.3503592783048583D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5281547442266309D+00
       B=0.3439863920645423D+00
       V=0.3541854792663162D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5645346989813992D+00
       B=0.3815229456121914D+00
       V=0.3565995517909428D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5988181252159848D+00
       B=0.4175752420966734D+00
       V=0.3578802078302898D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2850425424471603D+00
       B=0.3562149509862536D-01
       V=0.2958644592860982D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3324619433027876D+00
       B=0.7330318886871096D-01
       V=0.3119548129116835D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3785848333076282D+00
       B=0.1123226296008472D+00
       V=0.3250745225005984D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4232891028562115D+00
       B=0.1521084193337708D+00
       V=0.3355153415935208D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4664287050829722D+00
       B=0.1921844459223610D+00
       V=0.3435847568549328D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5078458493735726D+00
       B=0.2321360989678303D+00
       V=0.3495786831622488D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5473779816204180D+00
       B=0.2715886486360520D+00
       V=0.3537767805534621D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5848617133811376D+00
       B=0.3101924707571355D+00
       V=0.3564459815421428D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6201348281584888D+00
       B=0.3476121052890973D+00
       V=0.3578464061225468D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3852191185387871D+00
       B=0.3763224880035108D-01
       V=0.3239748762836212D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4325025061073423D+00
       B=0.7659581935637135D-01
       V=0.3345491784174287D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4778486229734490D+00
       B=0.1163381306083900D+00
       V=0.3429126177301782D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5211663693009000D+00
       B=0.1563890598752899D+00
       V=0.3492420343097421D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5623469504853703D+00
       B=0.1963320810149200D+00
       V=0.3537399050235257D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6012718188659246D+00
       B=0.2357847407258738D+00
       V=0.3566209152659172D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6378179206390117D+00
       B=0.2743846121244060D+00
       V=0.3581084321919782D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4836936460214534D+00
       B=0.3895902610739024D-01
       V=0.3426522117591512D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5293792562683797D+00
       B=0.7871246819312640D-01
       V=0.3491848770121379D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5726281253100033D+00
       B=0.1187963808202981D+00
       V=0.3539318235231476D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6133658776169068D+00
       B=0.1587914708061787D+00
       V=0.3570231438458694D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6515085491865307D+00
       B=0.1983058575227646D+00
       V=0.3586207335051714D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5778692716064976D+00
       B=0.3977209689791542D-01
       V=0.3541196205164025D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6207904288086192D+00
       B=0.7990157592981152D-01
       V=0.3574296911573953D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6608688171046802D+00
       B=0.1199671308754309D+00
       V=0.3591993279818963D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6656263089489130D+00
       B=0.4015955957805969D-01
       V=0.3595855034661997D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD3470(X,Y,Z,W,N)
       DOUBLE PRECISION X(3470)
       DOUBLE PRECISION Y(3470)
       DOUBLE PRECISION Z(3470)
       DOUBLE PRECISION W(3470)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV 3470-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] V.I. LEBEDEV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
CHVD        ORDER OF ACCURACY"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
CHVD
CHVD   [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
CHVD       "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
CHVD
CHVD   [4] V.I. LEBEDEV
CHVD       "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
CHVD       SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
CHVD
CHVD   [5] V.I. LEBEDEV
CHVD       "QUADRATURES ON A SPHERE"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
CHVD       1976, PP. 10-24.
CHVD
CHVD   [6] V.I. LEBEDEV
CHVD       "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
CHVD        ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
CHVD        OCTAHEDRON GROUP WITH INVERSION"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
CHVD       1975, PP. 44-51.
CHVD
       N=1
       V=0.2040382730826330D-04
       CALL LEB_GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3178149703889544D-03
       CALL LEB_GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1721420832906233D-01
       V=0.8288115128076110D-04
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4408875374981770D-01
       V=0.1360883192522954D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7594680813878681D-01
       V=0.1766854454542662D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1108335359204799D+00
       V=0.2083153161230153D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1476517054388567D+00
       V=0.2333279544657158D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1856731870860615D+00
       V=0.2532809539930247D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2243634099428821D+00
       V=0.2692472184211158D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2633006881662727D+00
       V=0.2819949946811885D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3021340904916283D+00
       V=0.2920953593973030D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3405594048030089D+00
       V=0.2999889782948352D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3783044434007372D+00
       V=0.3060292120496902D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4151194767407910D+00
       V=0.3105109167522192D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4507705766443257D+00
       V=0.3136902387550312D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4850346056573187D+00
       V=0.3157984652454632D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5176950817792470D+00
       V=0.3170516518425422D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5485384240820989D+00
       V=0.3176568425633755D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6039117238943308D+00
       V=0.3177198411207062D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6279956655573113D+00
       V=0.3175519492394733D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6493636169568952D+00
       V=0.3174654952634756D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6677644117704504D+00
       V=0.3175676415467654D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6829368572115624D+00
       V=0.3178923417835410D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6946195818184121D+00
       V=0.3183788287531909D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7025711542057026D+00
       V=0.3188755151918807D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7066004767140119D+00
       V=0.3191916889313849D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5132537689946062D-01
       V=0.1231779611744508D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1297994661331225D+00
       V=0.1924661373839880D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2188852049401307D+00
       V=0.2380881867403424D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3123174824903457D+00
       V=0.2693100663037885D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4064037620738195D+00
       V=0.2908673382834366D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4984958396944782D+00
       V=0.3053914619381535D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5864975046021365D+00
       V=0.3143916684147777D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6686711634580175D+00
       V=0.3187042244055363D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8715738780835950D-01
       B=0.2557175233367578D-01
       V=0.1635219535869790D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1248383123134007D+00
       B=0.5604823383376681D-01
       V=0.1968109917696070D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1638062693383378D+00
       B=0.8968568601900765D-01
       V=0.2236754342249974D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2035586203373176D+00
       B=0.1254086651976279D+00
       V=0.2453186687017181D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2436798975293774D+00
       B=0.1624780150162012D+00
       V=0.2627551791580541D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2838207507773806D+00
       B=0.2003422342683208D+00
       V=0.2767654860152220D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3236787502217692D+00
       B=0.2385628026255263D+00
       V=0.2879467027765895D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3629849554840691D+00
       B=0.2767731148783578D+00
       V=0.2967639918918702D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4014948081992087D+00
       B=0.3146542308245309D+00
       V=0.3035900684660351D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4389818379260225D+00
       B=0.3519196415895088D+00
       V=0.3087338237298308D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4752331143674377D+00
       B=0.3883050984023654D+00
       V=0.3124608838860167D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5100457318374018D+00
       B=0.4235613423908649D+00
       V=0.3150084294226743D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5432238388954868D+00
       B=0.4574484717196220D+00
       V=0.3165958398598402D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5745758685072442D+00
       B=0.4897311639255524D+00
       V=0.3174320440957372D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1723981437592809D+00
       B=0.3010630597881105D-01
       V=0.2182188909812599D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2149553257844597D+00
       B=0.6326031554204694D-01
       V=0.2399727933921445D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2573256081247422D+00
       B=0.9848566980258631D-01
       V=0.2579796133514652D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2993163751238106D+00
       B=0.1350835952384266D+00
       V=0.2727114052623535D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3407238005148000D+00
       B=0.1725184055442181D+00
       V=0.2846327656281355D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3813454978483264D+00
       B=0.2103559279730725D+00
       V=0.2941491102051334D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4209848104423343D+00
       B=0.2482278774554860D+00
       V=0.3016049492136107D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4594519699996300D+00
       B=0.2858099509982883D+00
       V=0.3072949726175648D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4965640166185930D+00
       B=0.3228075659915428D+00
       V=0.3114768142886460D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5321441655571562D+00
       B=0.3589459907204151D+00
       V=0.3143823673666223D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5660208438582166D+00
       B=0.3939630088864310D+00
       V=0.3162269764661535D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5980264315964364D+00
       B=0.4276029922949089D+00
       V=0.3172164663759821D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2644215852350733D+00
       B=0.3300939429072552D-01
       V=0.2554575398967435D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3090113743443063D+00
       B=0.6803887650078501D-01
       V=0.2701704069135677D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3525871079197808D+00
       B=0.1044326136206709D+00
       V=0.2823693413468940D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3950418005354029D+00
       B=0.1416751597517679D+00
       V=0.2922898463214289D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4362475663430163D+00
       B=0.1793408610504821D+00
       V=0.3001829062162428D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4760661812145854D+00
       B=0.2170630750175722D+00
       V=0.3062890864542953D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5143551042512103D+00
       B=0.2545145157815807D+00
       V=0.3108328279264746D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5509709026935597D+00
       B=0.2913940101706601D+00
       V=0.3140243146201245D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5857711030329428D+00
       B=0.3274169910910705D+00
       V=0.3160638030977130D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6186149917404392D+00
       B=0.3623081329317265D+00
       V=0.3171462882206275D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3586894569557064D+00
       B=0.3497354386450040D-01
       V=0.2812388416031796D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4035266610019441D+00
       B=0.7129736739757095D-01
       V=0.2912137500288045D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4467775312332510D+00
       B=0.1084758620193165D+00
       V=0.2993241256502206D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4883638346608543D+00
       B=0.1460915689241772D+00
       V=0.3057101738983822D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5281908348434601D+00
       B=0.1837790832369980D+00
       V=0.3105319326251432D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5661542687149311D+00
       B=0.2212075390874021D+00
       V=0.3139565514428167D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6021450102031452D+00
       B=0.2580682841160985D+00
       V=0.3161543006806366D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6360520783610050D+00
       B=0.2940656362094121D+00
       V=0.3172985960613294D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4521611065087196D+00
       B=0.3631055365867002D-01
       V=0.2989400336901431D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4959365651560963D+00
       B=0.7348318468484350D-01
       V=0.3054555883947677D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5376815804038283D+00
       B=0.1111087643812648D+00
       V=0.3104764960807702D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5773314480243768D+00
       B=0.1488226085145408D+00
       V=0.3141015825977616D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6148113245575056D+00
       B=0.1862892274135151D+00
       V=0.3164520621159896D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6500407462842380D+00
       B=0.2231909701714456D+00
       V=0.3176652305912204D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5425151448707213D+00
       B=0.3718201306118944D-01
       V=0.3105097161023939D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5841860556907931D+00
       B=0.7483616335067346D-01
       V=0.3143014117890550D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6234632186851500D+00
       B=0.1125990834266120D+00
       V=0.3168172866287200D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6602934551848843D+00
       B=0.1501303813157619D+00
       V=0.3181401865570968D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6278573968375105D+00
       B=0.3767559930245720D-01
       V=0.3170663659156037D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6665611711264577D+00
       B=0.7548443301360158D-01
       V=0.3185447944625510D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD3890(X,Y,Z,W,N)
       DOUBLE PRECISION X(3890)
       DOUBLE PRECISION Y(3890)
       DOUBLE PRECISION Z(3890)
       DOUBLE PRECISION W(3890)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV 3890-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] V.I. LEBEDEV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
CHVD        ORDER OF ACCURACY"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
CHVD
CHVD   [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
CHVD       "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
CHVD
CHVD   [4] V.I. LEBEDEV
CHVD       "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
CHVD       SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
CHVD
CHVD   [5] V.I. LEBEDEV
CHVD       "QUADRATURES ON A SPHERE"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
CHVD       1976, PP. 10-24.
CHVD
CHVD   [6] V.I. LEBEDEV
CHVD       "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
CHVD        ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
CHVD        OCTAHEDRON GROUP WITH INVERSION"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
CHVD       1975, PP. 44-51.
CHVD
       N=1
       V=0.1807395252196920D-04
       CALL LEB_GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2848008782238827D-03
       CALL LEB_GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2836065837530581D-03
       CALL LEB_GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1587876419858352D-01
       V=0.7013149266673816D-04
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4069193593751206D-01
       V=0.1162798021956766D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7025888115257997D-01
       V=0.1518728583972105D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1027495450028704D+00
       V=0.1798796108216934D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1371457730893426D+00
       V=0.2022593385972785D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1727758532671953D+00
       V=0.2203093105575464D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2091492038929037D+00
       V=0.2349294234299855D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2458813281751915D+00
       V=0.2467682058747003D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2826545859450066D+00
       V=0.2563092683572224D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3191957291799622D+00
       V=0.2639253896763318D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3552621469299578D+00
       V=0.2699137479265108D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3906329503406230D+00
       V=0.2745196420166739D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4251028614093031D+00
       V=0.2779529197397593D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4584777520111870D+00
       V=0.2803996086684265D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4905711358710193D+00
       V=0.2820302356715842D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5212011669847385D+00
       V=0.2830056747491068D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5501878488737995D+00
       V=0.2834808950776839D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6025037877479342D+00
       V=0.2835282339078929D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6254572689549016D+00
       V=0.2833819267065800D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6460107179528248D+00
       V=0.2832858336906784D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6639541138154251D+00
       V=0.2833268235451244D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6790688515667495D+00
       V=0.2835432677029253D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6911338580371512D+00
       V=0.2839091722743049D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6999385956126490D+00
       V=0.2843308178875841D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7053037748656896D+00
       V=0.2846703550533846D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4732224387180115D-01
       V=0.1051193406971900D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1202100529326803D+00
       V=0.1657871838796974D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2034304820664855D+00
       V=0.2064648113714232D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2912285643573002D+00
       V=0.2347942745819741D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3802361792726768D+00
       V=0.2547775326597726D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4680598511056146D+00
       V=0.2686876684847025D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5528151052155599D+00
       V=0.2778665755515867D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6329386307803041D+00
       V=0.2830996616782929D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8056516651369069D-01
       B=0.2363454684003124D-01
       V=0.1403063340168372D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1156476077139389D+00
       B=0.5191291632545936D-01
       V=0.1696504125939477D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1520473382760421D+00
       B=0.8322715736994519D-01
       V=0.1935787242745390D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1892986699745931D+00
       B=0.1165855667993712D+00
       V=0.2130614510521968D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2270194446777792D+00
       B=0.1513077167409504D+00
       V=0.2289381265931048D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2648908185093273D+00
       B=0.1868882025807859D+00
       V=0.2418630292816186D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3026389259574136D+00
       B=0.2229277629776224D+00
       V=0.2523400495631193D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3400220296151384D+00
       B=0.2590951840746235D+00
       V=0.2607623973449605D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3768217953335510D+00
       B=0.2951047291750847D+00
       V=0.2674441032689209D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4128372900921884D+00
       B=0.3307019714169930D+00
       V=0.2726432360343356D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4478807131815630D+00
       B=0.3656544101087634D+00
       V=0.2765787685924545D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4817742034089257D+00
       B=0.3997448951939695D+00
       V=0.2794428690642224D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5143472814653344D+00
       B=0.4327667110812024D+00
       V=0.2814099002062895D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5454346213905650D+00
       B=0.4645196123532293D+00
       V=0.2826429531578994D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5748739313170252D+00
       B=0.4948063555703345D+00
       V=0.2832983542550884D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1599598738286342D+00
       B=0.2792357590048985D-01
       V=0.1886695565284976D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1998097412500951D+00
       B=0.5877141038139065D-01
       V=0.2081867882748234D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2396228952566202D+00
       B=0.9164573914691377D-01
       V=0.2245148680600796D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2792228341097746D+00
       B=0.1259049641962687D+00
       V=0.2380370491511872D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3184251107546741D+00
       B=0.1610594823400863D+00
       V=0.2491398041852455D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3570481164426244D+00
       B=0.1967151653460898D+00
       V=0.2581632405881230D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3949164710492144D+00
       B=0.2325404606175168D+00
       V=0.2653965506227417D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4318617293970503D+00
       B=0.2682461141151439D+00
       V=0.2710857216747087D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4677221009931678D+00
       B=0.3035720116011973D+00
       V=0.2754434093903659D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5023417939270955D+00
       B=0.3382781859197439D+00
       V=0.2786579932519380D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5355701836636128D+00
       B=0.3721383065625942D+00
       V=0.2809011080679474D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5672608451328771D+00
       B=0.4049346360466055D+00
       V=0.2823336184560987D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5972704202540162D+00
       B=0.4364538098633802D+00
       V=0.2831101175806309D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2461687022333596D+00
       B=0.3070423166833368D-01
       V=0.2221679970354546D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2881774566286831D+00
       B=0.6338034669281885D-01
       V=0.2356185734270703D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3293963604116978D+00
       B=0.9742862487067941D-01
       V=0.2469228344805590D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3697303822241377D+00
       B=0.1323799532282290D+00
       V=0.2562726348642046D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4090663023135127D+00
       B=0.1678497018129336D+00
       V=0.2638756726753028D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4472819355411712D+00
       B=0.2035095105326114D+00
       V=0.2699311157390862D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4842513377231437D+00
       B=0.2390692566672091D+00
       V=0.2746233268403837D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5198477629962928D+00
       B=0.2742649818076149D+00
       V=0.2781225674454771D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5539453011883145D+00
       B=0.3088503806580094D+00
       V=0.2805881254045684D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5864196762401251D+00
       B=0.3425904245906614D+00
       V=0.2821719877004913D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6171484466668390D+00
       B=0.3752562294789468D+00
       V=0.2830222502333124D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3350337830565727D+00
       B=0.3261589934634747D-01
       V=0.2457995956744870D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3775773224758284D+00
       B=0.6658438928081572D-01
       V=0.2551474407503706D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4188155229848973D+00
       B=0.1014565797157954D+00
       V=0.2629065335195311D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4586805892009344D+00
       B=0.1368573320843822D+00
       V=0.2691900449925075D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4970895714224235D+00
       B=0.1724614851951608D+00
       V=0.2741275485754276D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5339505133960747D+00
       B=0.2079779381416412D+00
       V=0.2778530970122595D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5691665792531440D+00
       B=0.2431385788322288D+00
       V=0.2805010567646741D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6026387682680377D+00
       B=0.2776901883049853D+00
       V=0.2822055834031040D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6342676150163307D+00
       B=0.3113881356386632D+00
       V=0.2831016901243473D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4237951119537067D+00
       B=0.3394877848664351D-01
       V=0.2624474901131803D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4656918683234929D+00
       B=0.6880219556291447D-01
       V=0.2688034163039377D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5058857069185980D+00
       B=0.1041946859721635D+00
       V=0.2738932751287636D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5443204666713996D+00
       B=0.1398039738736393D+00
       V=0.2777944791242523D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5809298813759742D+00
       B=0.1753373381196155D+00
       V=0.2806011661660987D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6156416039447128D+00
       B=0.2105215793514010D+00
       V=0.2824181456597460D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6483801351066604D+00
       B=0.2450953312157051D+00
       V=0.2833585216577828D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5103616577251688D+00
       B=0.3485560643800719D-01
       V=0.2738165236962878D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5506738792580681D+00
       B=0.7026308631512033D-01
       V=0.2778365208203180D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5889573040995292D+00
       B=0.1059035061296403D+00
       V=0.2807852940418966D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6251641589516930D+00
       B=0.1414823925236026D+00
       V=0.2827245949674705D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6592414921570178D+00
       B=0.1767207908214530D+00
       V=0.2837342344829828D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5930314017533384D+00
       B=0.3542189339561672D-01
       V=0.2809233907610981D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6309812253390175D+00
       B=0.7109574040369549D-01
       V=0.2829930809742694D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6666296011353230D+00
       B=0.1067259792282730D+00
       V=0.2841097874111479D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6703715271049922D+00
       B=0.3569455268820809D-01
       V=0.2843455206008783D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD4334(X,Y,Z,W,N)
       DOUBLE PRECISION X(4334)
       DOUBLE PRECISION Y(4334)
       DOUBLE PRECISION Z(4334)
       DOUBLE PRECISION W(4334)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV 4334-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] V.I. LEBEDEV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
CHVD        ORDER OF ACCURACY"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
CHVD
CHVD   [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
CHVD       "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
CHVD
CHVD   [4] V.I. LEBEDEV
CHVD       "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
CHVD       SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
CHVD
CHVD   [5] V.I. LEBEDEV
CHVD       "QUADRATURES ON A SPHERE"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
CHVD       1976, PP. 10-24.
CHVD
CHVD   [6] V.I. LEBEDEV
CHVD       "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
CHVD        ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
CHVD        OCTAHEDRON GROUP WITH INVERSION"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
CHVD       1975, PP. 44-51.
CHVD
       N=1
       V=0.1449063022537883D-04
       CALL LEB_GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2546377329828424D-03
       CALL LEB_GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1462896151831013D-01
       V=0.6018432961087496D-04
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3769840812493139D-01
       V=0.1002286583263673D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6524701904096891D-01
       V=0.1315222931028093D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9560543416134648D-01
       V=0.1564213746876724D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1278335898929198D+00
       V=0.1765118841507736D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1613096104466031D+00
       V=0.1928737099311080D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1955806225745371D+00
       V=0.2062658534263270D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2302935218498028D+00
       V=0.2172395445953787D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2651584344113027D+00
       V=0.2262076188876047D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2999276825183209D+00
       V=0.2334885699462397D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3343828669718798D+00
       V=0.2393355273179203D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3683265013750518D+00
       V=0.2439559200468863D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4015763206518108D+00
       V=0.2475251866060002D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4339612026399770D+00
       V=0.2501965558158773D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4653180651114582D+00
       V=0.2521081407925925D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4954893331080803D+00
       V=0.2533881002388081D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5243207068924930D+00
       V=0.2541582900848261D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5516590479041704D+00
       V=0.2545365737525860D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6012371927804176D+00
       V=0.2545726993066799D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6231574466449819D+00
       V=0.2544456197465555D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6429416514181271D+00
       V=0.2543481596881064D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6604124272943595D+00
       V=0.2543506451429194D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6753851470408250D+00
       V=0.2544905675493763D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6876717970626160D+00
       V=0.2547611407344429D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6970895061319234D+00
       V=0.2551060375448869D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7034746912553310D+00
       V=0.2554291933816039D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7067017217542295D+00
       V=0.2556255710686343D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4382223501131123D-01
       V=0.9041339695118195D-04
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1117474077400006D+00
       V=0.1438426330079022D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1897153252911440D+00
       V=0.1802523089820518D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2724023009910331D+00
       V=0.2060052290565496D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3567163308709902D+00
       V=0.2245002248967466D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4404784483028087D+00
       V=0.2377059847731150D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5219833154161411D+00
       V=0.2468118955882525D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5998179868977553D+00
       V=0.2525410872966528D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6727803154548222D+00
       V=0.2553101409933397D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7476563943166086D-01
       B=0.2193168509461185D-01
       V=0.1212879733668632D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1075341482001416D+00
       B=0.4826419281533887D-01
       V=0.1472872881270931D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1416344885203259D+00
       B=0.7751191883575742D-01
       V=0.1686846601010828D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1766325315388586D+00
       B=0.1087558139247680D+00
       V=0.1862698414660208D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2121744174481514D+00
       B=0.1413661374253096D+00
       V=0.2007430956991861D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2479669443408145D+00
       B=0.1748768214258880D+00
       V=0.2126568125394796D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2837600452294113D+00
       B=0.2089216406612073D+00
       V=0.2224394603372113D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3193344933193984D+00
       B=0.2431987685545972D+00
       V=0.2304264522673135D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3544935442438745D+00
       B=0.2774497054377770D+00
       V=0.2368854288424087D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3890571932288154D+00
       B=0.3114460356156915D+00
       V=0.2420352089461772D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4228581214259090D+00
       B=0.3449806851913012D+00
       V=0.2460597113081295D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4557387211304052D+00
       B=0.3778618641248256D+00
       V=0.2491181912257687D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4875487950541643D+00
       B=0.4099086391698978D+00
       V=0.2513528194205857D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5181436529962997D+00
       B=0.4409474925853973D+00
       V=0.2528943096693220D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5473824095600661D+00
       B=0.4708094517711291D+00
       V=0.2538660368488136D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5751263398976174D+00
       B=0.4993275140354637D+00
       V=0.2543868648299022D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1489515746840028D+00
       B=0.2599381993267017D-01
       V=0.1642595537825183D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1863656444351767D+00
       B=0.5479286532462190D-01
       V=0.1818246659849308D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2238602880356348D+00
       B=0.8556763251425254D-01
       V=0.1966565649492420D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2612723375728160D+00
       B=0.1177257802267011D+00
       V=0.2090677905657991D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2984332990206190D+00
       B=0.1508168456192700D+00
       V=0.2193820409510504D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3351786584663333D+00
       B=0.1844801892177727D+00
       V=0.2278870827661928D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3713505522209120D+00
       B=0.2184145236087598D+00
       V=0.2348283192282090D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4067981098954663D+00
       B=0.2523590641486229D+00
       V=0.2404139755581477D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4413769993687534D+00
       B=0.2860812976901373D+00
       V=0.2448227407760734D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4749487182516394D+00
       B=0.3193686757808996D+00
       V=0.2482110455592573D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5073798105075426D+00
       B=0.3520226949547602D+00
       V=0.2507192397774103D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5385410448878654D+00
       B=0.3838544395667890D+00
       V=0.2524765968534880D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5683065353670530D+00
       B=0.4146810037640963D+00
       V=0.2536052388539425D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5965527620663510D+00
       B=0.4443224094681121D+00
       V=0.2542230588033068D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2299227700856157D+00
       B=0.2865757664057584D-01
       V=0.1944817013047896D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2695752998553267D+00
       B=0.5923421684485993D-01
       V=0.2067862362746635D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3086178716611389D+00
       B=0.9117817776057715D-01
       V=0.2172440734649114D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3469649871659077D+00
       B=0.1240593814082605D+00
       V=0.2260125991723423D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3845153566319655D+00
       B=0.1575272058259175D+00
       V=0.2332655008689523D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4211600033403215D+00
       B=0.1912845163525413D+00
       V=0.2391699681532458D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4567867834329882D+00
       B=0.2250710177858171D+00
       V=0.2438801528273928D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4912829319232061D+00
       B=0.2586521303440910D+00
       V=0.2475370504260665D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5245364793303812D+00
       B=0.2918112242865407D+00
       V=0.2502707235640574D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5564369788915756D+00
       B=0.3243439239067890D+00
       V=0.2522031701054241D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5868757697775287D+00
       B=0.3560536787835351D+00
       V=0.2534511269978784D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6157458853519617D+00
       B=0.3867480821242581D+00
       V=0.2541284914955151D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3138461110672113D+00
       B=0.3051374637507278D-01
       V=0.2161509250688394D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3542495872050569D+00
       B=0.6237111233730755D-01
       V=0.2248778513437852D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3935751553120181D+00
       B=0.9516223952401907D-01
       V=0.2322388803404617D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4317634668111147D+00
       B=0.1285467341508517D+00
       V=0.2383265471001355D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4687413842250821D+00
       B=0.1622318931656033D+00
       V=0.2432476675019525D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5044274237060283D+00
       B=0.1959581153836453D+00
       V=0.2471122223750674D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5387354077925727D+00
       B=0.2294888081183837D+00
       V=0.2500291752486870D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5715768898356105D+00
       B=0.2626031152713945D+00
       V=0.2521055942764682D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6028627200136111D+00
       B=0.2950904075286713D+00
       V=0.2534472785575503D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6325039812653463D+00
       B=0.3267458451113286D+00
       V=0.2541599713080121D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3981986708423407D+00
       B=0.3183291458749821D-01
       V=0.2317380975862936D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4382791182133300D+00
       B=0.6459548193880908D-01
       V=0.2378550733719775D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4769233057218166D+00
       B=0.9795757037087952D-01
       V=0.2428884456739118D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5140823911194238D+00
       B=0.1316307235126655D+00
       V=0.2469002655757292D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5496977833862983D+00
       B=0.1653556486358704D+00
       V=0.2499657574265851D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5837047306512727D+00
       B=0.1988931724126510D+00
       V=0.2521676168486082D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6160349566926879D+00
       B=0.2320174581438950D+00
       V=0.2535935662645334D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6466185353209440D+00
       B=0.2645106562168662D+00
       V=0.2543356743363214D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4810835158795404D+00
       B=0.3275917807743992D-01
       V=0.2427353285201535D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5199925041324341D+00
       B=0.6612546183967181D-01
       V=0.2468258039744386D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5571717692207494D+00
       B=0.9981498331474143D-01
       V=0.2500060956440310D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5925789250836378D+00
       B=0.1335687001410374D+00
       V=0.2523238365420979D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6261658523859670D+00
       B=0.1671444402896463D+00
       V=0.2538399260252846D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6578811126669331D+00
       B=0.2003106382156076D+00
       V=0.2546255927268069D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5609624612998100D+00
       B=0.3337500940231335D-01
       V=0.2500583360048449D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5979959659984670D+00
       B=0.6708750335901803D-01
       V=0.2524777638260203D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6330523711054002D+00
       B=0.1008792126424850D+00
       V=0.2540951193860656D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6660960998103972D+00
       B=0.1345050343171794D+00
       V=0.2549524085027472D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6365384364585819D+00
       B=0.3372799460737052D-01
       V=0.2542569507009158D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6710994302899275D+00
       B=0.6755249309678028D-01
       V=0.2552114127580376D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD4802(X,Y,Z,W,N)
       DOUBLE PRECISION X(4802)
       DOUBLE PRECISION Y(4802)
       DOUBLE PRECISION Z(4802)
       DOUBLE PRECISION W(4802)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV 4802-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] V.I. LEBEDEV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
CHVD        ORDER OF ACCURACY"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
CHVD
CHVD   [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
CHVD       "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
CHVD
CHVD   [4] V.I. LEBEDEV
CHVD       "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
CHVD       SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
CHVD
CHVD   [5] V.I. LEBEDEV
CHVD       "QUADRATURES ON A SPHERE"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
CHVD       1976, PP. 10-24.
CHVD
CHVD   [6] V.I. LEBEDEV
CHVD       "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
CHVD        ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
CHVD        OCTAHEDRON GROUP WITH INVERSION"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
CHVD       1975, PP. 44-51.
CHVD
       N=1
       V=0.9687521879420705D-04
       CALL LEB_GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2307897895367918D-03
       CALL LEB_GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2297310852498558D-03
       CALL LEB_GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2335728608887064D-01
       V=0.7386265944001919D-04
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4352987836550653D-01
       V=0.8257977698542210D-04
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6439200521088801D-01
       V=0.9706044762057630D-04
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9003943631993181D-01
       V=0.1302393847117003D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1196706615548473D+00
       V=0.1541957004600968D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1511715412838134D+00
       V=0.1704459770092199D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1835982828503801D+00
       V=0.1827374890942906D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2165081259155405D+00
       V=0.1926360817436107D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2496208720417563D+00
       V=0.2008010239494833D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2827200673567900D+00
       V=0.2075635983209175D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3156190823994346D+00
       V=0.2131306638690909D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3481476793749115D+00
       V=0.2176562329937335D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3801466086947226D+00
       V=0.2212682262991018D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4114652119634011D+00
       V=0.2240799515668565D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4419598786519751D+00
       V=0.2261959816187525D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4714925949329543D+00
       V=0.2277156368808855D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4999293972879466D+00
       V=0.2287351772128336D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5271387221431248D+00
       V=0.2293490814084085D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5529896780837761D+00
       V=0.2296505312376273D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6000856099481712D+00
       V=0.2296793832318756D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6210562192785175D+00
       V=0.2295785443842974D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6401165879934240D+00
       V=0.2295017931529102D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6571144029244334D+00
       V=0.2295059638184868D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6718910821718863D+00
       V=0.2296232343237362D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6842845591099010D+00
       V=0.2298530178740771D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6941353476269816D+00
       V=0.2301579790280501D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7012965242212991D+00
       V=0.2304690404996513D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7056471428242644D+00
       V=0.2307027995907102D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4595557643585895D-01
       V=0.9312274696671092D-04
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1049316742435023D+00
       V=0.1199919385876926D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1773548879549274D+00
       V=0.1598039138877690D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2559071411236127D+00
       V=0.1822253763574900D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3358156837985898D+00
       V=0.1988579593655040D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4155835743763893D+00
       V=0.2112620102533307D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4937894296167472D+00
       V=0.2201594887699007D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5691569694793316D+00
       V=0.2261622590895036D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6405840854894251D+00
       V=0.2296458453435705D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7345133894143348D-01
       B=0.2177844081486067D-01
       V=0.1006006990267000D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1009859834044931D+00
       B=0.4590362185775188D-01
       V=0.1227676689635876D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1324289619748758D+00
       B=0.7255063095690877D-01
       V=0.1467864280270117D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1654272109607127D+00
       B=0.1017825451960684D+00
       V=0.1644178912101232D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1990767186776461D+00
       B=0.1325652320980364D+00
       V=0.1777664890718961D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2330125945523278D+00
       B=0.1642765374496765D+00
       V=0.1884825664516690D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2670080611108287D+00
       B=0.1965360374337889D+00
       V=0.1973269246453848D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3008753376294316D+00
       B=0.2290726770542238D+00
       V=0.2046767775855328D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3344475596167860D+00
       B=0.2616645495370823D+00
       V=0.2107600125918040D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3675709724070786D+00
       B=0.2941150728843141D+00
       V=0.2157416362266829D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4001000887587812D+00
       B=0.3262440400919066D+00
       V=0.2197557816920721D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4318956350436028D+00
       B=0.3578835350611916D+00
       V=0.2229192611835437D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4628239056795531D+00
       B=0.3888751854043678D+00
       V=0.2253385110212775D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4927563229773636D+00
       B=0.4190678003222840D+00
       V=0.2271137107548774D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5215687136707969D+00
       B=0.4483151836883852D+00
       V=0.2283414092917525D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5491402346984905D+00
       B=0.4764740676087880D+00
       V=0.2291161673130077D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5753520160126075D+00
       B=0.5034021310998277D+00
       V=0.2295313908576598D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1388326356417754D+00
       B=0.2435436510372806D-01
       V=0.1438204721359031D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1743686900537244D+00
       B=0.5118897057342652D-01
       V=0.1607738025495257D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2099737037950268D+00
       B=0.8014695048539634D-01
       V=0.1741483853528379D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2454492590908548D+00
       B=0.1105117874155699D+00
       V=0.1851918467519151D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2807219257864278D+00
       B=0.1417950531570966D+00
       V=0.1944628638070613D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3156842271975842D+00
       B=0.1736604945719597D+00
       V=0.2022495446275152D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3502090945177752D+00
       B=0.2058466324693981D+00
       V=0.2087462382438514D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3841684849519686D+00
       B=0.2381284261195919D+00
       V=0.2141074754818308D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4174372367906016D+00
       B=0.2703031270422569D+00
       V=0.2184640913748162D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4498926465011892D+00
       B=0.3021845683091309D+00
       V=0.2219309165220329D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4814146229807701D+00
       B=0.3335993355165720D+00
       V=0.2246123118340624D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5118863625734701D+00
       B=0.3643833735518232D+00
       V=0.2266062766915125D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5411947455119144D+00
       B=0.3943789541958179D+00
       V=0.2280072952230796D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5692301500357246D+00
       B=0.4234320144403542D+00
       V=0.2289082025202583D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5958857204139576D+00
       B=0.4513897947419260D+00
       V=0.2294012695120025D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2156270284785766D+00
       B=0.2681225755444491D-01
       V=0.1722434488736947D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2532385054909710D+00
       B=0.5557495747805614D-01
       V=0.1830237421455091D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2902564617771537D+00
       B=0.8569368062950249D-01
       V=0.1923855349997633D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3266979823143256D+00
       B=0.1167367450324135D+00
       V=0.2004067861936271D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3625039627493614D+00
       B=0.1483861994003304D+00
       V=0.2071817297354263D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3975838937548699D+00
       B=0.1803821503011405D+00
       V=0.2128250834102103D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4318396099009774D+00
       B=0.2124962965666424D+00
       V=0.2174513719440102D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4651706555732742D+00
       B=0.2445221837805913D+00
       V=0.2211661839150214D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4974752649620969D+00
       B=0.2762701224322987D+00
       V=0.2240665257813102D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5286517579627517D+00
       B=0.3075627775211328D+00
       V=0.2262439516632620D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5586001195731895D+00
       B=0.3382311089826877D+00
       V=0.2277874557231869D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5872229902021319D+00
       B=0.3681108834741399D+00
       V=0.2287854314454994D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6144258616235123D+00
       B=0.3970397446872839D+00
       V=0.2293268499615575D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2951676508064861D+00
       B=0.2867499538750441D-01
       V=0.1912628201529828D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3335085485472725D+00
       B=0.5867879341903510D-01
       V=0.1992499672238701D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3709561760636381D+00
       B=0.8961099205022284D-01
       V=0.2061275533454027D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4074722861667498D+00
       B=0.1211627927626297D+00
       V=0.2119318215968572D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4429923648839117D+00
       B=0.1530748903554898D+00
       V=0.2167416581882652D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4774428052721736D+00
       B=0.1851176436721877D+00
       V=0.2206430730516600D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5107446539535904D+00
       B=0.2170829107658179D+00
       V=0.2237186938699523D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5428151370542935D+00
       B=0.2487786689026271D+00
       V=0.2260480075032884D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5735699292556964D+00
       B=0.2800239952795016D+00
       V=0.2277098884558542D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6029253794562866D+00
       B=0.3106445702878119D+00
       V=0.2287845715109671D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6307998987073145D+00
       B=0.3404689500841194D+00
       V=0.2293547268236294D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3752652273692719D+00
       B=0.2997145098184479D-01
       V=0.2056073839852528D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4135383879344028D+00
       B=0.6086725898678011D-01
       V=0.2114235865831876D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4506113885153907D+00
       B=0.9238849548435643D-01
       V=0.2163175629770551D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4864401554606072D+00
       B=0.1242786603851851D+00
       V=0.2203392158111650D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5209708076611709D+00
       B=0.1563086731483386D+00
       V=0.2235473176847839D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5541422135830122D+00
       B=0.1882696509388506D+00
       V=0.2260024141501235D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5858880915113817D+00
       B=0.2199672979126059D+00
       V=0.2277675929329182D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6161399390603444D+00
       B=0.2512165482924867D+00
       V=0.2289102112284834D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6448296482255090D+00
       B=0.2818368701871888D+00
       V=0.2295027954625118D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4544796274917948D+00
       B=0.3088970405060312D-01
       V=0.2161281589879992D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4919389072146628D+00
       B=0.6240947677636835D-01
       V=0.2201980477395102D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5279313026985183D+00
       B=0.9430706144280313D-01
       V=0.2234952066593166D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5624169925571135D+00
       B=0.1263547818770374D+00
       V=0.2260540098520838D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5953484627093287D+00
       B=0.1583430788822594D+00
       V=0.2279157981899988D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6266730715339185D+00
       B=0.1900748462555988D+00
       V=0.2291296918565571D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6563363204278871D+00
       B=0.2213599519592567D+00
       V=0.2297533752536649D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5314574716585696D+00
       B=0.3152508811515374D-01
       V=0.2234927356465995D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5674614932298185D+00
       B=0.6343865291465561D-01
       V=0.2261288012985219D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6017706004970264D+00
       B=0.9551503504223951D-01
       V=0.2280818160923688D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6343471270264178D+00
       B=0.1275440099801196D+00
       V=0.2293773295180159D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6651494599127802D+00
       B=0.1593252037671960D+00
       V=0.2300528767338634D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6050184986005704D+00
       B=0.3192538338496105D-01
       V=0.2281893855065666D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6390163550880400D+00
       B=0.6402824353962306D-01
       V=0.2295720444840727D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6711199107088448D+00
       B=0.9609805077002909D-01
       V=0.2303227649026753D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6741354429572275D+00
       B=0.3211853196273233D-01
       V=0.2304831913227114D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD5294(X,Y,Z,W,N)
       DOUBLE PRECISION X(5294)
       DOUBLE PRECISION Y(5294)
       DOUBLE PRECISION Z(5294)
       DOUBLE PRECISION W(5294)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV 5294-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] V.I. LEBEDEV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
CHVD        ORDER OF ACCURACY"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
CHVD
CHVD   [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
CHVD       "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
CHVD
CHVD   [4] V.I. LEBEDEV
CHVD       "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
CHVD       SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
CHVD
CHVD   [5] V.I. LEBEDEV
CHVD       "QUADRATURES ON A SPHERE"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
CHVD       1976, PP. 10-24.
CHVD
CHVD   [6] V.I. LEBEDEV
CHVD       "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
CHVD        ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
CHVD        OCTAHEDRON GROUP WITH INVERSION"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
CHVD       1975, PP. 44-51.
CHVD
       N=1
       V=0.9080510764308163D-04
       CALL LEB_GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2084824361987793D-03
       CALL LEB_GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2303261686261450D-01
       V=0.5011105657239616D-04
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3757208620162394D-01
       V=0.5942520409683854D-04
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5821912033821852D-01
       V=0.9564394826109721D-04
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8403127529194872D-01
       V=0.1185530657126338D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1122927798060578D+00
       V=0.1364510114230331D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1420125319192987D+00
       V=0.1505828825605415D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1726396437341978D+00
       V=0.1619298749867023D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2038170058115696D+00
       V=0.1712450504267789D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2352849892876508D+00
       V=0.1789891098164999D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2668363354312461D+00
       V=0.1854474955629795D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2982941279900452D+00
       V=0.1908148636673661D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3295002922087076D+00
       V=0.1952377405281833D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3603094918363593D+00
       V=0.1988349254282232D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3905857895173920D+00
       V=0.2017079807160050D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4202005758160837D+00
       V=0.2039473082709094D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4490310061597227D+00
       V=0.2056360279288953D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4769586160311491D+00
       V=0.2068525823066865D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5038679887049750D+00
       V=0.2076724877534488D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5296454286519961D+00
       V=0.2081694278237885D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5541776207164850D+00
       V=0.2084157631219326D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5990467321921213D+00
       V=0.2084381531128593D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6191467096294587D+00
       V=0.2083476277129307D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6375251212901849D+00
       V=0.2082686194459732D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6540514381131168D+00
       V=0.2082475686112415D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6685899064391510D+00
       V=0.2083139860289915D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6810013009681648D+00
       V=0.2084745561831237D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6911469578730340D+00
       V=0.2087091313375890D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6988956915141736D+00
       V=0.2089718413297697D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7041335794868720D+00
       V=0.2092003303479793D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7067754398018567D+00
       V=0.2093336148263241D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3840368707853623D-01
       V=0.7591708117365267D-04
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9835485954117399D-01
       V=0.1083383968169186D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1665774947612998D+00
       V=0.1403019395292510D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2405702335362910D+00
       V=0.1615970179286436D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3165270770189046D+00
       V=0.1771144187504911D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3927386145645443D+00
       V=0.1887760022988168D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4678825918374656D+00
       V=0.1973474670768214D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5408022024266935D+00
       V=0.2033787661234659D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6104967445752438D+00
       V=0.2072343626517331D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6760910702685738D+00
       V=0.2091177834226918D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6655644120217392D-01
       B=0.1936508874588424D-01
       V=0.9316684484675566D-04
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9446246161270182D-01
       B=0.4252442002115869D-01
       V=0.1116193688682976D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1242651925452509D+00
       B=0.6806529315354374D-01
       V=0.1298623551559414D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1553438064846751D+00
       B=0.9560957491205369D-01
       V=0.1450236832456426D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1871137110542670D+00
       B=0.1245931657452888D+00
       V=0.1572719958149914D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2192612628836257D+00
       B=0.1545385828778978D+00
       V=0.1673234785867195D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2515682807206955D+00
       B=0.1851004249723368D+00
       V=0.1756860118725188D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2838535866287290D+00
       B=0.2160182608272384D+00
       V=0.1826776290439367D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3159578817528521D+00
       B=0.2470799012277111D+00
       V=0.1885116347992865D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3477370882791392D+00
       B=0.2781014208986402D+00
       V=0.1933457860170574D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3790576960890540D+00
       B=0.3089172523515731D+00
       V=0.1973060671902064D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4097938317810200D+00
       B=0.3393750055472244D+00
       V=0.2004987099616311D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4398256572859637D+00
       B=0.3693322470987730D+00
       V=0.2030170909281499D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4690384114718480D+00
       B=0.3986541005609877D+00
       V=0.2049461460119080D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4973216048301053D+00
       B=0.4272112491408562D+00
       V=0.2063653565200186D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5245681526132446D+00
       B=0.4548781735309936D+00
       V=0.2073507927381027D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5506733911803888D+00
       B=0.4815315355023251D+00
       V=0.2079764593256122D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5755339829522475D+00
       B=0.5070486445801855D+00
       V=0.2083150534968778D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1305472386056362D+00
       B=0.2284970375722366D-01
       V=0.1262715121590664D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1637327908216477D+00
       B=0.4812254338288384D-01
       V=0.1414386128545972D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1972734634149637D+00
       B=0.7531734457511935D-01
       V=0.1538740401313898D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2308694653110130D+00
       B=0.1039043639882017D+00
       V=0.1642434942331432D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2643899218338160D+00
       B=0.1334526587117626D+00
       V=0.1729790609237496D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2977171599622171D+00
       B=0.1636414868936382D+00
       V=0.1803505190260828D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3307293903032310D+00
       B=0.1942195406166568D+00
       V=0.1865475350079657D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3633069198219073D+00
       B=0.2249752879943753D+00
       V=0.1917182669679069D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3953346955922727D+00
       B=0.2557218821820032D+00
       V=0.1959851709034382D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4267018394184914D+00
       B=0.2862897925213193D+00
       V=0.1994529548117882D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4573009622571704D+00
       B=0.3165224536636518D+00
       V=0.2022138911146548D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4870279559856109D+00
       B=0.3462730221636496D+00
       V=0.2043518024208592D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5157819581450322D+00
       B=0.3754016870282835D+00
       V=0.2059450313018110D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5434651666465393D+00
       B=0.4037733784993613D+00
       V=0.2070685715318472D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5699823887764627D+00
       B=0.4312557784139123D+00
       V=0.2077955310694373D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5952403350947741D+00
       B=0.4577175367122110D+00
       V=0.2081980387824712D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2025152599210369D+00
       B=0.2520253617719557D-01
       V=0.1521318610377956D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2381066653274425D+00
       B=0.5223254506119000D-01
       V=0.1622772720185755D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2732823383651612D+00
       B=0.8060669688588620D-01
       V=0.1710498139420709D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3080137692611118D+00
       B=0.1099335754081255D+00
       V=0.1785911149448736D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3422405614587601D+00
       B=0.1399120955959857D+00
       V=0.1850125313687736D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3758808773890420D+00
       B=0.1702977801651705D+00
       V=0.1904229703933298D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4088458383438932D+00
       B=0.2008799256601680D+00
       V=0.1949259956121987D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4410450550841152D+00
       B=0.2314703052180836D+00
       V=0.1986161545363960D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4723879420561312D+00
       B=0.2618972111375892D+00
       V=0.2015790585641370D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5027843561874343D+00
       B=0.2920013195600270D+00
       V=0.2038934198707418D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5321453674452458D+00
       B=0.3216322555190551D+00
       V=0.2056334060538251D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5603839113834030D+00
       B=0.3506456615934198D+00
       V=0.2068705959462289D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5874150706875146D+00
       B=0.3789007181306267D+00
       V=0.2076753906106002D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6131559381660038D+00
       B=0.4062580170572782D+00
       V=0.2081179391734803D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2778497016394506D+00
       B=0.2696271276876226D-01
       V=0.1700345216228943D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3143733562261912D+00
       B=0.5523469316960465D-01
       V=0.1774906779990410D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3501485810261827D+00
       B=0.8445193201626464D-01
       V=0.1839659377002642D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3851430322303653D+00
       B=0.1143263119336083D+00
       V=0.1894987462975169D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4193013979470415D+00
       B=0.1446177898344475D+00
       V=0.1941548809452595D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4525585960458567D+00
       B=0.1751165438438091D+00
       V=0.1980078427252384D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4848447779622947D+00
       B=0.2056338306745660D+00
       V=0.2011296284744488D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5160871208276894D+00
       B=0.2359965487229226D+00
       V=0.2035888456966776D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5462112185696926D+00
       B=0.2660430223139146D+00
       V=0.2054516325352142D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5751425068101757D+00
       B=0.2956193664498032D+00
       V=0.2067831033092635D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6028073872853596D+00
       B=0.3245763905312779D+00
       V=0.2076485320284876D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6291338275278409D+00
       B=0.3527670026206972D+00
       V=0.2081141439525255D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3541797528439391D+00
       B=0.2823853479435550D-01
       V=0.1834383015469222D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3908234972074657D+00
       B=0.5741296374713106D-01
       V=0.1889540591777677D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4264408450107590D+00
       B=0.8724646633650199D-01
       V=0.1936677023597375D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4609949666553286D+00
       B=0.1175034422915616D+00
       V=0.1976176495066504D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4944389496536006D+00
       B=0.1479755652628428D+00
       V=0.2008536004560983D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5267194884346086D+00
       B=0.1784740659484352D+00
       V=0.2034280351712291D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5577787810220990D+00
       B=0.2088245700431244D+00
       V=0.2053944466027758D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5875563763536670D+00
       B=0.2388628136570763D+00
       V=0.2068077642882360D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6159910016391269D+00
       B=0.2684308928769185D+00
       V=0.2077250949661599D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6430219602956268D+00
       B=0.2973740761960252D+00
       V=0.2082062440705320D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4300647036213646D+00
       B=0.2916399920493977D-01
       V=0.1934374486546626D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4661486308935531D+00
       B=0.5898803024755659D-01
       V=0.1974107010484300D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5009658555287261D+00
       B=0.8924162698525409D-01
       V=0.2007129290388658D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5344824270447704D+00
       B=0.1197185199637321D+00
       V=0.2033736947471293D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5666575997416371D+00
       B=0.1502300756161382D+00
       V=0.2054287125902493D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5974457471404752D+00
       B=0.1806004191913564D+00
       V=0.2069184936818894D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6267984444116886D+00
       B=0.2106621764786252D+00
       V=0.2078883689808782D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6546664713575417D+00
       B=0.2402526932671914D+00
       V=0.2083886366116359D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5042711004437253D+00
       B=0.2982529203607657D-01
       V=0.2006593275470817D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5392127456774380D+00
       B=0.6008728062339922D-01
       V=0.2033728426135397D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5726819437668618D+00
       B=0.9058227674571398D-01
       V=0.2055008781377608D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6046469254207278D+00
       B=0.1211219235803400D+00
       V=0.2070651783518502D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6350716157434952D+00
       B=0.1515286404791580D+00
       V=0.2080953335094320D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6639177679185454D+00
       B=0.1816314681255552D+00
       V=0.2086284998988521D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5757276040972253D+00
       B=0.3026991752575440D-01
       V=0.2055549387644668D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6090265823139755D+00
       B=0.6078402297870770D-01
       V=0.2071871850267654D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6406735344387661D+00
       B=0.9135459984176636D-01
       V=0.2082856600431965D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6706397927793709D+00
       B=0.1218024155966590D+00
       V=0.2088705858819358D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6435019674426665D+00
       B=0.3052608357660639D-01
       V=0.2083995867536322D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6747218676375681D+00
       B=0.6112185773983089D-01
       V=0.2090509712889637D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE LD5810(X,Y,Z,W,N)
       DOUBLE PRECISION X(5810)
       DOUBLE PRECISION Y(5810)
       DOUBLE PRECISION Z(5810)
       DOUBLE PRECISION W(5810)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    LEBEDEV 5810-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   THIS CODE WAS DISTRIBUTED THROUGH CCL (HTTP://WWW.CCL.NET/).
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] V.I. LEBEDEV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF 59TH ALGEBRAIC
CHVD        ORDER OF ACCURACY"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 50, 1995, PP. 283-286.
CHVD
CHVD   [3] V.I. LEBEDEV, AND A.L. SKOROKHODOV
CHVD       "QUADRATURE FORMULAS OF ORDERS 41, 47, AND 53 FOR THE SPHERE"
CHVD       RUSSIAN ACAD. SCI. DOKL. MATH., VOL. 45, 1992, PP. 587-592.
CHVD
CHVD   [4] V.I. LEBEDEV
CHVD       "SPHERICAL QUADRATURE FORMULAS EXACT TO ORDERS 25-29"
CHVD       SIBERIAN MATHEMATICAL JOURNAL, VOL. 18, 1977, PP. 99-107.
CHVD
CHVD   [5] V.I. LEBEDEV
CHVD       "QUADRATURES ON A SPHERE"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 16,
CHVD       1976, PP. 10-24.
CHVD
CHVD   [6] V.I. LEBEDEV
CHVD       "VALUES OF THE NODES AND WEIGHTS OF NINTH TO SEVENTEENTH
CHVD        ORDER GAUSS-MARKOV QUADRATURE FORMULAE INVARIANT UNDER THE
CHVD        OCTAHEDRON GROUP WITH INVERSION"
CHVD       COMPUTATIONAL MATHEMATICS AND MATHEMATICAL PHYSICS, VOL. 15,
CHVD       1975, PP. 44-51.
CHVD
       N=1
       V=0.9735347946175486D-05
       CALL LEB_GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.1907581241803167D-03
       CALL LEB_GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.1901059546737578D-03
       CALL LEB_GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1182361662400277D-01
       V=0.3926424538919212D-04
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3062145009138958D-01
       V=0.6667905467294382D-04
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5329794036834243D-01
       V=0.8868891315019135D-04
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7848165532862220D-01
       V=0.1066306000958872D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1054038157636201D+00
       V=0.1214506743336128D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1335577797766211D+00
       V=0.1338054681640871D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1625769955502252D+00
       V=0.1441677023628504D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1921787193412792D+00
       V=0.1528880200826557D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2221340534690548D+00
       V=0.1602330623773609D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2522504912791132D+00
       V=0.1664102653445244D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2823610860679697D+00
       V=0.1715845854011323D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3123173966267560D+00
       V=0.1758901000133069D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3419847036953789D+00
       V=0.1794382485256736D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3712386456999758D+00
       V=0.1823238106757407D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3999627649876828D+00
       V=0.1846293252959976D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4280466458648093D+00
       V=0.1864284079323098D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4553844360185711D+00
       V=0.1877882694626914D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4818736094437834D+00
       V=0.1887716321852025D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5074138709260629D+00
       V=0.1894381638175673D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5319061304570707D+00
       V=0.1898454899533629D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5552514978677286D+00
       V=0.1900497929577815D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5981009025246183D+00
       V=0.1900671501924092D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6173990192228116D+00
       V=0.1899837555533510D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6351365239411131D+00
       V=0.1899014113156229D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6512010228227200D+00
       V=0.1898581257705106D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6654758363948120D+00
       V=0.1898804756095753D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6778410414853370D+00
       V=0.1899793610426402D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6881760887484110D+00
       V=0.1901464554844117D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6963645267094598D+00
       V=0.1903533246259542D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7023010617153579D+00
       V=0.1905556158463228D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7059004636628753D+00
       V=0.1907037155663528D-03
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3552470312472575D-01
       V=0.5992997844249967D-04
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9151176620841283D-01
       V=0.9749059382456978D-04
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1566197930068980D+00
       V=0.1241680804599158D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2265467599271907D+00
       V=0.1437626154299360D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2988242318581361D+00
       V=0.1584200054793902D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3717482419703886D+00
       V=0.1694436550982744D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4440094491758889D+00
       V=0.1776617014018108D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5145337096756642D+00
       V=0.1836132434440077D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5824053672860230D+00
       V=0.1876494727075983D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6468283961043370D+00
       V=0.1899906535336482D-03
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6095964259104373D-01
       B=0.1787828275342931D-01
       V=0.8143252820767350D-04
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8811962270959388D-01
       B=0.3953888740792096D-01
       V=0.9998859890887728D-04
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1165936722428831D+00
       B=0.6378121797722990D-01
       V=0.1156199403068359D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1460232857031785D+00
       B=0.8985890813745037D-01
       V=0.1287632092635513D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1761197110181755D+00
       B=0.1172606510576162D+00
       V=0.1398378643365139D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2066471190463718D+00
       B=0.1456102876970995D+00
       V=0.1491876468417391D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2374076026328152D+00
       B=0.1746153823011775D+00
       V=0.1570855679175456D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2682305474337051D+00
       B=0.2040383070295584D+00
       V=0.1637483948103775D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2989653312142369D+00
       B=0.2336788634003698D+00
       V=0.1693500566632843D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3294762752772209D+00
       B=0.2633632752654219D+00
       V=0.1740322769393633D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3596390887276086D+00
       B=0.2929369098051601D+00
       V=0.1779126637278296D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3893383046398812D+00
       B=0.3222592785275512D+00
       V=0.1810908108835412D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4184653789358347D+00
       B=0.3512004791195743D+00
       V=0.1836529132600190D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4469172319076166D+00
       B=0.3796385677684537D+00
       V=0.1856752841777379D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4745950813276976D+00
       B=0.4074575378263879D+00
       V=0.1872270566606832D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5014034601410262D+00
       B=0.4345456906027828D+00
       V=0.1883722645591307D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5272493404551239D+00
       B=0.4607942515205134D+00
       V=0.1891714324525297D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5520413051846366D+00
       B=0.4860961284181720D+00
       V=0.1896827480450146D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5756887237503077D+00
       B=0.5103447395342790D+00
       V=0.1899628417059528D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1225039430588352D+00
       B=0.2136455922655793D-01
       V=0.1123301829001669D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1539113217321372D+00
       B=0.4520926166137188D-01
       V=0.1253698826711277D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1856213098637712D+00
       B=0.7086468177864818D-01
       V=0.1366266117678531D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2174998728035131D+00
       B=0.9785239488772918D-01
       V=0.1462736856106918D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2494128336938330D+00
       B=0.1258106396267210D+00
       V=0.1545076466685412D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2812321562143480D+00
       B=0.1544529125047001D+00
       V=0.1615096280814007D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3128372276456111D+00
       B=0.1835433512202753D+00
       V=0.1674366639741759D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3441145160177973D+00
       B=0.2128813258619585D+00
       V=0.1724225002437900D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3749567714853510D+00
       B=0.2422913734880829D+00
       V=0.1765810822987288D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4052621732015610D+00
       B=0.2716163748391453D+00
       V=0.1800104126010751D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4349335453522385D+00
       B=0.3007127671240280D+00
       V=0.1827960437331284D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4638776641524965D+00
       B=0.3294470677216479D+00
       V=0.1850140300716308D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4920046410462687D+00
       B=0.3576932543699155D+00
       V=0.1867333507394938D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5192273554861704D+00
       B=0.3853307059757764D+00
       V=0.1880178688638289D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5454609081136522D+00
       B=0.4122425044452694D+00
       V=0.1889278925654758D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5706220661424140D+00
       B=0.4383139587781027D+00
       V=0.1895213832507346D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5946286755181518D+00
       B=0.4634312536300553D+00
       V=0.1898548277397420D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1905370790924295D+00
       B=0.2371311537781979D-01
       V=0.1349105935937341D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2242518717748009D+00
       B=0.4917878059254806D-01
       V=0.1444060068369326D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2577190808025936D+00
       B=0.7595498960495142D-01
       V=0.1526797390930008D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2908724534927187D+00
       B=0.1036991083191100D+00
       V=0.1598208771406474D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3236354020056219D+00
       B=0.1321348584450234D+00
       V=0.1659354368615331D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3559267359304543D+00
       B=0.1610316571314789D+00
       V=0.1711279910946440D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3876637123676956D+00
       B=0.1901912080395707D+00
       V=0.1754952725601440D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4187636705218842D+00
       B=0.2194384950137950D+00
       V=0.1791247850802529D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4491449019883107D+00
       B=0.2486155334763858D+00
       V=0.1820954300877716D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4787270932425445D+00
       B=0.2775768931812335D+00
       V=0.1844788524548449D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5074315153055574D+00
       B=0.3061863786591120D+00
       V=0.1863409481706220D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5351810507738336D+00
       B=0.3343144718152556D+00
       V=0.1877433008795068D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5619001025975381D+00
       B=0.3618362729028427D+00
       V=0.1887444543705232D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5875144035268046D+00
       B=0.3886297583620408D+00
       V=0.1894009829375006D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6119507308734495D+00
       B=0.4145742277792031D+00
       V=0.1897683345035198D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2619733870119463D+00
       B=0.2540047186389353D-01
       V=0.1517327037467653D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2968149743237949D+00
       B=0.5208107018543989D-01
       V=0.1587740557483543D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3310451504860488D+00
       B=0.7971828470885599D-01
       V=0.1649093382274097D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3646215567376676D+00
       B=0.1080465999177927D+00
       V=0.1701915216193265D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3974916785279360D+00
       B=0.1368413849366629D+00
       V=0.1746847753144065D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4295967403772029D+00
       B=0.1659073184763559D+00
       V=0.1784555512007570D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4608742854473447D+00
       B=0.1950703730454614D+00
       V=0.1815687562112174D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4912598858949903D+00
       B=0.2241721144376724D+00
       V=0.1840864370663302D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5206882758945558D+00
       B=0.2530655255406489D+00
       V=0.1860676785390006D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5490940914019819D+00
       B=0.2816118409731066D+00
       V=0.1875690583743703D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5764123302025542D+00
       B=0.3096780504593238D+00
       V=0.1886453236347225D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6025786004213506D+00
       B=0.3371348366394987D+00
       V=0.1893501123329645D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6275291964794956D+00
       B=0.3638547827694396D+00
       V=0.1897366184519868D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3348189479861771D+00
       B=0.2664841935537443D-01
       V=0.1643908815152736D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3699515545855295D+00
       B=0.5424000066843495D-01
       V=0.1696300350907768D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4042003071474669D+00
       B=0.8251992715430854D-01
       V=0.1741553103844483D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4375320100182624D+00
       B=0.1112695182483710D+00
       V=0.1780015282386092D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4699054490335947D+00
       B=0.1402964116467816D+00
       V=0.1812116787077125D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5012739879431952D+00
       B=0.1694275117584291D+00
       V=0.1838323158085421D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5315874883754966D+00
       B=0.1985038235312689D+00
       V=0.1859113119837737D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5607937109622117D+00
       B=0.2273765660020893D+00
       V=0.1874969220221698D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5888393223495521D+00
       B=0.2559041492849764D+00
       V=0.1886375612681076D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6156705979160163D+00
       B=0.2839497251976899D+00
       V=0.1893819575809276D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6412338809078123D+00
       B=0.3113791060500690D+00
       V=0.1897794748256767D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4076051259257167D+00
       B=0.2757792290858463D-01
       V=0.1738963926584846D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4423788125791520D+00
       B=0.5584136834984293D-01
       V=0.1777442359873466D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4760480917328258D+00
       B=0.8457772087727143D-01
       V=0.1810010815068719D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5085838725946297D+00
       B=0.1135975846359248D+00
       V=0.1836920318248129D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5399513637391218D+00
       B=0.1427286904765053D+00
       V=0.1858489473214328D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5701118433636380D+00
       B=0.1718112740057635D+00
       V=0.1875079342496592D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5990240530606021D+00
       B=0.2006944855985351D+00
       V=0.1887080239102310D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6266452685139695D+00
       B=0.2292335090598907D+00
       V=0.1894905752176822D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6529320971415942D+00
       B=0.2572871512353714D+00
       V=0.1898991061200695D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4791583834610126D+00
       B=0.2826094197735932D-01
       V=0.1809065016458791D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5130373952796940D+00
       B=0.5699871359683649D-01
       V=0.1836297121596799D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5456252429628476D+00
       B=0.8602712528554394D-01
       V=0.1858426916241869D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5768956329682385D+00
       B=0.1151748137221281D+00
       V=0.1875654101134641D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6068186944699046D+00
       B=0.1442811654136362D+00
       V=0.1888240751833503D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6353622248024907D+00
       B=0.1731930321657680D+00
       V=0.1896497383866979D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6624927035731797D+00
       B=0.2017619958756061D+00
       V=0.1900775530219121D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5484933508028488D+00
       B=0.2874219755907391D-01
       V=0.1858525041478814D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5810207682142106D+00
       B=0.5778312123713695D-01
       V=0.1876248690077947D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6120955197181352D+00
       B=0.8695262371439526D-01
       V=0.1889404439064607D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6416944284294319D+00
       B=0.1160893767057166D+00
       V=0.1898168539265290D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6697926391731260D+00
       B=0.1450378826743251D+00
       V=0.1902779940661772D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6147594390585488D+00
       B=0.2904957622341456D-01
       V=0.1890125641731815D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6455390026356783D+00
       B=0.5823809152617197D-01
       V=0.1899434637795751D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6747258588365477D+00
       B=0.8740384899884715D-01
       V=0.1904520856831751D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6772135750395347D+00
       B=0.2919946135808105D-01
       V=0.1905534498734563D-03
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE OD0246(X,Y,Z,W,N)
       DOUBLE PRECISION X( 246)
       DOUBLE PRECISION Y( 246)
       DOUBLE PRECISION Z( 246)
       DOUBLE PRECISION W( 246)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    OCTAHEDRAL SYMMETRY (W/O INVERSION) 246-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] A.S. POPOV
CHVD       "THE SEARCH FOR THE SPHERE OF THE BEST CUBATURE FORMULAE
CHVD        INVARIANT UNDER OCTAHEDRAL GROUP OF ROTATIONS"
CHVD       SIB. ZH. VYCHISL. MAT., VOL. 5, NO. 4, 2002, PP. 367???372
CHVD
       N=1
       V=0.3897235138937088D-02
       CALL LEB_GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3493936461315880D-02
       A=0.9566807992266356D+00
       B=0.2833386636495198D+00
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3591650554371342D-02
       A=0.9757117308401445D+00
       B=0.1167317879067011D+00
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3967953538531415D-02
       A=0.9001427110863452D+00
       B=0.1004016871311435D+00
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.4044472616286787D-02
       A=0.9016081884868434D+00
       B=0.3046563962638665D+00
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.4158726826302798D-02
       A=0.8657185254531113D+00
       B=0.4789871418590519D+00
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.4223724523384502D-02
       A=0.7958577819235757D+00
       B=0.2817783845311169D+00
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.4287698313089224D-02
       A=0.7828854398832512D+00
       B=0.3632315411734175D-01
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.4295521563270202D-02
       A=0.7797755454487808D+00
       B=0.4876721073724419D+00
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.4306264958684786D-02
       A=0.6482368665122326D+00
       B=0.4557085192757556D+00
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.4322408526695458D-02
       A=0.7237785649544816D+00
       B=0.6543863740765325D+00
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE OD0264(X,Y,Z,W,N)
       DOUBLE PRECISION X( 264)
       DOUBLE PRECISION Y( 264)
       DOUBLE PRECISION Z( 264)
       DOUBLE PRECISION W( 264)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    OCTAHEDRAL SYMMETRY (W/O INVERSION) 264-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] A.S. POPOV
CHVD       "THE SEARCH FOR THE SPHERE OF THE BEST CUBATURE FORMULAE
CHVD        INVARIANT UNDER OCTAHEDRAL GROUP OF ROTATIONS"
CHVD       SIB. ZH. VYCHISL. MAT., VOL. 5, NO. 4, 2002, PP. 367???372
CHVD
       N=1
       V=3.404743138426950D-03
       A=7.501162447877730D-01
       B=6.569874882943270D-01
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=3.527346280748710D-03
       A=8.810738528375910D-01
       B=2.804539315326130D-01
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=3.560705552891180D-03
       A=6.985522168362600D-01
       B=6.638349164875560D-01
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=3.661015502447230D-03
       A=9.884416293559190D-01
       B=4.187763670001430D-02
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=3.734364518407020D-03
       A=8.703452942381050D-01
       B=4.831216465940420D-01
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=3.786218818475490D-03
       A=7.577370457280040D-01
       B=3.737914001474930D-01
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=3.812295213144400D-03
       A=9.340175312890680D-01
       B=7.308888363903200D-02
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=3.922282615281700D-03
       A=6.702365796250490D-01
       B=5.696274199623510D-01
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=3.934290436352150D-03
       A=8.218581364908940D-01
       B=4.788479450589300D-01
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=4.143508784124530D-03
       A=9.465732078368030D-01
       B=2.750687576636780D-01
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=4.179895806367260D-03
       A=8.172239405554010D-01
       B=1.335491891058270D-01
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE OD0342(X,Y,Z,W,N)
       DOUBLE PRECISION X( 342)
       DOUBLE PRECISION Y( 342)
       DOUBLE PRECISION Z( 342)
       DOUBLE PRECISION W( 342)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    OCTAHEDRAL SYMMETRY (W/O INVERSION) 264-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] A.S. POPOV
CHVD       "CUBATURE FORMULAE FOR THE SPHERE
CHVD        INVARIANT UNDER OCTAHEDRAL GROUP OF ROTATIONS" (IN RUSSIAN)
CHVD       COMPUT. MATH. MATH. PHYS., VOL. 38, NO. 1, 1998, PP. 30-37
CHVD
       N=1
       V=0.3023648748223408D-02
       CALL LEB_GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2640009601341421D-02
       A=0.8596638928344739D+00
       B=0.4062242558266217D+00
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2743272994096360D-02
       A=0.9048298313140403D+00
       B=0.2370253708753523D+00
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2892274286735084D-02
       A=0.7919900954039140D+00
       B=0.3635572640746994D+00
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2898956168072910D-02
       A=0.9246930024023610D+00
       B=0.5409950346576234D-01
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2922211185898755D-02
       A=0.7983845764909109D+00
       B=0.5645514094465199D+00
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2926667482335111D-02
       A=0.9629723259634715D+00
       B=0.2038117136355577D+00
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2931773085197849D-02
       A=0.9120669247729113D+00
       B=0.3877575423799307D+00
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2936869297379136D-02
       A=0.7377581016329626D+00
       B=0.5411494384156677D+00
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2975927626183587D-02
       A=0.9810590688132057D+00
       B=0.1311454620366076D-01
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2986111941772548D-02
       A=0.6522834517550465D+00
       B=0.4833339354284559D+00
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2987229421263409D-02
       A=0.6952847727766393D+00
       B=0.2919213249111464D+00
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2995506191909761D-02
       A=0.8278662627227689D+00
       B=0.1738700791494905D+00
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3035869254195839D-02
       A=0.8316507582159955D+00
       B=0.5548939979903705D+00
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3038075943229046D-02
       A=0.7126773308879356D+00
       B=0.9613108329184178D-01
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE OHD0350(X,Y,Z,W,N)
       DOUBLE PRECISION X( 350)
       DOUBLE PRECISION Y( 350)
       DOUBLE PRECISION Z( 350)
       DOUBLE PRECISION W( 350)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    Oh SYMMETRY 350-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] A.S. POPOV
CHVD       "THE SEARCH FOR THE SPHERE OF THE BEST CUBATURE FORMULAE
CHVD        INVARIANT UNDER OCTAHEDRAL GROUP OF ROTATIONS
CHVD        WITH INVERSION FOR A SPHERE"
CHVD       SIB. ZH. VYCHISL. MAT., VOL. 8, NO. 2, 2005, PP. 143???148
CHVD
       N=1
       V=0.3022655957073096D-02
       CALL LEB_GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3053589782677049D-02
       CALL LEB_GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.1606436947071299D-02
       A=0.7240037770684867D+00
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2894462580088764D-02
       A=0.9249389183479164D+00
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2992303278280904D-02
       A=0.9810676749667857D+00
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2705550334445959D-02
       A=0.3622473127982810D+00
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2966242835900409D-02
       A=0.1917690066919518D+00
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2999330613486812D-02
       A=0.4801685990565060D+00
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3041905028826839D-02
       A=0.6502682026896775D+00
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3044469577443018D-02
       A=0.6935461975932684D+00
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2846671108207593D-02
       A=0.9077761244847474D+00
       B=0.3765168541264181D+00
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2954001654360567D-02
       A=0.7932087174843227D+00
       B=0.5354121187940431D+00
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3020546347912858D-02
       A=0.8290201309029455D+00
       B=0.5505983319551926D+00
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)

       N=N-1
       RETURN
       END
       SUBROUTINE OHD0398(X,Y,Z,W,N)
       DOUBLE PRECISION X( 398)
       DOUBLE PRECISION Y( 398)
       DOUBLE PRECISION Z( 398)
       DOUBLE PRECISION W( 398)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    Oh SYMMETRY 398-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] A.S. POPOV
CHVD       "THE SEARCH FOR THE SPHERE OF THE BEST CUBATURE FORMULAE
CHVD        INVARIANT UNDER OCTAHEDRAL GROUP OF ROTATIONS
CHVD        WITH INVERSION FOR A SPHERE"
CHVD       SIB. ZH. VYCHISL. MAT., VOL. 8, NO. 2, 2005, PP. 143???148
CHVD
       N=1
       V=0.1822084579093247D-02
       CALL LEB_GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2827377032685751D-02
       CALL LEB_GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.1575484799105965D-02
       A=0.7062503174494366D+00
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2580087158579662D-02
       A=0.1133790435116091D+00
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2639846092802613D-02
       A=0.6901038654954956D+00
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2783844769874123D-02
       A=0.6457774767552608D+00
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2846001949687294D-02
       A=0.4796087844363489D+00
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2873971674037779D-02
       A=0.3632128775162115D+00
       CALL LEB_GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2553168380704661D-02
       A=0.9648952557917490D+00
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2654875161144892D-02
       A=0.9014928604201154D+00
       CALL LEB_GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2220113527157410D-02
       A=0.8797615447057731D+00
       B=0.4380336735755111D+00
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2326306559741858D-02
       A=0.9388803094389767D+00
       B=0.2946263848640117D+00
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2611549635094245D-02
       A=0.8091315291511537D+00
       B=0.5793558553216747D+00
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2722733540537044D-02
       A=0.7838817431841564D+00
       B=0.5439488406950955D+00
       CALL LEB_GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE OD0432(X,Y,Z,W,N)
       DOUBLE PRECISION X( 432)
       DOUBLE PRECISION Y( 432)
       DOUBLE PRECISION Z( 432)
       DOUBLE PRECISION W( 432)
       INTEGER N
       DOUBLE PRECISION A,B,V
CVW
CVW    OCTAHEDRAL SYMMETRY (W/O INVERSION) 264-POINT ANGULAR GRID
CVW
CHVD
CHVD   THIS ROUTINE IS PART OF A SET OF ROUTINES THAT GENERATE
CHVD   LEBEDEV GRIDS [1-6] FOR INTEGRATION ON A SPHERE. THE ORIGINAL
CHVD   C-CODE [1] WAS KINDLY PROVIDED BY DR. DMITRI N. LAIKOV AND
CHVD   TRANSLATED INTO FORTRAN BY DR. CHRISTOPH VAN WUELLEN.
CHVD   THIS ROUTINE WAS TRANSLATED USING A C TO FORTRAN77 CONVERSION
CHVD   TOOL WRITTEN BY DR. CHRISTOPH VAN WUELLEN.
CHVD
CHVD   USERS OF THIS CODE ARE ASKED TO INCLUDE REFERENCE [1] IN THEIR
CHVD   PUBLICATIONS, AND IN THE USER- AND PROGRAMMERS-MANUALS
CHVD   DESCRIBING THEIR CODES.
CHVD
CHVD   [1] V.I. LEBEDEV, AND D.N. LAIKOV
CHVD       "A QUADRATURE FORMULA FOR THE SPHERE OF THE 131ST
CHVD        ALGEBRAIC ORDER OF ACCURACY"
CHVD       DOKLADY MATHEMATICS, VOL. 59, NO. 3, 1999, PP. 477-481.
CHVD
CHVD   [2] A.S. POPOV
CHVD       "CUBATURE FORMULAE FOR THE SPHERE
CHVD        INVARIANT UNDER OCTAHEDRAL GROUP OF ROTATIONS" (IN RUSSIAN)
CHVD       COMPUT. MATH. MATH. PHYS., VOL. 38, NO. 1, 1998, PP. 30-37
CHVD
       N=1
       V=0.2073915331320519D-02
       A=0.7489806258088988D+00
       B=0.6108701068802251D+00
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2107255030397068D-02
       A=0.9161730444214968D+00
       B=0.3197412043673950D+00
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2115342322975654D-02
       A=0.9641217858844791D+00
       B=0.1775069067845282D+00
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2139217879967886D-02
       A=0.9932625280244080D+00
       B=0.4909487312548763D-01
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2213455223136293D-02
       A=0.9350867590143857D+00
       B=0.1021948630931257D+00
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2220017128107630D-02
       A=0.8328350332317131D+00
       B=0.5345156354592039D+00
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2255553046127313D-02
       A=0.9113314027451442D+00
       B=0.4001239125939150D+00
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2317129463567652D-02
       A=0.8726339023393756D+00
       B=0.5449493821656796D-01
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2319753441076136D-02
       A=0.9678242462766119D+00
       B=0.2485055252796990D+00
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2346020122932169D-02
       A=0.7078265977814910D+00
       B=0.1321992811820306D+00
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2365479437733548D-02
       A=0.7070053872543563D+00
       B=0.3277599108679579D+00
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2373641673707936D-02
       A=0.7815225804060471D+00
       B=0.6238577434353162D+00
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2434041150670229D-02
       A=0.7200100968665680D+00
       B=0.5549556437629704D+00
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2458995407627726D-02
       A=0.8068208852018608D+00
       B=0.1984078314629552D+00
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2465111275848897D-02
       A=0.8385089184456621D+00
       B=0.4494208832677718D+00
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2466377894866924D-02
       A=0.8801102003441889D+00
       B=0.2566114675058704D+00
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2496086422248489D-02
       A=0.7879155638776240D+00
       B=0.3902320295286575D+00
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2499274414354599D-02
       A=0.6504900372238729D+00
       B=0.4941000883702588D+00
       CALL LEB_GEN_OH( 7, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE ID0192(X,Y,Z,W,N)
       DOUBLE PRECISION X( 192)
       DOUBLE PRECISION Y( 192)
       DOUBLE PRECISION Z( 192)
       DOUBLE PRECISION W( 192)
       INTEGER N
       INTEGER I
       DOUBLE PRECISION A,B,V
       N=1
       V = 0.4164880157580243D-02
       CALL GEN_IH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V = 0.5077679306001083D-02
       A = 0.8032626540631647D+00
       B = 0.5443960115691368D+00
       CALL GEN_IH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V = 0.5349546802216597D-02
       A = 0.6978513452732168D+00
       B = 0.5282061903582911D+00
       CALL GEN_IH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V = 0.5406464526932938D-02
       A = 0.7323722630950026D+00
       B = 0.2934089823272757D+00
       CALL GEN_IH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE ID0212(X,Y,Z,W,N)
       DOUBLE PRECISION X( 212)
       DOUBLE PRECISION Y( 212)
       DOUBLE PRECISION Z( 212)
       DOUBLE PRECISION W( 212)
       INTEGER N
       INTEGER I
       DOUBLE PRECISION A,B,V
       N=1
       V = 0.3741936943619837D-02
       CALL GEN_IH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V = 0.4935771576102468D-02
       CALL GEN_IH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V = 0.4551001693496304D-02
       A = 0.8452548782475415D+00
       B = 0.4824136177694623D+00
       CALL GEN_IH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V = 0.4855145722266675D-02
       A = 0.6779115110470655D+00
       B = 0.6479557053061163D+00
       CALL GEN_IH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V = 0.4866874670145564D-02
       A = 0.7753898899839309D+00
       B = 0.4269810983100433D+00
       CALL GEN_IH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
       SUBROUTINE ID0242(X,Y,Z,W,N)
       DOUBLE PRECISION X( 242)
       DOUBLE PRECISION Y( 242)
       DOUBLE PRECISION Z( 242)
       DOUBLE PRECISION W( 242)
       INTEGER N
       INTEGER I
       DOUBLE PRECISION A,B,V
       N=1
       V = 0.3773233796543805D-02
       CALL GEN_IH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V = 0.4742703016240548D-02
       CALL GEN_IH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V = 0.2350527849974790D-02
       CALL GEN_IH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V = 0.3976391306805966D-02
       A = 0.7588635353153991D+00
       B = 0.4631578933426286D+00
       CALL GEN_IH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V = 0.4548524050513189D-02
       A = 0.8418750816582262D+00
       B = 0.4876800613655622D+00
       CALL GEN_IH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V = 0.4630939619637840D-02
       A = 0.6718027690441337D+00
       B = 0.6591034254439872D+00
       CALL GEN_IH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
       END
