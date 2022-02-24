C*MODULE DFTGRD  *DECK RADPT
C>    @brief Determine the radial points and weights
C>    @date : March 2019 - Vladimir Mironov
C>            Introduce Mura-Knowles and Treutler-Ahlrichs grids
C
      SUBROUTINE RADPT(PTRAD,WTRAD,NRAD)
      USE params, ONLY: rad_grid_type
      USE mod_grid_storage, ONLY: save_rad_grid
C
      IMPLICIT NONE
C
      DIMENSION PTRAD(NRAD),WTRAD(NRAD)
      DOUBLE PRECISION PTRAD,WTRAD
      INTEGER NRAD
C
      DOUBLE PRECISION ONE, TWO
      PARAMETER (ONE=1.0D+00)
      PARAMETER (TWO=2.0D+00)
C     Alpha0 is scaling factor needed for Log-3 type radial grid by Mura
C     and Knowles:
C     ri = -alpha*log(1-(xi)**3),
C     wi = 3*alpha * (xi*ri)**2 * i/((nrad+1)*(1-xi**3))
C     xi = i*(nrad+1)
C     In the original paper authors suggest alpha=7.0 for 1st and 2nd groups
C     of periodic table and alpha=5.0 otherwise. Here different value
C     is used:
C     alpha=alpha0*Rbs,
C     where Rbs is Bragg-Slater radius. Similar approach is used in NWChem.
C     Alpha0 value was selected using R-B3LYP energy calculations of few
C     small orgnic molecules. It looks reasonable, but more benchmarks are
C     needed to tune this parameter.
C     Note, roots and weights for Log-3 quadrature will be
C     scaled by Rbs and Rbs**3 respecitevely in subroutine 'OCT'.
      DOUBLE PRECISION ALPHA0
      PARAMETER (ALPHA0 = 3.95D+00)
C      PARAMETER (ALPHA0 = 5.8D+00)
C      PARAMETER (ALPHA0 = 7.9D+00)
C     Instead of mapping equidistant grid,
C     one may use Chebyshev roots and weights.
C     However, the latter does not improve the results significantly.
C     To use Chebyshev mapping uncomment CHMK blocks here
C     and in subroutine DFTGRD.
CHMK  PARAMETER (ALPHA0 = 1.0D+00)
CHMK  DOUBLE PRECISION ALPHA1
CHMK  PARAMETER (ALPHA1 = ALPHA0*0.5d0)
C
C     Next, parameters for Treutler-Ahlrichs grid. This grid is based on
C     the following variable transformation:
C     ri = R0/log(2) * (1+(xi))**ta_pow * log(2/(1-xi))
C     to map interval (-1, +1) to (0, +inf). Chebyshev 2nd kind grid is
C     used. R0 are per-atom scaling coefficients. They are specified
C     in TARADS array in GRDDFT and substitute Bragg-Slater radii.
C     Treutler-Ahlrichs grid parameters:
      DOUBLE PRECISION LOG2M1, TA_POW, PI
      PARAMETER (LOG2M1 = 1.0D0/log(2.0D0))
      PARAMETER (TA_POW = 0.6D0)
      PARAMETER (PI = 3.141592653589793238D+00)
      DOUBLE PRECISION xtmp, xnp1, xi, den
      DOUBLE PRECISION t0, t, y, dy, r, tpow, tlog, dr, w0, s
      INTEGER i, irad
C
      IF (rad_grid_type.EQ.0) THEN
C     Murray-Handy-Laming grid
      XTMP=NRAD
      XNP1=XTMP+ONE
      DO 10 IRAD=1,NRAD
         XI=IRAD
C
C     ----- DETERMINE RADIAL POINT -----
C
         PTRAD(IRAD)=(XI*XI)/((XNP1-XI)*(XNP1-XI))
C
C     ----- CALCULATE WEIGHT FOR RADIAL QUADRATURE -----
C
         DEN=(XNP1-XI)*(XNP1-XI)*(XNP1-XI)*(XNP1-XI)
     >        *(XNP1-XI)*(XNP1-XI)*(XNP1-XI)
         WTRAD(IRAD)=(TWO*XNP1*XI*XI*XI*XI*XI)/DEN
C
 10   CONTINUE
C
C     Mura-Knowles Log-3 grid (Chebyshev 2nd kind grid mapping)
CHMK  ELSE IF (rad_grid_type.EQ.1) THEN
CHMK    t0 = PI/(nrad+1.0D0)
CHMK    DO i = 1, nrad
CHMK      t  =    cos(t0*(nrad-i+1)) ! to ensure ascending root order
CHMK      w0 = t0*sin(t0*(nrad-i+1))*0.5d0
CHMK      t  = (t+1.0)*0.5d0
CHMK      y  = 1.0D0-t*t*t
CHMK      dy = 3.0D0 * t*t
CHMK      r  = -ALPHA1*log(y)
CHMK      dr =  ALPHA1*dy/y
CHMK      ptrad(i) = r
CHMK      wtrad(i) = r*r * w0*dr
CHMK    END DO
C
C     Mura-Knowles Log-3 grid (evently spaced grid mapping)
      ELSE IF (rad_grid_type.EQ.1) THEN
        t0 = 1.0D0/(nrad+1.0D0)
        w0 = t0
        DO i = 1, nrad
          t  = i*t0
          y  = 1.0D0-t*t*t
          dy = 3.0D0 * t*t
          r  = -ALPHA0*log(y)
          dr =  ALPHA0*dy/y
          ptrad(i) = r
          wtrad(i) = r*r * w0*dr
        END DO
C
C     Treutler-Ahlrichs grid (Chebyshev 2nd kind grid mapping)
      ELSE IF (rad_grid_type.EQ.2) THEN
        t0 = PI/(nrad+1.0D0)
        DO i = 1, nrad
          t  =    cos(t0*(nrad-i+1)) ! to ensure ascending root order
          w0 = t0*sin(t0*(nrad-i+1))
          tpow = (1.0d0+t)**TA_POW
          tlog = log(2.0d0/(1.0d0-t))
          r  = LOG2M1 * tpow * tlog
          dr = LOG2M1 * (TA_POW * tpow*tlog/(1.0d0+t) + tpow/(1.0d0-t))
          ptrad(i) = r
          wtrad(i) = r*r * w0*dr
        END DO
C
C     Treutler-Ahlrichs grid (evently spaced grid mapping)
C     ELSE IF (rad_grid_type.EQ.2) THEN
C       t0 = 2.0D0/(nrad+1.0D0)
C       w0 = t0
C       DO i = 1, nrad
C         t = i*t0 - 1.0
C         tpow = (1.0+t)**TA_POW
C         tlog = log(2.0d0/(1.0d0-t))
C         r  = LOG2M1 * tpow * tlog
C         dr = LOG2M1*(TA_POW* tpow/(1.0+t)*tlog + tpow/(1.0-t))
C         ptrad(i) = r
C         wtrad(i) = r*r * w0*dr
C       END DO
C
C     Unknown grid type
      ELSE
        WRITE(*,*) 'UNKNOWN RADIAL GRID TYPE =', rad_grid_type
c        CALL abrt
        CALL exit()
      END IF
C
C     Fill in radial grid data to use new DFT code
      CALL save_rad_grid(nrad, ptrad, wtrad)
      RETURN
      END
