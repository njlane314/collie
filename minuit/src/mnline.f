*
* $Id: mnline.f,v 1.2 2010/07/02 18:33:52 wfisher Exp $
*
* $Log: mnline.f,v $
* Revision 1.2  2010/07/02 18:33:52  wfisher
* Added rebinning algo to CollieIOfile
*
* Revision 1.1  2010/02/16 20:35:21  wfisher
*
* First commit for V00-04-00
*
* Revision 1.1.1.1  1996/03/07 14:31:30  mclareni
* Minuit
*
*
      SUBROUTINE MNLINE(FCN,START,FSTART,STEP,SLOPE,TOLER,FUTIL)
*
* $Id: mnline.f,v 1.2 2010/07/02 18:33:52 wfisher Exp $
*
* $Log: mnline.f,v $
* Revision 1.2  2010/07/02 18:33:52  wfisher
* Added rebinning algo to CollieIOfile
*
* Revision 1.1  2010/02/16 20:35:21  wfisher
*
* First commit for V00-04-00
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506dp.inc
*
C ************ DOUBLE PRECISION VERSION *************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CC        Perform a line search from position START
CC        along direction STEP, where the length of vector STEP
CC                   gives the expected position of minimum.
CC        FSTART is value of function at START
CC        SLOPE (if non-zero) is df/dx along STEP at START
CC        TOLER is initial tolerance of minimum in direction STEP
*
* $Id: mnline.f,v 1.2 2010/07/02 18:33:52 wfisher Exp $
*
* $Log: mnline.f,v $
* Revision 1.2  2010/07/02 18:33:52  wfisher
* Added rebinning algo to CollieIOfile
*
* Revision 1.1  2010/02/16 20:35:21  wfisher
*
* First commit for V00-04-00
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
*
*
* d506cm.inc
*
      PARAMETER (MNE=500 , MNI=475)
      PARAMETER (MNIHL=MNI*(MNI+1)/2)
      CHARACTER*10 CPNAM
      COMMON
     1/MN7NAM/ CPNAM(MNE)
     2/MN7EXT/ U(MNE)     ,ALIM(MNE)  ,BLIM(MNE)
     3/MN7ERR/ ERP(MNI)   ,ERN(MNI)   ,WERR(MNI)  ,GLOBCC(MNI)
     4/MN7INX/ NVARL(MNE) ,NIOFEX(MNE),NEXOFI(MNI)
     5/MN7INT/ X(MNI)     ,XT(MNI)    ,DIRIN(MNI)
     6/MN7FX2/ XS(MNI)    ,XTS(MNI)   ,DIRINS(MNI)
     7/MN7DER/ GRD(MNI)   ,G2(MNI)    ,GSTEP(MNI) ,GIN(MNE) ,DGRD(MNI)
     8/MN7FX3/ GRDS(MNI)  ,G2S(MNI)   ,GSTEPS(MNI)
     9/MN7FX1/ IPFIX(MNI) ,NPFIX
     A/MN7VAR/ VHMAT(MNIHL)
     B/MN7VAT/ VTHMAT(MNIHL)
     C/MN7SIM/ P(MNI,MNI+1),PSTAR(MNI),PSTST(MNI) ,PBAR(MNI),PRHO(MNI)
C
      PARAMETER (MAXDBG=10, MAXSTK=10, MAXCWD=20, MAXP=30, MAXCPT=101)
      PARAMETER (ZERO=0.0,  ONE=1.0,   HALF=0.5)
      COMMON
     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
     E/MN7IOU/ ISYSRD ,ISYSWR ,ISYSSA ,NPAGWD ,NPAGLN ,NEWPAG
     E/MN7IO2/ ISTKRD(MAXSTK) ,NSTKRD ,ISTKWR(MAXSTK) ,NSTKWR
     F/MN7TIT/ CFROM  ,CSTATU ,CTITL  ,CWORD  ,CUNDEF ,CVRSN ,COVMES
     G/MN7FLG/ ISW(7) ,IDBG(0:MAXDBG) ,NBLOCK ,ICOMND
     H/MN7MIN/ AMIN   ,UP     ,EDM    ,FVAL3  ,EPSI   ,APSI  ,DCOVAR
     I/MN7CNV/ NFCN   ,NFCNMX ,NFCNLC ,NFCNFR ,ITAUR,ISTRAT,NWRMES(2)
     J/MN7ARG/ WORD7(MAXP)
     K/MN7LOG/ LWARN  ,LREPOR ,LIMSET ,LNOLIM ,LNEWMN ,LPHEAD
     L/MN7CNS/ EPSMAC ,EPSMA2 ,VLIMLO ,VLIMHI ,UNDEFI ,BIGEDM,UPDFLT
     M/MN7RPT/ XPT(MAXCPT)    ,YPT(MAXCPT)
     N/MN7CPT/ CHPT(MAXCPT)
     o/MN7XCR/ XMIDCR ,YMIDCR ,XDIRCR ,YDIRCR ,KE1CR  ,KE2CR
      CHARACTER CTITL*50, CWORD*(MAXCWD), CUNDEF*10, CFROM*8,
     +          CVRSN*6,  COVMES(0:3)*22, CSTATU*10, CHPT*1
      LOGICAL   LWARN, LREPOR, LIMSET, LNOLIM, LNEWMN, LPHEAD
      EXTERNAL FCN,FUTIL
      DIMENSION START(*), STEP(*)
      PARAMETER (MAXPT=12)
      DIMENSION XPQ(MAXPT),YPQ(MAXPT)
      CHARACTER*1 CHPQ(MAXPT)
      DIMENSION XVALS(3),FVALS(3),COEFF(3)
      CHARACTER*26 CHARAL
      CHARACTER*60 CMESS
      PARAMETER (SLAMBG=5.,ALPHA=2.)
C SLAMBG and ALPHA control the maximum individual steps allowed.
C The first step is always =1. The max length of second step is SLAMBG.
C The max size of subsequent steps is the maximum previous successful
C   step multiplied by ALPHA + the size of most recent successful step,
C   but cannot be smaller than SLAMBG.
      LOGICAL LDEBUG
      DATA CHARAL / 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' /
      LDEBUG = (IDBG(1).GE.1)
C                  starting values for overall limits on total step SLAM
      OVERAL = 1000.
      UNDRAL = -100.
C                              debug check if start is ok
      IF (LDEBUG)  THEN
         CALL MNINEX(START)
         CALL FCN(NPARX,GIN,F1,U,4,FUTIL)
         NFCN=NFCN+1
         IF (F1 .NE. FSTART) THEN
             WRITE (ISYSWR,'(A/2E14.5/2X,10F10.5)')
     + ' MNLINE start point not consistent, F values, parameters=',
     +  (X(KK),KK=1,NPAR)
         ENDIF
      ENDIF
C                                      . set up linear search along STEP
      FVMIN = FSTART
      XVMIN = ZERO
      NXYPT = 1
      CHPQ(1) = CHARAL(1:1)
      XPQ(1) = 0.
      YPQ(1) = FSTART
C               SLAMIN = smallest possible value of ABS(SLAM)
      SLAMIN = ZERO
      DO 20 I= 1, NPAR
      IF (STEP(I) .EQ. ZERO)  GO TO 20
      RATIO = ABS(START(I)/STEP(I))
      IF (SLAMIN .EQ. ZERO)     SLAMIN = RATIO
      IF (RATIO .LT. SLAMIN)  SLAMIN = RATIO
   20 X(I) = START(I) + STEP(I)
      IF (SLAMIN .EQ. ZERO)  SLAMIN = EPSMAC
      SLAMIN = SLAMIN*EPSMA2
      NPARX = NPAR
C
      CALL MNINEX(X)
      CALL FCN(NPARX,GIN,F1,U,4,FUTIL)
      NFCN=NFCN+1
      NXYPT = NXYPT + 1
      CHPQ(NXYPT) = CHARAL(NXYPT:NXYPT)
      XPQ(NXYPT) = 1.
      YPQ(NXYPT) = F1
      IF (F1 .LT. FSTART) THEN
         FVMIN = F1
         XVMIN = 1.0
      ENDIF
C                         . quadr interp using slope GDEL and two points
      SLAM = 1.
      TOLER8 = TOLER
      SLAMAX = SLAMBG
      FLAST = F1
C                         can iterate on two-points (cut) if no imprvmnt
   25 CONTINUE
      DENOM = 2.0*(FLAST-FSTART-SLOPE*SLAM)/SLAM**2
C     IF (DENOM .EQ. ZERO)  DENOM = -0.1*SLOPE
                            SLAM  = 1.
      IF (DENOM .NE. ZERO)  SLAM = -SLOPE/DENOM
      IF (SLAM  .LT. ZERO)  SLAM = SLAMAX
      IF (SLAM .GT. SLAMAX)  SLAM = SLAMAX
      IF (SLAM .LT. TOLER8)  SLAM = TOLER8
      IF (SLAM .LT. SLAMIN)  GO TO 80
      IF (ABS(SLAM-1.0).LT.TOLER8 .AND. F1.LT.FSTART)  GO TO 70
      IF (ABS(SLAM-1.0).LT.TOLER8) SLAM = 1.0+TOLER8
      IF (NXYPT .GE. MAXPT) GO TO 65
      DO 30 I= 1, NPAR
   30 X(I) = START(I) + SLAM*STEP(I)
      CALL MNINEX(X)
      CALL FCN(NPAR,GIN,F2,U,4,FUTIL)
      NFCN = NFCN + 1
      NXYPT = NXYPT + 1
      CHPQ(NXYPT) = CHARAL(NXYPT:NXYPT)
      XPQ(NXYPT) = SLAM
      YPQ(NXYPT) = F2
      IF (F2 .LT. FVMIN)  THEN
         FVMIN = F2
         XVMIN = SLAM
      ENDIF
      IF (FSTART .EQ. FVMIN) THEN
         FLAST = F2
         TOLER8 = TOLER*SLAM
         OVERAL = SLAM-TOLER8
         SLAMAX = OVERAL
         GO TO 25
      ENDIF
C                                        . quadr interp using 3 points
      XVALS(1) = XPQ(1)
      FVALS(1) = YPQ(1)
      XVALS(2) = XPQ(NXYPT-1)
      FVALS(2) = YPQ(NXYPT-1)
      XVALS(3) = XPQ(NXYPT)
      FVALS(3) = YPQ(NXYPT)
C                             begin iteration, calculate desired step
   50 CONTINUE
      SLAMAX = MAX(SLAMAX,ALPHA*ABS(XVMIN))
      CALL MNPFIT(XVALS,FVALS,3,COEFF,SDEV)
      IF (COEFF(3) .LE. ZERO)  THEN
         SLOPEM = 2.0*COEFF(3)*XVMIN + COEFF(2)
         IF (SLOPEM .LE. ZERO) THEN
            SLAM = XVMIN + SLAMAX
         ELSE
            SLAM = XVMIN - SLAMAX
         ENDIF
      ELSE
         SLAM = -COEFF(2)/(2.0*COEFF(3))
         IF (SLAM .GT. XVMIN+SLAMAX)  SLAM = XVMIN+SLAMAX
         IF (SLAM .LT. XVMIN-SLAMAX)  SLAM = XVMIN-SLAMAX
      ENDIF
      IF (SLAM .GT. ZERO) THEN
          IF (SLAM .GT. OVERAL) SLAM = OVERAL
      ELSE
          IF (SLAM .LT. UNDRAL) SLAM = UNDRAL
      ENDIF
C               come here if step was cut below
   52 CONTINUE
      TOLER9 = MAX(TOLER8,ABS(TOLER8*SLAM))
      DO 55 IPT= 1, 3
      IF (ABS(SLAM-XVALS(IPT)) .LT. TOLER9)  GO TO 70
   55 CONTINUE
C                take the step
      IF (NXYPT .GE. MAXPT) GO TO 65
      DO 60 I= 1, NPAR
   60 X(I) = START(I)+SLAM*STEP(I)
      CALL MNINEX(X)
      CALL FCN(NPARX,GIN,F3,U,4,FUTIL)
      NFCN = NFCN + 1
      NXYPT = NXYPT + 1
      CHPQ(NXYPT) = CHARAL(NXYPT:NXYPT)
      XPQ(NXYPT) = SLAM
      YPQ(NXYPT) = F3
C             find worst previous point out of three
      FVMAX = FVALS(1)
      NVMAX = 1
      IF (FVALS(2) .GT. FVMAX) THEN
         FVMAX = FVALS(2)
         NVMAX = 2
      ENDIF
      IF (FVALS(3) .GT. FVMAX) THEN
         FVMAX = FVALS(3)
         NVMAX = 3
      ENDIF
C              if latest point worse than all three previous, cut step
      IF (F3 .GE. FVMAX)  THEN
          IF (NXYPT .GE. MAXPT) GO TO 65
          IF (SLAM .GT. XVMIN) OVERAL = MIN(OVERAL,SLAM-TOLER8)
          IF (SLAM .LT. XVMIN) UNDRAL = MAX(UNDRAL,SLAM+TOLER8)
          SLAM = 0.5*(SLAM+XVMIN)
          GO TO 52
      ENDIF
C              prepare another iteration, replace worst previous point
      XVALS(NVMAX) = SLAM
      FVALS(NVMAX) = F3
      IF (F3 .LT. FVMIN)  THEN
         FVMIN = F3
         XVMIN = SLAM
      ELSE
         IF (SLAM .GT. XVMIN) OVERAL = MIN(OVERAL,SLAM-TOLER8)
         IF (SLAM .LT. XVMIN) UNDRAL = MAX(UNDRAL,SLAM+TOLER8)
      ENDIF
      IF (NXYPT .LT. MAXPT)  GO TO 50
C                                            . . end of iteration . . .
C            stop because too many iterations
   65 CMESS = ' LINE SEARCH HAS EXHAUSTED THE LIMIT OF FUNCTION CALLS '
      IF (LDEBUG) THEN
        WRITE (ISYSWR,'(A/(2X,6G12.4))') ' MNLINE DEBUG: steps=',
     +    (STEP(KK),KK=1,NPAR)
      ENDIF
      GO TO 100
C            stop because within tolerance
   70 CONTINUE
      CMESS = ' LINE SEARCH HAS ATTAINED TOLERANCE '
      GO TO 100
   80 CONTINUE
      CMESS = ' STEP SIZE AT ARITHMETICALLY ALLOWED MINIMUM'
  100 CONTINUE
      AMIN = FVMIN
      DO 120 I= 1, NPAR
      DIRIN(I) = STEP(I)*XVMIN
  120 X(I) = START(I) + DIRIN(I)
      CALL MNINEX(X)
      IF (XVMIN .LT. 0.)      CALL MNWARN('D','MNLINE',
     +                   ' LINE MINIMUM IN BACKWARDS DIRECTION')
      IF (FVMIN .EQ. FSTART)  CALL MNWARN('D','MNLINE',
     +                     ' LINE SEARCH FINDS NO IMPROVEMENT ')
      IF (LDEBUG)  THEN
         WRITE (ISYSWR,'('' AFTER'',I3,'' POINTS,'',A)') NXYPT,CMESS
         CALL MNPLOT(XPQ,YPQ,CHPQ,NXYPT,ISYSWR,NPAGWD,NPAGLN)
      ENDIF
      RETURN
      END
