*
* $Id: mnderi.f,v 1.2 2010/07/02 18:33:52 wfisher Exp $
*
* $Log: mnderi.f,v $
* Revision 1.2  2010/07/02 18:33:52  wfisher
* Added rebinning algo to CollieIOfile
*
* Revision 1.1  2010/02/16 20:35:20  wfisher
*
* First commit for V00-04-00
*
* Revision 1.2  1996/03/15 18:02:43  james
*     Modified Files:
* mnderi.F eliminate possible division by zero
* mnexcm.F suppress print on STOP when print flag=-1
*          set FVAL3 to flag if FCN already called with IFLAG=3
* mninit.F set version 96.03
* mnlims.F remove arguments, not needed
* mnmigr.F VLEN -> LENV in debug print statement
* mnparm.F move call to MNRSET to after NPAR redefined, to zero all
* mnpsdf.F eliminate possible division by zero
* mnscan.F suppress printout when print flag =-1
* mnset.F  remove arguments in call to MNLIMS
* mnsimp.F fix CSTATU so status is PROGRESS only if new minimum
* mnvert.F eliminate possible division by zero
*
* Revision 1.1.1.1  1996/03/07 14:31:29  mclareni
* Minuit
*
*
      SUBROUTINE MNDERI(FCN,FUTIL)
*
* $Id: mnderi.f,v 1.2 2010/07/02 18:33:52 wfisher Exp $
*
* $Log: mnderi.f,v $
* Revision 1.2  2010/07/02 18:33:52  wfisher
* Added rebinning algo to CollieIOfile
*
* Revision 1.1  2010/02/16 20:35:20  wfisher
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
CC        Calculates the first derivatives of FCN (GRD),
CC        either by finite differences or by transforming the user-
CC        supplied derivatives to internal coordinates,
CC        according to whether ISW(3) is zero or one.
CC
*
* $Id: mnderi.f,v 1.2 2010/07/02 18:33:52 wfisher Exp $
*
* $Log: mnderi.f,v $
* Revision 1.2  2010/07/02 18:33:52  wfisher
* Added rebinning algo to CollieIOfile
*
* Revision 1.1  2010/02/16 20:35:20  wfisher
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
      LOGICAL LDEBUG
      CHARACTER CBF1*22
      NPARX = NPAR
      LDEBUG = (IDBG(2) .GE. 1)
      IF (AMIN .EQ. UNDEFI)  CALL MNAMIN(FCN,FUTIL)
      IF (ISW(3) .EQ. 1)  GO TO 100
      IF (LDEBUG) THEN
C                       make sure starting at the right place
        CALL MNINEX(X)
        NPARX = NPAR
        CALL FCN(NPARX,GIN,FS1,U,4,FUTIL)
        NFCN = NFCN + 1
        IF (FS1 .NE. AMIN) THEN
           DF = AMIN - FS1
           WRITE (CBF1(1:12),'(G12.3)') DF
           CALL MNWARN('D','MNDERI',
     +         'function value differs from AMIN by '//CBF1(1:12) )
           AMIN = FS1
        ENDIF
          WRITE
     +   (ISYSWR,'(/''  FIRST DERIVATIVE DEBUG PRINTOUT.  MNDERI''/
     +   '' PAR    DERIV     STEP      MINSTEP   OPTSTEP '',
     +   '' D1-D2    2ND DRV'')')
      ENDIF
      DFMIN = 8. * EPSMA2*(ABS(AMIN)+UP)
      VRYSML = 8.* EPSMAC**2
      IF (ISTRAT .LE. 0) THEN
         NCYC = 2
         TLRSTP = 0.5
         TLRGRD = 0.1
      ELSE IF (ISTRAT .EQ. 1) THEN
         NCYC = 3
         TLRSTP = 0.3
         TLRGRD = 0.05
      ELSE
         NCYC = 5
         TLRSTP = 0.1
         TLRGRD = 0.02
      ENDIF
C                                loop over variable parameters
      DO 60  I=1,NPAR
      EPSPRI = EPSMA2 + ABS(GRD(I)*EPSMA2)
C         two-point derivatives always assumed necessary
C         maximum number of cycles over step size depends on strategy
      XTF = X(I)
      STEPB4 = 0.
C                               loop as little as possible here!
      DO 45 ICYC= 1, NCYC
C                 ........ theoretically best step
      OPTSTP = SQRT(DFMIN/(ABS(G2(I))+EPSPRI))
C                     step cannot decrease by more than a factor of ten
      STEP = MAX(OPTSTP, ABS(0.1*GSTEP(I)))
C                 but if parameter has limits, max step size = 0.5
      IF (GSTEP(I).LT.ZERO .AND. STEP.GT.0.5)  STEP=0.5
C                 and not more than ten times the previous step
      STPMAX = 10.*ABS(GSTEP(I))
      IF (STEP .GT. STPMAX)  STEP = STPMAX
C                 minimum step size allowed by machine precision
      STPMIN = MAX(VRYSML, 8.*ABS(EPSMA2*X(I)))
      IF (STEP .LT. STPMIN)  STEP = STPMIN
C                 end of iterations if step change less than factor 2
      IF (ABS((STEP-STEPB4)/STEP) .LT. TLRSTP)  GO TO 50
C         take step positive
      GSTEP(I) = SIGN(STEP, GSTEP(I))
      STEPB4 = STEP
      X(I) = XTF + STEP
      CALL MNINEX(X)
      CALL FCN(NPARX,GIN,FS1,U,4,FUTIL)
      NFCN=NFCN+1
C         take step negative
      X(I) = XTF - STEP
      CALL MNINEX(X)
      CALL FCN(NPARX,GIN,FS2,U,4,FUTIL)
      NFCN=NFCN+1
      GRBFOR = GRD(I)
      GRD(I) = (FS1-FS2)/(2.0*STEP)
      G2(I) = (FS1+FS2-2.0*AMIN)/(STEP**2)
      X(I) = XTF
      IF (LDEBUG) THEN
         D1D2 = (FS1+FS2-2.0*AMIN)/STEP
         WRITE (ISYSWR,41) I,GRD(I),STEP,STPMIN,OPTSTP,D1D2,G2(I)
   41    FORMAT (I4,2G11.3,5G10.2)
      ENDIF
C         see if another iteration is necessary
      IF (ABS(GRBFOR-GRD(I))/(ABS(GRD(I))+DFMIN/STEP) .LT. TLRGRD)
     +        GO TO 50
   45 CONTINUE
C                           end of ICYC loop. too many iterations
      IF (NCYC .EQ. 1)  GO TO 50
         WRITE (CBF1,'(2E11.3)')  GRD(I),GRBFOR
         CALL MNWARN('D','MNDERI',
     +         'First derivative not converged. '//CBF1)
   50 CONTINUE
C
   60 CONTINUE
      CALL MNINEX(X)
      RETURN
C                                        .  derivatives calc by fcn
  100 DO 150 IINT= 1, NPAR
      IEXT = NEXOFI(IINT)
      IF (NVARL(IEXT) .GT. 1)  GO TO 120
      GRD(IINT) = GIN(IEXT)
      GO TO 150
  120 DD = (BLIM(IEXT)-ALIM(IEXT))*0.5 *COS(X(IINT))
      GRD(IINT) = GIN(IEXT)*DD
  150 CONTINUE
  200 RETURN
      END
