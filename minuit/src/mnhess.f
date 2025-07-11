*
* $Id: mnhess.f,v 1.2 2010/07/02 18:33:52 wfisher Exp $
*
* $Log: mnhess.f,v $
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
      SUBROUTINE MNHESS(FCN,FUTIL)
*
* $Id: mnhess.f,v 1.2 2010/07/02 18:33:52 wfisher Exp $
*
* $Log: mnhess.f,v $
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
CC        Calculates the full second-derivative matrix of FCN
CC        by taking finite differences. When calculating diagonal
CC        elements, it may iterate so that step size is nearly that
CC        which gives function change= UP/10. The first derivatives
CC        of course come as a free side effect, but with a smaller
CC        step size in order to obtain a known accuracy.
CC
*
* $Id: mnhess.f,v 1.2 2010/07/02 18:33:52 wfisher Exp $
*
* $Log: mnhess.f,v $
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
      DIMENSION YY(MNI)
      LOGICAL LDEBUG
      CHARACTER CBF1*22
C
      LDEBUG = (IDBG(3) .GE. 1)
      IF (AMIN .EQ. UNDEFI)  CALL MNAMIN(FCN,FUTIL)
      IF (ISTRAT .LE. 0) THEN
         NCYC = 3
         TLRSTP = 0.5
         TLRG2  = 0.1
      ELSE IF (ISTRAT .EQ. 1) THEN
         NCYC = 5
         TLRSTP = 0.3
         TLRG2  = 0.05
      ELSE
         NCYC = 7
         TLRSTP = 0.1
         TLRG2  = 0.02
      ENDIF
      IF (ISW(5).GE.2 .OR. LDEBUG)  WRITE (ISYSWR,'(A)')
     +   '   START COVARIANCE MATRIX CALCULATION.'
      CFROM = 'HESSE   '
      NFCNFR = NFCN
      CSTATU= 'OK        '
      NPARD = NPAR
C                 make sure starting at the right place
      CALL MNINEX(X)
      NPARX = NPAR
      CALL FCN(NPARX,GIN,FS1,U,4,FUTIL)
      NFCN = NFCN + 1
      IF (FS1 .NE. AMIN) THEN
         DF = AMIN - FS1
         WRITE (CBF1(1:12),'(G12.3)') DF
         CALL MNWARN('D','MNHESS',
     +       'function value differs from AMIN by '//CBF1(1:12) )
      ENDIF
      AMIN = FS1
      IF (LDEBUG) WRITE (ISYSWR,'(A,A)') ' PAR D   GSTEP          ',
     +' D          G2         GRD         SAG    '
C                                        . . . . . . diagonal elements .
C         ISW(2) = 1 if approx, 2 if not posdef, 3 if ok
C         AIMSAG is the sagitta we are aiming for in second deriv calc.
      AIMSAG = SQRT(EPSMA2)*(ABS(AMIN)+UP)
C         Zero the second derivative matrix
      NPAR2 = NPAR*(NPAR+1)/2
      DO 10 I= 1,NPAR2
   10 VHMAT(I) = 0.
C
C         Loop over variable parameters for second derivatives
      IDRV = 2
      DO 100 ID= 1, NPARD
      I = ID + NPAR - NPARD
      IEXT = NEXOFI(I)
      IF (G2(I) .EQ. ZERO) THEN
           WRITE (CBF1(1:4),'(I4)') IEXT
           CALL MNWARN('W','HESSE',
     +      'Second derivative enters zero, param '//CBF1(1:4) )
        WINT = WERR(I)
        IF (NVARL(IEXT) .GT. 1) THEN
           CALL MNDXDI(X(I),I,DXDI)
           IF (ABS(DXDI) .LT. .001) THEN
              WINT = .01
           ELSE
              WINT = WINT/ABS(DXDI)
           ENDIF
        ENDIF
        G2(I) = UP/WINT**2
      ENDIF
      XTF = X(I)
      DMIN = 8.*EPSMA2*ABS(XTF)
C
C                               find step which gives sagitta = AIMSAG
      D = ABS(GSTEP(I))
      DO 40 ICYC= 1, NCYC
C                               loop here only if SAG=0.
      DO 25 MULTPY= 1, 5
C           take two steps
         X(I) = XTF + D
         CALL MNINEX(X)
         NPARX = NPAR
         CALL FCN(NPARX,GIN,FS1,U,4,FUTIL)
         NFCN = NFCN + 1
         X(I) = XTF - D
         CALL MNINEX(X)
         CALL FCN(NPARX,GIN,FS2,U,4,FUTIL)
         NFCN = NFCN + 1
         X(I) = XTF
         SAG = 0.5*(FS1+FS2-2.0*AMIN)
         IF (SAG .NE. ZERO) GO TO 30
         IF (GSTEP(I) .LT. ZERO) THEN
           IF (D .GE. .5)  GO TO 26
           D = 10.*D
           IF (D .GT. 0.5)  D = 0.51
           GO TO 25
         ENDIF
         D = 10.*D
   25 CONTINUE
   26      WRITE (CBF1(1:4),'(I4)') IEXT
           CALL MNWARN('W','HESSE',
     +      'Second derivative zero for parameter'//CBF1(1:4) )
           GO TO 390
C                             SAG is not zero
   30 G2BFOR = G2(I)
      G2(I) = 2.*SAG/D**2
      GRD(I) = (FS1-FS2)/(2.*D)
      IF (LDEBUG) WRITE (ISYSWR,31) I,IDRV,GSTEP(I),D,G2(I),GRD(I),SAG
   31 FORMAT (I4,I2,6G12.5)
      GSTEP(I) = SIGN(D,GSTEP(I))
      DIRIN(I) = D
      YY(I) = FS1
      DLAST = D
      D = SQRT(2.0*AIMSAG/ABS(G2(I)))
C         if parameter has limits, max int step size = 0.5
      STPINM = 0.5
      IF (GSTEP(I) .LT. ZERO)  D = MIN(D,STPINM)
      IF (D .LT. DMIN)  D = DMIN
C           see if converged
      IF (ABS((D-DLAST)/D)          .LT. TLRSTP)  GO TO 50
      IF (ABS((G2(I)-G2BFOR)/G2(I)) .LT. TLRG2 )  GO TO 50
      D = MIN(D, 10.*DLAST)
      D = MAX(D, 0.1*DLAST)
   40 CONTINUE
C                       end of step size loop
      WRITE (CBF1,'(I2,2E10.2)') IEXT,SAG,AIMSAG
      CALL MNWARN('D','MNHESS','Second Deriv. SAG,AIM= '//CBF1)
C
   50 CONTINUE
      NDEX = I*(I+1)/2
      VHMAT(NDEX) = G2(I)
  100 CONTINUE
C                              end of diagonal second derivative loop
      CALL MNINEX(X)
C                                     refine the first derivatives
      IF (ISTRAT .GT. 0) CALL MNHES1(FCN,FUTIL)
      ISW(2) = 3
      DCOVAR = 0.
C                                        . . . .  off-diagonal elements
      IF (NPAR .EQ. 1)  GO TO 214
      DO 200 I= 1, NPAR
      DO 180 J= 1, I-1
      XTI = X(I)
      XTJ = X(J)
      X(I) = XTI + DIRIN(I)
      X(J) = XTJ + DIRIN(J)
      CALL MNINEX(X)
      CALL FCN(NPARX,GIN,FS1,U,4,FUTIL)
      NFCN = NFCN + 1
      X(I) = XTI
      X(J) = XTJ
      ELEM = (FS1+AMIN-YY(I)-YY(J)) / (DIRIN(I)*DIRIN(J))
      NDEX = I*(I-1)/2 + J
      VHMAT(NDEX) = ELEM
  180 CONTINUE
  200 CONTINUE
  214 CALL MNINEX(X)
C                  verify matrix positive-definite
      CALL MNPSDF
      DO 220 I= 1, NPAR
      DO 219 J= 1, I
      NDEX = I*(I-1)/2 + J
      P(I,J) = VHMAT(NDEX)
  219 P(J,I) = P(I,J)
  220 CONTINUE
      CALL MNVERT(P,MAXINT,MAXINT,NPAR,IFAIL)
      IF (IFAIL .GT. 0)  THEN
        CALL MNWARN('W','HESSE', 'Matrix inversion fails.')
        GO TO 390
      ENDIF
C                                        . . . . . . .  calculate  e d m
      EDM = 0.
        DO 230 I= 1, NPAR
C                              off-diagonal elements
        NDEX = I*(I-1)/2
          DO 225 J= 1, I-1
          NDEX = NDEX + 1
          ZTEMP = 2.0 * P(I,J)
          EDM = EDM + GRD(I)*ZTEMP*GRD(J)
  225     VHMAT(NDEX) = ZTEMP
C                              diagonal elements
        NDEX = NDEX + 1
        VHMAT(NDEX) = 2.0 * P(I,I)
        EDM = EDM  + P(I,I) * GRD(I)**2
  230   CONTINUE
      IF (ISW(5).GE.1 .AND. ISW(2).EQ.3 .AND. ITAUR.EQ.0)
     + WRITE(ISYSWR,'(A)')' COVARIANCE MATRIX CALCULATED SUCCESSFULLY'
      GO TO 900
C                              failure to invert 2nd deriv matrix
  390 ISW(2) = 1
      DCOVAR = 1.
      CSTATU = 'FAILED    '
      IF (ISW(5) .GE. 0) WRITE (ISYSWR,'(A)')
     +        '  MNHESS FAILS AND WILL RETURN DIAGONAL MATRIX. '
      DO 395 I= 1, NPAR
      NDEX = I*(I-1)/2
      DO 394 J= 1, I-1
      NDEX = NDEX + 1
  394 VHMAT(NDEX) = 0.0
      NDEX = NDEX +1
      G2I = G2(I)
      IF (G2I .LE. ZERO)  G2I = 1.0
  395 VHMAT(NDEX) = 2.0/G2I
  900 RETURN
      END
