*
* $Id: mnparm.f,v 1.2 2010/07/02 18:33:52 wfisher Exp $
*
* $Log: mnparm.f,v $
* Revision 1.2  2010/07/02 18:33:52  wfisher
* Added rebinning algo to CollieIOfile
*
* Revision 1.1  2010/02/16 20:35:21  wfisher
*
* First commit for V00-04-00
*
* Revision 1.2  1996/03/15 18:02:50  james
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
* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni
* Minuit
*
*
      SUBROUTINE MNPARM(K,CNAMJ,UK,WK,A,B,IERFLG)
*
* $Id: mnparm.f,v 1.2 2010/07/02 18:33:52 wfisher Exp $
*
* $Log: mnparm.f,v $
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
CC        Called from MNPARS and user-callable
CC    Implements one parameter definition, that is:
CC          K     (external) parameter number
CC          CNAMK parameter name
CC          UK    starting value
CC          WK    starting step size or uncertainty
CC          A, B  lower and upper physical parameter limits
CC    and sets up (updates) the parameter lists.
CC    Output: IERFLG=0 if no problems
CC                  >0 if MNPARM unable to implement definition
CC
*
* $Id: mnparm.f,v 1.2 2010/07/02 18:33:52 wfisher Exp $
*
* $Log: mnparm.f,v $
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
      CHARACTER*(*) CNAMJ
      CHARACTER  CNAMK*10, CHBUFI*4
C
      CNAMK = CNAMJ
      KINT = NPAR
      IF (K.LT.1 .OR. K.GT.MAXEXT) THEN
C                     parameter number exceeds allowed maximum value
        WRITE (ISYSWR,9)  K,MAXEXT
    9   FORMAT (/' MINUIT USER ERROR.  PARAMETER NUMBER IS',I11/
     +         ',  ALLOWED RANGE IS ONE TO',I4/)
        GO TO 800
      ENDIF
C                     normal parameter request
      KTOFIX = 0
      IF (NVARL(K) .LT. 0) GO TO 50
C         previously defined parameter is being redefined
C                                     find if parameter was fixed
      DO 40 IX= 1, NPFIX
      IF (IPFIX(IX) .EQ. K)  KTOFIX = K
   40 CONTINUE
      IF (KTOFIX .GT. 0)  THEN
         CALL MNWARN('W','PARAM DEF','REDEFINING A FIXED PARAMETER.')
         IF (KINT .GE. MAXINT)  THEN
            WRITE (ISYSWR,'(A)') ' CANNOT RELEASE. MAX NPAR EXCEEDED.'
            GO TO 800
            ENDIF
         CALL MNFREE(-K)
         ENDIF
C                       if redefining previously variable parameter
      IF(NIOFEX(K) .GT. 0) KINT = NPAR-1
   50 CONTINUE
C
C                                      . . .print heading
      IF (LPHEAD .AND. ISW(5).GE.0) THEN
        WRITE (ISYSWR,61)
        LPHEAD = .FALSE.
      ENDIF
   61 FORMAT(/' PARAMETER DEFINITIONS:'/
     +        '    NO.   NAME         VALUE      STEP SIZE      LIMITS')
      IF (WK .GT. ZERO)  GO TO 122
C                                        . . .constant parameter . . . .
      IF (ISW(5) .GE. 0)  WRITE (ISYSWR, 82)  K,CNAMK,UK
   82 FORMAT (1X,I5,1X,1H',A10,1H',1X,G13.5, '  constant')
      NVL = 0
      GO TO 200
  122 IF (A.EQ.ZERO .AND. B.EQ.ZERO) THEN
C                                      variable parameter without limits
      NVL = 1
      IF (ISW(5) .GE. 0)  WRITE (ISYSWR, 127)  K,CNAMK,UK,WK
  127 FORMAT (1X,I5,1X,1H',A10,1H',1X,2G13.5, '     no limits')
      ELSE
C                                         variable parameter with limits
      NVL = 4
      LNOLIM = .FALSE.
      IF (ISW(5) .GE. 0)  WRITE (ISYSWR, 132)  K,CNAMK,UK,WK,A,B
  132 FORMAT(1X,I5,1X,1H',A10,1H',1X,2G13.5,2X,2G13.5)
      ENDIF
C                             . . request for another variable parameter
      KINT = KINT + 1
      IF (KINT .GT. MAXINT)  THEN
         WRITE (ISYSWR,135)  MAXINT
  135    FORMAT (/' MINUIT USER ERROR.   TOO MANY VARIABLE PARAMETERS.'/
     +   ' THIS VERSION OF MINUIT DIMENSIONED FOR',I4//)
         GO TO 800
         ENDIF
      IF (NVL .EQ. 1)  GO TO 200
      IF (A .EQ. B)  THEN
        WRITE (ISYSWR,'(/A,A/A/)') ' USER ERROR IN MINUIT PARAMETER',
     +   ' DEFINITION',' UPPER AND LOWER LIMITS EQUAL.'
        GO TO 800
        ENDIF
      IF (B .LT. A) THEN
         SAV = B
         B = A
         A = SAV
         CALL MNWARN('W','PARAM DEF','PARAMETER LIMITS WERE REVERSED.')
         IF (LWARN) LPHEAD=.TRUE.
         ENDIF
      IF ((B-A) .GT. 1.0E7)  THEN
         WRITE (CHBUFI,'(I4)') K
         CALL MNWARN('W','PARAM DEF',
     +               'LIMITS ON PARAM'//CHBUFI//' TOO FAR APART.')
         IF (LWARN) LPHEAD=.TRUE.
      ENDIF
      DANGER = (B-UK)*(UK-A)
      IF (DANGER .LT. 0.)
     +     CALL MNWARN('W','PARAM DEF','STARTING VALUE OUTSIDE LIMITS.')
      IF (DANGER .EQ. 0.)
     +     CALL MNWARN('W','PARAM DEF','STARTING VALUE IS AT LIMIT.')
  200 CONTINUE
C                           . . . input OK, set values, arrange lists,
C                                    calculate step sizes GSTEP, DIRIN
      CFROM = 'PARAMETR'
      NFCNFR = NFCN
      CSTATU= 'NEW VALUES'
      NU = MAX(NU,K)
      CPNAM(K) = CNAMK
      U(K) = UK
      ALIM(K) = A
      BLIM(K) = B
      NVARL(K) = NVL
C                             K is external number of new parameter
C           LASTIN is the number of var. params with ext. param. no.< K
      LASTIN = 0
      DO 240 IX= 1, K-1
      IF (NIOFEX(IX) .GT. 0)  LASTIN=LASTIN+1
  240 CONTINUE
C                 KINT is new number of variable params, NPAR is old
      IF (KINT .EQ. NPAR)  GO TO 280
      IF (KINT .GT. NPAR) THEN
C                          insert new variable parameter in list
         DO 260 IN= NPAR,LASTIN+1,-1
         IX = NEXOFI(IN)
         NIOFEX(IX) = IN+1
         NEXOFI(IN+1)= IX
         X    (IN+1) = X    (IN)
         XT   (IN+1) = XT   (IN)
         DIRIN(IN+1) = DIRIN(IN)
         G2   (IN+1) = G2   (IN)
         GSTEP(IN+1) = GSTEP(IN)
  260    CONTINUE
      ELSE
C                          remove variable parameter from list
         DO 270 IN= LASTIN+1,KINT
         IX = NEXOFI(IN+1)
         NIOFEX(IX) = IN
         NEXOFI(IN)= IX
         X     (IN)= X    (IN+1)
         XT    (IN)= XT   (IN+1)
         DIRIN (IN)= DIRIN(IN+1)
         G2    (IN)= G2   (IN+1)
         GSTEP (IN)= GSTEP(IN+1)
  270    CONTINUE
      ENDIF
  280 CONTINUE
      IX = K
      NIOFEX(IX) = 0
      NPAR = KINT
      CALL MNRSET(1)
C                                       lists are now arranged . . . .
      IF (NVL .GT. 0)  THEN
         IN = LASTIN+1
         NEXOFI(IN) = IX
         NIOFEX(IX) = IN
         SAV = U(IX)
         CALL MNPINT(SAV,IX,PINTI)
         X(IN) = PINTI
         XT(IN) = X(IN)
         WERR(IN) = WK
         SAV2 = SAV + WK
         CALL MNPINT(SAV2,IX,PINTI)
         VPLU = PINTI - X(IN)
         SAV2 = SAV - WK
         CALL MNPINT(SAV2,IX,PINTI)
         VMINU = PINTI - X(IN)
         DIRIN(IN) = 0.5 * (ABS(VPLU) +ABS(VMINU))
         G2(IN) = 2.0*UP / DIRIN(IN)**2
         GSMIN = 8.*EPSMA2*ABS(X(IN))
         GSTEP(IN) = MAX (GSMIN, 0.1*DIRIN(IN))
         IF (AMIN .NE. UNDEFI) THEN
             SMALL = SQRT(EPSMA2*(AMIN+UP)/UP)
             GSTEP(IN) = MAX(GSMIN, SMALL*DIRIN(IN))
         ENDIF
         GRD  (IN) = G2(IN)*DIRIN(IN)
C                   if parameter has limits
         IF (NVARL(K) .GT. 1) THEN
            IF (GSTEP(IN).GT. 0.5)  GSTEP(IN)=0.5
            GSTEP(IN) = -GSTEP(IN)
         ENDIF
      ENDIF
      IF (KTOFIX .GT. 0)  THEN
         KINFIX = NIOFEX(KTOFIX)
         IF (KINFIX .GT. 0)  CALL MNFIXP(KINFIX,IERR)
         IF (IERR .GT. 0)  GO TO 800
      ENDIF
      IERFLG = 0
      RETURN
C                   error on input, unable to implement request  . . . .
  800 CONTINUE
      IERFLG = 1
      RETURN
      END
