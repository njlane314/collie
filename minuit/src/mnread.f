*
* $Id: mnread.f,v 1.2 2010/07/02 18:33:52 wfisher Exp $
*
* $Log: mnread.f,v $
* Revision 1.2  2010/07/02 18:33:52  wfisher
* Added rebinning algo to CollieIOfile
*
* Revision 1.1  2010/02/16 20:35:21  wfisher
*
* First commit for V00-04-00
*
* Revision 1.1.1.1  1996/03/07 14:31:31  mclareni
* Minuit
*
*
      SUBROUTINE MNREAD(FCN,IFLGIN,IFLGUT,FUTIL)
*
* $Id: mnread.f,v 1.2 2010/07/02 18:33:52 wfisher Exp $
*
* $Log: mnread.f,v $
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
CC        Called from MINUIT.  Reads all user input to MINUIT.
CC     This routine is highly unstructured and defies normal logic.
CC
CC     IFLGIN indicates the function originally requested:
CC           = 1: read one-line title
CC             2: read parameter definitions
CC             3: read MINUIT commands
CC
CC     IFLGUT= 1: reading terminated normally
CC             2: end-of-data on input
CC             3: unrecoverable read error
CC             4: unable to process parameter requests
CC             5: more than 100 incomprehensible commands
CC internally,
CC     IFLGDO indicates the subfunction to be performed on the next
CC         input record: 1: read a one-line title
CC                       2: read a parameter definition
CC                       3: read a command
CC                       4: read in covariance matrix
CC     for example, when IFLGIN=3, but IFLGDO=1, then it should read
CC       a title, but this was requested by a command, not by MINUIT.
CC
*
* $Id: mnread.f,v 1.2 2010/07/02 18:33:52 wfisher Exp $
*
* $Log: mnread.f,v $
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
      CHARACTER  CRDBUF*80, CUPBUF*10
      CHARACTER CPROMT(3)*40, CLOWER*26, CUPPER*26
      LOGICAL LEOF
      DATA CPROMT/' ENTER MINUIT TITLE, or "SET INPUT n" : ',
     +            ' ENTER MINUIT PARAMETER DEFINITION:     ',
     +            ' ENTER MINUIT COMMAND:                  '/
C
      DATA CLOWER/'abcdefghijklmnopqrstuvwxyz'/
      DATA CUPPER/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
C
      IFLGUT = 1
      IFLGDO = IFLGIN
      LEOF = .FALSE.
      INCOMP = 0
C                                           . . . . read next record
   10 CONTINUE
      IF (ISW(6) .EQ. 1) THEN
           WRITE (ISYSWR,'(A)') CPROMT(IFLGDO)
           IF (IFLGDO .EQ. 2)  LPHEAD = .FALSE.
      ENDIF
      CRDBUF = '   '
      READ (ISYSRD,'(A)',ERR=500,END=45)  CRDBUF
C
C                 CUPBUF is the first few characters in upper case
      CUPBUF(1:10) = CRDBUF(1:10)
      DO 12 I= 1, 10
      IF (CRDBUF(I:I) .EQ. '''') GO TO 13
         DO 11 IC= 1, 26
         IF (CRDBUF(I:I) .EQ. CLOWER(IC:IC)) CUPBUF(I:I)=CUPPER(IC:IC)
   11    CONTINUE
   12 CONTINUE
   13 CONTINUE
C                                           . .   preemptive commands
      LEOF = .FALSE.
      IF (INDEX(CUPBUF,'*EOF') .EQ. 1)    THEN
         WRITE (ISYSWR,'(A,I3)') ' *EOF ENCOUNTERED ON UNIT NO.',ISYSRD
         LPHEAD = .TRUE.
         GO TO 50
         ENDIF
      IF (INDEX(CUPBUF,'SET INP') .EQ. 1)    THEN
         ICOMND = ICOMND + 1
         WRITE (ISYSWR, 21) ICOMND,CRDBUF(1:50)
   21    FORMAT (' **********'/' **',I5,' **',A/' **********')
         LPHEAD = .TRUE.
         GO TO 50
         ENDIF
      GO TO 80
C                                    . . hardware EOF on current ISYSRD
   45 CRDBUF = '*EOF '
      WRITE (ISYSWR,'(A,I3)') ' END OF DATA ON UNIT NO.',ISYSRD
C                                     or SET INPUT command
   50 CONTINUE
         CALL MNSTIN(CRDBUF,IERR)
         IF (IERR .EQ. 0)  GO TO 10
         IF (IERR .EQ. 2)  THEN
            IF (.NOT. LEOF) THEN
               WRITE (ISYSWR,'(A,A/)') ' TWO CONSECUTIVE EOFs ON ',
     +              'PRIMARY INPUT FILE WILL TERMINATE EXECUTION.'
               LEOF = .TRUE.
               GO TO 10
            ENDIF
         ENDIF
         IFLGUT = IERR
         GO TO 900
   80 IF (IFLGDO .GT. 1) GO TO 100
C                            read title        . . . . .   IFLGDO = 1
C              if title is 'SET TITLE', skip and read again
      IF (INDEX(CUPBUF,'SET TIT') .EQ. 1)  GO TO 10
      CALL MNSETI(CRDBUF(1:50))
      WRITE (ISYSWR,'(1X,A50)')  CTITL
      WRITE (ISYSWR,'(1X,78(1H*))')
         LPHEAD = .TRUE.
      IF (IFLGIN .EQ. IFLGDO)  GO TO 900
      IFLGDO = IFLGIN
      GO TO 10
C                            data record is not a title.
  100 CONTINUE
      IF (IFLGDO .GT. 2)  GO TO 300
C                          expect parameter definitions.   IFLGDO = 2
C              if parameter def is 'PARAMETER', skip and read again
      IF (INDEX(CUPBUF,'PAR') .EQ. 1)  GO TO 10
C              if line starts with SET TITLE, read a title first
      IF (INDEX(CUPBUF,'SET TIT') .EQ. 1)  THEN
         IFLGDO = 1
         GO TO 10
         ENDIF
C                      we really have parameter definitions now
      CALL MNPARS(CRDBUF,ICONDP)
      IF (ICONDP .EQ. 0)  GO TO 10
C          format error
      IF (ICONDP .EQ. 1)  THEN
         IF (ISW(6) .EQ. 1)  THEN
           WRITE (ISYSWR,'(A)') ' FORMAT ERROR.  IGNORED.  ENTER AGAIN.'
           GO TO 10
         ELSE
           WRITE (ISYSWR,'(A)') ' ERROR IN PARAMETER DEFINITION'
           IFLGUT = 4
           GO TO 900
         ENDIF
      ENDIF
C                     ICONDP = 2            . . . end parameter requests
      IF (ISW(5).GE.0 .AND. ISW(6).LT.1) WRITE (ISYSWR,'(4X,75(1H*))')
      LPHEAD = .TRUE.
      IF (IFLGIN .EQ. IFLGDO)  GO TO 900
      IFLGDO = IFLGIN
      GO TO 10
C                                              . . . . .   IFLGDO = 3
C                                           read commands
  300 CONTINUE
      CALL MNCOMD(FCN,CRDBUF,ICONDN,FUTIL)
CC     ICONDN = 0: command executed normally
CC              1: command is blank, ignored
CC              2: command line unreadable, ignored
CC              3: unknown command, ignored
CC              4: abnormal termination (e.g., MIGRAD not converged)
CC              5: command is a request to read PARAMETER definitions
CC              6: 'SET INPUT' command
CC              7: 'SET TITLE' command
CC              8: 'SET COVAR' command
CC              9: reserved
CC             10: END command
CC             11: EXIT or STOP command
CC             12: RETURN command
      IF (ICONDN .EQ. 2 .OR. ICONDN .EQ. 3) THEN
         INCOMP = INCOMP + 1
         IF (INCOMP .GT. 100) THEN
            IFLGUT = 5
            GO TO 900
            ENDIF
         ENDIF
C                         parameter
      IF (ICONDN .EQ. 5)  IFLGDO = 2
C                         SET INPUT
      IF (ICONDN .EQ. 6)  GO TO 50
C                         SET TITLE
      IF (ICONDN .EQ. 7)  IFLGDO = 1
C                                        . . . . . . . . . . set covar
      IF (ICONDN .EQ. 8) THEN
         ICOMND = ICOMND + 1
         WRITE (ISYSWR,405) ICOMND,CRDBUF(1:50)
  405    FORMAT (1H ,10(1H*)/' **',I5,' **',A)
         WRITE (ISYSWR, '(1H ,10(1H*))' )
         NPAR2 = NPAR*(NPAR+1)/2
         READ (ISYSRD,420,ERR=500,END=45)  (VHMAT(I),I=1,NPAR2)
  420    FORMAT (BN,7E11.4,3X)
         ISW(2) = 3
         DCOVAR = 0.0
         IF (ISW(5) .GE. 0)  CALL MNMATU(1)
         IF (ISW(5) .GE. 1)  CALL MNPRIN(2,AMIN)
         GO TO 10
         ENDIF
      IF (ICONDN .LT. 10) GO TO 10
      GO TO 900
C                                              . . . . error conditions
  500 IFLGUT = 3
  900 RETURN
      END
