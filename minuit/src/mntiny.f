*
* $Id: mntiny.f,v 1.1 2010/02/16 20:35:21 wfisher Exp $
*
* $Log: mntiny.f,v $
* Revision 1.1  2010/02/16 20:35:21  wfisher
*
* First commit for V00-04-00
*
* Revision 1.1.1.1  1996/03/07 14:31:32  mclareni
* Minuit
*
*
      SUBROUTINE MNTINY(EPSP1,EPSBAK)
*
* $Id: mntiny.f,v 1.1 2010/02/16 20:35:21 wfisher Exp $
*
* $Log: mntiny.f,v $
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
CC        Compares its argument with the value 1.0, and returns
CC        the value .TRUE. if they are equal.  To find EPSMAC
CC        safely by foiling the Fortran optimizer
CC
      PARAMETER (ONE=1.0)
      EPSBAK =  EPSP1  - ONE
      RETURN
      END
