C ***********************************************************************************************************
C GCG3AL
C ALLOCATE STORAGE IN THE X AND IX ARRAYS FOR GCG ARRAYS
C ***********************************************************************************************************
      subroutine gcg_alloc

      use rt_global
      use mf_rt_link
      use GLOBAL,     only:NCOL,NROW,NLAY

      IMPLICIT NONE
      INTEGER  I,J,K

C     PRINT A MESSAGE IDENTIFYING GCG PACKAGE
      WRITE(IOUT,1) INGCG
    1 FORMAT(1X,'GCG3 -- GENERALIZED CONJUGATE GRADIENT SOLVER PACKAGE',
     & ', VER DoD_3.0, JUNE 1998',' INPUT READ FROM UNIT',I3)
C
C     READ AND PRINT RT_MXITER AND ISOLVE
      READ(INGCG,*) RT_ITER1,RT_MXITER,ISOLVE,NCRS
      WRITE(IOUT,3) RT_MXITER,RT_ITER1
    3 FORMAT(1X,'MAXIMUM OF',I5,' OUTER ITERATIONS',
     &      /1X,'       AND',I5,' INNER ITERATIONS ALLOWED FOR CLOSURE')
      IF(RT_MXITER.LE.0) THEN
        WRITE(*,5)
        STOP
      ELSEIF(RT_ITER1.LE.0) THEN
        WRITE(*,7)
        STOP
      ENDIF
      IF(ISOLVE.EQ.1) THEN
        WRITE(IOUT,13)
      ELSEIF(ISOLVE.EQ.2) THEN
        WRITE(IOUT,23)
      ELSEIF(ISOLVE.EQ.3) THEN
        WRITE(IOUT,33)
      ELSE
        WRITE(IOUT,43)
        STOP
      ENDIF
    5 FORMAT(/1X,'ERROR: OUTER ITERATION NUMBER MUST BE > 0.')
    7 FORMAT(/1X,'ERROR: INNER ITERATION NUMBER MUST BE > 0.')
   13 FORMAT(1X,'THE PRECONDITIONING TYPE SELECTED IS JACOBI.')
   23 FORMAT(1X,'THE PRECONDITIONING TYPE SELECTED IS SSOR.')
   33 FORMAT(1X,'THE PRECONDITIONING TYPE SELECTED IS ',
     &          'MODIFIED INCOMPLETE CHOLESKY (MIC).')
   43 FORMAT(1X,'ERROR: INVALID PRECONDITIONING TYPE.')
C
      IF(NCRS.GT.0) THEN
        WRITE(IOUT,50)
      ELSE
        WRITE(IOUT,52)
      ENDIF
   50 FORMAT(1X,'FULL DISPERSION TENSOR INCLUDED IN IMPLICIT SOLUTION')
   52 FORMAT(1X,'DISPERSION CROSS TERMS LUMPED INTO RIGHT-HAND-SIDE')

C     SET NCRS TO 0 FOR 1D PROBLEMS
      IF(NCOL*NROW.EQ.1 .OR. NCOL*NLAY.EQ.1 .OR. NROW*NLAY.EQ.1) NCRS=0

C     The number of cells in the 3D grid
      NODES=NLAY*NROW*NCOL

C     Initialize the array pointers for allocating memory for the integer
C     and real arrays
      ISUM2 = 1
      ISUM = 1

C     INTEGER ARRAYS (IX array)
      LCLRCH = ISUM2
      ISUM2 = ISUM2+ (3*RT_MXITER*RT_ITER1) !Size of LRCH array
      LCIBND = ISUM2
      ISUM2 = ISUM2 + NODES

C     REAL ARRAYS (X array)
C     A array      
      LCA = ISUM
      IF(NCRS.GT.0) THEN
        ISUM = ISUM + (NODES*19)
      ELSE
        ISUM = ISUM + (NODES*7)
      ENDIF
C     Q array
      LCQ = ISUM
      IF(ISOLVE.EQ.3) THEN
        IF(NCRS.GT.0) THEN
          ISUM = ISUM + (NODES*19)
        ELSE
          ISUM = ISUM + (NODES*7)
        ENDIF
      ENDIF
C     WK array
      LCWK = ISUM
      ISUM = ISUM + (NODES*7)
C     CNCG array
      LCCNCG = ISUM
      ISUM = ISUM + (RT_MXITER*RT_ITER1)
C     RHS array
      LCRHS = ISUM
      ISUM = ISUM + NODES
C     CNEW array
      LCCNEWN = ISUM
      ISUM = ISUM + NODES

C     CHECK HOW MANY ELEMENTS OF THE X AND IX ARRAYS ARE USED FOR THE GCG ARRAY
      WRITE(IOUT,1090) ISUM,ISUM2
 1090 FORMAT(1X,I10,' ELEMENTS OF THE  X ARRAY USED BY THE GCG PACKAGE'
     & /1X,I10,' ELEMENTS OF THE IX ARRAY USED BY THE GCG PACKAGE'/)

*     Set values for implicit solver
      IMPSOL = 0
      if(TRNOP(5)) then
        IMPSOL = 1
        ISPD = 1
        if(MIXELM.EQ.0) ISPD = 0
      endif
      allocate(X(ISUM))
      allocate(IX(ISUM2))
      do I=1,ISUM
        X(I) = 0.
      enddo
      do I=1,ISUM2
        IX(I) = 0
      enddo

C     READ RT_ACCL,CCLOSE,IPRGCG
      READ(INGCG,*) RT_ACCL,CCLOSE,IPRGCG
      IF(RT_ACCL.EQ.0.) RT_ACCL = 1.

C     PRINT DATA VALUES JUST READ
      WRITE(IOUT,100)
  100 FORMAT(///47X,'SOLUTION BY THE GENERALIZED CONJUGATE GRADIENT',
     & ' ISOLVE'/47X,53('-'))
      WRITE(IOUT,115) RT_MXITER
  115 FORMAT(37X,'MAXIMUM OUTER ITERATIONS ALLOWED FOR CLOSURE =',I9)
      WRITE(IOUT,116) RT_ITER1
  116 FORMAT(37X,'MAXIMUM INNER ITERATIONS ALLOWED FOR CLOSURE =',I9)
      WRITE(IOUT,117) ISOLVE
  117 FORMAT(52X,'PRECONDITIONING TYPE SELECTED =',I5)
      WRITE(IOUT,120) RT_ACCL
  120 FORMAT(59X,'ACCELERATION PARAMETER =',G15.5)
      WRITE(IOUT,125) CCLOSE
  125 FORMAT(39X,'CONCENTRATION CHANGE CRITERION FOR CLOSURE =',E15.5)
      IF(IPRGCG.LE.0) IPRGCG=999
      WRITE(IOUT,130) IPRGCG
  130 FORMAT(39X,'GCG CONCENTRATION CHANGE PRINTOUT INTERVAL =',I9)

      RETURN
      END


C************************************************************************************************************
C SGCG3P
C PRINT MAXIMUM CONCENTRATION CHANGES FOR EACH ITERATION DURING A TRANSPORT TIME STEP
C************************************************************************************************************
      SUBROUTINE SGCG3P(CNCG,LRCH)

      use rt_global
      use mf_rt_link
      use GLOBAL,     only:NCOL,NROW,NLAY

      IMPLICIT  NONE
      INTEGER   LRCH,I,J
      REAL      CNCG
      DIMENSION CNCG(RT_MXITER*RT_ITER1),LRCH(3,RT_MXITER*RT_ITER1)
C
      WRITE(IOUT,5)
    5 FORMAT(1X,' MAXIMUM CONCENTRATION CHANGES FOR EACH ITERATION:'
     &      /1X, 5(' MAX. CHANGE LAYER,ROW,COL')/1X,132('-'))
      WRITE(IOUT,10) (CNCG(J),(LRCH(I,J),I=1,3),J=1,ITP)
   10 FORMAT((1X,5(G12.4,' (',I3,',',I3,',',I3,')')))
C
      RETURN
      END


C ***********************************************************************************************************
C GCG3AP
C SOLUTION BY THE GENERALIZED CONJUGATE GRADIENT ISOLVES, USING ONE OF THE THREE PRECONDITIONERS (JACOBI, 
C SSOR, AND MIC) WITH LANCZOS/ORTHOMIN ACCELERATION, UP TO ITER1 ITERATIONS. HOWEVER, IF IT IS A FLOW PROBLEM 
C OR IT IS A SYMMETRIC CASE, (I.E., ISPD>0), THE ORDINARY CG IS USED.
C ***********************************************************************************************************
C
C                           PARAMETER LIST
C 
C  IOUT    : INTEGER, OUTPUT UNIT NUMBER.(INPUT)
C  RT_MXITER  : INTEGER, OUTER LOOP MAX NUMBER.(INPUT)
C  RT_ITER1   : INTEGER, MAXIMUM NUMBER OF ITERATIONS ALLOWED.(INPUT)
C  ITO     : INTEGER, OUTER LOOP COUNTER.(INPUT)
C  ITP     : INTEGER, ACTUAL NUMBER OF ITERATION REQUIRED.(OUTPUT)
C  ISOLVE  : INTEGER, BASIC PRECONDITIONER SELECTION.(INPUT)
C                 1=JACOBI, 2=SSOR, 3=MIC
C  ACCL   : REAL, RELAXATION FACTOR FOR THE SSOR ISOLVE.(INPUT)
C  CCLOSE  : REAL, STOPPING CRITERION.(INPUT)
C  ICNVG   : INTEGER, GLOBAL CONVERGENCE FLAG.(OUTPUT)
C                 1=GLOBAL CONVERGENCE, 0=NOT CONVERGED
C  LICNVG  : INTEGER, LOCAL CONVERGENCE FLAG.(OUTPUT)
C                 1=LOCAL CONVERGENCE, 0=NO CONVERGED
C  LITP    : INTEGER, LOCAL ITERATION COUNTER.
C  CNCG    : REAL ARRAY, CONTAINS MAXIMUM CHANGE OF CONCENTRATION
C            AT EVERY ITERATION.(OUTPUT)
C  LRCH    : INTEGER 2-D ARRAY, CONTAINS THE LOCATIONS OF THE
C            MAXIMUM CHANGE OF CONCENTRATION AT EVERY ITERATION.
C  NCOL,NROW,NLAY:
C            INTEGERS, NUMBER OF COLUMNS, ROWS, AND LAYERS.(INPUT)
C  NODES   : INTEGER, DIMENSION OF THE MATRIX.(INPUT)
C  NTRANS  : INTEGER, TRANSPORT STEP INDEX.(INPUT)
C  KSTP    : INTEGER, FLOW TIME STEP INDEX.(INPUT)
C  KPER    : INTEGER, STRESS PERIOD INDEX.(INPUT)
C  TIME2   : REAL, TOTAL ELAPSED TIME.(INPUT)
C  T2     : REAL, ACCUMULATED TIME IN CURRENT STRESS PERIOD.(INPUT)
C  UPDLHS  : LOGICAL, FLAG FOR UPDATING COEFF. MATRIX.(INPUT)
C  IPRGCG  : INTEGER, INTERVAL FOR PRINTING MAX. CHANGES.(INPUT)
C  ICBUND  : INTEGER, BOUNDART TYPE INDICATOR.(INPUT)
C  A       : COEFF. MATRIX.(INPUT)
C  CNEW    : REAL ARRAY, CONTAINS CONCENTRATION.(OUTPUT)
C  RHS     : REAL ARRAY, CONTAINS RIGHT HAND SIDE OF THE SYSTEM.(INPUT)
C  Q,DQ    : BASIC ITERATIVE ISOLVE PRECONDITIONING MATRIX STORAGE.
C  WK      : REAL ARRAY OF LENGTH 7*NCOL*NROW*NLAY, WORK SPACES.
C  NCRS    : INTEGER, 7 OR 19 DIAGONALS INDICATOR.
C  ISPD    : INPUT INTEGER, SYMMETRIC CASE INDICATOR.
C ***********************************************************************************************************
      subroutine gcg_solve(CNCG,LRCH,A,CNEWN,RHS,ICBUNDN,Q,WK)

      use rt_global
      use mf_rt_link
      use GLOBAL,     only:NCOL,NROW,NLAY

      IMPLICIT  NONE
      INTEGER   I,J,K,N,IDELTA,IRDEL,IPN,LRCH,
     &          IPA,IAPN,IATPN,IWK,II,LITP,LICNVG,IPLACE,
     &          ICBUNDN,KLAYER,IROW,JCOLMN,IJ
      REAL      CNEWN,A,RHS,Q,WK,CNCG,
     &          DELCC,DCOLD,ALN,CHANGE,CTEMP,CHTMP,
     &          RLN,DELPN,PNAPN,SCALE,
     &          DELAPN,APNPA,GSTOP,RHSNORM,TINY
      DIMENSION A(19*NODES),CNEWN(NODES),
     &          RHS(NODES),Q(19*NODES),WK(7*NODES),
     &          CNCG(RT_MXITER*RT_ITER1),LRCH(3,RT_MXITER*RT_ITER1),
     &          ICBUNDN(NODES)
      PARAMETER (TINY=1.E-30)

C
      IDELTA = 1
      IPN    = IDELTA+NODES
      IAPN   = IPN   +NODES
      IWK    = IAPN  +NODES
C
C--THE FOLLOWING WORK SPACE ALLOCATIONS ARE NEEDED FOR
C--NONSYMMETRIC CASES
      IF (ISPD.EQ.0) THEN
      IRDEL  = IWK   +NODES
      IPA    = IRDEL +NODES
      IATPN  = IPA   +NODES
      ENDIF
C
      IF(ITO.EQ.1) ITP=0
      LITP = 0
      ICNVG=0
      LICNVG=0
C
C--NORMALIZE CNEW AND RHS TO AVOID TOO BIG OR TOO SMALL QUANTITIES
C--NORMALIZATION FACTOR IS MAX (CNEW)
      SCALE = 0.
      DO I = 1,NODES
         IF (ICBUNDN(I).NE.0) THEN
           IF (ABS(CNEWN(I)).GT.SCALE) SCALE = CNEWN(I)
         ENDIF
      ENDDO
      IF (SCALE .GT. TINY) THEN
         DO I = 1,NODES
           IF (ICBUNDN(I).NE.0) THEN
             RHS(I) = RHS(I)/SCALE
             CNEWN(I) = CNEWN(I)/SCALE
           ENDIF
         ENDDO
      ENDIF
C
C--COMPUTE RESIDUAL VECTOR R=RHS-A*CNEW AND TEST FOR SOLUTION
      CALL MVPRD(ICBUNDN,A,CNEWN,WK(IWK))
      GSTOP=0.
      RHSNORM = 0.
      DO I = 1,NODES
         IF (ICBUNDN(I).NE.0) THEN
            WK(IWK-1+I) = RHS(I)-WK(IWK-1+I)
            RHSNORM = RHSNORM+RHS(I)*RHS(I)
            GSTOP =GSTOP+WK(IWK-1+I)*WK(IWK-1+I)
         ENDIF
      ENDDO
      RHSNORM = SQRT(RHSNORM)
      GSTOP = SQRT(GSTOP)
      IF (RHSNORM.NE.0) GSTOP = GSTOP / RHSNORM
      IF(GSTOP .LE. 1.0E-6) THEN
          CHANGE = 0.
          IPLACE = 1
          ITP = ITP + 1
          LITP= LITP+ 1
          IF(UPDLHS.AND.ISOLVE.EQ.3) CALL MIC(A,Q)
          GO TO 300
      ENDIF
C
C--COMPUTE PSEUDO-RESIDUAL DELTA = Q^-1*(RHS-A*CNEW)
      IF (ISOLVE.EQ.3) THEN
        IF(UPDLHS) CALL MIC(A,Q)
      ENDIF
C...... JACOBI OR SSOR ISOLVE .....
      IF(ISOLVE.EQ.1.OR.ISOLVE.EQ.2) THEN
        CALL QSOLVE(A,WK(IWK),WK(IDELTA))
      ELSE
C...... MIC ISOLVE .....
        CALL QSOLVE(Q,WK(IWK),WK(IDELTA))
      ENDIF
C
C--BRANCH FOR THE NONSYMMETRIC CASE
C--LANCZOS/ORTHOMIN ACCELERATION IS USED
      IF (ISPD.NE.0) GO TO 200
      DO II = 1,NODES
         WK(IPN-1+II) = WK(IDELTA-1+II)
         WK(IRDEL-1+II) = WK(IDELTA-1+II)
         WK(IPA-1+II) = WK(IDELTA-1+II)
      ENDDO
      DELCC = 0.
      DO II = 1,NODES
         IF (ICBUNDN(II).NE.0)
     &       DELCC = DELCC + WK(IDELTA-1+II)*WK(IRDEL-1+II)
      ENDDO
      GO TO 100
C
C--COMPUTE DIRECTION VECTORS
 90   CONTINUE
      DCOLD = DELCC
      DELCC = 0.
      DO II = 1, NODES
         IF (ICBUNDN(II).NE.0)
     &       DELCC = DELCC + WK(IDELTA-1+II)*WK(IRDEL-1+II)
      ENDDO
      ALN = DELCC / DCOLD
      DO I = 1,NODES
         WK(IPN-1+I) = WK(IDELTA-1+I) + ALN * WK(IPN-1+I)
         WK(IPA-1+I) = WK(IRDEL-1+I) + ALN * WK(IPA-1+I)
      ENDDO
C
C--COMPUTE NEW ITERATES
 100  CONTINUE
      ITP = ITP + 1
      LITP = LITP + 1
      CALL MVPRD(ICBUNDN,A,WK(IPN),WK(IWK))
      IF(ISOLVE.EQ.1.OR.ISOLVE.EQ.2) THEN
        CALL QSOLVE(A,WK(IWK),WK(IAPN))
        CALL QTSLVE(A,WK(IPA),WK(IWK)) 
      ELSE
        CALL QSOLVE(Q,WK(IWK),WK(IAPN))
        CALL QTSLVE(Q,WK(IPA),WK(IWK)) 
      ENDIF
C
      CALL MTVPRD(ICBUNDN,A,WK(IWK),WK(IATPN))
      APNPA = 0.
      DO II = 1, NODES
         IF(ICBUNDN(II).NE.0)
     &      APNPA = APNPA + WK(IAPN-1+II)*WK(IPA-1+II)
      ENDDO
      IF(APNPA.NE.0) RLN = DELCC / APNPA
      CHANGE = 0.0
      IPLACE = 1
      DO I = 1,NODES
         CTEMP = CNEWN(I)
         CNEWN(I) = CNEWN(I) + RLN * WK(IPN-1+I)
         WK(IDELTA-1+I) = WK(IDELTA-1+I) - RLN * WK(IAPN-1+I)
         WK(IRDEL-1+I) = WK(IRDEL-1+I) - RLN * WK(IATPN-1+I)
         CHTMP = ABS(CNEWN(I)-CTEMP)
         IF (CHTMP.GT.CHANGE) THEN
            CHANGE = CHTMP
            IPLACE = I
         ENDIF
      ENDDO
      GO TO 300
 200  CONTINUE
C
C--THE FOLLOWING IS FOR THE SYMMETRIC CASE
C--ORDINARY CG ACCELERATION IS USED
      DO II = 1, NODES
         WK(IPN-1+II) = WK(IDELTA-1+II)
      ENDDO
      GO TO 220
C
C--COMPUTE DIRECTION VECTORS
 210  CONTINUE
      DELAPN = 0.
      DO II = 1, NODES
         IF(ICBUNDN(II).NE.0)
     &      DELAPN = DELAPN + WK(IDELTA-1+II)*WK(IAPN-1+II)
      ENDDO
      ALN = -DELAPN / PNAPN
      DO I = 1,NODES
         WK(IPN-1+I) = WK(IDELTA-1+I) + ALN * WK(IPN-1+I)
      ENDDO
C
C--COMPUTE NEW ITERATES
 220  CONTINUE
      ITP = ITP + 1
      LITP = LITP + 1
      DELPN = 0.
      DO II = 1, NODES
         IF(ICBUNDN(II).NE.0)
     &      DELPN = DELPN + WK(IDELTA-1+II)*WK(IPN-1+II)
      ENDDO
      CALL MVPRD(ICBUNDN,A,WK(IPN),WK(IWK))
      IF(ISOLVE.EQ.1.OR.ISOLVE.EQ.2) THEN
        CALL QSOLVE(A,WK(IWK),WK(IAPN))
      ELSE
        CALL QSOLVE(Q,WK(IWK),WK(IAPN))
      ENDIF
      PNAPN = 0.
      DO II = 1, NODES
         IF(ICBUNDN(II).NE.0)
     &      PNAPN = PNAPN + WK(IAPN-1+II)*WK(IPN-1+II)
      ENDDO
      IF(PNAPN.NE.0) RLN = DELPN / PNAPN
      CHANGE = 0.0
      IPLACE = 1
      DO I = 1,NODES
         CTEMP = CNEWN(I)
         CNEWN(I) = CNEWN(I) + RLN * WK(IPN-1+I)
         WK(IDELTA-1+I) = WK(IDELTA-1+I) - RLN * WK(IAPN-1+I)
         CHTMP = ABS(CNEWN(I)-CTEMP)
         IF (CHTMP.GT.CHANGE) THEN
            CHANGE = CHTMP
            IPLACE = I
         ENDIF
      ENDDO
C
 300  CONTINUE
C
C--STORE MAXIMUM CHANGE VALUE AND LOCATION
      CNCG(ITP) = CHANGE
      KLAYER = (IPLACE-1) / (NCOL*NROW) + 1
      IJ = IPLACE - (KLAYER-1)*NCOL*NROW
      IROW = (IJ-1)/NCOL + 1
      JCOLMN = IJ - (IROW-1)*NCOL
        LRCH(1,ITP) = KLAYER
        LRCH(2,ITP) = IROW
        LRCH(3,ITP) = JCOLMN
C
      !WRITE(*,1111) ITO,LITP,CHANGE,KLAYER,IROW,JCOLMN
 1111 FORMAT(1X,'Outer Iter.',I3,'  Inner Iter.',I3,
     &  ':  Max. DC =',G12.4,'  [K,I,J]',3I5)
C
C--CHECK CONVERGENCE ......
      IF(CHANGE.LE.CCLOSE) THEN
         LICNVG=1
      ENDIF
      IF(RT_MXITER.EQ.1) THEN
        IF(LICNVG.EQ.1) ICNVG=1
      ELSEIF(ITO.GT.1) THEN
        IF(LICNVG.EQ.1.AND.LITP.EQ.1) ICNVG=1
      ENDIF
C
C--LOCAL CONVERGENCE NOT MET, LOOP BACK
      IF(LICNVG.EQ.0 .AND. LITP.LT.RT_ITER1) THEN
        IF (ISPD.EQ.0) THEN
           GOTO 90
        ELSE
           GOTO 210
        ENDIF
      ENDIF
      IF(ICNVG.EQ.0 .AND. ITO.NE.RT_MXITER) GOTO 600
      IF(NTRANS.EQ.1) WRITE(IOUT,1000)
 1000 FORMAT(/1X)
      WRITE(IOUT,1010) ITO,NTRANS,rt_kstp,rt_kper,ITP
 1010 FORMAT(1X,I5,' CALLS TO GCG PACKAGE FOR TRANSPORT TIME STEP',I4,
     & ' IN FLOW TIME STEP',I4,' STRESS PERIOD',I4,
     & /1X,I5,' TOTAL ITERATIONS')
C
      IF(ICNVG.EQ.0 .OR. TIME2.GE.T2 .OR. MOD(NTRANS,IPRGCG).EQ.0)
     & CALL SGCG3P(CNCG,LRCH)
C
  600 CONTINUE
C
C--BEFORE RETURN UNSCALE CNEW AND RHS
         IF (SCALE .GT. TINY) THEN
            DO I = 1,NODES
              IF (ICBUNDN(I).NE.0) THEN
                RHS(I) = RHS(I)*SCALE
                CNEWN(I) = CNEWN(I)*SCALE
              ELSE
                CNEWN(I)=CINACT
              ENDIF
            ENDDO
         ENDIF
      
C     Put CNEW back into the J,I,K indices array (for use in the RCT3SV subroutine)
      DO K=1,NLAY
        DO I=1,NROW
          DO J=1,NCOL
            N=(K-1)*NCOL*NROW+(I-1)*NCOL+J
            CNEW(J,I,K,ICOMP) = CNEWN(N)
          ENDDO
        ENDDO
      ENDDO      
      
      RETURN
      END


C************************************************************************************************************
C MVPRD
C THIS SUBROUTINE, MVPRD, PERFORMS AX=Y WHERE THE MATRIX A IS STORED IN DIAGONAL FORM
C************************************************************************************************************
      SUBROUTINE MVPRD(ICBUNDN,A,XX,Y)

      use rt_global
      use GLOBAL,     only:NCOL,NROW,NLAY

      IMPLICIT  NONE
      INTEGER   L,IDIAG,I,JCOL,K,ICBUNDN
      REAL      A,XX,Y
      DIMENSION A(NODES,*),XX(NODES),Y(NODES),ICBUNDN(NODES)
      COMMON   /GCGIDX/L(19)
C
      IDIAG = 7
      IF(NCRS.GT.0) IDIAG = 19
      DO I=1,NODES
         Y(I) = 0.
         DO K = 1,IDIAG
            JCOL = I + L(K)
            IF (JCOL.GE.1.AND.JCOL.LE.NODES) THEN
             IF(ICBUNDN(JCOL).NE.0) Y(I) = Y(I)+A(I,K)*XX(JCOL)
            ENDIF
         ENDDO
      ENDDO
C
      RETURN
      END


C************************************************************************************************************
C MTVPRD
C THIS SUBROUTINE, MTVPRD, PERFORMS A^TX=Y
C************************************************************************************************************
      subroutine MTVPRD(ICBUNDN,A,XX,Y)

      use rt_global
      use GLOBAL,     only:NCOL,NROW,NLAY

      IMPLICIT  NONE
      INTEGER   L,IDIAG,I,JCOL,ICBUNDN,J
      REAL      A,XX,Y
      DIMENSION A(NODES,*),XX(NODES),Y(NODES),ICBUNDN(NODES)
      COMMON   /GCGIDX/L(19)
C
      IDIAG = 7
      IF(NCRS.GT.0) IDIAG = 19
      DO I=1,NODES
         Y(I) = 0.
      ENDDO
      DO I=1,NODES
         DO J = 1,IDIAG
            JCOL = I + L(J)
            IF(JCOL.GE.1.AND.JCOL.LE.NODES) THEN
             IF (ICBUNDN(JCOL).NE.0) Y(JCOL) = Y(JCOL)+A(I,J)*XX(I)
            ENDIF
         ENDDO
      ENDDO
C
      RETURN
      END


C************************************************************************************************************
C QSOLVE
C*... FUNCTION:
C*
C*         THIS SUBROUTINE, QSOLVE, PERFORMS THE FORWARD AND BACKWARD
C*         SUBSTITUTION.  WHICH SOLVES   Q * Y = SY  FOR THE JACOBI
C*         SSOR, AND MIC ISOLVES.
C*
C*... PARAMETER LIST
C*
C*         N       : INPUT INTEGER, THE DIMENSION OF THE SYSTEM
C*         ISOLVE  : INPUT INTEGER, DEFINES BASIC ITERATIVE ISOLVE
C*         ACCL   : INPUT REAL, RELAXATION FACTOR FOR SSOR
C*         NCRS    : INPUT INTEGER, 7 OR 19 DIAGONALS INDICATOR.
C*         A       :  INPUT REAL ARRAY. CONTAINS THE NONZERO ELEMENTS
C*                    OF THE COEFFICIENT MATRIX A.
C*         SY      :  INPUT REAL ARRAY. THE RIGHT HAND SIDE OF THE
C*                    SYSTEM.
C*         Y       :  OUTPUT REAL ARRAY. CONTAINS THE SOLUTION.
C************************************************************************************************************
      subroutine QSOLVE(A,SY,Y)

      use rt_global
      use GLOBAL,     only:NCOL,NROW,NLAY

      IMPLICIT  NONE
      INTEGER   LL,LU,I,J,K,II,IDIAG,JCOL,L
      REAL      Y,SY,A
      DIMENSION Y(NODES),SY(NODES),A(NODES,*),LL(9),LU(9)
      COMMON   /GCGIDX/ L(19)
C
      DO II = 1,NODES
         Y(II) = SY(II)
      ENDDO
C
      IF (ISOLVE .NE. 1) GO TO 5
C
C ... JACOBI ISOLVE
C ... NOTE: THE FIRST ELEMENT IN EACH ROW OF THE MATRIX A IS THE 
C ....      DIAGONAL ELEMENT OF THE ROW.
C
      DO I=1,NODES
        Y(I)=Y(I)/A(I,1)
      ENDDO
      RETURN
C
C     SOLVE LDUY=Y
C
 5    CONTINUE
      IF (ISOLVE .NE. 2) RT_ACCL = 1.0
      IF (NCRS.GT.0) THEN
         IDIAG = 9
      ELSE
         IDIAG = 3
      ENDIF
C
      LL(1) = 2
      LL(2) = 4
      LL(3) = 6
      LU(1) = 3
      LU(2) = 5
      LU(3) = 7
      IF (NCRS.GT.0) THEN
         LL(4) = 8
         LL(5) = 9
         LL(6) = 10
         LL(7) = 11
         LL(8) = 16
         LL(9) = 17
         LU(4) = 12
         LU(5) = 13
         LU(6) = 14
         LU(7) = 15
         LU(8) = 18
         LU(9) = 19
      ENDIF
C
C ... SOLVE LOWER TRIANGULAR SYSTEM FOR THE SSOR ISOLVE
C
      DO I=1,NODES
         DO J = 1,IDIAG
            JCOL = I + L(LL(J)) 
            IF (JCOL.GT.0) Y(I) = Y(I) - A(I,LL(J))*Y(JCOL)
         ENDDO
         Y(I) = RT_ACCL*Y(I) / A(I,1)
      ENDDO
C
C     SOLVE DY=Y
C
      IF (ISOLVE .EQ. 2) THEN
         DO I = 1,NODES
            Y(I) = (2.0-RT_ACCL)/RT_ACCL * Y(I) * A(I,1)
         ENDDO
      ELSE
         DO I=1,NODES
            Y(I)=Y(I)*A(I,1)
         ENDDO
      ENDIF
C
C     SOLVE UY=Y
C
      DO I=NODES,1,-1
         DO K=1,IDIAG
            JCOL = I + L(LU(K))
            IF (JCOL.LE.NODES) Y(I)=Y(I)-A(I,LU(K))*Y(JCOL)
         ENDDO
         Y(I)=RT_ACCL*Y(I)/A(I,1)
      ENDDO
C
      RETURN
      END





C************************************************************************************************************
C QTSLVE
C*... FUNCTION:
C*
C*         THIS SUBROUTINE, QTSLVE, PERFORMS THE DRIVER OF SOLVING
C*         THE SYSTEM  Q ** (T) * Y = SY  WHERE THE MATRIX Q IS IN THE
C*         L * D * U FORM.
C*
C*... PARAMETER LIST:
C*
C*         N          : INPUT INTEGER, THE DIMENSION OF THE SYSTEM
C*         ISOLVE     : INPUT INTEGER, DEFINES BASIC ITERATIVE ISOLVE
C*         ACCL      : INPUT REAL, RELAXATION FACTOR FOR SSOR
C*         NCRS       : INPUT INTEGER, 7 OR 19 DIAGONALS INDICATOR
C*         A          : INPUT REAL ARRAY, THE COEFFICIENT MATRIX
C*         SY         : INPUT REAL ARRAY. IT CONTAINS THE RIGHT HAND
C*                      SIDE OF THIS SYSTEM.
C*         Y          : OUTPUT REAL ARRAY. IT CONTAINS THE SOLUTION OF
C*                      THIS SYSTEM.
C************************************************************************************************************
      subroutine QTSLVE(A,SY,Y)

      use rt_global
      use GLOBAL,     only:NCOL,NROW,NLAY

      IMPLICIT  NONE
      INTEGER   LL,LU,I,J,K,II,IDIAG,JCOL,L
      REAL      Y,SY,A
      DIMENSION Y(NODES),SY(NODES),A(NODES,*),LL(9),LU(9)
      COMMON   /GCGIDX/ L(19)
C
      DO II = 1,NODES
        Y(II) = SY(II)
      ENDDO
C
      IF (ISOLVE .NE. 1) GOTO 5
C
C ... JACOBI ISOLVE 
C
      DO I = 1,NODES
        Y(I) = Y(I) / A(I,1)
      ENDDO
      RETURN
C
C ... SOLVE LDU ** (T) * Y = Y
C
 5    CONTINUE
      IF (ISOLVE .NE. 2) RT_ACCL = 1.0
      IF (NCRS.GT.0) THEN
         IDIAG = 9
      ELSE
         IDIAG = 3
      ENDIF
C
      LL(1) = 2
      LL(2) = 4
      LL(3) = 6
      LU(1) = 3
      LU(2) = 5
      LU(3) = 7
      IF (NCRS.GT.0) THEN
         LL(4) = 8
         LL(5) = 9
         LL(6) = 10
         LL(7) = 11
         LL(8) = 16
         LL(9) = 17
         LU(4) = 12
         LU(5) = 13
         LU(6) = 14
         LU(7) = 15
         LU(8) = 18
         LU(9) = 19
      ENDIF
C
C ... SOLVE (UT)Y = Y
C
      DO I=1,NODES
         Y(I) = RT_ACCL * Y(I) / A(I,1)
         DO J = 1,IDIAG
            JCOL = I + L(LU(J))
            IF (JCOL.LE.NODES) Y(JCOL) = Y(JCOL) - A(I,LU(J))*Y(I)
         ENDDO
      ENDDO
C
C     SOLVE DY=Y
C
      IF (ISOLVE .EQ. 2) THEN
         DO I = 1,NODES
            Y(I) = (2.0-RT_ACCL)/RT_ACCL * Y(I) * A(I,1)
         ENDDO
      ELSE
         DO I=1,NODES
            Y(I)=Y(I)*A(I,1)
         ENDDO
      ENDIF
C
C     SOLVE (LT)Y=Y
C
      DO I=NODES,1,-1
         Y(I) = RT_ACCL * Y(I) / A(I,1)
         DO K = 1,IDIAG
            JCOL = I + L(LL(K))
            IF (JCOL.GT.0) Y(JCOL) = Y(JCOL) - A(I,LL(K))*Y(I) 
         ENDDO
      ENDDO
C
      RETURN
      END





C************************************************************************************************************
C MIC
C*... FUNCTION:
C*
C*         THIS SUBROUTINE, MIC, PERFORMS
C*         THE MODIFIED INCOMPLETE CHOLESKY FACTORIZATION. 
C*
C*... PARAMETER USED IN THIS SUBROUTINE:
C*
C*         NODESN  INPUT INTEGER. DIMENSION OF THE MATRIX A. (= N)
C*         NCRS    INPUT INTEGER, 7 OR 19 DIAGONALS INDICATOR
C*         A       INPUT REAL VECTOR. CONTAINS THE NONZERO ELEMENTS
C*                 OF THE MATRIX.
C*         Q       OUTPUT, CONTAINS THE COMPACT LDU FORM FOR THE MIC
C************************************************************************************************************
      subroutine MIC(A,Q)

      use rt_global
      use GLOBAL,     only:NCOL,NROW,NLAY

      IMPLICIT  NONE
      INTEGER   K,II,IJ,IPVT,IICOL,ILST,IROW,JP,JCOL,KK,LL,IERR,
     &          LU,IDIAG,L,J
      REAL      A,Q,QMULT,TINY
      PARAMETER (TINY=1.E-30)
      DIMENSION A(NODES,*),Q(NODES,*),LU(9)
      COMMON   /GCGIDX/ L(19)
C
C ... SET INITIAL PARAMETERS:
C
      LU(1) = 3
      LU(2) = 5
      LU(3) = 7
      IF (NCRS.GT.0) THEN
         LU(4) = 12
         LU(5) = 13
         LU(6) = 14
         LU(7) = 15
         LU(8) = 18
         LU(9) = 19
      ENDIF
      IF (NCRS.GT.0) THEN
          IDIAG = 19
          ILST  = 9
      ELSE
          IDIAG = 7
          ILST  = 3
      ENDIF
      DO K = 1,NODES
         DO J = 1, IDIAG
            Q(K,J) = A(K,J)
         ENDDO
      ENDDO
      IF(ABS(Q(1,1)).LT.TINY) Q(1,1)=1.
C
      DO IPVT = 1,NODES-1
        DO KK = 1,ILST
          IROW = IPVT+L(LU(KK))
          IF(IROW.GT.NODES) CYCLE
          DO II = 1,IDIAG
            IICOL = IROW + L(II)
            IF (IICOL .EQ. IPVT) GO TO 95
          ENDDO
          CYCLE
 95       QMULT = Q(IROW,II)
          DO LL = 1,ILST
            JCOL = IPVT + L(LU(LL))
            IERR = 1
            DO IJ = 1,IDIAG
               JP = IROW + L(IJ)
               IF(JP .EQ. JCOL) THEN
                 IERR = 0
                 EXIT
               ENDIF
            ENDDO
            IF(IERR.EQ.0) THEN
              Q(IROW,IJ)=Q(IROW,IJ)-QMULT*Q(IPVT,LU(LL))/Q(IPVT,1)
            ELSEIF(ABS(Q(IROW,1)).LT.TINY) THEN
              Q(IROW,1)=1.
            ENDIF
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END
