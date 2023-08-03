C ************************************************************************************************************
C ADV (Advection) Package
C Handles subroutines related to Advection
C
C Included Subroutines:
C ADV3AL: Allocate computer memory
C ADV3RP: Read and Prepare
C ADV3SV: Calculate change in cell concentration due to advection
C ADV3FM: Formulate
C ADV3BD: Budget calculations
C ************************************************************************************************************ 


C ************************************************************************************************************
C  ADV3AL
C  ALLOCATES SPACE FOR ARRAYS NaEEDED BY THE ADVECTION (ADV) PACKAGE.
C ************************************************************************************************************
      subroutine adv_alloc

      use rt_global
      use mf_rt_link
      use GLOBAL,     only:NCOL,NROW,NLAY


      print *, 'Reading ADV file...'

*     Read advection solution option (MIXELM,PERCEL,MXPART,NADVFD)
      !MIXELM should = 0
      !NADVFD should = 1 or 2 (upstream or central finite-difference scheme)
      read(INADV,*) MIXELM,PERCEL,MXPART
      IF(MIXELM.EQ.0 .AND. NADVFD.NE.1.AND.NADVFD.NE.2) NADVFD=1
      
*     Echo the input data
      IF(MIXELM.EQ.0.AND.NADVFD.EQ.1) WRITE(IOUT,1006)
      IF(MIXELM.EQ.0.AND.NADVFD.EQ.2) WRITE(IOUT,1007)
      
*     Check for errors
      IF(.NOT.TRNOP(5).AND.MIXELM.EQ.0.AND.PERCEL.GT.1.0) THEN
        WRITE(*,1042)
        PERCEL=1.
      ENDIF

 1006 FORMAT(1X,'ADVECTION IS SOLVED WITH',
     & ' THE UPSTREAM FINITE DIFFERENCE SCHEME')
 1007 FORMAT(1X,'ADVECTION IS SOLVED WITH',
     & ' THE CENTRAL FINITE DIFFERENCE SCHEME')
 1042 FORMAT(/1X,'WARNING: COURANR NUMBER [PERCEL] MUST NOT EXCEED 1.0',
     & /1X,'FOR THE EXPLICIT FINITE-DIFFERENCE SCHEME; ',
     & 'RESET TO DEFAULT OF 1.0.')

C     Allocate (dynamically) space for ADV arrays
C     Integer Arrays
      IF(MIXELM.GT.0) allocate(NPCHEK(NCOL,NROW,NLAY,MCOMP))
      IF(NCOL.GT.1.AND.MIXELM.GT.0) allocate(INDEXX(MXPART,MCOMP))
      IF(NROW.GT.1.AND.MIXELM.GT.0) allocate(INDEXY(MXPART,MCOMP))
      IF(NLAY.GT.1.AND.MIXELM.GT.0) allocate(INDEXZ(MXPART,MCOMP))     

C     Real Arrays
      IF(NCOL.GT.1.AND.MIXELM.GT.0) allocate(XP(MXPART,MCOMP))
      IF(NROW.GT.1.AND.MIXELM.GT.0) allocate(YP(MXPART,MCOMP))
      IF(NLAY.GT.1.AND.MIXELM.GT.0) allocate(ZP(MXPART,MCOMP))
      IF(MIXELM.GT.0) allocate(CNPT(MXPART,2,MCOMP))

      return
      end


C ***********************************************************************************************************
C ADV3FM
C THIS subroutine FORMULATES COEFFICIENT MATRICES FOR THE ADVECTION TERM WITH THE OPTIONS OF UPSTREAM 
C (NADVFD=1) AND CENTRAL (NADVFD=2) WEIGHTING.
C USED IF IMPLICIT SCHEME IS USED.
C ***********************************************************************************************************
      subroutine adv_imp(A)

      use rt_global
      use VST_ARRAYS
      use mf_rt_link
      use GLOBAL,     only:NCOL,NROW,NLAY,DELR,DELC

      IMPLICIT  NONE
      INTEGER   J,DUM,I,K,N,NCR,IUPS,ICTRL
      REAL      WW,THKSAT,AREA,ALPHA,FLUXZ,DFTEMP,DH2,QX2,QY2,QZ2,A
      DIMENSION DH2(NODES),QX2(NODES),QY2(NODES),QZ2(NODES),A(NODES,19)
      

      PARAMETER (IUPS=1,ICTRL=2)

      NCR=NROW*NCOL
      DO K=1,NLAY
        DO I=1,NROW
          DO J=1,NCOL
            N=(K-1)*NCR + (I-1)*NCOL + J
            DH2(N) = DH(J,I,K)
            QX2(N) = QX(J,I,K)
            QY2(N) = QY(J,I,K)
            QZ2(N) = QZ(J,I,K)
          ENDDO
        ENDDO
      ENDDO

C     RETURN IF COEFF MATRICES ARE NOT TO BE UPDATED
      IF(.NOT.UPDLHS) GOTO 999

C     LOOP THROUGH ALL ACTIVE CELLS
      NCR=NROW*NCOL
      DO K=1,NLAY
        DO I=1,NROW
          DO J=1,NCOL
            N=(K-1)*NCR + (I-1)*NCOL + J

C           SKIP IF INACTIVE OR CONSTANT CELL
            IF(ICBUND(J,I,K,ICOMP).LE.0) CYCLE

C           CALCULATE IN THE Z DIRECTION ------------------------------------------------------------------------------
            IF(NLAY.LT.2) GOTO 410
            AREA=DELR(J)*DELC(I)
C           TOP FACE ----------------------------------------------------------
            IF(K.GT.1) THEN
              ALPHA = 0.
              IF(NADVFD.EQ.ICTRL) ALPHA=DH2(N-NCR)/(DH2(N-NCR)+DH2(N))
              IF(NADVFD.EQ.IUPS.AND.QZ2(N-NCR).LT.0.) ALPHA=1.0
              A(N,1)=A(N,1)+ALPHA*QZ2(N-NCR)*AREA
              A(N,2)=A(N,2)+(1.-ALPHA)*QZ2(N-NCR)*AREA
            ENDIF
C           BOTTOM FACE -------------------------------------------------------
            IF(K.LT.NLAY) THEN
              ALPHA = 0.
              IF(NADVFD.EQ.ICTRL) ALPHA=DH2(N)/(DH2(N)+DH2(N+NCR))
              IF(NADVFD.EQ.IUPS.AND.QZ2(N).LT.0.) ALPHA=1.0
              A(N,1)=A(N,1)-(1.-ALPHA)*QZ2(N)*AREA
              A(N,3)=A(N,3)-ALPHA*QZ2(N)*AREA
            ENDIF

C           CALCULATE IN THE Y DIRECTION ------------------------------------------------------------------------------
  410       IF(NROW.LT.2) GOTO 420    
C           BACK FACE ---------------------------------------------------------
            IF(I.GT.1) THEN
                WW=DELC(I)/(DELC(I)+DELC(I-1))
                THKSAT=DH2(N-NCOL)*WW+DH2(N)*(1.-WW)
                AREA=DELR(J)*THKSAT
                ALPHA = 0.
                IF(NADVFD.EQ.ICTRL) ALPHA=DELC(I-1)/(DELC(I-1)+DELC(I))
                IF(NADVFD.EQ.IUPS.AND.QY2(N-NCOL).LT.0.) ALPHA=1.0
                A(N,1)=A(N,1)+ALPHA*QY2(N-NCOL)*AREA
                A(N,4)=A(N,4)+(1.-ALPHA)*QY2(N-NCOL)*AREA   
            ENDIF
C           FRONT FACE --------------------------------------------------------
            IF(I.LT.NROW) THEN
                WW=DELC(I+1)/(DELC(I+1)+DELC(I))
                THKSAT=DH2(N)*WW+DH2(N+NCOL)*(1.-WW)
                AREA=DELR(J)*THKSAT
                ALPHA = 0.
                IF(NADVFD.EQ.ICTRL) ALPHA=DELC(I)/(DELC(I)+DELC(I+1))
                IF(NADVFD.EQ.IUPS.AND.QY2(N).LT.0.) ALPHA=1.0
                A(N,1)=A(N,1)-(1.-ALPHA)*QY2(N)*AREA
                A(N,5)=A(N,5)-ALPHA*QY2(N)*AREA
            ENDIF

C           CALCULATE IN THE X DIRECTION ------------------------------------------------------------------------------
  420       IF(NCOL.LT.2) GOTO 430
C           LEFT FACE ---------------------------------------------------------
            IF(J.GT.1) THEN
                WW=DELR(J)/(DELR(J)+DELR(J-1))
                THKSAT=DH2(N-1)*WW+DH2(N)*(1.-WW)
                AREA=DELC(I)*THKSAT
                ALPHA = 0.
                IF(NADVFD.EQ.ICTRL) ALPHA=DELR(J-1)/(DELR(J-1)+DELR(J))
                IF(NADVFD.EQ.IUPS.AND.QX2(N-1).LT.0.) ALPHA=1.0
                A(N,1)=A(N,1)+ALPHA*QX2(N-1)*AREA
                A(N,6)=A(N,6)+(1.-ALPHA)*QX2(N-1)*AREA
            ENDIF
C           RIGHT FACE --------------------------------------------------------  
            IF(J.LT.NCOL) THEN
                WW=DELR(J+1)/(DELR(J+1)+DELR(J))
                THKSAT=DH2(N)*WW+DH2(N+1)*(1.-WW)
                AREA=DELC(I)*THKSAT
                ALPHA = 0.
                IF(NADVFD.EQ.ICTRL) ALPHA=DELR(J)/(DELR(J)+DELR(J+1))
                IF(NADVFD.EQ.IUPS.AND.QX2(N).LT.0.) ALPHA=1.0
                A(N,1)=A(N,1)-(1-ALPHA)*QX2(N)*AREA
                A(N,7)=A(N,7)-ALPHA*QX2(N)*AREA   
            ENDIF
C
  430       CONTINUE
          ENDDO
        ENDDO
      ENDDO

  999 RETURN
      END


C ***********************************************************************************************************
C ADV3BD
C THIS subroutine CALCULATES MASS BUDGET OF CONSTANT-CONCENTRATION NODES
C DUE TO ADVECTION.
C ***********************************************************************************************************
      subroutine adv_imp_budget

      use rt_global
      use mf_rt_link
      use GLOBAL,     only:NCOL,NROW,NLAY

      IMPLICIT  NONE
      INTEGER   J,I,K
      REAL      SADV3Q,QCTMP

*     LOOP OVER ALL MODEL CELLS
      DO K=1,NLAY
        DO I=1,NROW
          DO J=1,NCOL
            !skip if not a constant concentration cell
            IF(ICBUND(J,I,K,ICOMP).LT.0) THEN
              QCTMP = SADV3Q(J,I,K)
              IF(QCTMP.GT.0) THEN
                RMASIO(6,1,ICOMP)=RMASIO(6,1,ICOMP)+QCTMP
              ELSE
                RMASIO(6,2,ICOMP)=RMASIO(6,2,ICOMP)+QCTMP
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      RETURN
      END


C ***********************************************************************************************************
C SADV3Q
C THIS FUNCTION COMPUTES ADVECTIVE MASS FLUX BETWEEN CELL (JJ,II,KK) AND THE SURROUNDING CELLS DURING TIME 
C INCREMENT DT.  MASS ISMOVING OUT OF THE CELL IF SADV3Q > 0, INTO THE CELL IF SADV3Q < 0.
C NADVFD=1 IS FOR THE UPSTREAM SCHEME; NADVFD=2 IS FOR THE CENTRAL WEIGHTING SCHEME.
C ***********************************************************************************************************
      FUNCTION SADV3Q(JJ,II,KK)

      use rt_global
      use mf_rt_link
      use GLOBAL,     only:NCOL,NROW,NLAY,DELR,DELC

      IMPLICIT  NONE
      INTEGER   JJ,II,KK
      REAL      SADV3Q,AREA,QCTMP,WW,THKSAT,ALPHA,CTMP

C
C--SET QCTMP = 0 FOR ACCUMULATING Q*C*DTRANS IN ALL FACES
      QCTMP=0.
C
C--CALCULATE IN THE Z DIRECTION
      IF(NLAY.LT.2) GOTO 410
      AREA=DELR(JJ)*DELC(II)
C--TOP FACE
      IF(KK.GT.1) THEN
        IF(ICBUND(JJ,II,KK-1,ICOMP).NE.0) THEN
        WW=DH(JJ,II,KK)/(DH(JJ,II,KK-1)+DH(JJ,II,KK))
        ALPHA=0.
        IF(QZ(JJ,II,KK-1).GT.0) ALPHA=1.
        IF(NADVFD.EQ.2) ALPHA=WW
        CTMP=COLD(JJ,II,KK-1,ICOMP)*ALPHA + COLD(JJ,II,KK,ICOMP)*
     &      (1.-ALPHA)
        QCTMP=QCTMP-QZ(JJ,II,KK-1)*CTMP*AREA*DT
        ENDIF
      ENDIF
C--BOTTOM FACE
      IF(KK.LT.NLAY) THEN
        IF(ICBUND(JJ,II,KK+1,ICOMP).NE.0) THEN
        WW=DH(JJ,II,KK+1)/(DH(JJ,II,KK)+DH(JJ,II,KK+1))
        ALPHA=0.
        IF(QZ(JJ,II,KK).GT.0) ALPHA=1.
        IF(NADVFD.EQ.2) ALPHA=WW
        CTMP=COLD(JJ,II,KK,ICOMP)*ALPHA + COLD(JJ,II,KK+1,ICOMP)*
     &      (1.-ALPHA)
        QCTMP=QCTMP+QZ(JJ,II,KK)*CTMP*AREA*DT
        ENDIF
      ENDIF
C
C--CALCULATE IN THE Y DIRECTION
  410 IF(NROW.LT.2) GOTO 420
C--BACK FACE
      IF(II.GT.1) THEN
        IF(ICBUND(JJ,II-1,KK,ICOMP).NE.0) THEN
        WW=DELC(II)/(DELC(II)+DELC(II-1))
        THKSAT=DH(JJ,II-1,KK)*WW+DH(JJ,II,KK)*(1.-WW)
        AREA=DELR(JJ)*THKSAT
        ALPHA=0.
        IF(QY(JJ,II-1,KK).GT.0) ALPHA=1.
        IF(NADVFD.EQ.2) ALPHA=WW
        CTMP=COLD(JJ,II-1,KK,ICOMP)*ALPHA + COLD(JJ,II,KK,ICOMP)*
     &       (1.-ALPHA)
        QCTMP=QCTMP-QY(JJ,II-1,KK)*CTMP*AREA*DT
        ENDIF
      ENDIF
C--FRONT FACE
      IF(II.LT.NROW) THEN
        IF(ICBUND(JJ,II+1,KK,ICOMP).NE.0) THEN
        WW=DELC(II+1)/(DELC(II+1)+DELC(II))
        THKSAT=DH(JJ,II,KK)*WW+DH(JJ,II+1,KK)*(1.-WW)
        AREA=DELR(JJ)*THKSAT
        ALPHA=0.
        IF(QY(JJ,II,KK).GT.0) ALPHA=1.
        IF(NADVFD.EQ.2) ALPHA=WW
        CTMP=COLD(JJ,II,KK,ICOMP)*ALPHA + COLD(JJ,II+1,KK,ICOMP)*
     &      (1.-ALPHA)
        QCTMP=QCTMP+QY(JJ,II,KK)*CTMP*AREA*DT
        ENDIF
      ENDIF
C
C--CALCULATE IN THE X DIRECTION
  420 IF(NCOL.LT.2) GOTO 430
C--LEFT FACE
      IF(JJ.GT.1) THEN
        IF(ICBUND(JJ-1,II,KK,ICOMP).NE.0) THEN
        WW=DELR(JJ)/(DELR(JJ)+DELR(JJ-1))
        THKSAT=DH(JJ-1,II,KK)*WW+DH(JJ,II,KK)*(1.-WW)
        AREA=DELC(II)*THKSAT
        ALPHA=0.
        IF(QX(JJ-1,II,KK).GT.0) ALPHA=1.
        IF(NADVFD.EQ.2) ALPHA=WW
        CTMP=COLD(JJ-1,II,KK,ICOMP)*ALPHA + COLD(JJ,II,KK,ICOMP)*
     &      (1.-ALPHA)
        QCTMP=QCTMP-QX(JJ-1,II,KK)*CTMP*AREA*DT
        ENDIF
      ENDIF
C--RIGHT FACE
      IF(JJ.LT.NCOL) THEN
        IF(ICBUND(JJ+1,II,KK,ICOMP).NE.0) THEN
        WW=DELR(JJ+1)/(DELR(JJ+1)+DELR(JJ))
        THKSAT=DH(JJ,II,KK)*WW+DH(JJ+1,II,KK)*(1.-WW)
        AREA=DELC(II)*THKSAT
        ALPHA=0.
        IF(QX(JJ,II,KK).GT.0) ALPHA=1.
        IF(NADVFD.EQ.2) ALPHA=WW
        CTMP=COLD(JJ,II,KK,ICOMP)*ALPHA + COLD(JJ+1,II,KK,ICOMP)*
     &      (1.-ALPHA)
        QCTMP=QCTMP+QX(JJ,II,KK)*CTMP*AREA*DT
        ENDIF
      ENDIF
C
C--ASSIGN QCTMP TO THE FUNCTION AND RETURN
  430 SADV3Q=QCTMP

      return
      end
