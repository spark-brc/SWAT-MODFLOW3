C ************************************************************************************************************
C FMI (Flow Model Interface) Package
C Handles subroutines related to interfacing with the flow model results (hydraulic head and fluxes)
C
C July 12 2010: Changed READ from unformatted to formatted linker file (RTB)
C
C Included Subroutines:
C fmi_alloc: Allocate computer memory
C fmi_prepare:Prepare saturated cell thickness and fluxes across cell interfaces
C ************************************************************************************************************ 


C ************************************************************************************************************
C fmi_alloc
C ************************************************************************************************************
      subroutine fmi_alloc

      use rt_global
      use mf_rt_link
      use GLOBAL,     only:NCOL,NROW,NLAY

      implicit none

      INTEGER	MTWEL,MTDRN,MTRCH,MTEVT,MTRIV,MTGHB,MTCHD,MTAPPL
      CHARACTER VERSION*11
      CHARACTER HEADER*80

*     Initialize variables
      ISS = 0
      IVER = 2

      !these are already set to FALSE
      !FWEL = .FALSE.
      !FDRN = .FALSE.
      !FRCH = .FALSE.
      !FEVT = .FALSE.
      !FRIV = .FALSE.
      !FGHB = .FALSE.
      !FCNH = .FALSE.

*     Initialize for variably-saturated flow
      if (TRNOP(7)) then
        FAPP= .FALSE.
      endif

      WRITE(IOUT,*)
      WRITE(IOUT,309)
      IF(FWEL) WRITE(IOUT,340)
      IF(FDRN) WRITE(IOUT,341)
      IF(FRCH) WRITE(IOUT,342)
      IF(FEVT) WRITE(IOUT,343)
      IF(FRIV) WRITE(IOUT,344)
      IF(FGHB) WRITE(IOUT,345)
      IF(FCNH) WRITE(IOUT,346)
      WRITE(IOUT,*)

  309 FORMAT(5X,'Flow Source/Sink Components:')
  340 FORMAT(5X,'                                WELL')
  341 FORMAT(5X,'                                DRAIN')
  342 FORMAT(5X,'                                RECHARGE')
  343 FORMAT(5X,'                                EVAPOTRANSPIRATION')
  344 FORMAT(5X,'                                RIVER OR STREAM')
  345 FORMAT(5X,'                                GENERAL-HEAD')
  346 FORMAT(5X,'                                CONSTANT-HEAD')  

      return
      end


C ***********************************************************************************************************
C rtmf_prepare
C This subroutine prepares data from MODFLOW simulation (flows and head information)
C ***********************************************************************************************************
      subroutine rtmf_prepare
      
      use rt_global
      use VST_ARRAYS
      use mf_rt_link
      use GLOBAL,     only:NCOL,NROW,NLAY,DELR,DELC

      IMPLICIT	NONE
      INTEGER	J,I,K,JTRACK,KTRACK,INDEX,
     &          KKPER,KKSTP,NC,NR,NL,KSSM,ISSM,JSSM,NUM
      REAL	    WW,WTBL,THKFLUID,TK,CTMP,CREWET,THKMIN0,VOLAQU,
     &          TM
      CHARACTER TEXT*16,LABEL*16,HEADER*80


C     Use values for further calculations -----------------------------------------------
C     Assemble arrays for variably-saturated flow and transport
C     Assemble: APPL,WT,WATVOL,SAT,QZVAR
      if (TRNOP(7)) then
        call vst_prep
      endif

C     Calculate the "thickness" of fluid in each cell [DF]. If the conditions
C     are fully saturated, than this is simply the saturated thickness
C     DH. If VST package is used, then DF is calculated using the saturation 
C     index of each cell.
      do I=1,NROW
        do J=1,NCOL
          do K=1,NLAY
            IF (TRNOP(7)) THEN !VST package is used
              IF (DH(J,I,K).LT.DZ(J,I,K) .AND. 
     &            DH(J,I,K).GT.0.00) THEN !Water Table cell
C               Use weighted value of saturation for the cell portion above
C               the water table.                
                DF(J,I,K) = DH(J,I,K)+((DZ(J,I,K)-DH(J,I,K))*
     &                      SATNEW(J,I,K))
              ELSE
C               For all other cells (fully-unsaturated or fully-saturated)
                DF(J,I,K) = DZ(J,I,K)*SATNEW(J,I,K)
              ENDIF
            ELSE !Fully-saturated conditions
              DF(J,I,K) = DH(J,I,K)
            ENDIF 
          enddo
        enddo
      enddo

C--SET ICBUND=0 IF CELL IS DRY OR INACTIVE (H=1.E30)
C--AND REACTIVATE DRY CELL IF REWET AND ASSIGN CONC AT REWET CELL
      DO K=1,NLAY
	  DO I=1,NROW
	    DO J=1,NCOL
	      IF(ABS(DF(J,I,K)-1.E30).LT.1.E-5) THEN
	        ICBUND(J,I,K,1)=0
	      ELSEIF(ICBUND(J,I,K,1).EQ.0.AND.PRSITY(J,I,K).GT.0) THEN
	        ICBUND(J,I,K,1)=30000
	        DO INDEX=1,NCOMP
		        CTMP=CREWET(J,I,K,INDEX)
		        CTMP=(COLD(J,I,K,INDEX)*(RETA(J,I,K,INDEX)-1.0)+CTMP)
     &		    /RETA(J,I,K,INDEX)
		        CNEW(J,I,K,INDEX)=CTMP
	        ENDDO
	      ENDIF
	    ENDDO
	  ENDDO
      ENDDO

C--GETDF FOR CONFINED LAYERS FOR FILE SAVED BY LKMT PACKAGE VERSION 2
      IF(IVER.EQ.2) THEN
	DO K=1,NLAY
	  DO I=1,NROW
	    DO J=1,NCOL
	      IF(ICBUND(J,I,K,1).NE.0) THEN
		      IF(LAYCON(K).EQ.0 .OR. INT(DF(J,I,K)).EQ.-111) THEN
		        DF(J,I,K)=DZ(J,I,K)
		      ENDIF
              THKMIN0=THKMIN*DZ(J,I,K)
              IF((DF(J,I,K).LT.THKMIN0).or.(DF(J,I,K).eq.0.0)) THEN
                WRITE(IOUT,105) DF(J,I,K),K,I,J,THKMIN0
		        ICBUND(J,I,K,1)=0
		      ENDIF
	      ENDIF
	    ENDDO
	  ENDDO
	ENDDO

C--FOR FILE SAVED BY LKMT PACKAGE VERSION 1
      ELSEIF(IVER.EQ.1) THEN
C
	DO K=1,NLAY
	  DO I=1,NROW
	    DO J=1,NCOL
	      IF(ICBUND(J,I,K,1).NE.0) THEN
		      IF(LAYCON(K).EQ.0) THEN
		        THKFLUID=DZ(J,I,K)
		      ELSE
		        WTBL=HORIGN-DF(J,I,K)
		        THKFLUID=ZBC(J,I,K)+0.5*DZ(J,I,K)-WTBL
		        THKFLUID=MIN(THKFLUID,DZ(J,I,K))
		      ENDIF
              THKMIN0=THKMIN*DZ(J,I,K)
              IF((DF(J,I,K).LT.THKMIN0).or.(DF(J,I,K).eq.0.0)) THEN
                WRITE(IOUT,105) DF(J,I,K),K,I,J,THKMIN0
		        ICBUND(J,I,K,1)=0
		      ENDIF
	      ENDIF
	    ENDDO
	  ENDDO
	ENDDO
C
      ENDIF

  105 FORMAT(/1X,'SATURATED THICKNESS =',G13.5,'AT CELL (K,I,J):',3I5,
     & /1X,'BELOW SPECIFIED MINIMUM =',G13.5,'CELL RESET AS INACTIVE')
  355 FORMAT(/1X,'ERROR: SATURATED THICKNESS CALCULATED TO BE ZERO',
     &	' OR NEGATIVE FOR ACTIVE CELL'/1X,'       AT LAYER=',I2,
     &	'  ROW=',I4,'  COLUMN=',I4,
     &	/1X,'CHECK [HTOP] AND [DZ] ARRAYS IN [BTN] INPUT FILE.')
C
C--DETERMINE MAXIMUM TIME INCREMENT DURING WHICH ANY PARTICLE
C--CANNOT MOVE MORE THAN ONE CELL IN ANY DIRECTION.
      DTRACK=1.E30
C
      open(222222,file='dtrack')

      IF(NCOL.LT.2) GOTO 410
      DO K=1,NLAY
	  DO I=1,NROW
	    DO J=2,NCOL
	      IF(ICBUND(J,I,K,1).NE.0) THEN
	        TK=0.5*(QX(J-1,I,K)+QX(J,I,K))
	        IF(TK.EQ.0) CYCLE
	        TK=DELR(J)*DELC(I)*DF(J,I,K)*PRSITY(J,I,K)/TK
	        IF(ABS(TK).LT.DTRACK) THEN
		        DTRACK=ABS(TK)
                write(222222,*) j,i,k,dtrack,prsity(j,i,k)
		        JTRACK=J
		        ITRACK=I
		         KTRACK=K
	        ENDIF
	      ENDIF
	    ENDDO
	  ENDDO
      ENDDO
C
  410 IF(NROW.LT.2) GOTO 420
      DO K=1,NLAY
	  DO J=1,NCOL
	    DO I=2,NROW
	      IF(ICBUND(J,I,K,1).NE.0) THEN
	        TK=0.5*(QY(J,I-1,K)+QY(J,I,K))
	        IF(TK.EQ.0) CYCLE
	        TK=DELR(J)*DELC(I)*DF(J,I,K)*PRSITY(J,I,K)/TK
	        IF(ABS(TK).LT.DTRACK) THEN
		        DTRACK=ABS(TK)
                write(222222,*) j,i,k,dtrack,prsity(j,i,k)
		        JTRACK=J
		        ITRACK=I
		        KTRACK=K
	        ENDIF
	      ENDIF
	    ENDDO
	  ENDDO
      ENDDO
C
  420 IF(NLAY.LT.2) GOTO 430
      DO J=1,NCOL
	  DO I=1,NROW
	    DO K=2,NLAY
	      IF(ICBUND(J,I,K,1).NE.0) THEN
	        TK=0.5*(QZ(J,I,K-1)+QZ(J,I,K))
	        IF(TK.EQ.0) CYCLE
	        TK=DELR(J)*DELC(I)*DF(J,I,K)*PRSITY(J,I,K)/TK
	        IF(ABS(TK).LT.DTRACK) THEN
		        DTRACK=ABS(TK)
                write(222222,*) j,i,k,dtrack,prsity(j,i,k)
		        JTRACK=J
		        ITRACK=I
		        KTRACK=K
	        ENDIF
	      ENDIF
	    ENDDO
	  ENDDO
      ENDDO
C
C--PRINT INFORMATION ON DTRACK
  430 CONTINUE
C      WRITE(IOUT,500) DTRACK,KTRACK,ITRACK,JTRACK
  500 FORMAT(/1X,'MAXIMUM STEPSIZE DURING WHICH ANY PARTICLE CANNOT',
     & ' MOVE MORE THAN ONE CELL'/1X,'=',G11.4,
     & '(WHEN MIN. R.F.=1)  AT K=',I4,', I=',I4,
     & ', J=',I4)
C
C--DETERMINE STABILITY CRITERION ASSOCIATED WITH EXPLICIT FINITE
C--DIFFERENCE SOLUTION OF THE ADVECTION TERM
      DTRACK2=1.E30
C
      DO K=1,NLAY
	  DO I=1,NROW
	    DO J=1,NCOL
	      IF(ICBUND(J,I,K,1).EQ.0) CYCLE
	      TK=0.
	      IF(J.GT.1) TK=TK+MAX(ABS(QX(J-1,I,K)),ABS(QX(J,I,K)))
	      IF(I.GT.1) TK=TK+MAX(ABS(QY(J,I-1,K)),ABS(QY(J,I,K)))
	      IF(K.GT.1) TK=TK+MAX(ABS(QZ(J,I,K-1)),ABS(QZ(J,I,K)))
	      IF(TK.EQ.0) CYCLE
	      TK=DELR(J)*DELC(I)*DZ(J,I,K)*PRSITY(J,I,K)/TK
	      IF(TK.LT.DTRACK2) THEN
	        DTRACK2=TK
	        JTRACK=J
	        ITRACK=I
	        KTRACK=K
	      ENDIF
	    ENDDO
	  ENDDO
      ENDDO
C
C--PRINT INFORMATION ON DTRACK2
C      WRITE(IOUT,550) DTRACK2,KTRACK,ITRACK,JTRACK
  550 FORMAT(/1X,'MAXIMUM STEPSIZE WHICH MEETS STABILITY CRITERION',
     & ' OF THE ADVECTION TERM'/1X,
     & '(FOR PURE FINITE-DIFFERENCE OPTION, MIXELM=0) '/1X,'=',G11.4,
     & '(WHEN MIN. R.F.=1)  AT K=',I4,', I=',I4,', J=',I4)

C     DIVIDE VOLUMETRIC QX, QY AND QZ BY AREAS TO GET SPECIFIC DISCHARGE 
C     (L/T) ACROSS EACH CELL INTERFACE
C     QX
C     For VST package, DF handles cases of variably-saturated content.          
      IF(NCOL.LT.2) GOTO 910
      DO K=1,NLAY
	  DO I=1,NROW
	    DO J=1,NCOL-1
	      WW=DELR(J+1)/(DELR(J+1)+DELR(J)) !Weighted, using adjacent cells
	      THKFLUID=DF(J,I,K)*WW+DF(J+1,I,K)*(1.-WW)
	      IF(THKFLUID.LE.0.OR.ICBUND(J,I,K,1).EQ.0) THEN
	        QX(J,I,K)=0
              IF(J.GT.1) QX(J-1,I,K)=0.
            ELSE
	        QX(J,I,K)=QX(J,I,K)/(DELC(I)*THKFLUID)
	      ENDIF
	    ENDDO
          IF(ICBUND(NCOL,I,K,1).EQ.0) QX(NCOL-1,I,K)=0.
        ENDDO
      ENDDO
      
C     QY
C     For VST package, DF handles cases of variably-saturated content.
  910 IF(NROW.LT.2) GOTO 920
      DO K=1,NLAY
	  DO J=1,NCOL
	    DO I=1,NROW-1
	      WW=DELC(I+1)/(DELC(I+1)+DELC(I))
	      THKFLUID=DF(J,I,K)*WW+DF(J,I+1,K)*(1.-WW)
	      IF(THKFLUID.LE.0.OR.ICBUND(J,I,K,1).EQ.0) THEN
	        QY(J,I,K)=0
              IF(I.GT.1) QY(J,I-1,K)=0.
            ELSE
	        QY(J,I,K)=QY(J,I,K)/(DELR(J)*THKFLUID)
	      ENDIF
	    ENDDO
          IF(ICBUND(J,NROW,K,1).EQ.0) QY(J,NROW-1,K)=0.
        ENDDO
      ENDDO

C     QZ
  920 IF(NLAY.LT.2) GOTO 980
      DO J=1,NCOL
	  DO I=1,NROW
          DO K=1,NLAY
	      THKFLUID=DF(J,I,K)
	      IF(THKFLUID.LE.0.OR.ICBUND(J,I,K,1).EQ.0) THEN
	        QZ(J,I,K)=0
              IF(K.GT.1) QZ(J,I,K-1)=0.
            ELSE
C             If VST package is used, use SAT index for the
C             bottom cross-section area.
	        QZ(J,I,K)=QZ(J,I,K)/(DELR(J)*DELC(I))
	      ENDIF
	    ENDDO
	  ENDDO
      ENDDO

C     DIVIDE STORAGE BY CELL VOLUME TO GET DIMENSION (1/TIME)
  980 CONTINUE
      DO K=1,NLAY
	  DO I=1,NROW
	    DO J=1,NCOL
	      THKFLUID=DF(J,I,K)
	      IF(THKFLUID.LE.0.OR.ICBUND(J,I,K,1).EQ.0) THEN
	        QSTO(J,I,K)=0
	      ELSE
	        QSTO(J,I,K)=QSTO(J,I,K)/(THKFLUID*DELR(J)*DELC(I))
	      ENDIF
	    ENDDO
	  ENDDO
      ENDDO
      
C     DIVIDE VOLUMETRIC QZVAR BY BOTTOM AREA TO GET SPECIFIC DISCHARGE 
C     (L/T) ACROSS EACH CELL INTERFACE
      IF (TRNOP(7)) THEN
        DO J=1,NCOL
          DO I=1,NROW
            DO K=1,NLAY
              THKFLUID=DF(J,I,K)
              IF(THKFLUID.LE.0.OR.ICBUND(J,I,K,1).EQ.0) THEN
                QZVAR(J,I,K)=0
                IF(K.GT.1) QZVAR(J,I,K-1)=0.
              ELSE            
                QZVAR(J,I,K)=QZVAR(J,I,K)/(DELR(J)*DELC(I))
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF

C     DIVIDE APPL BY AQUIFER VOLUME (VOLUME THAT RECEIVES SOURCE WATER)
C     TO GET (1/T)
      IF (TRNOP(7)) THEN
        DO I=1,NROW
	    DO J=1,NCOL
	      K=1 !APPL occurs only for top cell layer
	      IF(K.EQ.0) CYCLE
	      VOLAQU = DELR(J)*DELC(I)*DZ(J,I,K)
	      IF(ICBUND(J,I,K,1).EQ.0.OR.VOLAQU.LE.0) THEN
	        APPL(J,I) = 0.0
	      ELSE
	        APPL(J,I) = APPL(J,I) / VOLAQU
	      ENDIF
	    ENDDO
        ENDDO
      ENDIF

C     DIVIDE UZQSTO BY AQUIFER VOLUME (VOLUME WHICH RECEIVES SOURCE WATER)
C     TO GET DIMENSION (1/T)
      IF (TRNOP(7)) THEN
        DO K=1,NLAY
	    DO I=1,NROW
	      DO J=1,NCOL
c	        THKFLUID=DF(J,I,K)
	        THKFLUID=DZ(J,I,K)
	        IF(THKFLUID.LE.0.OR.ICBUND(J,I,K,1).EQ.0) THEN
	          UZQSTO(J,I,K)=0.0
	        ELSE
	          UZQSTO(J,I,K)=UZQSTO(J,I,K)/(THKFLUID*DELR(J)*DELC(I))
	        ENDIF
	      ENDDO
	    ENDDO
        ENDDO
      ENDIF

C     SYNCHRONIZE ICBUND CONDITIONS OF ALL SPECIES
      IF(NCOMP.EQ.1) GOTO 999
      DO K=1,NLAY
	  DO I=1,NROW
	    DO J=1,NCOL
	      DO INDEX=2,NCOMP
	        IF(ICBUND(J,I,K,INDEX).GE.0) THEN
                ICBUND(J,I,K,INDEX)=IABS(ICBUND(J,I,K,1))
	        ELSEIF(ICBUND(J,I,K,1).EQ.0) THEN
		        ICBUND(J,I,K,INDEX)=0
	        ENDIF
	      ENDDO
	    ENDDO
	  ENDDO
      ENDDO

C     RESET NUMBER OF POINT SOURCES OF SPECIFIED CONCENTRATIONS [NSS]
C     AND MAKE [ICBUND] VALUE AT ACTIVE CELLS EQUAL TO 1
 999  NTSS=NSS
      DO K=1,NLAY
	  DO I=1,NROW
	    DO J=1,NCOL
	      IF(ICBUND(J,I,K,1).GT.0) ICBUND(J,I,K,1)=1
	    ENDDO
	  ENDDO
      ENDDO

C     DIVIDE RECH, EVTR, Q_SS BY AQUIFER VOLUME TO GET FLUXES OF SINKS/SOURCES
C     PER UNIT AQUIFER VOLUME. ALSO DETERMINE STEPSIZE WHICH MEETS STABILITY 
C     CRITERION FOR SOLVING THE SINK/SOURCE TERM

C     For general conditions, the variable DF (fluid thickness) has replaced DH (saturated
C     thickness) as the value to use when calculating effective aquifer volumes for each 
C     cell.

C     Keep track of transport step size required by sources/sinks (to be used in BTN3AD when
C     calculating next transport stez size)      
      DTSSM=1.E30

C     RECHARGE
      IF(FRCH) THEN
        DO I=1,NROW
	    DO J=1,NCOL
	      K=RT_IRCH(J,I)
	      IF(K.EQ.0) CYCLE
C           VOLAQU is the bulk volume of aquifer (in the cell) that can receive
C           the source fluid and accompanying solute mass.
	      VOLAQU = DELR(J)*DELC(I)*DF(J,I,K)
	      IF(ICBUND(J,I,K,1).EQ.0.OR.VOLAQU.LE.0) THEN
	        RT_RECH(J,I)=0.
	      ELSE
	        RT_RECH(J,I) = RT_RECH(J,I)/VOLAQU
	      ENDIF
            IF(RT_RECH(J,I).LE.0 .OR. ICBUND(J,I,K,1).EQ.0) CYCLE
	      TM = PRSITY(J,I,K)/RT_RECH(J,I)
	      IF(ABS(TM).LT.DTSSM) THEN
	        DTSSM=ABS(TM)
	        KSSM=K
	        ISSM=I
	        JSSM=J
	      ENDIF
	    ENDDO
        ENDDO
      ENDIF

C     POINT SOURCES/SINKS (WEL,DRN,RIV,GHB,CHD)
      
      !Wells
      do num=1,rt_nwel
        k = WEL(num,1)
        i = WEL(num,2)
        j = WEL(num,3)
        !VOLAQU is the bulk volume of aquifer (in the cell) that can receive
        !the source fluid and accompanying solute mass.
        VOLAQU = DELR(J)*DELC(I)*DZ(J,I,K)
        if(ICBUND(J,I,K,1).EQ.0.OR.VOLAQU.LE.0) then
	    WEL(num,4)=0
	  else
	    WEL(num,4) = WEL(num,4) / VOLAQU
	  endif
        if(WEL(num,4).LE.0 .OR. ICBUND(J,I,K,1).EQ.0) cycle
        TM = PRSITY(J,I,K) / WEL(num,4)
	  if(ABS(TM).LT.DTSSM) then
	    DTSSM=ABS(TM)
	    KSSM=K
	    ISSM=I
	    JSSM=J
	  endif
      enddo

      !Drains
      do num=1,rt_ndrn
        k = DRN(num,1)
        i = DRN(num,2)
        j = DRN(num,3)
        !VOLAQU is the bulk volume of aquifer (in the cell) that can receive
        !the source fluid and accompanying solute mass.
        VOLAQU = DELR(J)*DELC(I)*DZ(J,I,K)
        if(ICBUND(J,I,K,1).EQ.0.OR.VOLAQU.LE.0) then
	    DRN(num,4)=0
	  else
	    DRN(num,4) = DRN(num,4) / VOLAQU
	  endif
        if(DRN(num,4).LE.0 .OR. ICBUND(J,I,K,1).EQ.0) cycle
        TM = PRSITY(J,I,K) / DRN(num,4)
	  if(ABS(TM).LT.DTSSM) then
	    DTSSM=ABS(TM)
	    KSSM=K
	    ISSM=I
	    JSSM=J
	  endif
      enddo

      !River
      do num=1,rt_nriv
        k = RIV(num,1)
        i = RIV(num,2)
        j = RIV(num,3)
        !VOLAQU is the bulk volume of aquifer (in the cell) that can receive
        !the source fluid and accompanying solute mass.
        VOLAQU = DELR(J)*DELC(I)*DZ(J,I,K)
        if(ICBUND(J,I,K,1).EQ.0.OR.VOLAQU.LE.0) then
	    RIV(num,4)=0
	  else
	    RIV(num,4) = RIV(num,4) / VOLAQU
	  endif
        if(RIV(num,4).LE.0 .OR. ICBUND(J,I,K,1).EQ.0) cycle
        TM = PRSITY(J,I,K) / RIV(num,4)
	  if(ABS(TM).LT.DTSSM) then
	    DTSSM=ABS(TM)
	    KSSM=K
	    ISSM=I
	    JSSM=J
	  endif
      enddo

      !Constant-head cells
      do num=1,rt_ncnh
        k = CNH(num,1)
        i = CNH(num,2)
        j = CNH(num,3)
        !VOLAQU is the bulk volume of aquifer (in the cell) that can receive
        !the source fluid and accompanying solute mass.
        VOLAQU = DELR(J)*DELC(I)*DZ(J,I,K)
        if(ICBUND(J,I,K,1).EQ.0.OR.VOLAQU.LE.0) then
	    CNH(num,4)=0
	  else
	    CNH(num,4) = CNH(num,4) / VOLAQU
	  endif
        if(CNH(num,4).LE.0 .OR. ICBUND(J,I,K,1).EQ.0) cycle
        TM = PRSITY(J,I,K) / CNH(num,4)
	  if(ABS(TM).LT.DTSSM) then
	    DTSSM=ABS(TM)
	    KSSM=K
	    ISSM=I
	    JSSM=J
	  endif
      enddo

      !General-head boundary cells
      do num=1,rt_nghb
        k = GHB(num,1)
        i = GHB(num,2)
        j = GHB(num,3)
        !VOLAQU is the bulk volume of aquifer (in the cell) that can receive
        !the source fluid and accompanying solute mass.
        VOLAQU = DELR(J)*DELC(I)*DZ(J,I,K)
        if(ICBUND(J,I,K,1).EQ.0.OR.VOLAQU.LE.0) then
	    GHB(num,4)=0
	  else
	    GHB(num,4) = GHB(num,4) / VOLAQU
	  endif
        if(GHB(num,4).LE.0 .OR. ICBUND(J,I,K,1).EQ.0) cycle
        TM = PRSITY(J,I,K) / GHB(num,4)
	  if(ABS(TM).LT.DTSSM) then
	    DTSSM=ABS(TM)
	    KSSM=K
	    ISSM=I
	    JSSM=J
	  endif
      enddo 

C--PRINT INFORMATION ON DTSSM
  990 CONTINUE
C      WRITE(IOUT,1000) DTSSM,KSSM,ISSM,JSSM
 1000 FORMAT(/1X,'MAXIMUM STEPSIZE WHICH MEETS STABILITY CRITERION',
     & ' OF THE SINK & SOURCE TERM'/1X,'=',G11.4,
     & '(WHEN MIN. R.F.=1)  AT K=',I4,', I=',I4,
     & ', J=',I4)

C--SYNCHRONIZE ICBUND CONDITIONS OF ALL SPECIES
      IF(NCOMP.EQ.1) GOTO 1999
      DO K=1,NLAY
	  DO I=1,NROW
	    DO J=1,NCOL
	      DO INDEX=2,NCOMP
	        IF(ICBUND(J,I,K,INDEX).GE.0) THEN
                ICBUND(J,I,K,INDEX)=IABS(ICBUND(J,I,K,1))
	        ELSEIF(ICBUND(J,I,K,1).EQ.0) THEN
		        ICBUND(J,I,K,INDEX)=0
	        ENDIF
	      ENDDO
	    ENDDO
	  ENDDO
      ENDDO

 1999 return
      end


C ***********************************************************************************************************
C CREWET
C THIS FUNCTION OBTAINS CONCENTRATION AT A REWET CELL (JJ,II,KK) FROM CONCENTRATIONS AT NEIGHBORING NODES 
C WITH INVERSE DISTANCE(POWER 2) WEIGHTING .
C ***********************************************************************************************************
      FUNCTION CREWET(JJ,II,KK,INDEX)

      use rt_global
      use GLOBAL,     only:NCOL,NROW,NLAY

      IMPLICIT	NONE
      INTEGER	JJ,II,KK,INDEX
      REAL	    CTMP,CREWET,D2,D2SUM
C
C--INITIALIZE
      D2SUM=0
      CTMP=0
C
C--ACCUMULATE CONCENTRATIONS AT NEIGHBORING NODELS
C--IN THE LAYER DIRECTION
      IF(NLAY.EQ.1) GOTO 10
      IF(KK-1.GT.0) THEN
	IF(ICBUND(JJ,II,KK-1,INDEX).NE.0) THEN
	  D2=(ZBC(JJ,II,KK)-ZBC(JJ,II,KK-1))**2
	  IF(D2.NE.0) THEN
	    D2SUM=D2SUM+1./D2
	    CTMP=CTMP+CNEW(JJ,II,KK-1,INDEX)/D2
	  ELSE
	    CTMP=CNEW(JJ,II,KK-1,INDEX)
	    GOTO 100
	  ENDIF
	ENDIF
      ENDIF
      IF(KK+1.LE.NLAY) THEN
	IF(ICBUND(JJ,II,KK+1,INDEX).NE.0) THEN
	  D2=(ZBC(JJ,II,KK)-ZBC(JJ,II,KK+1))**2
	  IF(D2.NE.0) THEN
	    D2SUM=D2SUM+1./D2
	    CTMP=CTMP+CNEW(JJ,II,KK+1,INDEX)/D2
	  ELSE
	    CTMP=CNEW(JJ,II,KK+1,INDEX)
	    GOTO 100
	  ENDIF
	ENDIF
      ENDIF
C
C--IN THE ROW DIRECTION
   10 IF(NROW.EQ.1) GOTO 20
      IF(II-1.GT.0) THEN
	IF(ICBUND(JJ,II-1,KK,INDEX).NE.0) THEN
	  D2=(YBC(II)-YBC(II-1))**2
	  IF(D2.NE.0) THEN
	    D2SUM=D2SUM+1./D2
	    CTMP=CTMP+CNEW(JJ,II-1,KK,INDEX)/D2
	  ELSE
	    CTMP=CNEW(JJ,II-1,KK,INDEX)
	    GOTO 100
	  ENDIF
	ENDIF
      ENDIF
      IF(II+1.LE.NROW) THEN
	IF(ICBUND(JJ,II+1,KK,INDEX).NE.0) THEN
	  D2=(YBC(II)-YBC(II+1))**2
	  IF(D2.NE.0) THEN
	    D2SUM=D2SUM+1./D2
	    CTMP=CTMP+CNEW(JJ,II+1,KK,INDEX)/D2
	  ELSE
	    CTMP=CNEW(JJ,II+1,KK,INDEX)
	    GOTO 100
	  ENDIF
	ENDIF
      ENDIF
C
C--IN THE COLUMN DIRECTION
   20 IF(NCOL.EQ.1) GOTO 30
      IF(JJ-1.GT.0) THEN
	IF(ICBUND(JJ-1,II,KK,INDEX).NE.0) THEN
	  D2=(XBC(JJ)-XBC(JJ-1))**2
	  IF(D2.NE.0) THEN
	    D2SUM=D2SUM+1./D2
	    CTMP=CTMP+CNEW(JJ-1,II,KK,INDEX)/D2
	  ELSE
	    CTMP=CNEW(JJ-1,II,KK,INDEX)
	    GOTO 100
	  ENDIF
	ENDIF
      ENDIF
      IF(JJ+1.LE.NCOL) THEN
	IF(ICBUND(JJ+1,II,KK,INDEX).NE.0) THEN
	  D2=(XBC(JJ)-XBC(JJ+1))**2
	  IF(D2.NE.0) THEN
	    D2SUM=D2SUM+1./D2
	    CTMP=CTMP+CNEW(JJ+1,II,KK,INDEX)/D2
	  ELSE
	    CTMP=CNEW(JJ+1,II,KK,INDEX)
	    GOTO 100
	  ENDIF
	ENDIF
      ENDIF
C
C--OBTAIN WEIGHTED CONCENTRATION
   30 IF(D2SUM.EQ.0) THEN
	ICBUND(JJ,II,KK,INDEX)=0
      ELSE
	CTMP=CTMP/D2SUM
      ENDIF
C
C--ASSIGN WEIGHTED CONCENTRATION TO CREWET
  100 CREWET=CTMP

      return
      end
