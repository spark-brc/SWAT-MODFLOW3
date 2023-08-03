C ************************************************************************************************************
C btn_define
C THIS subroutine READS PROBLEM DIMENSIONS AND TRANSPORT OPTIONS.
C ************************************************************************************************************
      subroutine btn_define

      use rt_global
      use VST_ARRAYS
      use mf_rt_link
      use GLOBAL,     only:NCOL,NROW,NLAY

      implicit none

*     local variables
      INTEGER   I,J,COMPTYPE,INDEX,NUMNUTRIENTS,STPERLEN,LENFLAG
      CHARACTER HEADNG(2)*80,HEADER*80,COMP*80
	DIMENSION STPERLEN(1000)

*     Set maximum values for certain variables
      LENX=999999999
      LENIX=999999999
      !maximum size of variables      
      MXSTP=1000000
      MXOBS=1000
      MXCOMP=100
      MXCOL=200
      MXROW=200
      MXLAY=20
      MXCROP=15
      MXPTS=1000
      MXFTIMES=30
      MXYR=20
      MXCYCLE=10
      MXPROC=28

*     Read number of species (and number of mobile species)
      READ(INBTN,*) HEADER
      READ(INBTN,*) NCOMP,MCOMP

*     Allocate arrays that track mass and mass balance errors
      allocate(TMASIO(28,2,NCOMP))
      allocate(RMASIO(28,2,NCOMP))
      allocate(TMASS(6,NCOMP))
      allocate(TMASIN(NCOMP))
      allocate(TMASOUT(NCOMP))
      allocate(ERROR(NCOMP))
      allocate(ERROR2(NCOMP))
      allocate(ERROR2_TS(NCOMP))
      allocate(SPECIES(NCOMP)) !names of the species
      allocate(TSLNGH(MXSTP))

      !set total mass arrays to 0
      TMASIO = 0.
      
C     Read in flag for concentration output. 0=ASCII, 1=BINARY, 2=BOTH
      READ(INBTN,*) HEADER
      READ(INBTN,*) WRITEFLAG

C     Open up output files for each species
      READ(INBTN,*) HEADER
      INDEX = 1
      DO I=1,NCOMP
C       Read in species name and index number
        READ(INBTN,*) SPECIES(I),PHASE(I),WRITE_OUT(I)
        COMP = SPECIES(I)

C       Open up file to output XY (surface contour) concentration at each output time
        OPEN(ICONC+I,FILE='swatmf_out_RT_CONC'//COMP,
     &     STATUS='UNKNOWN')
      ENDDO   

C     Allocate arrays for leaching calculations
      allocate(TLEACH(NCOL,NROW,NLAY,NCOMP))
      allocate(LEACH(NCOL,NROW,NLAY,NCOMP))
      TLEACH = 0.0 !Initialize
      LEACH = 0.0 !Initialize   

C     Total number of cells in the 3D grid
      NODES = NCOL*NROW*NLAY

C     Write out information to main output file
      WRITE(IOUT,1000) NLAY,NROW,NCOL
      WRITE(IOUT,1001) NODES
      WRITE(IOUT,*)
      WRITE(IOUT,1002) NPER
      WRITE(IOUT,*)
      WRITE(IOUT,1004) NCOMP
      WRITE(IOUT,1005) MCOMP
            
C     Read and print units of time, length, and mass
      READ(INBTN,*) TUNIT,LUNIT,MUNIT
      WRITE(IOUT,1006) TUNIT,LUNIT,MUNIT
 
C     Write out the packages used in the simulation
      WRITE(IOUT,*)
      WRITE(IOUT,1007)
      IF (TRNOP(1)) THEN
        WRITE(IOUT,1008)
      ENDIF
      IF (TRNOP(2)) THEN
        WRITE(IOUT,1009)
      ENDIF
      IF (TRNOP(3)) THEN
        WRITE(IOUT,1010)
      ENDIF
      IF (TRNOP(4)) THEN
        WRITE(IOUT,1011)
      ENDIF
      IF (TRNOP(5)) THEN
        WRITE(IOUT,1012)
      ENDIF
      IF (TRNOP(7)) THEN
        WRITE(IOUT,1014)
      ENDIF
 
 1000 FORMAT(5X,'                    3D Grid:',1X,I5,' Layers',I5,
     & ' Rows',I5,' Columns')
 1001 FORMAT(5X,'         Total No. of Cells:',1X,I5)
 1002 FORMAT(5X,'      No. of Stress Periods:',3X,I5)
 1003 FORMAT(5X,' Total length of simulation:',1X,F10.2)    
 1004 FORMAT(5X,'             No. of Species:',1X,I5)
 1005 FORMAT(5X,'      No. of Mobile Species:',I5)   
 1006 FORMAT(5X,'  Units of Time,Length,Mass:'4X,A4,2X,A4,2X,A4)
 1007 FORMAT(5X,'              Packages used:')
 1008 FORMAT(5X,'                                ADV')
 1009 FORMAT(5X,'                                DSP')
 1010 FORMAT(5X,'                                SSM')
 1011 FORMAT(5X,'                                RCT')
 1012 FORMAT(5X,'                                GCG')
 1013 FORMAT(5X,'                                NTR')
 1014 FORMAT(5X,'                                VST')
 1015 FORMAT(5X,'                                IRG')

      return
      end


C ***********************************************************************************************************
C btn_alloc
C THIS subroutine ALLOCATES SPACE FOR ARRAYS NEEDED BY THE ENTIRE MODEL
C ***********************************************************************************************************
      subroutine btn_alloc

      use rt_global
      use mf_rt_link
      use GLOBAL,     only:NCOL,NROW,NLAY

      implicit none

*     Allocate general arrays for the simulation

      !Allocate integer arrays
      allocate(LAYCON(NLAY))
      allocate(ICBUND(NCOL,NROW,NLAY,NCOMP))

      !Allocate arrays for grid formulation and coordinates
      allocate(HTOP(NCOL,NROW))
      allocate(CELLTOP(NCOL,NROW,NLAY))
      allocate(DZ(NCOL,NROW,NLAY))
      allocate(DF(NCOL,NROW,NLAY))            
      allocate(XBC(NCOL))
      allocate(YBC(NROW))
      allocate(ZBC(NCOL,NROW,NLAY))
      
      !Allocate arrays for concentration requirements           
      allocate(CNEW(NCOL,NROW,NLAY,NCOMP))
      allocate(COLD(NCOL,NROW,NLAY,NCOMP))
      allocate(CWGT(NCOL,NROW,NLAY,NCOMP))
      allocate(PRSITY(NCOL,NROW,NLAY))
      allocate(RETA(NCOL,NROW,NLAY,NCOMP))
      allocate(CADV(NCOL,NROW,NLAY,NCOMP))
      allocate(CSR(NCOL,NROW,NLAY,NCOMP))

      !Allocate arrays for flow terms
      allocate(QX(NCOL,NROW,NLAY))
      allocate(QY(NCOL,NROW,NLAY))
      allocate(QZ(NCOL,NROW,NLAY))
      allocate(QSTO(NCOL,NROW,NLAY))
      allocate(DH(NCOL,NROW,NLAY))   

      return
      end


C ***********************************************************************************************************
C btn_read
C THIS subroutine READS AND PREPARES INPUT DATA RELEVANT TO THE ENTIRE SIMULATION
C************************************************************************************************************
      subroutine btn_read

      use rt_global
      USE VST_ARRAYS
      use mf_rt_link
      use GLOBAL,     only:NCOL,NROW,NLAY,DELR,DELC,BOTM,IBOUND

      implicit none

*     local variables
      INTEGER   N,J,I,K,INDEX,READNUMBER,DUM,II,JJ,KK
      INTEGER   IVERSION,IOBTY,N6,ISFLT,JSFLT,ISFLG,JSFLG,ISCL,INODE,
     &          IELEM,NNP,INAME,IERR,IP1
      INTEGER   icbund_values(NCOL,NROW)
      REAL      TEMP,TOP,VOLUME,CONSTANT,SOLID,ZZ
      CHARACTER ANAME*24
	CHARACTER EXT*80, HEADER*80
	INTEGER*4 II3
	CHARACTER*80 FLNAME2*40
      
      print *, 'Reading BTN file...'
      print *, ''

C     STORE TOP ELEVATION OF 1st LAYER (from MODFLOW grid)
      do I=1,NROW
        do J=1,NCOL
          HTOP(J,I) = BOTM(J,I,0)
        enddo
      enddo
      
C     STORE THICKNESS OF EACH LAYER
      DO K=1,NLAY
        DO I=1,NROW
          DO J=1,NCOL
            DZ(J,I,K) = BOTM(J,I,K-1) - BOTM(J,I,K)
          ENDDO
        ENDDO
      ENDDO      

C     Calculate the top elevation of every cell in the grid
      DO I=1,NROW
        DO J=1,NCOL
          TOP = HTOP(J,I)
          DO K=1,NLAY
            CELLTOP(J,I,K) = TOP
            TOP = TOP - DZ(J,I,K)
          ENDDO
        ENDDO
      ENDDO

C     Read porosity
      READ(INBTN,*) HEADER
      DO K=1,NLAY
        !READ(INBTN,*) READNUMBER,CONSTANT
        !IF (READNUMBER.EQ.0) THEN
        !  DO I=1,NROW
        !    DO J=1,NCOL
        !      PRSITY(J,I,K) = CONSTANT
        !    ENDDO
        !  ENDDO
        !ELSE
          READ(INBTN,*) ((PRSITY(J,I,K),J=1,NCOL),I=1,NROW)
        !ENDIF
      ENDDO      

C     Read ICBUND (active cell) array
      read(inbtn,*)
      do i=1,nrow
        read(inbtn,*) (icbund_values(j,i),j=1,ncol)
      enddo
      do i=1,nrow
        do j=1,ncol
          do K=1,NLAY
            ICBUND(j,i,k,1) = icbund_values(j,i)
          enddo
        enddo
      enddo

C     Read in initial concentrations for each species
      READ(INBTN,*) HEADER
      DO INDEX=1,NCOMP
        DO K=1,NLAY
          READ(INBTN,*) READNUMBER,CONSTANT
          IF (READNUMBER.EQ.0) THEN
            DO I=1,NROW
              DO J=1,NCOL
                COLD(J,I,K,INDEX) = CONSTANT
              ENDDO
            ENDDO
          ELSE
            READ(INBTN,*) ((COLD(J,I,K,INDEX),J=1,NCOL),I=1,NROW)
          ENDIF           
        ENDDO
      ENDDO        

C     If Unsaturated conditions, read in initial moisture content for each cell. 
C     This is required at this point in the data input so that the appropriate 
C     Retardation Coefficients (using moisture content instead of porosity) can 
C     be calculated in the RCT3RP subroutine. 
C     If fully saturated conditions are simulated, copy PRSITY into MC array
C     (Moisture Content = Porosity).
C     READ MOISTURE CONTENT VALUES FOR EACH CELL IN THE GRID
      IF (TRNOP(7)) THEN
        READ(INBTN,*) HEADER
        READ(INBTN,*) (((MC(J,I,K),J=1,NCOL),I=1,NROW),K=1,NLAY)
C       If there are cells with MC = 0, change to residual MC (otherwise,
C       calculations such as for Retardation Coefficient will = Infinity)
        DO I=1,NROW
          DO J=1,NCOL
            DO K=1,NLAY
              IF (MC(J,I,K).LE.0) THEN                             
                MC(J,I,K) = 0.100
              ENDIF
              IF (MC(J,I,K).GT.PRSITY(J,I,K)) THEN
                MC(J,I,K) = PRSITY(J,I,K)
              ENDIF
              MCOLD(J,I,K) = MC(J,I,K)
            ENDDO
          ENDDO
        ENDDO         
C       Calculate cell saturation
        DO I=1,NROW
          DO J=1,NCOL
            DO K=1,NLAY
              SATOLD(J,I,K) = MC(J,I,K) / PRSITY(J,I,K)
            ENDDO
          ENDDO
        ENDDO      
      ENDIF

C     Calculate MOLD (Mass Old, for each solute in each cell) using COLD,
C     and the volume and moisture content of the cell. 
C     Only if the VST package is used.
      IF (TRNOP(7)) THEN
        DO INDEX=1,NCOMP
          DO J=1,NCOL
            DO I=1,NROW
              DO K=1,NLAY
                VOLUME = DELR(J)*DELC(I)*DZ(J,I,K)
                IF (PHASE(INDEX).EQ.1) THEN !Specie is mobile (fluid-phase)
C                 M = [M/L3fluid] * [L3fluid/L3bulk] * [L3bulk]                  
                  MOLD(J,I,K,INDEX)=COLD(J,I,K,INDEX)*MC(J,I,K)*VOLUME
                  MTEMP(J,I,K,INDEX)=MOLD(J,I,K,INDEX)
                ELSE !Specie is immobile (solid-phase)
C                 M = [M/L3solid] * [L3fluid/L3bulk] * [L3bulk]
                  SOLID = 1 - PRSITY(J,I,K)
                  MOLD(J,I,K,INDEX)=COLD(J,I,K,INDEX)*SOLID*VOLUME
                  MTEMP(J,I,K,INDEX)=MOLD(J,I,K,INDEX)
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO  
      ENDIF

C     READ AND ECHO CINACT,THKMIN
      READ(INBTN,*) HEADER
      READ(INBTN,'(2F10.0)',ERR=50,IOSTAT=IERR) CINACT,THKMIN
      IF(THKMIN.LE.0.0) THKMIN=0.0
   50 IF(IERR.NE.0) THEN
        BACKSPACE (INBTN)
        READ(INBTN,'(F10.0)') CINACT
        THKMIN=0.0
      ENDIF
C      WRITE(IOUT,1020) CINACT,THKMIN
      IF(THKMIN .GT. 0.1) THEN
        WRITE(IOUT,1022)
        THKMIN=0.01
      ENDIF
 1020 FORMAT(/1X,'VALUE INDICATING INACTIVE CONCENTRATION CELLS = ',
     & G15.7/1X,'MINIMUM SATURATED THICKNESS [THKMIN] ',
     & 'ALLOWED =',F8.4,' OF TOTAL CELL THICKNESS')
 1022 FORMAT(1X,'[THKMIN] MUST BE >= 0 AND <= 0.1;  ',
     & ' RESET TO DEFAULT VALUE OF 0.01 [1% OF TOTAL CELL THICKNESS]')

C     READ AND ECHO OUTPUT CONTROL OPTIONS
      READ(INBTN,*) HEADER
      READ(INBTN,*) IFMTCN,IFMTNP,IFMTRF,IFMTDP,SAVUCN
      SAVCBM=.FALSE.
C      WRITE(IOUT,1025)

C      IF(IFMTCN.NE.0) WRITE(IOUT,1030) IFMTCN
C      IF(IFMTCN.EQ.0) WRITE(IOUT,1032)
 1025 FORMAT(//1X,'OUTPUT CONTROL OPTIONS'/1X,22('-'))
 1030 FORMAT(/1X,'PRINT CELL CONCENTRATION USING FORMAT CODE:',I5)
 1032 FORMAT(/1X,'DO NOT PRINT CELL CONCENTRATION')
C
C      IF(IFMTNP.NE.0) WRITE(IOUT,1034) IFMTNP
C      IF(IFMTNP.EQ.0) WRITE(IOUT,1036)
 1034 FORMAT(1X,'PRINT PARTICLE NUMBER IN EACH CELL',
     &          ' USING FORMAT CODE:',I5)
 1036 FORMAT(1X,'DO NOT PRINT PARTICLE NUMBER IN EACH CELL')
C
C      IF(IFMTNP.NE.0) WRITE(IOUT,1038) IFMTRF
C      IF(IFMTNP.EQ.0) WRITE(IOUT,1040)
 1038 FORMAT(1X,'PRINT RETARDATION FACTOR USING FORMAT CODE:',I5)
 1040 FORMAT(1X,'DO NOT PRINT RETARDATION FACTOR')
C
C      IF(IFMTDP.NE.0) WRITE(IOUT,1042) IFMTDP
C      IF(IFMTDP.EQ.0) WRITE(IOUT,1044)
 1042 FORMAT(1X,'PRINT DISPERSION COEFFICIENT USING FORMAT CODE:',I5)
 1044 FORMAT(1X,'DO NOT PRINT DISPERSION COEFFICIENT')

 1046 FORMAT(1X,'SAVE CONCENTRATION IN UNFORMATTED FILE ',
     & '[RT3Dnnn.UCN] ON UNIT ',I3)
 1047 FORMAT(1X,'DO NOT SAVE CONCENTRATION IN UNFORMATTED FILE')

C--READ NUMBER OF TIMES AT WHICH SIMULATION RESULTS SHOULD BE SAVED
C--IN STANDARD OUTPUT FILE OR RECORDED IN UNFORMATTED FILE
      READ(INBTN,*) HEADER
      READ(INBTN,'(I10)') NPRS
      allocate(TIMPRS(NPRS))
      IF(NPRS.GT.0) THEN
        READ(INBTN,*) HEADER
        READ(INBTN,*) (TIMPRS(I),I=1,NPRS)
C
C--MAKE SURE ELEMENTS IN ARRAY [TIMPRS] ARE MONOTONICALLY INCREASING
        DO I=1,NPRS-1
          DO IP1=I+1,NPRS
            IF(TIMPRS(I).GT.TIMPRS(IP1)) THEN
              TEMP=TIMPRS(I)
              TIMPRS(I)=TIMPRS(IP1)
              TIMPRS(IP1)=TEMP
            ENDIF
          ENDDO
        ENDDO
        WRITE(IOUT,*)
        WRITE(IOUT,1055) (TIMPRS(I),I=1,NPRS)
      ENDIF
 1050 FORMAT(/1X,'SIMULATION RESULTS ARE SAVED EVERY ',I3,
     & ' TRANSPORT STEP(S)')
 1052 FORMAT(/1X,'NUMBER OF TIMES AT WHICH SIMULATION RESULTS',
     &' ARE SAVED =',I5)
 1055 FORMAT(1X,'TOTAL ELAPSED TIMES AT WHICH SIMULATION RESULTS ',
     & 'ARE SAVED: ',100(/1X,8G13.5))
C
C     Read number of observation cells
      READ(INBTN,*) HEADER
      READ(INBTN,*) NOBS,NPROBS
      IF(NPROBS.LT.1) NPROBS=1
C      WRITE(IOUT,1056) NOBS
C	  Check to see if the number of cells exceeds the maximum
	  IF(NOBS.GT.MXOBS) THEN
        WRITE(*,1058) MXOBS
        STOP
      ENDIF
C	  If there are observation cells, record the i,j,k indices      
	  IF(NOBS.GT.0) THEN
C        WRITE(IOUT,1062) NOBS,NPROBS
C        WRITE(IOUT,1060)
        allocate(LOCOBS(3,NOBS))
        DO N=1,NOBS
          READ(INBTN,*) (LOCOBS(I,N),I=1,3)
C          WRITE(IOUT,'(I4,4X,3(I5,2X))') N,(LOCOBS(I,N),I=1,3)
        ENDDO
C		Create and open up files for observation cell concentration output        
		DO INDEX=1,NCOMP
		  FLNAME = 'swatmf_out_RT_OBS'//SPECIES(INDEX)
		  OPEN(IOBS+INDEX, FILE=FLNAME, STATUS='UNKNOWN')
		  WRITE(IOBS+INDEX,1080) ((LOCOBS(I,N),I=1,3),N=1,NOBS)
		ENDDO

      ENDIF
 1056 FORMAT(/1X,'NUMBER OF OBSERVATION POINTS =',I5)
 1058 FORMAT(/1X,'ERROR: MAXIMUM NUMBER OF OBSERVATION POINTS IS',I5,
     & /1X,'INCREASE DIMENSION OF [MXOBS] IN THE MAIN PROGRAM')
 1060 FORMAT(1X,'LOCATION OF OBSERVATION POINTS'/1X,30('.')
     &      /1X,'NUMBER  LAYER   ROW   COLUMN')
 1062 FORMAT(1X,'CONCENTRATION AT OBSERVATION POINTS SAVED IN FILE',
     & ' [RT3Dnnn.OBS] ON UNIT ',I3,' EVERY',I4,' TRANSPORT STEPS')
 1063 FORMAT(1X,' STEP   TOTAL TIME',
     & '             LOCATION OF OBSERVATION POINTS (I,J,K)'
     & /1X,17X,16(1X,3I4,1X)/(1X,17X,16(1X,3I4,1X)))
 1080 FORMAT(1X,' STEP   TOTAL TIME',
     & '             LOCATION OF OBSERVATION POINTS (I,J,K)'
     & /400(I5))
C
C--READ AND ECHO LOGICAL FLAG CHKMAS
      READ(INBTN,*) HEADER
      READ(INBTN,*) CHKMAS
      NPRMAS=1
      IF(CHKMAS) THEN
C        WRITE(IOUT,1064) IMAS,NPRMAS
        EXT='mas'
      ELSE
        WRITE(IOUT,1065)
      ENDIF
 1064 FORMAT(/1X,'A ONE-LINE SUMMRY OF MASS BALANCE SAVED IN FILE',
     & ' [RT3Dnnn.MAS] ON UNIT ',I3,' EVERY',I4,' TRANSPORT STEPS')
 1065 FORMAT(/1X,'A ONE-LINE SUMMRY OF MASS BALANCE',
     & ' WILL NOT BE SAVED')
 1066 FORMAT(1X,'     TIME       TOTAL IN      TOTAL OUT      ',
     & 'SOURCES        SINKS       TOTAL MASS         DISCREPANCY(%)')
 1067 FORMAT(20(A15))
 1068 FORMAT(1X,5(4X,'(',A4,')',4X),
     & '  IN AQUIFER   (TOTAL IN-OUT)  (ALTERNATIVE)')

C--ASSIGN SHARED ICBUND ARRAY TO ALL SPECIES
      DO INDEX=2,NCOMP
        DO K=1,NLAY
          DO I=1,NROW
            DO J=1,NCOL
              ICBUND(J,I,K,INDEX) = ICBUND(J,I,K,1)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C--ASSIGN CINACT TO INACTIVE CELLS AND COPY COLD TO CNEW
      DO INDEX=1,NCOMP
        DO K=1,NLAY
          DO I=1,NROW
            DO J=1,NCOL
              IF(ICBUND(J,I,K,INDEX).EQ.0) COLD(J,I,K,INDEX)=CINACT
              CNEW(J,I,K,INDEX)=COLD(J,I,K,INDEX)
              CADV(J,I,K,INDEX)=COLD(J,I,K,INDEX)
C             Do the same for MOLD and MNEW
              IF (TRNOP(7)) THEN !Only if VST package is used   
                MNEW(J,I,K,INDEX)=MOLD(J,I,K,INDEX)
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C--SET PRSITY=0 IN CELLS WHERE ICBUND=0
C--AND ENSURE ICBUND=0 IF POROSITY IS ZERO
      DO K=1,NLAY
        DO I=1,NROW
          DO J=1,NCOL
            IF(ICBUND(J,I,K,1).EQ.0) THEN
              PRSITY(J,I,K)=0
            ELSEIF(PRSITY(J,I,K).EQ.0) THEN
              DO INDEX=1,NCOMP
                ICBUND(J,I,K,INDEX)=0
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO
C
C--CALCULATE COORDINATE ARRAYS XBC, YBC AND ZBC
C--AND MAXIMUN WIDTHS ALONG ROWS, COLUMNS AND LAYERS
      HORIGN=HTOP(1,1)
      XBC(1)=DELR(1)/2.
C     Calculate the upper-most coordinate (Y) of the grid
      TOP = 0
      DO I=1,NROW
        TOP = TOP + DELC(I)
      ENDDO      
      YBC(1)= TOP - (DELC(1)/2)
C     Calculate the Cell X coordinates (starting from upper left of grid)
      DO J=2,NCOL
        XBC(J)=XBC(J-1)+(DELR(J-1)+DELR(J))/2.
      ENDDO
C     Calculate the Cell Y coordaintes (starting from upper left of grid)      
      DO I=2,NROW
        YBC(I) = YBC(I-1)-(DELC(I-1)+DELC(I))/2. 
C        YBC(I)=YBC(I-1)+(DELC(I-1)+DELC(I))/2.
      ENDDO
      XMAX=XBC(NCOL)+DELR(NCOL)/2.
      YMAX=YBC(1)+DELC(NROW)/2.
      ZMAX=0
C     Calculate cell elevations (Z)      
      DO I=1,NROW
        DO J=1,NCOL
          ZBC(J,I,1) = HTOP(J,I) - DZ(J,I,1)/2
          ZZ=DZ(J,I,1)
          IF (NLAY.GT.1) THEN
            DO K=2,NLAY
              ZBC(J,I,K) = ZBC(J,I,K-1)-((DZ(J,I,K-1)+DZ(J,I,K))/2)
            ENDDO
          ENDIF
          ZMAX=MAX(ZMAX,ZZ)
        ENDDO
      ENDDO
C      WRITE(IOUT,1300) XMAX,YMAX,ZMAX
 1300 FORMAT(/1X,'MAXIMUM LENGTH ALONG THE X (J) AXIS =',G15.7,
     &       /1X,'MAXIMUM LENGTH ALONG THE Y (I) AXIS =',G15.7,
     &       /1X,'MAXIMUM LENGTH ALONG THE Z (K) AXIS =',G15.7)

C     INITIALIZE RETARDATION FACTOR ARRAY AND THE MINIMUM
      DO INDEX=1,NCOMP
        DO K=1,NLAY
          DO I=1,NROW
            DO J=1,NCOL
              RETA(J,I,K,INDEX)=1.
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      RFMIN=1.

 1900 FORMAT(30(A15))

      return
      end


C ***********************************************************************************************************
C btn_ts
C THIS subroutine ADVANCES THE TRANSPORT SIMULATION ONE STEP, DETERMINING THE STEPSIZE TO BE USED AND 
C WHETHER PRINTOUT IS REQUIRED FOR NEXT TRANSPORT STEP.
C ************************************************************************************************************
      subroutine btn_ts

      use rt_global
      use VST_ARRAYS
      use mf_rt_link
      use GLOBAL,     only:NCOL,NROW,NLAY,DELR,DELC,mf_interval

      IMPLICIT  NONE
      INTEGER   INDEX,K,I,J,N,S,MOBILE_INDEX,dum
      REAL      DTOLD,VOLUME,EPSILON,xx,SOLID,O2_SAT
      PARAMETER (EPSILON=0.5E-6)

C     SAVE PREVIOUS TRANSPORT STEPSIZE
      DTOLD = DT

C     DETERMINE STEPSIZE FOR NEXT TRANSPORT STEP
      IF(IMPSOL.EQ.0) THEN
        IF(MIXELM.GT.0) THEN
          IF(TRNOP(1)) DT=MIN(DT,DTRACK*PERCEL*RFMIN)
          IF(TRNOP(2)) DT=MIN(DT,DTDISP*RFMIN)
          IF(TRNOP(3)) DT=MIN(DT,DTSSM*RFMIN)
c          IF(TRNOP(4)) DT=MIN(DT,DTRCT) !no rxn constraint for RT3D
        ELSE
          DT=0.
          IF(TRNOP(1)) DT=DT+1./(DTRACK2*PERCEL*RFMIN)
          IF(TRNOP(2)) DT=DT+1./(DTDISP*RFMIN)
          IF(TRNOP(3)) DT=DT+1./(DTSSM*RFMIN)
c          IF(TRNOP(4)) DT=DT+1./DTRCT !no rxn constraint for RT3D
          DT=1./DT
        ENDIF
        IF(DT0.GT.0.AND.DT0.LT.DT) DT=DT0
      ELSEIF(IMPSOL.EQ.1) THEN
        IF(MIXELM.EQ.0) THEN
          IF(TRNOP(1)) DT=MIN(DT,DTRACK*PERCEL*RFMIN)
          IF(DT0.GT.0) DT=DT0
          IF(TTSMULT.GT.1) DT=DT*TTSMULT**(NTRANS-1)
          IF(TTSMAX.GT.0.AND.DT.GT.TTSMAX) DT=TTSMAX
        ELSE
          IF(TRNOP(1).AND.MIXELM.GT.0) THEN
            DT=MIN(DT,DTRACK*PERCEL*RFMIN)
          ELSEIF(TRNOP(1).AND.MIXELM.LT.0) THEN
            DT=MIN(DT,DTRACK2*PERCEL*RFMIN)
          ENDIF
          IF(DT0.GT.0.AND.DT0.LT.DT) DT=DT0
          IF(TTSMAX.GT.0.AND.DT.GT.TTSMAX) DT=TTSMAX
        ENDIF
      ENDIF

C     IF DT0 NEGATIVE, USE |DT0| AS DEFAULT TRANSPORT STEPSIZE
      IF(DT0.LT.0) DT=ABS(DT0)

      !DT = rt_delt / mf_interval !rtb MODFLOW

C     UPDATE TOTAL ELASPED TIME
      TIME1=TIME2
      TIME2=TIME1+DT

C     DETERMINE IF PRINTOUT OF SIMULATION RESULTS IS NEEDED FOR NEXT STEP
      PRTOUT=.FALSE.
      IF(NTRANS.EQ.MXSTRN) THEN
        PRTOUT=.TRUE.
      ENDIF
      IF(NPRS.LT.0) THEN
        IF(MOD(NTRANS,-NPRS).EQ.0) PRTOUT=.TRUE.
      ENDIF

C     IF TOTAL ELAPSED TIME AT NEXT STEP EXCEEDS TIME AT WHICH
C     PRINTOUT IS REQUESTED, CUT DOWN STEPSIZE
      IF(NPRS.GT.0.AND.NPS.LE.NPRS) THEN
        IF(TIME2.GE.TIMPRS(NPS)) THEN
          IF(TIME2.GT.TIMPRS(NPS)) THEN
            TIME2=TIMPRS(NPS)
            DT=TIME2-TIME1
          ENDIF
          PRTOUT=.TRUE.
          NPS=NPS+1
        ENDIF
      ENDIF

C     IF TOTAL ELAPSED TIME AT NEXT TRANSPORT STEP EXCEEDS
C     THE LIMIT OF CURRENT TIME STEP, CUT DOWN STEPSIZE
      IF(TIME2.GE.T2) THEN
        TIME2=T2
        DT=TIME2-TIME1
      ENDIF
      UPDLHS=.TRUE.
      IF(ABS(DT-DTOLD).LT.ABS(DT+DTOLD)*EPSILON
     & .AND. NTRANS.GT.1) UPDLHS=.FALSE.

C     PRINT OUT AN IDENTIFYING MESSAGE
      !WRITE(*,70) NTRANS,DT,TIME2,SIMLEN
      !WRITE(IOUT,70) NTRANS,DT,TIME2,SIMLEN
   70 FORMAT(1X,'Transp. Step:',I4,2X,'Step Size:',G9.3,
     & ' Elapsed Time: ',F10.2,' of',F10.2)

C     COPY ARRAY [CNEW] TO [COLD]
      DO INDEX=1,NCOMP
        DO K=1,NLAY
          DO I=1,NROW
            DO J=1,NCOL
              
              if(isnan(cnew(j,i,k,index))) then
                dum = 10
              endif
              
              IF(ICBUND(J,I,K,INDEX).EQ.0) THEN
                CNEW(J,I,K,INDEX) = CINACT 
              ELSE
                COLD(J,I,K,INDEX)=CNEW(J,I,K,INDEX)
              ENDIF
             ENDDO
          ENDDO
        ENDDO
      ENDDO
      
C     Copy array MNEW to MOLD
      IF (TRNOP(7)) THEN !Only if VST package is used
        DO INDEX=1,NCOMP
          DO J=1,NCOL
            DO I=1,NROW
              DO K=1,NLAY
                MOLD(J,I,K,INDEX) = MNEW(J,I,K,INDEX)
                MTEMP(J,I,K,INDEX) = MNEW(J,I,K,INDEX)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF

C     CLEAR RMASIO ARRAY FOR ACCUMULATING MASS IN/OUT AT NEXT TRANSPORT STEP
      DO INDEX=1,NCOMP
        DO I=1,28
          RMASIO(I,1,INDEX)=0.
          RMASIO(I,2,INDEX)=0.
        ENDDO
      ENDDO
C     Clear LEACH array for next transport step
      LEACH = 0.0

C     If first transport step, calculate the initial mass (dissolved, sorbed, solid)
C     in the aquifer system.
      IF(rt_kper*rt_kstp*NTRANS.GT.1) GOTO 9999
      MOBILE_INDEX = 1
      DO INDEX=1,NCOMP

        TMASS(1,INDEX) = 0.0 !Dissolved phase
        TMASS(2,INDEX) = 0.0 !Sorbed phase
        TMASS(3,INDEX) = 0.0 !Solid phase

C       Loop through each cell in the grid        
        DO K=1,NLAY
          DO I=1,NROW
            DO J=1,NCOL
              VOLUME = DELR(J)*DELC(I)*DZ(J,I,K)
C             Only proceed if cell is active
              IF(ICBUND(J,I,K,INDEX).GT.0) THEN
                IF (TRNOP(7)) THEN !VST package is used
                  IF(PHASE(INDEX).EQ.1) THEN !Mass in dissolved phase
                    TMASS(1,INDEX) = TMASS(1,INDEX) + MOLD(J,I,K,INDEX)
C                   Calculate initial sorbed mass if the species has sorption
C                   parameters specified (in RCT file)             
                    IF(TRNOP(4)) THEN !If the RCT package is used
                      IF (ISOTHM.GT.0) THEN
                        TMASS(2,INDEX) = TMASS(2,INDEX) + VOLUME*
     &                        COLD(J,I,K,INDEX)*(RHOB(J,I,K)/MC(J,I,K))
                      ENDIF
                    ENDIF
                  ELSE !Mass in solid phase
                    TMASS(3,INDEX) = TMASS(3,INDEX) + MOLD(J,I,K,INDEX)
                  ENDIF
C                 Store 1,2,3 into 4,5,6 (so that the current total mass
C                 in the aquifer = the initial total mass)
                  TMASS(4,INDEX) = TMASS(1,INDEX)
                  TMASS(5,INDEX) = TMASS(2,INDEX)
                  TMASS(6,INDEX) = TMASS(3,INDEX)                  
                ELSE
                  VOLUME=DELR(J)*DELC(I)*DF(J,I,K)*PRSITY(J,I,K)
C                 SUM UP MASS IN EACH CELL (DISSOLVED PHASE)
                  TMASS(1,INDEX)=TMASS(1,INDEX)+VOLUME*COLD(J,I,K,INDEX)
C                 SUM UP MASS IN EACH CELL (SOLID PHASE)
                  IF(ISOTHM.GT.0)
     &            TMASS(2,INDEX)=TMASS(2,INDEX)+VOLUME*COLD(J,I,K,INDEX)
     &             *RHOB(J,I,K)/PRSITY(J,I,K)
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO
C     Update the index of the mobile species (from 1 to MCOMP)
      IF(PHASE(INDEX).EQ.1) THEN
        MOBILE_INDEX = MOBILE_INDEX + 1
      ENDIF
      ENDDO
      
C     Update Retardation Factor if nonlinear sorption is simulated or if
C     unsaturated conditions are simulated (thus, moisture content changes with
C     each time step, and need to recalculate R).
      IF(TRNOP(7) .OR. ISOTHM.EQ.2 .OR. ISOTHM.EQ.3) THEN
        index = 1        
        DO ICOMP=1,NCOMP
          IF (PHASE(ICOMP).EQ.1) THEN !Specie is mobile
            call reta_calc(index)          
            index = index + 1
          ENDIF
        ENDDO
      ENDIF           
      
 9999 RETURN
      END


C ***********************************************************************************************************
C btn_solve
C THIS subroutine UPDATES CELL CONCENTRATION AND MASS IN/OUT ACCUMULATING ARRAY TO PREPARE FOR SIMULATION 
C AT NEXT STEP.
C ***********************************************************************************************************
      subroutine btn_solve

      use rt_global
      use mf_rt_link
      use GLOBAL,     only:NCOL,NROW,NLAY

      IMPLICIT  NONE
      INTEGER   J,I,K,DUM


C     COPY CNEW TO CWGT
      DO K=1,NLAY
        DO I=1,NROW
          DO J=1,NCOL
            IF(ICBUND(J,I,K,ICOMP).EQ.0) CYCLE
            CWGT(J,I,K,ICOMP)=CNEW(J,I,K,ICOMP)
          ENDDO
        ENDDO
      ENDDO

C     CLEAR RMASIO ARRAY FOR ACCUMULATING MASS IN/OUT AT NEXT TRANSPORT STEP
      DO I=1,28
        RMASIO(I,1,ICOMP)=0.
        RMASIO(I,2,ICOMP)=0.
      ENDDO

      RETURN
      END


C ***********************************************************************************************************
C btn_imp_update
C THIS SUBROUTINE UPDATES CELL CONCENTRATION AND MASS IN/OUT ACCUMULATING ARRAY TO PREPARE FOR SIMULATION 
C AT NEXT STEP, IF THE IMPLICIT SCHEME (GCG SOLVER) IS USED
C ***********************************************************************************************************
      subroutine btn_imp_update(CNEWN,ICBUNDN)

      use rt_global
      use mf_rt_link
      use GLOBAL,     only:NCOL,NROW,NLAY

      IMPLICIT  NONE
      INTEGER   J,I,K,N,ICBUNDN
      REAL      CNEWN
      DIMENSION CNEWN(NODES),ICBUNDN(NODES)

C     COPY ICBUND TO ICBUNDN (array for the implicit scheme)
      DO K=1,NLAY
        DO I=1,NROW
          DO J=1,NCOL
            N=(K-1)*NCOL*NROW+(I-1)*NCOL+J
            ICBUNDN(N) = ICBUND(J,I,K,ICOMP)
          ENDDO
        ENDDO
      ENDDO

C     COPY CNEW TO CNEWN
      DO K=1,NLAY
        DO I=1,NROW
          DO J=1,NCOL
            IF(ICBUND(J,I,K,ICOMP).EQ.0) CYCLE
            N=(K-1)*NCOL*NROW+(I-1)*NCOL+J
            CNEWN(N) = CNEW(J,I,K,ICOMP)
          ENDDO
        ENDDO
      ENDDO

      RETURN
      END


C ***********************************************************************************************************
C btn_budget
C THIS subroutine SUMMARIZES VOLUMETRIC MASS BUDGETS AND CALCULATES
C MASS BALANCE DISCREPANCY SINCE THE BEGINNING OF THE SIMULATION.
C ***********************************************************************************************************
      subroutine btn_budget

      use rt_global
      USE VST_ARRAYS
      use mf_rt_link
      use GLOBAL,     only:NCOL,NROW,NLAY,DELR,DELC

      IMPLICIT  NONE
      INTEGER   K,I,J,IQ,INDEX,N,DUM
      REAL      DMSTRG,TM1,TM2,FLOWRATE,VOLUME,mass_old,mass_new,
     &          DIFF1,LEFT,RIGHT,OLDMASS,SOURCE_TS,SINK_TS,DM_TS
     

C     MASS BUDGET CALCULATIONS

      !Mass balance equation: NEW MASS = INITIAL MASS + (MASS GAINED - MASS LOST)

      !The new mass is calculated using the total change in species mass since the
      !beginning of the simulation.

      !The Mass Gained - Mass Lost is calculated using the sources and sinks of the
      !species since the beginning of the simulation.

      !If the numerical algorithms in the code are accurate, then the equation should
      !balance out within 1 percent (approximately).


C     1. Calculate the INITIAL MASS in the aquifer (from beginning of simulation) -----------------
C        This is a combination of dissolved + sorbed + solid phases
      INITMASS = TMASS(1,ICOMP) + TMASS(2,ICOMP) + TMASS(3,ICOMP)

C     2. Calculate the OLD MASS in the aquifer (from the last time step) --------------------------
      OLDMASS = TMASS(4,ICOMP) + TMASS(5,ICOMP) + TMASS(6,ICOMP) 
      
C     3. Calculate the NEW MASS in the aquifer (from current time step) ---------------------------
C     First, calculate the change in mass from the previous transport step.
      DO K=1,NLAY
        DO I=1,NROW
          DO J=1,NCOL

            !Only proceed if the cell is active
            IF(ICBUND(J,I,K,ICOMP).GT.0) THEN

              !Calculate the change in mass storage since the last transport step
              IF(TRNOP(7)) THEN
                volume = DELR(J)*DELC(I)*DZ(J,I,K) 
                mass_old = cold(j,i,k,icomp)*mcold(j,i,k)*volume
                mass_new = cnew(j,i,k,icomp)*mc(j,i,k)*volume
                DMSTRG = mass_new - mass_old 
              ELSE
                DMSTRG=(CNEW(J,I,K,ICOMP)-COLD(J,I,K,ICOMP))
     &                     *DELR(J)*DELC(I)*DH(J,I,K)*PRSITY(J,I,K)
              ENDIF

              !Store the mass storage for the dissolved, sorbed, and solid phases
              IF(PHASE(ICOMP).EQ.1) THEN !Dissolved and sorbed phases
                IF(DMSTRG.GE.0) THEN !Increase in storage
                  RMASIO(19,1,ICOMP)=RMASIO(19,1,ICOMP)+DMSTRG
                  RMASIO(20,1,ICOMP)=RMASIO(20,1,ICOMP)+
     &                             ((RETA(J,I,K,ICOMP)-1)*DMSTRG)
                ELSE !Decrease in storage
                  RMASIO(19,2,ICOMP)=RMASIO(19,2,ICOMP)+DMSTRG
                  RMASIO(20,2,ICOMP)=RMASIO(20,2,ICOMP)+
     &                             ((RETA(J,I,K,ICOMP)-1)*DMSTRG)
                ENDIF
              ELSEIF(PHASE(ICOMP).EQ.2) THEN !Sorbed species
                IF(DMSTRG.GE.0) THEN
                  RMASIO(20,1,ICOMP) = RMASIO(20,1,ICOMP) + DMSTRG
                ELSE
                  RMASIO(20,2,ICOMP) = RMASIO(20,2,ICOMP) + DMSTRG
                ENDIF  
              ELSE !Solid phase
                IF(DMSTRG.GE.0) THEN !Increase in storage
                  RMASIO(21,1,ICOMP)=RMASIO(21,1,ICOMP)+DMSTRG
                ELSE !Decrease in storage
                  RMASIO(21,2,ICOMP)=RMASIO(21,2,ICOMP)+DMSTRG
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      !Calculate the cumulative change in mass from the beginning of the 
      !simulation (19=dissolved, 20=sorbed, 21=solid)     
      DO IQ=19,21
        TMASIO(IQ,1,ICOMP) = TMASIO(IQ,1,ICOMP) + RMASIO(IQ,1,ICOMP)
        TMASIO(IQ,2,ICOMP) = TMASIO(IQ,2,ICOMP) + RMASIO(IQ,2,ICOMP)
      ENDDO
      
      !Calculate the new mass in the aquifer
      !NEW MASS = INITIAL MASS + (INCREASE IN MASS - DECREASE IN MASS)
      !dissolved phase
      TMASS(4,ICOMP) = TMASS(1,ICOMP) + 
     &                  (TMASIO(19,1,ICOMP) + TMASIO(19,2,ICOMP))
      !sorbed phase
      TMASS(5,ICOMP) = TMASS(2,ICOMP) + 
     &                  (TMASIO(20,1,ICOMP) + TMASIO(20,2,ICOMP))
      !solid phase
      TMASS(6,ICOMP) = TMASS(3,ICOMP) + 
     &                  (TMASIO(21,1,ICOMP) + TMASIO(21,2,ICOMP))
      NEWMASS = TMASS(4,ICOMP) + TMASS(5,ICOMP) + TMASS(6,ICOMP)


C     4. Calculate the change in mass (i) from the beginning of the simulation, and
C     (ii) during the last time step
      DM_TOTAL = NEWMASS - INITMASS !Total change in mass
      DM_TS = NEWMASS - OLDMASS !Change in mass from last time step


C     5. Calculate the Sources and Sinks of mass since the last time step
C     SOURCE_TS = Source of mass during the previous time step
C     SINK_TS = Sink of mass during the previous time step      
      SOURCE_TS = 0.0
      SINK_TS = 0.0
      DO IQ=1,18
        SOURCE_TS = SOURCE_TS + RMASIO(IQ,1,ICOMP)
        SINK_TS = SINK_TS + RMASIO(IQ,2,ICOMP)
      ENDDO
      
            
C     6. Calculate the Sources and Sinks of mass since the beginning of the simulation
      !First, tally the cumulative mass for each source/sink/process
      DO IQ=1,18       
        TMASIO(IQ,1,ICOMP) = TMASIO(IQ,1,ICOMP) + RMASIO(IQ,1,ICOMP)
        TMASIO(IQ,2,ICOMP) = TMASIO(IQ,2,ICOMP) + RMASIO(IQ,2,ICOMP)
      ENDDO    
      !Second, add all sources together and sinks together      
      SOURCE_TOT = 0.0
      SINK_TOT = 0.0
      DO IQ=1,18
        SOURCE_TOT = SOURCE_TOT + TMASIO(IQ,1,ICOMP)
        SINK_TOT = SINK_TOT + TMASIO(IQ,2,ICOMP)
      ENDDO


C     7. Calculate the mass budget discrepancy for the entire simulation
      !MASS GAINED + INITIAL MASS = MASS LOST + NEW MASS
      LEFT = ABS(SOURCE_TOT) + INITMASS
      RIGHT = ABS(SINK_TOT) + NEWMASS
      ERROR2(ICOMP) = 0.0
      IF((LEFT+RIGHT).NE.0) THEN
        ERROR2(ICOMP) = (RIGHT-LEFT)/(0.5*(RIGHT+LEFT))*100.0
      ENDIF
      

C     8. Calculate the mass budget discrepancy for the current time step
      LEFT = ABS(SOURCE_TS) + OLDMASS
      RIGHT = ABS(SINK_TS) + NEWMASS
      ERROR2_TS(ICOMP) = 0.0
      IF((LEFT+RIGHT).NE.0) THEN
        ERROR2_TS(ICOMP) = (RIGHT-LEFT)/(0.5*(RIGHT+LEFT))*100.0
      ENDIF
      
C     9. Write a one-line mass budget summary for the specie
      DM_SS = ABS(SOURCE_TOT) - ABS(SINK_TOT)
 1012 FORMAT(F15.1,G15.5,G15.5,G15.5,G15.5,G15.5,G15.5,G15.5,G15.5,
     &       F15.4,F15.4)       



C     LEACHING CALCULATIONS
C     Calculate leaching rates for the current species, for all grid cells.
C     Only do for mobile species
      INDEX = 1
      IF (PHASE(ICOMP).EQ.1) THEN !Only if the species is mobile
      DO I=1,NROW
        DO J=1,NCOL
          DO K=1,NLAY
C           Calculate leaching for current time step. Leaching is calculated using
C           cell concentration and flow rate from above cell, and the time step
C           LEACH = [M/L3]*[L3/T]*T = [M]
            IF (K.GT.1) THEN !Do not calculate for first layer
              IF (TRNOP(7)) THEN !Variably-Saturated Flow
                FLOWRATE = QZVAR(J,I,K-1)*(DELR(J)*DELC(I))
              ELSE
                FLOWRATE = QZ(J,I,K-1)*(DELR(J)*DELC(I))
              ENDIF
              LEACH(J,I,K,ICOMP) = CNEW(J,I,K-1,ICOMP)*FLOWRATE*DT
            ENDIF
C           Calculate cumulative leaching
            TLEACH(J,I,K,ICOMP)=TLEACH(J,I,K,ICOMP)+LEACH(J,I,K,ICOMP)
          ENDDO
        ENDDO
      ENDDO
      ENDIF
 
      return
      end


C ***********************************************************************************************************
C btn_output
C THIS subroutine SAVES SIMULATION RESULTS IN THE STANDARD OUTPUT FILE AND VARIOUS OPTIONAL OUTPUT FILES, 
C ACCORDING TO THE OUTPUT CONTROL OPTIONS SPECIFIED IN THE BASIC TRANSPORT INPUT FILE.
C ***********************************************************************************************************
      subroutine btn_output

      use rt_global
      use VST_ARRAYS
      use mf_rt_link
      use GLOBAL,     only:NCOL,NROW,NLAY,DELR,DELC

      implicit none
      integer   K,I,J,N
      real      VOLUME,MASS,MASS_SAT,MASS_UNSAT,SUM_SAT,SUM_UNSAT,
     &          SAT_PORTION,UNSAT_PORTION
      CHARACTER TEXT*16

      IF (WRITE_OUT(ICOMP).EQ.1) THEN
C     Write out concentration and leaching at observation cells -----------------------------------
      IF(NOBS.GT.0.AND.(MOD(NTRANS-1,NPROBS).EQ.0.OR.PRTOUT)) THEN
C       Concentrations at observation cells          
        WRITE(IOBS+ICOMP,1081) TIME2,
     &  (CNEW(LOCOBS(2,N),LOCOBS(1,N),LOCOBS(3,N),ICOMP),N=1,NOBS)
      ENDIF
      ENDIF


C     Write out concentrations and leaching -------------------------------------------------------
C     Either ASCII or BINARY, depending on input flag in BTN file      
      IF (WRITE_OUT(ICOMP).EQ.1) THEN
      IF (PRTOUT) THEN

C       Layer by Layer (XY plane)
	  DO K=1,NLAY
C		  For the current layer, write out information
		  IF (WRITEFLAG.EQ.0 .OR. WRITEFLAG.EQ.2) THEN
		    WRITE(ICONC+ICOMP,1022) TIME2,'LAYER',K
		  ENDIF
		  IF (WRITEFLAG.EQ.1 .OR. WRITEFLAG.EQ.2) THEN
		    WRITE(ICBIN+ICOMP) TIME2,'LAYER',K
		  ENDIF
C         BINARY output
          IF (WRITEFLAG.EQ.1 .OR. WRITEFLAG.EQ.2) THEN	  
		    WRITE(ICBIN+ICOMP) ((CNEW(J,I,K,ICOMP),J=1,NCOL),I=1,NROW)
		  ENDIF
C         ASCII output
		  IF (WRITEFLAG.EQ.0 .OR. WRITEFLAG.EQ.2) THEN  
		  DO I=1, NROW
		    DO J=1, NCOL
C			  For current cell, write out species concentration 			    
			  if(trnop(7)) then
              WRITE(ICONC+ICOMP,1021) XBC(J),YBC(I),CNEW(J,I,K,ICOMP),
     &                                MC(J,I,K) !ASCII
              else
              WRITE(ICONC+ICOMP,1021) XBC(J),YBC(I),CNEW(J,I,K,ICOMP)
              endif
			ENDDO
		  ENDDO
		  ENDIF
	  ENDDO

      ENDIF
      ENDIF

 1016 FORMAT(F12.4,F12.4,F12.4,F12.4,F12.4)
 1017 FORMAT(F12.4,F12.4,F12.4,F12.4)
 1018 FORMAT(F10.4,A15,A15,I4)
 1019 FORMAT(I4,I4,I4,F15.4,A15,A15,I4)
 1021 FORMAT(F12.4,F12.4,E18.8,F12.4)
 1022 FORMAT(F15.4,A15,I8)
 1081 format(1000(F15.8))

C     Print out total mass in unsaturated and saturated zones, for the current species
      IF (TRNOP(7)) THEN !VST package is used
      IF (TIME2.EQ.INT(TIME2)) THEN
      MASS_UNSAT = 0
      MASS_SAT = 0
      DO I=1,NROW
	  DO J=1,NCOL
          DO K=1,NLAY
C		    For current cell, calculate the mass
C           Mass = Concentration * Water Content * Bulk Volume
            VOLUME = DELR(J)*DELC(I)*DZ(J,I,K)
            MASS = CNEW(J,I,K,ICOMP)*MC(J,I,K)*VOLUME 		  
C           Determine if the cell is:
C               1. fully unsaturated
C               2. a water table cell
C               3. fully saturated
C           the unsaturated and saturated zones.
            IF(SATNEW(J,I,K).LT.1) THEN !Partially-saturated
              IF (DH(J,I,K).LT.DZ(J,I,K) .AND. !Water Table cell
     &            DH(J,I,K).GT.0.0) THEN
                SAT_PORTION = DH(J,I,K) / DZ(J,I,K)
                UNSAT_PORTION = 1 - SAT_PORTION
              ELSE
                SAT_PORTION = 0
                UNSAT_PORTION = 1
              ENDIF
            ELSE !Fully-saturated
              SAT_PORTION = 1
              UNSAT_PORTION = 0
            ENDIF
C           Partition the mass between the unsaturated and saturated zones
C           Report in kg (to avoid large numbers)            
            MASS_UNSAT = MASS_UNSAT + ((MASS*UNSAT_PORTION)/1000)
            MASS_SAT = MASS_SAT + ((MASS*SAT_PORTION)/1000)
		  ENDDO
		ENDDO
	ENDDO  
	ENDIF
	ENDIF     

C     Write out mass budget terms and calculations ------------------------------------------------
      IF(.NOT.PRTOUT) GOTO 9999

C     Print out the species index and name
      WRITE(IOUT,*)
      WRITE(IOUT,1020) ICOMP,SPECIES(ICOMP)
 1020 FORMAT('FOR SPECIES NO.:',I5,2X,A15)
          
C     Print out the total mass in the aquifer for the initial and 
C     current time steps. Report the cumulative change in mass.
      WRITE(IOUT,1100)
      WRITE(IOUT,1101)
      WRITE(IOUT,1102) TMASS(1,ICOMP),TMASS(2,ICOMP),TMASS(3,ICOMP),
     &                 INITMASS
      WRITE(IOUT,1103) TMASS(4,ICOMP),TMASS(5,ICOMP),TMASS(6,ICOMP),
     &                 NEWMASS
      WRITE(IOUT,1104) DM_TOTAL
      
C     Print out the total mass in and mass out for each BC/process      
      WRITE(IOUT,1105)
      WRITE(IOUT,1106)
      WRITE(IOUT,1108) TMASIO(6,1,ICOMP),TMASIO(6,2,ICOMP)
      IF(TRNOP(4)) WRITE(IOUT,1116) TMASIO(9,1,ICOMP),TMASIO(9,2,ICOMP)
C     Include new mass balance information (for IRG package)     
      IF (PHASE(ICOMP).EQ.1) THEN !Dissolved species
        WRITE(IOUT,1109) TMASIO(1,1,ICOMP),TMASIO(1,2,ICOMP)
        IF(FWEL) WRITE(IOUT,1110) TMASIO(2,1,ICOMP),TMASIO(2,2,ICOMP)
        IF(FDRN) WRITE(IOUT,1111) TMASIO(3,1,ICOMP),TMASIO(3,2,ICOMP)
        IF(FRIV) WRITE(IOUT,1112) TMASIO(4,1,ICOMP),TMASIO(4,2,ICOMP)
        IF(FGHB) WRITE(IOUT,1113) TMASIO(5,1,ICOMP),TMASIO(5,2,ICOMP)
        IF(FRCH) WRITE(IOUT,1114) TMASIO(7,1,ICOMP),TMASIO(7,2,ICOMP)
        IF(FEVT) WRITE(IOUT,1115) TMASIO(8,1,ICOMP),TMASIO(8,2,ICOMP)
      ELSE !Solid species
 
      ENDIF
      WRITE(IOUT,1107)
      WRITE(IOUT,1124) SOURCE_TOT,SINK_TOT
      WRITE(IOUT,1125) DM_SS
      WRITE(IOUT,*)
      
C     Write out the mass balance error (comparing DMASS to DMASS_SS)      
      WRITE(IOUT,1126) ERROR2(ICOMP)
      WRITE(IOUT,*)

C     Formats for DMASS1
 1100 FORMAT(/15X,'TOTAL CHANGE IN MASS USING CONCENTRATION (DMASS1)')
 1101 FORMAT(30X,'TOTAL MASS =   DISSOLVED   +   SORBED   +   SOLID')
 1102 FORMAT(30X,'INITIAL =',2X,G15.8,'+',G15.8,'+',G15.8,'=',G15.8)
 1103 FORMAT(30X,'CURRENT =',2X,G15.8,'+',G15.8,'+',G15.8,'=',G15.8)
 1104 FORMAT(30X,' DMASS1 =',2X,G15.8)
 
C     Formats for DMASS2 
 1105 FORMAT(/15X,'TOTAL CHANGE IN MASS USING SOURCES/SINKS (DMASS2)')
 1106 FORMAT(63X,'MASS IN',12X,'MASS OUT')
 1107 FORMAT(60X,30('-'))
 1108 FORMAT(30X,'     CONSTANT CONCENTRATION: ',G15.7,5X,G15.7)
 1109 FORMAT(30X,'              CONSTANT HEAD: ',G15.7,5X,G15.7)
 1110 FORMAT(30X,'                      WELLS: ',G15.7,5X,G15.7)
 1111 FORMAT(30X,'                     DRAINS: ',G15.7,5X,G15.7)
 1112 FORMAT(30X,'                     RIVERS: ',G15.7,5X,G15.7)
 1113 FORMAT(30X,'    HEAD-DEPENDENT BOUNDARY: ',G15.7,5X,G15.7)
 1114 FORMAT(30X,'                   RECHARGE: ',G15.7,5X,G15.7)
 1115 FORMAT(30X,'         EVAPOTRANSPIRATION: ',G15.7,5X,G15.7)
 1116 FORMAT(30X,'         CHEMICAL REACTIONS: ',G15.7,5X,G15.7)
 1117 FORMAT(30X,'           CANAL IRRIGATION: ',G15.7,5X,G15.7)
 1118 FORMAT(30X,'         AQUIFER IRRIGATION: ',G15.7,5X,G15.7)
 1119 FORMAT(30X,'                 FERTILIZER: ',G15.7,5X,G15.7)
 1120 FORMAT(30X,'                CROP UPTAKE: ',G15.7,5X,G15.7)
 1121 FORMAT(30X,'                     MANURE: ',G15.7,5X,G15.7)
 1122 FORMAT(30X,'                  ROOT MASS: ',G15.7,5X,G15.7)
 1123 FORMAT(30X,'                STOVER MASS: ',G15.7,5X,G15.7)
 1124 FORMAT(52X,'TOTAL: ',G15.7,5X,G15.7)
 1125 FORMAT(30X,'DMASS2 = ',G15.7)
C     Format for mass balance error output 
 1126 FORMAT(15X,'% DISCREPANCY (DMASS2 TO DMASS1) = ',F8.4)

 9999 RETURN
      END
      

C ***********************************************************************************************************
C BTN3FM
C This subroutine initializes all matrices for the implicit scheme
C ***********************************************************************************************************
      subroutine btn_imp(A,RHS)
      
      use rt_global
      USE VST_ARRAYS
      use mf_rt_link
      use GLOBAL,     only:NCOL,NROW,NLAY,DELR,DELC

      IMPLICIT  NONE
      INTEGER   J,I,K,N,NSIZE,ICBUND2,L
      REAL      TEMP,COLD2,CADV2,RETA2,PRSITY2,DH2,SATOLD2,DZ2,MC2,
     &          A,RHS
      DIMENSION ICBUND2(NODES),COLD2(NODES),
     &          CADV2(NODES),RETA2(NODES),PRSITY2(NODES),DH2(NODES),
     &          SATOLD2(NODES),DZ2(NODES),MC2(NODES),
     &          A(19*NODES),RHS(NODES)
      COMMON   /GCGIDX/L(19)

C     COPY NECESSARY ARRAYS TO USE IN THIS SUBROUTINE
      DO K=1,NLAY
        DO I=1,NROW
          DO J=1,NCOL
            N=(K-1)*NCOL*NROW + (I-1)*NCOL + J
            ICBUND2(N) = ICBUND(J,I,K,ICOMP)
            COLD2(N) = COLD(J,I,K,ICOMP)
            CADV2(N) = CADV(J,I,K,ICOMP)
            RETA2(N) = RETA(J,I,K,ICOMP)
            PRSITY2(N) = PRSITY(J,I,K)
            DH2(N) = DH(J,I,K)
            DZ2(N) = DZ(J,I,K)
            if(trnop(7)) then
              MC2(N) = MC(J,I,K)
              SATOLD2(N) = SATOLD(J,I,K)
            endif
          ENDDO
        ENDDO
      ENDDO 

C     GET RIGHT-HAND-SIDE ARRAY [RHS]
      DO K=1,NLAY
        DO I=1,NROW
          DO J=1,NCOL
            N=(K-1)*NCOL*NROW + (I-1)*NCOL + J
            IF(MIXELM.EQ.0) THEN
              TEMP=COLD2(N)
            ELSE
              TEMP=CADV2(N)
            ENDIF
            IF(ICBUND2(N).LE.0) THEN
              RHS(N)=-TEMP
            ELSE
C             UNITS of RHS: M/T     
              IF (TRNOP(7)) THEN
                RHS(N)=-TEMP*RETA2(N)/DT*PRSITY2(N)
     &                *SATOLD2(N)*DELR(J)*DELC(I)*DZ2(N)
              ELSE
                RHS(N)=-TEMP*RETA2(N)/DT*PRSITY2(N)
     &               *DELR(J)*DELC(I)*DH2(N)
              ENDIF  
            ENDIF
          ENDDO
        ENDDO
      ENDDO

C     RETURN IF COEFF MATRIX [A] NOT TO BE UPDATED
      IF(.NOT.UPDLHS) GOTO 999

C     RESET COEFF MATRIX [A]
      IF(NCRS.EQ.1) THEN
         NSIZE=19*NODES
      ELSE
         NSIZE=7*NODES
      ENDIF
      DO I=1,NSIZE
        A(I)=0.
      ENDDO

C     LOOP THROUGH ALL CELLS AND RESET A
      N=0
      DO K=1,NLAY
        DO I=1,NROW
          DO J=1,NCOL
            N=N+1
C           IF INACTIVE OR CONSTANT CELL
            IF(ICBUND2(N).LE.0) THEN
               A(N)=-1.
            ELSE
C              Units are L3(fluid)/T
               IF (TRNOP(7)) THEN !If VST package is used
                 A(N)=-RETA2(N)/DT*MC2(N)*DELR(J)*DELC(I)*DZ2(N)
               ELSE  
                 A(N)=-RETA2(N)/DT*PRSITY2(N)
     &              *DELR(J)*DELC(I)*DH2(N)
               ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO

C     CALCULATE MATRIX INDICES FOR THE GCG SOLVER
      NRC = NROW*NCOL
      L(1) =  0
      L(2) = -NRC
      L(3) =  NRC
      L(4) = -NCOL
      L(5) =  NCOL
      L(6) = -1
      L(7) =  1
      L(8) = -NCOL-NRC
      L(9) = -1-NRC
      L(10)= 1-NRC
      L(11)= NCOL-NRC
      L(12)=-NCOL+NRC
      L(13)=-1+NRC
      L(14)= 1+NRC
      L(15)= NCOL+NRC
      L(16)=-1-NCOL
      L(17)= 1-NCOL
      L(18)=-1+NCOL
      L(19)= 1+NCOL

  999 RETURN
      END


C ***********************************************************************************************************
C BTNRES
C Prints out final species concentrations and final cell moisture content (for
C re-starting the simulation)
C ***********************************************************************************************************
      subroutine btn_restart
      
      use rt_global
      USE VST_ARRAYS
      use mf_rt_link
      use GLOBAL,     only:NCOL,NROW,NLAY

      IMPLICIT  NONE
      INTEGER   I,J,K,N

C     Print out final species concentrations      
      DO N=1,NCOMP
        WRITE(OUTRES,2999) 2,0.0,SPECIES(N)
        DO K=1,NLAY
          WRITE(OUTRES,*) '#LAYER',K
          DO I=1,NROW
            WRITE(OUTRES,3000) (CNEW(J,I,K,N),J=1,NCOL)
          ENDDO
        ENDDO
      ENDDO
      
C     If VST package is used, then print out final moisture contents      
      IF (TRNOP(7)) THEN
        WRITE(OUTRES,*) 'INITIAL CONDITIONS OF SOIL MOISTURE CONTENT'
        DO K=1,NLAY
          WRITE(OUTRES,3000) ((MC(J,I,K),J=1,NCOL),I=1,NROW)
        ENDDO
      ENDIF
      
 2999 FORMAT(I6,F15.4,3X,A15)
 3000 FORMAT(30000(F15.4))
      
      return
      end
     