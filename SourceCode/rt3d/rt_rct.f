C ************************************************************************************************************
C RCT (Chemical Reaction) Package
C Handles subroutines related to chemical reactions
C
C Included Subroutines:
C RCT3AL: Allocate dynamic arrays and read input data
C RCT3SV: Calculate change in cell concentration due to chemical reactions
C ************************************************************************************************************ 


C ************************************************************************************************************
C RCT3AL
C THIS subroutine ALLOCATES SPACE FOR ARRAYS NEEDED BY THE CHEMICAL REACTION (RCT) PACKAGE.
C ************************************************************************************************************
      subroutine rct_alloc

      use rt_global
      use VST_ARRAYS
      use mf_rt_link
      use GLOBAL,     only:NCOL,NROW,NLAY

      IMPLICIT  NONE
      INTEGER   I,J,K,INDEX,READNUMBER,II,JJ,KK
      REAL*4    CONSTANT,RXN_VALUES
      CHARACTER ANAME*24
      CHARACTER HEADER*80
      DIMENSION RXN_VALUES(NCOL,NROW)


      print *, 'Reading RCT file...'

*     Read the first line of RCT data
	READ(INRCT,*) HEADER
	READ(INRCT,*) ISOTHM,IREACT,NCRXNDATA,NVRXNDATA,ISOLVER,IRCTOP

*     Write out the type of sorption for the simulation
      IF(ISOTHM.EQ.1) THEN
        WRITE(IOUT,1022)
      ELSEIF(ISOTHM.EQ.2) THEN
        WRITE(IOUT,1024)
      ELSEIF(ISOTHM.EQ.3) THEN
        WRITE(IOUT,1026)
      ELSE
        WRITE(IOUT,1028)
      ENDIF
      IF(IREACT.EQ.0) THEN
        WRITE(IOUT,1030)
      ELSE
        WRITE(IOUT,1032)
      ENDIF

*     Allocate (dynamically) the arrays for the RCT package
      allocate(RHOB(NCOL,NROW,NLAY))
      if(ISOTHM.NE.0) allocate(SP1(NCOL,NROW,NLAY,MCOMP))
      if(ISOTHM.NE.0) allocate(SP2(NCOL,NROW,NLAY,MCOMP))
      allocate(RC(NCRXNDATA))
      allocate(VRC(NCOL,NROW,NLAY,NVRXNDATA))

*     DTRCT is not calculated in RT3D; a dummy value is assigned (this is used
*     to calculate the transport time step size in BTN3AD)
	DTRCT = 9999.0

*     Initialize the Retardation Factor to 1.0 for all species
	RETA = 1.0
      WRITE(IOUT,1000)

*     Read Bulk Density values --------------------------------------------------------------------
      if(ISOTHM.GT.0) then
      ANAME='BULK DENSITY (RHOB)    '
      READ(INRCT,*) HEADER
      IF (IRCTOP.LT.1) THEN  !If 1: Read layer by layer
        READ(INRCT,*) READNUMBER,CONSTANT
        IF (READNUMBER.EQ.0) THEN
          DO K=1,NLAY
            DO I=1,NROW
              DO J=1,NCOL
                RHOB(J,I,K) = CONSTANT
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ELSE !read in values layer by layer
        DO K=1,NLAY  
          DO I=1,NROW
            READ(INRCT,*) (RHOB(J,I,K),J=1,NCOL)
          ENDDO
        ENDDO
      ENDIF
      endif

*     Read Sorption parameters --------------------------------------------------------------------
      IF(ISOTHM.GT.0) THEN      
        READ(INRCT,*) HEADER
        IF (IRCTOP.LT.1) THEN
          !First Sorption Parameter (SP1)          
          DO INDEX=1,MCOMP
            READ(INRCT,*) READNUMBER,CONSTANT
            IF (READNUMBER.EQ.0) THEN
              DO K=1,NLAY
                DO I=1,NROW
                  DO J=1,NCOL
                    SP1(J,I,K,INDEX) = CONSTANT
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDDO
          !Second Sorption Parameter (SP2)
          DO INDEX=1,MCOMP
            READ(INRCT,*) READNUMBER,CONSTANT
            IF (READNUMBER.EQ.0) THEN
              DO K=1,NLAY
                DO I=1,NROW
                  DO J=1,NCOL
                    SP2(J,I,K,INDEX) = CONSTANT
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDDO           
        ELSE !read layer-by-layer
          !First Sorption Parameter (SP1)
          DO INDEX=1,MCOMP
            DO K=1,NLAY  
              DO I=1,NROW
                READ(INRCT,*) (SP1(J,I,K,INDEX),J=1,NCOL)
              ENDDO
            ENDDO
          ENDDO      
          !Second Sorption Parameter (SP2)          
          DO INDEX=1,MCOMP  
            DO K=1,NLAY  
              DO I=1,NROW
                READ(INRCT,*) (SP2(J,I,K,INDEX),J=1,NCOL)
              ENDDO
            ENDDO
          ENDDO           
        ENDIF
        !Calculate retardation factor (RETA) for all mobile components,
        rfmin = 1.0e30
        INDEX = 1
        DO ICOMP=1,NCOMP
          IF (PHASE(ICOMP).EQ.1) THEN !Specie is mobile
            call reta_calc(index)
            index = index + 1
          ENDIF
        ENDDO
      ENDIF

*     Read parameters for chemical reactions ------------------------------------------------------    
*     Read spatially-constant reaction parameters: put in VRC    
      IF (NCRXNDATA.GT.0) THEN
        READ(INRCT,*) HEADER
        DO I=1,NCRXNDATA
          READ (INRCT,*) RC(I)
        ENDDO
      ENDIF
C     Read spatially-variable reaction parameters: put in VRC
      ANAME = 'Rxn Constants'
      IF (NVRXNDATA.GT.0) THEN   
          READ(INRCT,*) HEADER
          DO I=1,NVRXNDATA
            READ(INRCT,*) READNUMBER,CONSTANT
            !Read the value and copy to all grid cells
            IF (READNUMBER.EQ.0) THEN
              DO KK=1,NLAY
                DO II=1,NROW
                  DO JJ=1,NCOL
                    VRC(JJ,II,KK,I) = CONSTANT
                  ENDDO
                ENDDO
              ENDDO
            !Read an array of values for the first layer, copy to other layers
            ELSEIF(READNUMBER.EQ.1) THEN
              READ(INRCT,*) ((RXN_VALUES(JJ,II),JJ=1,NCOL),II=1,NROW)
              DO KK=1,NLAY
                DO II=1,NROW
                  DO JJ=1,NCOL
                    VRC(JJ,II,KK,I) = RXN_VALUES(JJ,II)
                  ENDDO
                ENDDO
              ENDDO
            !Read array of values for entire 3D grid
            ELSE
              DO II=1,NROW
                READ(INRCT,*) (VRC(JJ,II,KK,I),JJ=1,NCOL)
              ENDDO
            ENDIF               
          ENDDO
      ENDIF 	   

 1000 FORMAT(//1X,'SORPTION AND 1ST ORDER RATE REACTION PARAMETERS',
     & /1X,47('-')/)
 1022 FORMAT(1X,'TYPE OF SORPTION SELECTED IS LINEAR')
 1024 FORMAT(1X,'TYPE OF SORPTION SELECTED IS FREUNDLICH')
 1026 FORMAT(1X,'TYPE OF SORPTION SELECTED IS LANGMUIR')
 1028 FORMAT(1X,'NO SORPTION ISOTHERM IS SIMULATED')
 1030 FORMAT(1X,'NO REACTION IS SIMULATED')
 1032 FORMAT(1X,'KINETICS ARE SIMULATED')

  500 return
      end


C ************************************************************************************************************
C RCT3SV
C THIS subroutine CALCULATES THE CHANGES INRCT CELL CONCENTRATIONS DUE TO CHEMICAL REACTIONS--ONLY SORPTION AND 
C FIRST-ORDER DECAY (OR BIODEGRADATION) ARE CONSIDERED INRCT THE CURRENT VERSION.
C ************************************************************************************************************
      subroutine rct_solve

      use rt_global
      USE VST_ARRAYS
      use mf_rt_link
      use GLOBAL,     only:NCOL,NROW,NLAY,DELR,DELC

      IMPLICIT  NONE
      INTEGER   K,I,J,L,N,INDEX,DAYCOUNT,YRINDEX,dum
      REAL      SOLID,DCRCT,O2_SAT,DCONC,VOLUME,CELL_TOP,
     &          FERTDAY,PLANTDAY,PLOWDAY,RTDEPTH,RTDPELEV
      DOUBLE PRECISION CONC_OLD,CONC_NEW,CONC_RG,K_RG,T_RG,phi
      DIMENSION CONC_OLD(NCOMP),CONC_NEW(NCOMP),CONC_RG(NCOMP),
     &          K_RG(4,NCOMP),T_RG(4),phi(NCOMP)
      

*     Calculate change in concentration for each component for each cell    
      do k=1,nlay
        do I=1,nrow
          do J=1,ncol

            if(i.eq.112 .and. j.eq.80) then
              dum = 10
            endif
            
            !only proceed if the cell is active
            if(ICBUND(J,I,K,1).EQ.0) goto 120

            !get current concentration of all species
            do INDEX=1,NCOMP
              if(CNEW(J,I,K,INDEX) .LT. 0.0) then
                CNEW(J,I,K,INDEX) = 0.0
              endif
              CONC_OLD(INDEX) = dble(CNEW(J,I,K,INDEX)) !double prec.
            enddo                     
            
            !solve for new concentration. If KINETICS is equal to 10, use user-defined reaction
            !module
            if(IREACT.EQ.10) then
              
              !calculate the change in species concentrations using the 4th-order Runge-Kutta
              !scheme.for each slope, the runge-kutta slopes will be calculated using the
              !runge-kutta concentrations
              
              !K1 (first slope)
              T_RG(1) = TIME2
              do n=1,ncomp
                CONC_RG(n) = CONC_OLD(n)
              enddo
              call rxns(J,I,K,CONC_RG,K_RG,1) !k_rg will be filled up

              !K2 (second slope)
              T_RG(2) = TIME2 + (0.5*DT)
              do n=1,ncomp
                CONC_RG(n) = CONC_OLD(n) + (0.5*DT*K_RG(1,n))
              enddo
              call rxns(J,I,K,CONC_RG,K_RG,2)
              
              !K3 (third slope)
              T_RG(3) = TIME2 + (0.5*DT)
              do n=1,ncomp
                CONC_RG(n) = CONC_OLD(n) + (0.5*DT*K_RG(2,n))
              enddo
              call rxns(J,I,K,CONC_RG,K_RG,3) !k_rg will be filled up

              !K4 (fourth slope)
              T_RG(4) = TIME2 + DT
              do n=1,ncomp
                CONC_RG(n) = CONC_OLD(n) + (DT*K_RG(3,n))
              enddo
              call rxns(J,I,K,CONC_RG,K_RG,4) !k_rg will be filled up

              !calculate new concentration
              do n=1,ncomp
                !calculate the increment, then the new concentration
                phi(n) = (1./6.) * (k_rg(1,n) + (2*k_rg(2,n)) + 
     &                           (2*k_rg(3,n)) +   k_rg(4,n))
                CONC_NEW(n) = CONC_OLD(n) + (phi(n)*DT)
              enddo
            endif


C           Re-assert boundary conditions and calculate mass in/out
            do INDEX=1,NCOMP
              if (CONC_NEW(INDEX).LT.0.0) then
                CONC_NEW(INDEX) = 0.0
	        endif

              !If constant concentration cell, re-assert old concentration. 
	        if (ICBUND(J,I,K,INDEX).LT.0) then
	          CONC_NEW(INDEX) = CONC_OLD(INDEX)
	        endif
              
              !Copy new concentration into rt_global concentration array
              CNEW(J,I,K,INDEX) = CONC_NEW(INDEX)

              !Calculate change in concentration, and add to mass in/out array
              DCRCT = CONC_NEW(INDEX) - CONC_OLD(INDEX)

              !Calculate change in mass and store in RMASIO              
              IF(TRNOP(7)) THEN
                !Calculate new mass for the cell and store mass change in RMASIO
                VOLUME = DELR(J)*DELC(I)*DZ(J,I,K) !Bulk Volume
                IF (PHASE(INDEX).EQ.1) THEN !Specie is mobile
                  MNEW(J,I,K,INDEX)=CNEW(J,I,K,INDEX)*MC(J,I,K)*VOLUME
                  IF (DCRCT.LT.0) THEN !Decrease in concentration
                    RMASIO(9,2,INDEX) = RMASIO(9,2,INDEX) + 
     &               (DCRCT*VOLUME*MC(J,I,K)*RETA(J,I,K,INDEX))
                  ELSE !Increase in concentration
                    RMASIO(9,1,INDEX) = RMASIO(9,1,INDEX) + 
     &               (DCRCT*VOLUME*MC(J,I,K)*RETA(J,I,K,INDEX))
                  ENDIF                 
                ELSE !Specie is immobile (solid-phase)
                  SOLID = 1 - PRSITY(J,I,K)
                  MNEW(J,I,K,INDEX)=CNEW(J,I,K,INDEX)*SOLID*VOLUME
                  IF (DCRCT.LT.0) THEN !Decrease in concentration
                    RMASIO(9,2,INDEX) = RMASIO(9,2,INDEX) + 
     &               (DCRCT*VOLUME*SOLID)
                  ELSE !Increase in concentration
                    RMASIO(9,1,INDEX) = RMASIO(9,1,INDEX) + 
     &               (DCRCT*VOLUME*SOLID)
                  ENDIF                    
                ENDIF
              ELSE
                !Store mass change in RMASIO array.                         
                IF(DCRCT.LT.0) THEN
                  RMASIO(9,2,INDEX)=RMASIO(9,2,INDEX)+DCRCT*
     &            RETA(J,I,K,INDEX)*DELR(J)*DELC(I)*
     &            DF(J,I,K)*PRSITY(J,I,K)
                ELSE
                  RMASIO(9,1,INDEX)=RMASIO(9,1,INDEX)+DCRCT*
     &            RETA(J,I,K,INDEX)*DELR(J)*DELC(I)*
     &            DF(J,I,K)*PRSITY(J,I,K)
                ENDIF
              ENDIF
            enddo
                          
120       enddo
        enddo            
      enddo 


1000  FORMAT(400(F15.4))  
1001  FORMAT(12(F15.4))
1002  FORMAT(F10.4,I5,I5,I5,F15.4)
      
200   RETURN
      END


C ************************************************************************************************************
C reta_calc (Retardation Factor Coefficient)
C THIS subroutine EVALUATES RETARDATION FACTOR FOR THREE TYPES OF SORPTION ISOTHERMS: LINEAR, FREUNDLICH AND 
C LANGMUIR.
C The array MC holds the moisture content for each grid cell. If fully saturated conditions are simulated,
C then this is simply POROSITY.
C ************************************************************************************************************
      subroutine reta_calc(INDEX)

      use rt_global
      USE VST_ARRAYS
      use mf_rt_link
      use GLOBAL,     only:NCOL,NROW,NLAY

      IMPLICIT  NONE
      INTEGER   J,I,K,L,INDEX
      REAL*4    C,OM,VOLUME,OMENHANCE,SO4,PO4,INHIB,
     &          KD,KD_SED,TOTAL

C     EVALUATE RETARDATION FACTOR DEPENDING ON TYPES OF SORPTION SELECTED
      IF(ISOTHM.EQ.1) THEN
        DO 10 K=1,NLAY
          DO 20 I=1,NROW
            DO 30 J=1,NCOL
              IF(ICBUND(J,I,K,ICOMP).EQ.0) GOTO 30
              if(trnop(7)) then
                RETA(J,I,K,ICOMP)=1.+RHOB(J,I,K)/MC(J,I,K) *
     &                            SP1(J,I,K,INDEX)
              else
                RETA(J,I,K,ICOMP)=1.+RHOB(J,I,K)/PRSITY(J,I,K) *
     &                            SP1(J,I,K,INDEX)
              endif
              RFMIN=MIN(RFMIN,RETA(J,I,K,ICOMP))
   30       CONTINUE
   20     CONTINUE
   10   CONTINUE
      ELSEIF(ISOTHM.EQ.2) THEN
        DO 40 K=1,NLAY
          DO 50 I=1,NROW
            DO 60 J=1,NCOL
              IF(ICBUND(J,I,K,ICOMP).EQ.0) GOTO 60
              IF(COLD(J,I,K,ICOMP).LE.0) THEN
                RETA(J,I,K,ICOMP)=1.
              ELSE
                C = COLD(J,I,K,ICOMP)
                if(trnop(7)) then
                RETA(J,I,K,ICOMP)=1.+RHOB(J,I,K)/MC(J,I,K)*
     &          SP1(J,I,K,INDEX)*SP2(J,I,K,INDEX)*COLD(J,I,K,ICOMP) **
     &         (SP2(J,I,K,INDEX)-1.)
                else
                RETA(J,I,K,ICOMP)=1.+RHOB(J,I,K)/PRSITY(J,I,K)*
     &          SP1(J,I,K,INDEX)*SP2(J,I,K,INDEX)*COLD(J,I,K,ICOMP) **
     &         (SP2(J,I,K,INDEX)-1.)
                endif
              ENDIF
              RFMIN=MIN(RFMIN,RETA(J,I,K,ICOMP))
   60       CONTINUE
   50     CONTINUE
   40   CONTINUE
      ELSEIF(ISOTHM.EQ.3) THEN
        DO 70 K=1,NLAY
          DO 80 I=1,NROW
            DO 90 J=1,NCOL
              IF(ICBUND(J,I,K,ICOMP).EQ.0) GOTO 90
              if(trnop(7)) then
              RETA(J,I,K,ICOMP)=1.+RHOB(J,I,K)/MC(J,I,K)*
     &        SP1(J,I,K,INDEX)*SP2(J,I,K,INDEX)/(1.+SP1(J,I,K,INDEX) *
     &        COLD(J,I,K,ICOMP))**2
              else
              RETA(J,I,K,ICOMP)=1.+RHOB(J,I,K)/PRSITY(J,I,K)*
     &        SP1(J,I,K,INDEX)*SP2(J,I,K,INDEX)/(1.+SP1(J,I,K,INDEX) *
     &        COLD(J,I,K,ICOMP))**2
              endif
              RFMIN=MIN(RFMIN,RETA(J,I,K,ICOMP))
   90       CONTINUE
   80     CONTINUE
   70   CONTINUE
      ENDIF


      return 
      end
