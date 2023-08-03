C ************************************************************************************************************
C VST (Variably-Saturated Transport) Package
C Handles subroutines related to variably-saturated transport
C ************************************************************************************************************


C ***********************************************************************************************************
C vst_alloc
C Allocate arrays needed for variably-saturated transport
C ***********************************************************************************************************
C Written: 09-24-2010
      subroutine vst_alloc

      use rt_global
      use VST_ARRAYS
      use mf_rt_link
      use GLOBAL,     only:NCOL,NROW,NLAY

      implicit none

      allocate(MC(NCOL,NROW,NLAY))
      allocate(MCOLD(NCOL,NROW,NLAY))
      allocate(SATOLD(NCOL,NROW,NLAY))
      allocate(SATNEW(NCOL,NROW,NLAY))

      allocate(WATVOL(NCOL,NROW,NLAY))
      allocate(WT(NCOL,NROW))
      allocate(APPL(NCOL,NROW))
      allocate(UZQ(NCOL,NROW,NLAY))
      allocate(UZQSTO(NCOL,NROW,NLAY))
      allocate(GWQOUT(NCOL,NROW,NLAY))
      allocate(RT_GWET(NCOL,NROW,NLAY))
      allocate(QZVAR(NCOL,NROW,NLAY))
      allocate(CAPPL(NCOL,NROW,NCOMP))
      allocate(CUZET(NCOL,NROW,NLAY,NCOMP))
      allocate(CGWET(NCOL,NROW,NLAY,NCOMP))
      
      allocate(MNEW(NCOL,NROW,NLAY,NCOMP))
      allocate(MOLD(NCOL,NROW,NLAY,NCOMP))
      allocate(MTEMP(NCOL,NROW,NLAY,NCOMP))
      MTEMP = 0

      RETURN
      END


C ***********************************************************************************************************
C vst_prep
C Assemble arrays for the VST package
C Assemble:
C     SAT:      Saturation of each cell
C     APPL:     Applied amounts (irrigation water) at surface cells
C     WT:       Water Table elevation for each cell column
C     WATVOL:   Volume of water in each cell
C     QZVAR:    Flux in Z direction for variably-saturated conditions (combines UZQ and QZZ)
C This subroutine is called for every flow time step, since these variables change with the calculated
C water content of the cells.
C ***********************************************************************************************************
C Written: 09-21-2010
      subroutine vst_prep

      use rt_global
      use VST_ARRAYS
      use mf_rt_link
      use GLOBAL,     only:NCOL,NROW,NLAY

      IMPLICIT	NONE
      INTEGER	J,I,K,DUM,N
      REAL	    TOPELEV,SATTHICK,THICK,TOP,BOT,MOIST,W,D,
     &          UZQ_FLUX,
     &          PART1,PART2,
     &          SUM,RESIDUAL_WC

C     Update the Saturated Thickness values (DH) ------------------------------------------------------------
      DO I=1,NROW
        DO J=1,NCOL
          DO K=1,NLAY
            IF(DH(J,I,K).GT.1E10) THEN
              DH(J,I,K) = 0.
            ENDIF
          ENDDO
        ENDDO
      ENDDO

C     Water Content values read in from the HFF file must be corrected in some
C     instances. These instances deal with the water table and the saturated zone.
C     For cells in which either (i) the water table is present or (ii) full saturation
C     exists, the WC in the HFF file is set to 0. Hence, WC must be adjusted in
C     these situations.
C     Loop through all cells in the grid
      RESIDUAL_WC = 0.100
      DO I=1,NROW
        DO J=1,NCOL
          DO K=1,NLAY                      
C           1. Make sure that water contents are not more than porosity            
            IF (MC(J,I,K).GT.PRSITY(J,I,K)) THEN
              MC(J,I,K) = PRSITY(J,I,K)
C           2. Change MC value if it is equal to 0.           
            ElSEIF (MC(J,I,K).LE.0) THEN
C             Condition 1: Water Table is below the cell
              IF(DH(J,I,K).EQ.0) THEN
                MC(J,I,K) = RESIDUAL_WC
              ELSEIF (DH(J,I,K).GT.(DZ(J,I,K)*10)) THEN !DH has high value for unsat cells
                MC(J,I,K) = RESIDUAL_WC !Residual moisture content
              ELSEIF (DH(J,I,K).GT.DZ(J,I,K)) THEN
                MC(J,I,K) = PRSITY(J,I,K) !If DH just barely more than DZ, it is in saturated zone
C             Condition 2: Water Table resides in the cell. Only if DH > 0
C             (because if DH = 0, then cell is unsaturated cell).
              ELSEIF (DH(J,I,K).LT.DZ(J,I,K) .AND. DH(J,I,K).GT.0) THEN
                PART1 = (DZ(J,I,K)-DH(J,I,K)) / DZ(J,I,K)
                PART2 = DH(J,I,K) / DZ(J,I,K)
                MC(J,I,K) = (PART1*RESIDUAL_WC) + (PART2*PRSITY(J,I,K))
C             Condition 3: Water Table is above the cell (fully saturated)              
              ELSEIF (DH(J,I,K).EQ.DZ(J,I,K)) THEN
                MC(J,I,K) = PRSITY(J,I,K)
              ELSE !Default
                MC(J,I,K) = PRSITY(J,I,K)
              ENDIF
C           3. Change MC value to the residual water content if it is lower than
C              the residual water content.           
            ELSEIF(MC(J,I,K).LT.RESIDUAL_WC) THEN
              MC(J,I,K) = RESIDUAL_WC          
            ENDIF
          ENDDO
        ENDDO
      ENDDO

C     Calculate the saturation for each cell ----------------------------------------------------------------
C     If it is not the first flow time step, then copy SATNEW to SATOLD (so that
C     SATOLD, which corresponds to the SATURATION from the previuos flow time step,
C     can be used in BTN3FM in the Implicit Scheme)
      IF(rt_kper.NE.1 .OR. rt_kstp.NE.1) THEN
        DO K=1,NLAY
          DO I=1,NROW
            DO J=1,NCOL
              SATOLD(J,I,K) = SATNEW(J,I,K)
            ENDDO
          ENDDO
        ENDDO
      ENDIF  

C     Now, calculate the new saturation for each cell SATNEW, and re-calculate the Water Content
C     using the SATNEW value (to keep consistency betweeen SATNEW and MC)
      DO K=1,NLAY
        DO I=1,NROW
          DO J=1,NCOL
            IF(DH(J,I,K).GE.DZ(J,I,K)) THEN
              SATNEW(J,I,K)=1
            ELSEIF(.NOT.ICBUND(J,I,K,1).EQ.0) THEN
              SATNEW(J,I,K)=((DZ(J,I,K)-DH(J,I,K))/DZ(J,I,K))*
     &                        MC(J,I,K)/PRSITY(J,I,K)+
     &                        DH(J,I,K)/DZ(J,I,K)*1
            ENDIF
            MC(J,I,K) = SATNEW(J,I,K)*PRSITY(J,I,K)
          ENDDO
        ENDDO
      ENDDO

C     If is is the first flow time step of the first stress period, then copy SATNEW into
C     SATOLD, so that SATOLD incorporates DH and MC from the first time step.
      IF(rt_kper.EQ.1 .AND. rt_kstp.EQ.1) THEN
        DO K=1,NLAY
          DO I=1,NROW
            DO J=1,NCOL
              SATOLD(J,I,K) = SATNEW(J,I,K)
            ENDDO
          ENDDO
        ENDDO
      ENDIF   

C     Calculate Applied Amounts for each surface cell
C     Copy UZ FLUX values from layer 1 into APPL (applied amounts array),
C     for use in the SSM package.
      DO I=1,NROW
        DO J=1,NCOL
          APPL(J,I) = UZQ(J,I,1) !From K=1
        ENDDO
      ENDDO

C     Calculate Water Table elevation for each cell column. Also, keep track of
C     layer in which water table resides.
      DO I=1,NROW
        DO J=1,NCOL          
          SUM = 0
          DO K=1,NLAY
C           If the saturated thickness is less than the cell thickness,
C           calculate the placement of the water table within the cell.            
            IF (DH(J,I,K).NE.0 .AND. DH(J,I,K).LT.DZ(J,I,K)) THEN
              WT(J,I) = HTOP(J,I) - SUM - (DZ(J,I,K)-DH(J,I,K))
            ELSE
              SUM = SUM + DZ(J,I,K)
            ENDIF  
          ENDDO
        ENDDO
      ENDDO     

C     Calculate Flux in Z direction for each cell. Flux in  unsaturated zone is
C     contained in UZQ (and QZ=0). Flux in saturated zone is contained in QZ
C     (and UZQ = 0). Adding them together provides an array that has all of the 
C     Z flux information in both the unsaturated and saturated zones.

C     Furthermore: the first value in QZ refers to the flux across the BOTTOM of the
C     first cell, whereas the first value in UZQ refers to the flux across the TOP
C     of the cell. Thus, change indices for UZQ (shift up by 1). 
      DO I=1,NROW
        DO J=1,NCOL
          DO K=1,NLAY         
C           Get value of UZ FLUX            
            IF (K.EQ.NLAY) THEN
              UZQ_FLUX = 0.0
            ELSE
              UZQ_FLUX = UZQ(J,I,K+1)
            ENDIF
            QZVAR(J,I,K) = QZ(J,I,K) + UZQ_FLUX
          ENDDO
        ENDDO
      ENDDO 

C     Meld UZF and QZZ arrays so that UZQ doesn't need to be dealt
C     with later in the code. This is especially convenient for the ADV
C     subroutines.
      DO K=1,(NLAY-1)
        DO I=1,NROW
          DO J=1,NCOL
            IF(ICBUND(J,I,K,1).NE.0) THEN
              IF(DH(J,I,K).LT.1E-5) THEN
                QZ(J,I,K) = UZQ(J,I,K+1)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      RETURN
      END