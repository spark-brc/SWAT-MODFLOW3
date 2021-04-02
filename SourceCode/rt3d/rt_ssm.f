C ***********************************************************************************************************
C ssm_alloc
C THIS subroutine ALLOCATES SPACE FOR ARRAYS NEEDED INSSM THE SINK & SOURCE MIXING (SSM) PACKAGE.
C ***********************************************************************************************************
      subroutine ssm_alloc

      use rt_global
      use mf_rt_link
      use GLOBAL,     only:NCOL,NROW,NLAY

      implicit none
      character header*80

*     Read whether or not there are constant concentration cells (to be associated
*     with constant-head cells from MODFLOW or general-head boundary cells from
*     MODFLOW)
      read(INSSM,*) header
      read(INSSM,*) fcnc

*     If VST package is used, read in additional line of flags
      if (TRNOP(7)) then
        read(INSSM,*) header
        read(INSSM,*) FAPP
      endif     

*     Allocate (dynamically) arrays for the SSM package
      
      !point source-sink arrays
      allocate(CNH(maxcnh,4))
      allocate(WEL(maxwel,4))
      allocate(DRN(maxdrn,4))
      allocate(RIV(maxriv,4))
      allocate(GHB(maxghb,4))
      rt_ncnh = 0
      rt_nwel = 0
      rt_ndrn = 0
      rt_nriv = 0
      rt_nghb = 0

      !concentration arrays (only for CNH and GHB) - these will be filled up using
      !data from the SSM file
      if(fcnc) then
        allocate(ccnc(NCOL,NROW,NLAY,NCOMP))
        ccnc = 0
      endif

      !river concentration array (will be filled up by SWAT)
      allocate(rt_criv(maxriv,ncomp))
      rt_criv = 0.

      !river mass file (to be given back to SWAT)
      allocate(rt_rivmass(maxriv,ncomp))
      rt_rivmass = 0.

      !drain mass file (to be given back to SWAT)
      allocate(rt_drnmass(maxdrn,ncomp))
      rt_drnmass = 0.

      !diffuse arrays - integer
      IF(FRCH) allocate(RT_IRCH(NCOL,NROW))
      !IF(FEVT) allocate(RT_IEVT(NCOL,NROW))
      !diffuse arrays - real
      IF(FRCH) allocate(RT_RECH(NCOL,NROW))
      IF(FRCH) allocate(CRCH(NCOL,NROW,NCOMP))
      !IF(FEVT) allocate(RT_EVTR(NCOL,NROW))
      !IF(FEVT) allocate(CEVT(NCOL,NROW,NCOMP))

      allocate(SSCUMUL(NCOMP,10000)) !cumulative source-sink mass
      SSCUMUL = 0.

      return
      end


C ***********************************************************************************************************
C ssm_read
C THIS subroutine READS CONCENTRATIONS OF SOURCES OR SINKS NEEDED BY THE SINK AND SOURCE MIXING (SSM) PACKAGE
C Reads values for current transport stress period
C ***********************************************************************************************************
      subroutine ssm_read

      use rt_global
      use vst_arrays
      use mf_rt_link
      use GLOBAL,     only:NCOL,NROW,NLAY

      implicit none
      INTEGER   I,J,K,JJ,II,KK,NUM,IQ,INCRCH,NTMP,INDEX,ierr,
     &          INCAPPL,INCGWET,APPLREAD,IREAD
      REAL      CSS,CONSTANT,APPLVAL,dum_array
      CHARACTER ANAME*24,HEADER*80
      dimension dum_array(100)


C     Solute concentration of diffusive sources/sinks are supplied by SWAT (concentration
C       in deep percolation water)


C     This subroutine is only needed if constant-head or general-head boundary condition
C       solute concentrations are specified.

C     Point source/sinks of specified concentration
      
      !Constant-concentration cells
      if(fcnc) then
        read(INSSM,*) header
        read(INSSM,*) ntmp !number of constant concentration cells
	  do num=1,ntmp
          !Read in the K,I,J location and Concentration for each species
          read(INSSM,*) kk,ii,jj,(dum_array(index),index=1,mcomp)
          index = 1
          do k=1,ncomp
            if (phase(k).EQ.1) then !Specie is mobile
              ccnc(jj,ii,kk,k) = dum_array(index)
              index = index + 1
            else
              ccnc(jj,ii,kk,k) = 0
            endif
          enddo
        enddo
      endif

      return
      end


C ***********************************************************************************************************
C ssm_imp
C THIS subroutine FORMULATES MATRIX COEFFICIENTS FOR THE SINK/SOURCE TERMS IF THE IMPLICIT SCHEME IS USED
C ***********************************************************************************************************
      subroutine ssm_imp(A,RHS)

      use rt_global
      use vst_arrays
      use mf_rt_link
      use GLOBAL,     only:NCOL,NROW,NLAY,DELR,DELC

      IMPLICIT  NONE
      INTEGER   I,J,K,N,dum,NUM,IQ
      REAL      ctmp,QSS,A,RHS
      DIMENSION A(NODES),RHS(NODES)

      
C     Formulate [A] and [RHS] matrices

      !Areal sink/source terms --------------------------------------------------------------------
      
      !Applied Amounts (irrigation water)
      if(trnop(7)) then
        if(fapp) then !Only if flows and concentrations are defined
          do I=1,NROW
            do J=1,NCOL          
              k = 2 !layer below the soil profile
              if (ICBUND(J,I,K,ICOMP).LE.0) cycle
              n = (K-1)*NCOL*NROW+(I-1)*NCOL+J
              ctmp = cappl(j,i,icomp)
              if(APPL(j,i).lt.0) then
                if(UPDLHS) then
                  A(N) = A(N)+(APPL(J,I)*DELR(J)*DELC(I)*DZ(J,I,K))
                endif
              else
                RHS(N) = RHS(N)
     &           -APPL(J,I)*cappl(J,I,ICOMP)*DELR(J)*DELC(I)*DZ(J,I,K)
              endif
            enddo
          enddo
        endif
      endif

      !Recharge water
      !(flow and concentration are supplied by SWAT)
      if(frch) then  
        do I=1,NROW
          do J=1,NCOL
            K = RT_IRCH(J,I)
            if(K.GT.0 .AND. ICBUND(J,I,K,ICOMP).GT.0) then
              n = (K-1)*NCOL*NROW+(I-1)*NCOL+J
              if(RT_RECH(J,I).LT.0) then
                if(UPDLHS) then
                  A(N) = A(N)+RT_RECH(J,I)*DELR(J)*DELC(I)*DH(J,I,K)
                endif
              else
                RHS(N) = RHS(N)
     &           -RT_RECH(J,I)*CRCH(J,I,ICOMP)*DELR(J)*DELC(I)*DH(J,I,K)
              endif
            endif
          enddo
        enddo
      endif


C     Point Sink/Source Terms (CONST-HEAD, WELL, DRAIN, RIVER & GHB)  -----------------------------
      
      !Wells
      do num=1,rt_nwel !loop through the wells
        k = WEL(num,1)
        i = WEL(num,2)
        j = WEL(num,3)
        qss = WEL(num,4) !flow rate (pumping rate)
        ctmp = cold(j,i,k,icomp) !concentration
        if(ICBUND(J,I,K,ICOMP).GT.0) then
          N=(K-1)*NCOL*NROW+(I-1)*NCOL+J
          if(qss.lt.0) then !flow is leaving model domain
            if(UPDLHS) A(N)=A(N)+qss*DELR(J)*DELC(I)*DH(J,I,K)
          else !flow is entering model domain
            RHS(N)=RHS(N)-qss*ctmp*DELR(J)*DELC(I)*DH(J,I,K)
          endif
        endif
      enddo

      !Drains
      do num=1,rt_ndrn
        k = DRN(num,1)
        i = DRN(num,2)
        j = DRN(num,3)
        qss = DRN(num,4) !flow rate (pumping rate)
        ctmp = cold(j,i,k,icomp) !concentration
        if(ICBUND(J,I,K,ICOMP).GT.0) then
          N=(K-1)*NCOL*NROW+(I-1)*NCOL+J
          if(qss.lt.0) then !flow is leaving model domain
            if(UPDLHS) A(N)=A(N)+qss*DELR(J)*DELC(I)*DH(J,I,K)
          else !flow is entering model domain
            RHS(N)=RHS(N)-qss*ctmp*DELR(J)*DELC(I)*DH(J,I,K)
          endif
        endif
      enddo

      !River
      !If river seepage, than flow and concentration come from SWAT
      !If discharge to river, flow comes from MODFLOW and concentration comes from RT3D
      do num=1,rt_nriv
        k = RIV(num,1)
        i = RIV(num,2)
        j = RIV(num,3)
        qss = RIV(num,4) !flow rate
        if(qss.gt.0) then
          ctmp = rt_criv(num,icomp) !from SWAT
        else
          ctmp = cold(j,i,k,icomp) !from RT3D
        endif
        if(ICBUND(J,I,K,ICOMP).GT.0) then
          N=(K-1)*NCOL*NROW+(I-1)*NCOL+J
          if(qss.lt.0) then !flow is leaving model domain
            if(UPDLHS) A(N)=A(N)+qss*DELR(J)*DELC(I)*DH(J,I,K)
          else !flow is entering model domain
            RHS(N)=RHS(N)-qss*ctmp*DELR(J)*DELC(I)*DH(J,I,K)
          endif
        endif
      enddo

      !Constant-Head Boundary
      do num=1,rt_ncnh
        k = CNH(num,1)
        i = CNH(num,2)
        j = CNH(num,3)
        qss = CNH(num,4) !flow rate (pumping rate)
        if(qss.gt.0) then
          if(fcnc) then
            ctmp = ccnc(j,i,k,icomp) !from SSM input file
          else
            ctmp = 0.0
          endif
        else
          ctmp = cold(j,i,k,icomp) !from RT3D
        endif
        if(ICBUND(J,I,K,ICOMP).GT.0) then
          N=(K-1)*NCOL*NROW+(I-1)*NCOL+J
          if(qss.lt.0) then !flow is leaving model domain
            if(UPDLHS) A(N)=A(N)+qss*DELR(J)*DELC(I)*DH(J,I,K)
          else !flow is entering model domain
            RHS(N)=RHS(N)-qss*ctmp*DELR(J)*DELC(I)*DH(J,I,K)
          endif
        endif
      enddo

      !General-Head Boundary
      do num=1,rt_nghb
        k = GHB(num,1)
        i = GHB(num,2)
        j = GHB(num,3)
        qss = GHB(num,4) !flow rate (pumping rate)
        if(qss.gt.0) then
          if(fcnc) then
            ctmp = ccnc(j,i,k,icomp) !from SSM input file
          else
            ctmp = 0.0
          endif
        else
          ctmp = cold(j,i,k,icomp) !from RT3D
        endif
        if(ICBUND(J,I,K,ICOMP).GT.0) then
          N=(K-1)*NCOL*NROW+(I-1)*NCOL+J
          if(qss.lt.0) then !flow is leaving model domain
            if(UPDLHS) A(N)=A(N)+qss*DELR(J)*DELC(I)*DH(J,I,K)
          else !flow is entering model domain
            RHS(N)=RHS(N)-qss*ctmp*DELR(J)*DELC(I)*DH(J,I,K)
          endif
        endif
      enddo

      return
      end


C ***********************************************************************************************************
C ssm_imp_budget
C THIS subroutine CALCULATES MASS BUDGETS ASSOCIATED WITH ALL SINK/SOURCE TERMS.
C ***********************************************************************************************************
      subroutine ssm_imp_budget

      use rt_global
      use VST_ARRAYS
      use mf_rt_link
      use GLOBAL,     only:NCOL,NROW,NLAY,DELR,DELC

      IMPLICIT  NONE
      INTEGER   NUM,IQ,K,I,J,dum,N
      REAL      CTMP,QSS,mass_in,volume,mass


C     APPLIED AMOUNTS (associated with deep percolation)
      if(TRNOP(7)) then
        if(FAPP) then !Only if flows and concentrations are defined
          do I=1,NROW
            do J=1,NCOL                                   
              k = 2 !deep percolation layer
              if(ICBUND(J,I,K,ICOMP).LE.0) cycle
              volume = delr(j)*delc(i)*dz(j,i,k)
              mass_in = appl(j,i)*cappl(j,i,icomp)*volume*dt
              RMASIO(6,1,ICOMP) = RMASIO(6,1,ICOMP) + mass_in
            enddo
          enddo
        endif
      endif

C      IQ
C      1 = constant-head boundary
C      2 = well
C      3 = drain
C      4 = river
C      5 = general-head boundary

      !Constant-Head Boundary
      iq = 1
      do num=1,rt_ncnh
        k = CNH(num,1)
        i = CNH(num,2)
        j = CNH(num,3)
        qss = CNH(num,4) !flow rate
        if(qss.gt.0) then
          if(fcnc) then
            ctmp = ccnc(j,i,k,icomp) !from SSM input file
          else
            ctmp = 0.0
          endif
        else
          ctmp = cold(j,i,k,icomp) !from RT3D
        endif
        if(ICBUND(J,I,K,ICOMP).GT.0) then
          N=(K-1)*NCOL*NROW+(I-1)*NCOL+J
          IF(QSS.GT.0) THEN
            RMASIO(IQ,1,ICOMP)=RMASIO(IQ,1,ICOMP)+QSS*CTMP*DT*
     &       DELR(J)*DELC(I)*DH(J,I,K)
          ELSE
            RMASIO(IQ,2,ICOMP)=RMASIO(IQ,2,ICOMP)+QSS*CTMP*DT*
     &       DELR(J)*DELC(I)*DH(J,I,K)
          ENDIF
        endif
      enddo

C     Point Sink/Source Terms (CONST-HEAD, WELL, DRAIN, RIVER & G-H-B) 
      !Wells
      iq = 2
      do num=1,rt_nwel
        k = WEL(num,1)
        i = WEL(num,2)
        j = WEL(num,3)
        qss = WEL(num,4) !flow rate (pumping rate)
        ctmp = cnew(j,i,k,icomp) !concentration
        if(ICBUND(J,I,K,ICOMP).GT.0) then
          N=(K-1)*NCOL*NROW+(I-1)*NCOL+J
          IF(QSS.GT.0) THEN
            RMASIO(IQ,1,ICOMP)=RMASIO(IQ,1,ICOMP)+QSS*CTMP*DT*
     &       DELR(J)*DELC(I)*DH(J,I,K)
          ELSE
            RMASIO(IQ,2,ICOMP)=RMASIO(IQ,2,ICOMP)+QSS*CTMP*DT*
     &       DELR(J)*DELC(I)*DH(J,I,K)
          ENDIF
        endif
      enddo

      !Drains
      iq = 3
      do num=1,rt_ndrn
        k = DRN(num,1)
        i = DRN(num,2)
        j = DRN(num,3)
        qss = DRN(num,4) !flow rate
        ctmp = cold(j,i,k,icomp) !concentration
        volume = delr(j)*delc(i)*dh(j,i,k)
        qss = qss*volume !now units are L3/T
        mass = qss*ctmp*dt
        if(isnan(mass)) then
         dum = 10
        endif
        rt_drnmass(num,icomp) = rt_drnmass(num,icomp) + mass !store for mapping to SWAT
        if(ICBUND(J,I,K,ICOMP).GT.0) then
          N=(K-1)*NCOL*NROW+(I-1)*NCOL+J
          IF(QSS.GT.0) THEN
            RMASIO(IQ,1,ICOMP)=RMASIO(IQ,1,ICOMP)+QSS*CTMP*DT*
     &       DELR(J)*DELC(I)*DH(J,I,K)
          ELSE
            RMASIO(IQ,2,ICOMP)=RMASIO(IQ,2,ICOMP)+QSS*CTMP*DT*
     &       DELR(J)*DELC(I)*DH(J,I,K)
          ENDIF
        endif
      enddo

      !River
      !If river seepage, than flow and concentration come from SWAT
      !If discharge to river, flow comes and from MODFLOW and concentration from RT3D
      iq = 4
      do num=1,rt_nriv
        
        !calculate mass of discharge/seepage
        k = RIV(num,1)
        i = RIV(num,2)
        j = RIV(num,3)
        if(k.eq.1 .and. i.eq.334 .and. j.eq.235) then
          dum = 10
        endif
        qss = RIV(num,4) !flow rate
        if(qss.gt.0) then !flow is entering aquifer
          ctmp = rt_criv(num,icomp) !from SWAT
        else !flow is leaving aquifer
          ctmp = cold(j,i,k,icomp) !from RT3D
        endif
        volume = delr(j)*delc(i)*dh(j,i,k)
        qss = qss*volume !now units are L3/T
        mass = qss*ctmp*dt
        rt_rivmass(num,icomp) = rt_rivmass(num,icomp) + mass !store for mapping to SWAT
        !store mass for mass budget
        if(ICBUND(J,I,K,ICOMP).GT.0) then
          N=(K-1)*NCOL*NROW+(I-1)*NCOL+J
          if(QSS.GT.0) then
            RMASIO(IQ,1,ICOMP) = RMASIO(IQ,1,ICOMP) + mass
          else
            RMASIO(IQ,2,ICOMP) = RMASIO(IQ,2,ICOMP) + mass
          endif
        endif
      enddo

      !General-Head Boundary
      iq = 5
      do num=1,rt_nghb
        k = GHB(num,1)
        i = GHB(num,2)
        j = GHB(num,3)
        qss = GHB(num,4) !flow rate
        if(qss.gt.0) then
          if(fcnc) then
            ctmp = ccnc(j,i,k,icomp) !from SSM input file
          else
            ctmp = 0.0
          endif
        else
          ctmp = cold(j,i,k,icomp) !from RT3D
        endif
        if(ICBUND(J,I,K,ICOMP).GT.0) then
          N=(K-1)*NCOL*NROW+(I-1)*NCOL+J
          IF(QSS.GT.0) THEN
            RMASIO(IQ,1,ICOMP)=RMASIO(IQ,1,ICOMP)+QSS*CTMP*DT*
     &       DELR(J)*DELC(I)*DH(J,I,K)
          ELSE
            RMASIO(IQ,2,ICOMP)=RMASIO(IQ,2,ICOMP)+QSS*CTMP*DT*
     &       DELR(J)*DELC(I)*DH(J,I,K)
          ENDIF
        endif
      enddo


  400 return
      end
