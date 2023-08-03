*     This subroutine runs the MODFLOW simulator for the time period indicated. Typically,
*     the stress period and flow time step are the same as for SWAT (1 day).

*     This subroutine is called for each day of the SWAT simulation. 

      subroutine mf_run
C     Specify which modules are to be used
      USE GLOBAL
      USE GWFBASMODULE
      USE GWFHUFMODULE, ONLY:IOHUFHDS,IOHUFFLWS
      USE GWFEVTMODULE, ONLY:NEVTOP
      USE GWFRCHMODULE, ONLY:NRCHOP, RECH
      USE GWFLAKMODULE, ONLY:NLAKESAR,THETA,STGOLD,STGNEW,VOL
      USE GWFSFRMODULE, ONLY:NUMTAB
      USE GWFNWTMODULE, ONLY:LINMETH,ICNVGFLG,ITREAL
      USE PCGMODULE
      USE SIPMODULE
      USE DE4MODULE
      USE GMGMODULE
      USE LMTMODULE
      use mf_rt_link !tcw
      use smrt_parm
      
      INCLUDE 'mf_openspec.inc'

*     Declare variables
      integer   num_ts,t,kkstp,ii,jj,kk
      real, dimension (:), allocatable :: obs_heads
      real perlen_day(1)
      real mfDayLength(1)

      print *, '   Running MODFLOW  ','Period:',kkper,'  Day:',daytotal
      
      MF_leapyr = 1
      
      IGRID = 1
      if (read_stress.eq.1) then
        kstp = 0 !initialize flow time step for the stress period
            
        IF(IUNIT(62).GT.0 ) CALL GWF2UPWUPDATE(1,Igrid)
        CALL GWF2BAS7ST(kkper,IGRID)
        IF(IUNIT(19).GT.0) CALL GWF2IBS7ST(KKPER,IGRID)
        IF(IUNIT(54).GT.0) CALL GWF2SUB7ST(KKPER,IGRID)
        IF(IUNIT(57).GT.0) CALL GWF2SWT7ST(KKPER,IGRID)
        IF(IUNIT(2).GT.0) CALL GWF2WEL7RP(IUNIT(2),IGRID) !Well
        IF(IUNIT(3).GT.0) CALL GWF2DRN7RP(IUNIT(3),IGRID) !Drain
        IF(IUNIT(4).GT.0) CALL GWF2RIV7RP(IUNIT(4),IGRID) !River stage is provided by SWAT
        IF(IUNIT(5).GT.0) CALL GWF2EVT7RP(IUNIT(5),IGRID) !ET is provided by SWAT
        IF(IUNIT(7).GT.0) CALL GWF2GHB7RP(IUNIT(7),IGRID) !General Head
        IF(IUNIT(8).GT.0) CALL GWF2RCH7RP(IUNIT(8),IGRID) !recharge by is provided by SWAT
        IF(IUNIT(17).GT.0) CALL GWF2RES7RP(IUNIT(17),IGRID) !Reservoir
        IF(IUNIT(18).GT.0) CALL GWF2STR7RP(IUNIT(18),IGRID) !Stream stage is provided by SWAT
        IF(IUNIT(43).GT.0 .AND. IUNIT(18).GT.0)
     &                   CALL GWF2HYD7STR7RP(IUNIT(43),KKPER,IGRID)
        IF(IUNIT(20).GT.0) CALL GWF2CHD7RP(IUNIT(20),IGRID) !CHD
        IF(IUNIT(44).GT.0) CALL GWF2SFR7RP(IUNIT(44),IUNIT(15), !Stream stage is provided by SWAT
     &                                   IUNIT(22),KKPER,KKSTP,
     &                                   NSOL,IOUTS,IUNIT(1),
     &                                   IUNIT(23),IUNIT(37),
     &                                   IUNIT(62), IUNIT(55), IGRID)
        IF(IUNIT(43).GT.0 .AND. IUNIT(44).GT.0)
     &                   CALL GWF2HYD7SFR7RP(IUNIT(43),KKPER,IGRID)
        IF(IUNIT(55).GT.0)CALL GWF2UZF1RP(IUNIT(55),KKPER,!Infiltration and ET is provided by SWAT
     &                                   IUNIT(44),IGRID)
        IF(IUNIT(22).GT.0) CALL GWF2LAK7RP(IUNIT(22),IUNIT(1), !LAK
     &               IUNIT(15),IUNIT(23),IUNIT(37),IUNIT(44),IUNIT(55),
     &               IUNIT(62),KKPER,NSOL,IOUTS,IGRID)
        IF(IUNIT(46).GT.0.AND.KKPER.EQ.1) CALL GWF2GAG7RP(IUNIT(15),
     &               IUNIT(22),IUNIT(55),NSOL,IGRID)
        IF(IUNIT(39).GT.0) CALL GWF2ETS7RP(IUNIT(39),IGRID)
        IF(IUNIT(40).GT.0) CALL GWF2DRT7RP(IUNIT(40),IGRID)
        IF(IUNIT(50).GT.0) CALL GWF2MNW27RP(IUNIT(50),kper,IUNIT(9),
     &                       IUNIT(10),0,IUNIT(13),
     &                       IUNIT(15),IUNIT(63),IGRID)
        IF(IUNIT(51).GT.0.AND.KKPER.EQ.1) CALL GWF2MNW2I7RP(IUNIT(51),
     &                     0,IGRID)
        IF(IUNIT(52).GT.0) CALL GWF2MNW17RP(IUNIT(52),IUNIT(1),
     &                            IUNIT(23),IUNIT(37),IUNIT(62),KKPER,
     &                            IGRID)    
        IF(IUNIT(64).GT.0) CALL GWF2SWR7RP(IUNIT(64),KKPER,IGRID) 
      endif


      !rtb MODFLOW
      !pass SWAT information to MODFLOW -----------------------------------------------------------
      if(iunit(4).gt.0) call smrt_mfRiver !River package
      if(iunit(8).gt.0) call smrt_recharge !Recharge package
      if(iunit(5).gt.0) call smrt_evt !EVT package
      if(iunit(2).gt.0 .and. mf_irrigation_swat.eq.1) call smrt_well !pumping-irrigation linkage
      if(rt_active.eq.1) call smrt_rt3d_stress !RT3D
      

      !determine the time step duration -----------------------------------------------------------

      !convert the mf_interval value (in days) to MODFLOW's current units
      mfDayLength(1) = mf_interval
      call units(mfDayLength, 4, ITMUNI, 1, 1, MF_leapyr)

      !Initialize a variable for MODFLOW's current stress period length in units of days, used at the end to determine looping within stress periods, tcw
      perlen_day(1) = PERLEN(kper) !tcw
      call units(perlen_day, ITMUNI, 4, 1, 1, MF_leapyr) !tcw

      !if mf_interval is greater than the MODFLOW stress period, then set it
      !to the stress period length
      if(mfDayLength(1).gt.PERLEN(kper)) then
        mfDayLength(1) = PERLEN(kper)
      endif
      
      !calculate time step length (as specified in MODFLOW .dis file)
      DELT = PERLEN(kper) / NSTP(kper)

      !If time step is greater than mf_interval, decrease the time step to match mf_interval. This
      !results in only one time step for the current run.
      if(DELT.ge.mfDayLength(1)) then
        DELT = mfDayLength(1)
        num_ts = 1
      elseif(mfDayLength(1).gt.DELT) then
        !calculate required number of time steps
        num_ts = nint(mfDayLength(1)/DELT) + 1
        DELT = mfDayLength(1) / num_ts
      endif

      !copy the flow time step for use by UZF-RT3D
      rt_delt = DELT !rtb


*     Loop through the flow time steps ------------------------------------------------------------
      do t=1,num_ts
        
        kstp = kstp + 1
        kkstp = kstp
        rt_kstp = rt_kstp + 1

        IF (IRESTART.EQ.1) THEN
            IF (KPER.EQ.KPERSTART .AND. KSTP.EQ.KSTPSTART) THEN
              CALL RESTARTHEADS(IOUT)
            END IF
        END IF

        !calculate time step length. set HOLD=HNEW
        IF(IUNIT(62).GT.0 ) CALL GWF2UPWUPDATE(1,Igrid)
        CALL GWF2BAS7AD(KKPER,KKSTP,IGRID)
        IF(IUNIT(62).GT.0) CALL GWF2UPW1AD(IGRID)
        IF(IUNIT(20).GT.0) CALL GWF2CHD7AD(KKPER,IGRID)
        IF(IUNIT(1).GT.0) CALL GWF2BCF7AD(KKPER,IGRID)
        IF(IUNIT(17).GT.0) CALL GWF2RES7AD(KKSTP,KKPER,IGRID)
        IF(IUNIT(23).GT.0) CALL GWF2LPF7AD(KKPER,IGRID)
        IF(IUNIT(37).GT.0) CALL GWF2HUF7AD(KKPER,IGRID)
        IF(IUNIT(16).GT.0) CALL GWF2FHB7AD(IGRID)
        IF(IUNIT(22).GT.0) CALL GWF2LAK7AD(KKPER,KKSTP,IUNIT(15),IGRID)
        IF(IUNIT(44).GT.0 .AND. NUMTAB.GT.0 ) 
     &                             CALL GWF2SFR7AD(IUNIT(22))
          IF(IUNIT(50).GT.0) THEN
            IF (IUNIT(1).GT.0) THEN
              CALL GWF2MNW27BCF(KPER,IGRID)
            ELSE IF (IUNIT(23).GT.0) THEN
              CALL GWF2MNW27LPF(KPER,IGRID)
            ELSE IF(IUNIT(37).GT.0) THEN
              CALL GWF2MNW27HUF(KPER,IGRID)
            ELSE IF(IUNIT(62).GT.0) THEN
              CALL GWF2MNW27UPW(KPER,IGRID)
            ELSE
              WRITE(IOUT,1000)
 1000         FORMAT(/1X,
     &      '***ERROR: MNW2 PACKAGE DOES NOT SUPPORT',/,
     &      ' SELECTED FLOW PACKAGE',/,
     &  ' (MNW2 DOES FULLY SUPPORT BCF, LPF, HUF, AND UPW PACKAGES)',/,
     &      ' -- STOP EXECUTION')
              CALL USTOP('MNW2 error-flow package')
            END IF
            CALL GWF2MNW27AD(KKSTP,KKPER,IUNIT(62),IGRID)
          END IF
          IF(IUNIT(52).GT.0) CALL GWF2MNW17AD(IUNIT(1),IUNIT(23),
     1                                  IUNIT(37),IUNIT(62),IGRID)
!          IF(IUNIT(61).GT.0) THEN                                       !FMP2AD CALL ADDED BY SCHMID
!             IF(IRTFL.EQ.3.OR.ICUFL.EQ.3.OR.IPFL.EQ.3
!     1                           .OR.IEBFL.EQ.1.OR.IEBFL.EQ.3)
!     2       CALL FMP2AD(ISTARTFL,KKPER,IGRID)
!          ENDIF     
          IF(IUNIT(64).GT.0) CALL GWF2SWR7AD(KKPER,IGRID,IUNIT(54))
C
C         INDICATE IN PRINTOUT THAT SOLUTION IS FOR HEADS
          IF ( IRESTART.GT.0 ) THEN
            IF ( KPER.LT.KPERSTART ) THEN
              CALL UMESPR('SKIPPING STEP TO RESTART',' ',IOUT)
              WRITE(*,26)KPER,KSTP
            ELSE IF ( KPER.EQ.KPERSTART .AND. KSTP.LT.KSTPSTART ) THEN
              CALL UMESPR('SKIPPING STEP TO RESTART',' ',IOUT)
              WRITE(*,26)KPER,KSTP
            ELSE
              CALL UMESPR('SOLVING FOR HEAD',' ',IOUT)
              !WRITE(*,25)KPER,KSTP
            END IF
          ELSE
            CALL UMESPR('SOLVING FOR HEAD',' ',IOUT)
            !WRITE(*,25)KPER,KSTP
          END IF
   25     FORMAT(' Solving:  Stress period: ',i5,4x,
     &       'Time step: ',i5,4x,'Groundwater-Flow Eqn.')
   26     FORMAT('Skipping:  Stress period: ',i5,4x,
     &       'Time step: ',i5)

        !formulate and solve the flow equations. this is performed iteratively
        KITER = 0
        ITREAL2 = 0
        NOITER = 1
        if(IUNIT(63).GT.0) ITREAL = 0
        IF(IRESTART.GT.0) THEN
          NOITER = 0
          IF(KPER.GT.KPERSTART) THEN
            NOITER = 1
          ELSEIF (KPER.EQ.KPERSTART .AND. KSTP.GE.KSTPSTART) THEN
            NOITER = 1 
          ENDIF
        ENDIF

        !go through iterations. stop when convergence is achieved
        do while (ITREAL2.LT.MXITER .AND. NOITER.EQ.1)
          KITER = KITER + 1
          KKITER = KITER
          IF(IUNIT(63).EQ.0) ITREAL2 = KITER
          IF(IUNIT(62).GT.0) CALL GWF2UPWUPDATE(2,Igrid)

          !formulate the finite difference equations ----------------------------------------------
          CALL GWF2BAS7FM(IGRID) !set HCOF and RHS to 0
          IF(IUNIT(1).GT.0) CALL GWF2BCF7FM(KKITER,KKSTP,KKPER,IGRID)
          IF(IUNIT(62).GT.0) CALL GWF2UPWFMS(KKITER,KKSTP,KKPER,IGRID)
          IF(IUNIT(23).GT.0) CALL GWF2LPF7FM(KKITER,KKSTP,KKPER,IGRID)
          IF(IUNIT(37).GT.0) CALL GWF2HUF7FM(KKITER,
     &                             KKSTP,KKPER,IUNIT(47),IGRID)
          IF(IUNIT(62).EQ.0 ) THEN
            IF(IUNIT(21).GT.0) CALL GWF2HFB7FM(IGRID)
          END IF
          IF(IUNIT(2).GT.0) CALL GWF2WEL7FM(IUNIT(63),IGRID)
          IF(IUNIT(3).GT.0) CALL GWF2DRN7FM(IGRID)
          IF(IUNIT(4).GT.0) CALL GWF2RIV7FM(IGRID)
          
            !if SWAT tile drainage is active, then remove groundwater here
            !rtb tile
            !if(mf_tile_active.eq.5) then
            !do ii=1,NROW
            !  do jj=1,NCOL
            !    do kk=1,NLAY
            !      if(IBOUND(jj,ii,kk).gt.0) then
            !        RHS(jj,ii,kk) = RHS(jj,ii,kk) - qtile_cells(jj,ii)
            !      endif
            !    enddo
            !  enddo
            !enddo
            !endif
          
          IF(IUNIT(5).GT.0) THEN
            IF(IUNIT(22).GT.0.AND.NEVTOP.EQ.3) THEN
              CALL GWF2LAK7ST(0,IGRID)
            ENDIF
            CALL GWF2EVT7FM(IGRID)
            IF(IUNIT(22).GT.0.AND.NEVTOP.EQ.3) THEN
              CALL GWF2LAK7ST(1,IGRID)
            END IF
          ENDIF
          IF(IUNIT(7).GT.0) CALL GWF2GHB7FM(IGRID)
          IF(IUNIT(8).GT.0) THEN
            IF(IUNIT(22).GT.0.AND.NRCHOP.EQ.3) CALL GWF2LAK7ST(0,IGRID)
            CALL GWF2RCH7FM(IGRID)
            IF(IUNIT(22).GT.0.AND.NRCHOP.EQ.3) CALL GWF2LAK7ST(1,IGRID)
          END IF
          IF(IUNIT(16).GT.0) CALL GWF2FHB7FM(IGRID)
          IF(IUNIT(17).GT.0) CALL GWF2RES7FM(IGRID)
          IF(IUNIT(18).GT.0) CALL GWF2STR7FM(IGRID)
          IF(IUNIT(19).GT.0) CALL GWF2IBS7FM(KKPER,IGRID)
          IF(IUNIT(39).GT.0) CALL GWF2ETS7FM(IGRID)
          IF(IUNIT(40).GT.0) CALL GWF2DRT7FM(IGRID)
          IF(IUNIT(55).GT.0) CALL GWF2UZF1FM(KKPER,KKSTP,KKITER,
     &                           IUNIT(44),IUNIT(22),IUNIT(58),
     &                           IUNIT(63),IGRID)
          IF(IUNIT(44).GT.0) CALL GWF2SFR7FM(KKITER,KKPER,KKSTP,
     &                           IUNIT(22),IUNIT(63),IUNIT(8), 
     &                           IUNIT(55),IGRID)
          IF(IUNIT(22).GT.0) CALL GWF2LAK7FM(KKITER,KKPER,KKSTP,
     &                           IUNIT(44),IUNIT(55),IGRID)
          IF(IUNIT(50).GT.0) THEN
            IF (IUNIT(1).GT.0) THEN
              CALL GWF2MNW27BCF(KPER,IGRID)
            ELSE IF (IUNIT(23).GT.0) THEN
              CALL GWF2MNW27LPF(KPER,IGRID)
            ELSE IF(IUNIT(37).GT.0) THEN
              CALL GWF2MNW27HUF(KPER,IGRID)
            ELSE IF(IUNIT(62).GT.0) THEN
              CALL GWF2MNW27UPW(KPER,IGRID)
            END IF
            CALL GWF2MNW27FM(KKITER,IUNIT(62),kkstp,kkper,IGRID)
          END IF
          IF(IUNIT(52).GT.0) CALL GWF2MNW17FM(KKITER,IUNIT(1),
     &                               IUNIT(23),IUNIT(37),IUNIT(62),
     &                               IGRID)
          IF(IUNIT(54).GT.0) CALL GWF2SUB7FM(KKPER,KKITER,
     &                                         IUNIT(9),IGRID)
            IF(IUNIT(57).GT.0) CALL GWF2SWT7FM(KKPER,IGRID)
            IF(IUNIT(64).GT.0) CALL GWF2SWR7FM(KKITER,KKPER,KKSTP,IGRID)
      
          !calculate new head values
            IF (IUNIT(9).GT.0) THEN
                   CALL SIP7PNT(IGRID)
                   CALL SIP7AP(HNEW,IBOUND,CR,CC,CV,HCOF,RHS,EL,FL,GL,
     1               V,W,HDCG,LRCH,NPARM,KKITER,HCLOSE,ACCL,ICNVG,
     2               KKSTP,KKPER,IPCALC,IPRSIP,MXITER,NSTP(KKPER),
     3               NCOL,NROW,NLAY,NODES,IOUT,0,IERR)
            END IF
            IF (IUNIT(10).GT.0) THEN
                   CALL DE47PNT(IGRID)
                   CALL DE47AP(HNEW,IBOUND,AU,AL,IUPPNT,IEQPNT,D4B,MXUP,
     1               MXLOW,MXEQ,MXBW,CR,CC,CV,HCOF,RHS,ACCLDE4,KITER,
     2               ITMX,MXITER,NITERDE4,HCLOSEDE4,IPRD4,ICNVG,NCOL,
     3               NROW,NLAY,IOUT,LRCHDE4,HDCGDE4,IFREQ,KKSTP,KKPER,
     4               DELT,NSTP(KKPER),ID4DIR,ID4DIM,MUTD4,
     5               DELTL,NBWL,NUPL,NLOWL,NLOW,NEQ,NUP,NBW,IERR)
            END IF
            IF (IUNIT(13).GT.0) THEN
                   CALL PCG7PNT(IGRID)
                   CALL PCG7AP(HNEW,IBOUND,CR,CC,CV,HCOF,RHS,VPCG,SS,
     1               P,CD,HCHG,LHCH,RCHG,LRCHPCG,KKITER,NITER,
     2               HCLOSEPCG,RCLOSEPCG,ICNVG,KKSTP,KKPER,IPRPCG,
     3               MXITER,ITER1,NPCOND,NBPOL,NSTP(KKPER),NCOL,NROW,
     4               NLAY,NODES,RELAXPCG,IOUT,MUTPCG,IT1,DAMPPCG,BUFF,
     5               HCSV,IERR,HPCG,DAMPPCGT,ISSFLG(KKPER),HDRY,
     6               IHCOFADD)
            END IF
            IF (IUNIT(42).GT.0) THEN
                   CALL GMG7PNT(IGRID)
                   CALL GMG7AP(HNEW,RHS,CR,CC,CV,HCOF,HNOFLO,IBOUND,
     1                         IITER,MXITER,RCLOSEGMG,HCLOSEGMG,
     2                         KKITER,KKSTP,KKPER,NCOL,NROW,NLAY,ICNVG,
     3                         SITER,TSITER,DAMPGMG,IADAMPGMG,IOUTGMG,
     4                         IOUT,GMGID,
     5                         IUNITMHC,DUP,DLOW,CHGLIMIT,
     6                         BIGHEADCHG,HNEWLAST)
            ENDIF
          !Newton Solver
          IF(IUNIT(63).GT.0 ) 
     1          CALL GWF2NWT1FM(KKITER,ICNVG,KSTP,KPER,Mxiter,
     2                          IUNIT(22),IGRID)
          IF ( IUNIT(63).GT.0 )ITREAL2 = ITREAL
C
C-------ENSURE CONVERGENCE OF SWR - BASEFLOW CHANGES LESS THAN TOLF - JDH
            IF(IUNIT(64).GT.0) THEN
              CALL GWF2SWR7CV(KKITER,IGRID,ICNVG,MXITER)
            END IF

          !if convergence criterion has been met, stop iterating
          IF (ICNVG.EQ.1) GOTO 33
        
        enddo !go to the next iteration
        KITER = MXITER


 33     CONTINUE
        IF (NOITER.EQ.1) THEN
          IF(IUNIT(62).GT.0 ) CALL GWF2UPWUPDATE(2,Igrid)
        
        !determine if output is required for the current day
        CALL GWF2BAS7OC(KKSTP,KKPER,ICNVG,IUNIT(12),IGRID)

        !calculate budget terms. save cell-by-cell flow terms -----------------------------------------------
        MSUM = 1
        IF (IUNIT(1).GT.0) THEN
          CALL GWF2BCF7BDS(KKSTP,KKPER,IGRID)
          CALL GWF2BCF7BDCH(KKSTP,KKPER,IGRID)
          IBDRET=0
          IC1=1
          IC2=NCOL
          IR1=1
          IR2=NROW
          IL1=1
          IL2=NLAY
          DO IDIR=1,3
            CALL GWF2BCF7BDADJ(KKSTP,KKPER,IDIR,IBDRET,
     &                         IC1,IC2,IR1,IR2,IL1,IL2,IGRID)
          ENDDO
        ENDIF
        IF(IUNIT(23).GT.0) THEN
          CALL GWF2LPF7BDS(KKSTP,KKPER,IGRID)
          CALL GWF2LPF7BDCH(KKSTP,KKPER,IGRID)
          IBDRET=0
          IC1=1
          IC2=NCOL
          IR1=1
          IR2=NROW
          IL1=1
          IL2=NLAY
          DO IDIR=1,3
            CALL GWF2LPF7BDADJ(KKSTP,KKPER,IDIR,IBDRET,
     &                         IC1,IC2,IR1,IR2,IL1,IL2,IGRID)
          ENDDO
        ENDIF

        !ADDED CALLS FOR FLOW BUDGETS TO UPW PACKAGE
        IF(IUNIT(62).GT.0) THEN
          CALL GWF2UPWBDS(KKSTP,KKPER,IGRID)
          CALL GWF2UPWBDCH(KKSTP,KKPER,IGRID)
          IBDRET=0
          IC1=1
          IC2=NCOL
          IR1=1
          IR2=NROW
          IL1=1
          IL2=NLAY
          DO IDIR=1,3
            CALL GWF2UPWBDADJ(KKSTP,KKPER,IDIR,IBDRET,
     &                        IC1,IC2,IR1,IR2,IL1,IL2,IGRID)
          ENDDO
        ENDIF
        IF(IUNIT(37).GT.0) THEN
          CALL GWF2HUF7BDS(KKSTP,KKPER,IGRID)
          CALL GWF2HUF7BDCH(KKSTP,KKPER,IUNIT(47),IGRID)
          IBDRET=0
          IC1=1
          IC2=NCOL
          IR1=1
          IR2=NROW
          IL1=1
          IL2=NLAY
          DO IDIR=1,3
            CALL GWF2HUF7BDADJ(KKSTP,KKPER,IDIR,IBDRET,
     &                         IC1,IC2,IR1,IR2,IL1,IL2,IUNIT(47),IGRID)
          ENDDO
        ENDIF
        IF(IUNIT(2).GT.0) CALL GWF2WEL7BD(KKSTP,KKPER,IUNIT(63),IGRID)
        IF(IUNIT(3).GT.0) CALL GWF2DRN7BD(KKSTP,KKPER,IGRID)
        IF(IUNIT(4).GT.0) CALL GWF2RIV7BD(KKSTP,KKPER,IGRID)
        IF(IUNIT(5).GT.0) THEN
          IF(IUNIT(22).GT.0.AND.NEVTOP.EQ.3) CALL GWF2LAK7ST(0,IGRID)
          CALL GWF2EVT7BD(KKSTP,KKPER,IGRID)
          IF(IUNIT(22).GT.0.AND.NEVTOP.EQ.3) CALL GWF2LAK7ST(1,IGRID)
        END IF
        IF(IUNIT(7).GT.0) CALL GWF2GHB7BD(KKSTP,KKPER,IGRID)
        IF(IUNIT(8).GT.0) THEN
          IF(IUNIT(22).GT.0.AND.NRCHOP.EQ.3) CALL GWF2LAK7ST(0,IGRID)
          CALL GWF2RCH7BD(KKSTP,KKPER,IGRID) 
          IF(IUNIT(22).GT.0.AND.NRCHOP.EQ.3) CALL GWF2LAK7ST(1,IGRID)
        END IF
        IF(IUNIT(16).GT.0) CALL GWF2FHB7BD(KKSTP,KKPER,IGRID)
        IF(IUNIT(17).GT.0) CALL GWF2RES7BD(KKSTP,KKPER,IGRID)
        IF(IUNIT(18).GT.0) CALL GWF2STR7BD(KKSTP,KKPER,IGRID)
        IF(IUNIT(19).GT.0) CALL GWF2IBS7BD(KKSTP,KKPER,IGRID)
        IF(IUNIT(39).GT.0) CALL GWF2ETS7BD(KKSTP,KKPER,IGRID)
        IF(IUNIT(40).GT.0) CALL GWF2DRT7BD(KKSTP,KKPER,IGRID)  
        IF(IUNIT(44).GT.0) CALL GWF2SFR7BD(KKSTP,KKPER,IUNIT(15),
     &                        IUNIT(22),IUNIT(46),IUNIT(55),NSOL,
     &                        IUNIT(8),IGRID)
        !Moved call to UZF1BD to follow SFR7BD for printing net recharge in UZF.
        IF(IUNIT(55).GT.0) CALL GWF2UZF1BD(KKSTP,KKPER,IUNIT(22),
     &                             IUNIT(44),IGRID)
        IF(IUNIT(22).GT.0) CALL GWF2LAK7BD(KKSTP,KKPER,IUNIT(15),
     &                       IUNIT(46),IUNIT(44),IUNIT(55),NSOL,IGRID)
        IF(IUNIT(50).GT.0) CALL GWF2MNW27BD(KKSTP,KKPER,IUNIT(62),
     &                                        IGRID)
        IF(IUNIT(52).GT.0) CALL GWF2MNW17BD(NSTP(KPER),KKSTP,KKPER,
     &                      IGRID)
        IF(IUNIT(54).GT.0) CALL GWF2SUB7BD(KKSTP,KKPER,IGRID)
        IF(IUNIT(57).GT.0) CALL GWF2SWT7BD(KKSTP,KKPER,IGRID)     
        IF(IUNIT(64).GT.0) CALL GWF2SWR7BD(KKSTP,KKPER,IGRID)  !SWR - JDH
  
C--SWM: SWAP POINTERS FOR LMT DATA TO CURRENT GRID
        !populate the variables - if SWAT-MODFLOW is used, this is handled by smrt_read_river2grid subroutine, !tcw
        if(readMFinput)then !tcw
          CALL SLMT7PNT(IGRID)
        end if !tcw

*       Write out head and flow data for transport model, if requested
        !if(ILMTFMT.eq.0 .or. ILMTFMT.eq.1) then
        
        
        !average Head data for easy output (rtb avg)
        !calculate average Head for the current month
        do k=1,nlay
          do j=1,ncol
            do i=1,nrow
              if(ibound(j,i,k).eq.1) then
                smrt_HNEW_tot_mo(j,i,k) = smrt_HNEW_tot_mo(j,i,k) + 
     &                                    HNEW(j,i,k)
                smrt_HNEW_avg_mo(j,i,k) = smrt_HNEW_tot_mo(j,i,k) / 
     &                                    day_count_mo
              endif
            enddo
          enddo
        enddo

        !calculate average Head for the current year
        do k=1,nlay
          do j=1,ncol
            do i=1,nrow
              if(ibound(j,i,k).eq.1) then
                smrt_HNEW_tot_yr(j,i,k) = smrt_HNEW_tot_yr(j,i,k) + 
     &                                    HNEW(j,i,k)
                smrt_HNEW_avg_yr(j,i,k) = smrt_HNEW_tot_yr(j,i,k) / 
     &                                    day_count_yr
              endif
            enddo
          enddo
        enddo


        if(rt_active.eq.1) then

        !Write a notification line to MODFLOW output file
        WRITE(IOUT,9876) IUMT3D,KKSTP,KKPER
 9876   FORMAT(/1X,'SAVING SATURATED THICKNESS AND FLOW TERMS ON UNIT',
     &   I5,' FOR MT3DMS',/1X,'BY THE LINK-MT3DMS PACKAGE V7',
     &   ' AT TIME STEP',I5,', STRESS PERIOD',I5/)

        !Collect and save all relevant flow model information
        IF(IUNIT(55).GT.0) 
     &    CALL LMT7UZF1(ILMTFMT,ISSMT3D,IUMT3D,KSTP,KPER,IGRID)
        IF(IUNIT(1) .GT.0) 
     &    CALL LMT7BCF7(ILMTFMT,ISSMT3D,IUMT3D,KKSTP,KKPER,IGRID)
        IF(IUNIT(23).GT.0) 
     &    CALL LMT7LPF7(ILMTFMT,ISSMT3D,IUMT3D,KKSTP,KKPER,IGRID)
        IF(IUNIT(62).GT.0) 
     &    CALL LMT7UPW1(ILMTFMT,ISSMT3D,IUMT3D,KKSTP,KKPER,IGRID)
        IF(IUNIT(37).GT.0) 
     &    CALL LMT7HUF7(ILMTFMT,ISSMT3D,IUMT3D,
     &    KKSTP,KKPER,IUNIT(47),IGRID)
        IF(IUNIT(2) .GT.0) 
     &    CALL LMT7WEL7(IUNIT(62),ILMTFMT,IUMT3D,KKSTP,KKPER,IGRID)
        IF(IUNIT(3) .GT.0) 
     &    CALL LMT7DRN7(ILMTFMT,IUMT3D,KKSTP,KKPER,IGRID)
        IF(IUNIT(8) .GT.0) 
     &    CALL LMT7RCH7(ILMTFMT,IUMT3D,KKSTP,KKPER,IGRID)
        IF(IUNIT(5) .GT.0) 
     &    CALL LMT7EVT7(ILMTFMT,IUMT3D,KKSTP,KKPER,IGRID)
        IF(IUNIT(39).GT.0) 
     &    CALL LMT7ETS7(ILMTFMT,IUMT3D,KKSTP,KKPER,IGRID) 
        IF(IUNIT(55).GT.0) 
     &    CALL LMT7UZFET(ILMTFMT,ISSMT3D,IUMT3D,KSTP,KPER,IGRID)
        IF(IUNIT(4) .GT.0)
     &    CALL LMT7RIV7(ILMTFMT,IUMT3D,KKSTP,KKPER,IGRID)
        IF(IUNIT(7) .GT.0) 
     &    CALL LMT7GHB7(ILMTFMT,IUMT3D,KKSTP,KKPER,IGRID)
        IF(IUNIT(18).GT.0) 
     &    CALL LMT7STR7(ILMTFMT,IUMT3D,KKSTP,KKPER,IGRID)
        IF(IUNIT(17).GT.0) 
     &    CALL LMT7RES7(ILMTFMT,IUMT3D,KKSTP,KKPER,IGRID)
        IF(IUNIT(16).GT.0) 
     &    CALL LMT7FHB7(ILMTFMT,IUMT3D,KKSTP,KKPER,IGRID)
        IF(IUNIT(52).GT.0) 
     &    CALL LMT7MNW17(ILMTFMT,IUMT3D,KKSTP,KKPER,IGRID)
        IF(IUNIT(50).GT.0) 
     &    CALL LMT7MNW27(ILMTFMT,IUMT3D,KKSTP,KKPER,IGRID)
        IF(IUNIT(40).GT.0) 
     &    CALL LMT7DRT7(ILMTFMT,IUMT3D,KKSTP,KKPER,IGRID)

        endif


        !Set NWT heads to Hdry when head is below bottom
        IF(IUNIT(63).GT.0) CALL GWF2NWT1BD()

        !Observation simulated equivalents
        CALL OBS2BAS7SE(IUNIT(28),IGRID)
          IF(IUNIT(33).GT.0) CALL OBS2DRN7SE(IGRID)
          IF(IUNIT(34).GT.0) CALL OBS2RIV7SE(IGRID)
          IF(IUNIT(35).GT.0) CALL OBS2GHB7SE(IGRID)
          IF(IUNIT(36).GT.0) CALL OBS2STR7SE(IGRID)
          IF(IUNIT(38).GT.0) CALL OBS2CHD7SE(KKPER,IGRID)
          IF(IUNIT(43).GT.0) CALL GWF2HYD7BAS7SE(1,IGRID)
          IF(IUNIT(43).GT.0 .AND. IUNIT(19).GT.0)
     1                              CALL GWF2HYD7IBS7SE(1,IGRID)
          IF(IUNIT(43).GT.0 .AND. IUNIT(54).GT.0)
     1                              CALL GWF2HYD7SUB7SE(1,IGRID)
          IF(IUNIT(43).GT.0 .AND. IUNIT(18).GT.0)
     1                              CALL GWF2HYD7STR7SE(1,IGRID)
          IF(IUNIT(43).GT.0 .AND. IUNIT(44).GT.0)
     1                              CALL GWF2HYD7SFR7SE(1,IGRID)

          !Print and/or save data
          CALL GWF2BAS7OT(KKSTP,KKPER,ICNVG,1,IGRID,BUDPERC)
          IF(IUNIT(19).GT.0) CALL GWF2IBS7OT(KKSTP,KKPER,IUNIT(19),
     1                                       IGRID)
          IF(IUNIT(37).GT.0)THEN
            IF(IOHUFHDS .NE.0 .OR.IOHUFFLWS .NE.0)
     1         CALL GWF2HUF7OT(KKSTP,KKPER,ICNVG,1,IGRID)
          ENDIF
          IF(IUNIT(51).NE.0) CALL GWF2MNW2I7OT(NSTP(KKPER),KKSTP,
     1                       KKPER,IGRID)
          IF(IUNIT(54).GT.0) CALL GWF2SUB7OT(KKSTP,KKPER,IUNIT(54),
     1                                       IGRID)
          IF(IUNIT(57).GT.0) CALL GWF2SWT7OT(KKSTP,KKPER,IGRID)
          IF(IUNIT(43).GT.0) CALL GWF2HYD7BAS7OT(KKSTP,KKPER,IGRID)

          !Jump to end of program if convergence was not achieved
          IF(IUNIT(63).GT.0) THEN
            IF(ICNVGFLG.EQ.0) THEN
              IF(ICNVG.EQ.0) GO TO 110
            ENDIF
          ELSE
            IF(ICNVG.EQ.0) THEN
              NCVGERR=NCVGERR+1
              WRITE(IOUT,87) BUDPERC
   87         FORMAT(1X,'FAILURE TO MEET SOLVER CONVERGENCE CRITERIA',/
     &       1X,'BUDGET PERCENT DISCREPANCY IS',F10.4)
              IF(ABS(BUDPERC).GT.STOPER) THEN
                WRITE(IOUT,*) 'STOPPING SIMULATION'
                GO TO 110
              ELSE
                WRITE(IOUT,*) 'CONTINUING EXECUTION'
              ENDIF
            ENDIF
          ENDIF
        ENDIF

        !If solute transport is simulated, call RT3D subroutines to simulate reactive transport
        !for the current flow time step. !rtb
        if(rt_active.eq.1) then
          call smrt_rt3d_run
        endif

      enddo !go to the next flow time step

*     If solution did not converge, write out statement and stop the program
  110 continue
      if(ICNVG.EQ.0) THEN
        write(*,*) 'FAILED TO MEET SOLVER CONVERGENCE CRITERIA'
        if(rt_active.eq.1) then
          call rt_close
        endif
        stop
      endif

*     Increment 'daycount'. This keeps track of how many days have been simulated
*     for the current stress period.
      daycount = daycount + mf_interval
      daytotal = daytotal + mf_interval

      !write out head data for observation cells
      if(mf_obs_flag.eq.1) then
        if(num_MF_obs.gt.0) then
          allocate(obs_heads(num_MF_obs))
          do i=1,num_MF_obs
            obs_heads(i) = HNEW(MF_obs(i,2),MF_obs(i,1),MF_obs(i,3))
          enddo
          write(30051,100) (obs_heads(i),i=1,num_MF_obs)
        endif
      endif

*     Determine if stress period information should be read before the next flow
*     time step.
      if(daycount.gt.perlen_day(1)) then
        read_stress = 1
        kper = kper + 1
        kkper = kper
        rt_kper = kper
        daycount = mf_interval
      else
        read_stress = 0
      endif
         
  100 format(1000(1pe13.5)) !spark modified format for high altitude
           
      return
      end
