*     read MODFLOW data for the groundwater flow simulation
      subroutine mf_read
C     Specify which modules are to be used
      USE GLOBAL
      USE GWFBASMODULE
      USE GWFHUFMODULE, ONLY:IOHUFHDS,IOHUFFLWS
      USE GWFEVTMODULE, ONLY:NEVTOP
      USE GWFRCHMODULE, ONLY:NRCHOP
      USE GWFLAKMODULE, ONLY:NLAKESAR,THETA,STGOLD,STGNEW,VOL
      USE GWFSFRMODULE, ONLY:NUMTAB
      USE GWFNWTMODULE, ONLY:LINMETH,ICNVGFLG,ITREAL
      !USE GWFLPFMODULE, ONLY:SC2
      USE GWFBCFMODULE, ONLY:SC2
      USE PCGMODULE
      USE SIPMODULE
      USE DE4MODULE
      USE GMGMODULE
      use GWFUPWMODULE, only: LAYTYPUPW !tcw
      use LMTMODULE, only: ISSMT3D,IUMT3D,ILMTFMT !tcw
      use smrt_parm
      INCLUDE 'mf_openspec.inc'

C     Assign version number and date
      CHARACTER*40 VERSION,VERSION2
      CHARACTER*10 MFVNAM
      PARAMETER (VERSION='1.0.5 05/14/2012')
      PARAMETER (VERSION2='1.8.00 12/18/2009')
      PARAMETER (MFVNAM='-NWT')
      CHARACTER*80 HEADNG(2)
      CHARACTER*200 FNAME

      ALLOCATE(CUNIT(NIUNIT))
      ALLOCATE(IBDT(8))

      CUNIT(1) = 'BCF6'
      CUNIT(2) = 'WEL '
      CUNIT(3) = 'DRN '
      CUNIT(4) = 'RIV '
      CUNIT(5) = 'EVT '
      CUNIT(6) = '    '
      CUNIT(7) = 'GHB '
      CUNIT(8) = 'RCH '
      CUNIT(9) = 'SIP '
      CUNIT(10) = 'DE4 '
      CUNIT(11) = '    '
      CUNIT(12) = 'OC  '
      CUNIT(13) = 'PCG '
      CUNIT(14) = 'lmg '
      CUNIT(15) = 'gwt '
      CUNIT(16) = 'FHB '
      CUNIT(17) = 'RES '
      CUNIT(18) = 'STR '
      CUNIT(19) = 'IBS '
      CUNIT(20) = 'CHD '
      CUNIT(21) = 'HFB6'
      CUNIT(22) = 'LAK '
      CUNIT(23) = 'LPF '
      CUNIT(24) = 'DIS '
      CUNIT(25) = '    '
      CUNIT(26) = 'PVAL'
      CUNIT(27) = '    '
      CUNIT(28) = 'HOB '
      CUNIT(29) = '    '
      CUNIT(30) = '    '
      CUNIT(31) = 'ZONE'
      CUNIT(32) = 'MULT'
      CUNIT(33) = 'DROB'
      CUNIT(34) = 'RVOB'
      CUNIT(35) = 'GBOB'
      CUNIT(36) = '    '
      CUNIT(37) = 'HUF2'
      CUNIT(38) = 'CHOB'
      CUNIT(39) = 'ETS '
      CUNIT(40) = 'DRT '
      CUNIT(41) = '    '
      CUNIT(42) = 'GMG '
      CUNIT(43) = 'HYD '
      CUNIT(44) = 'SFR '
      CUNIT(45) = '    '
      CUNIT(46) = 'GAGE'
      CUNIT(47) = 'LVDA'
      CUNIT(48) = '    '
      CUNIT(49) = 'LMT6'
      CUNIT(50) = 'MNW2'
      CUNIT(51) = 'MNWI'
      CUNIT(52) = 'MNW1'
      CUNIT(53) = 'KDEP'
      CUNIT(54) = 'SUB '
      CUNIT(55) = 'UZF '
      CUNIT(56) = 'gwm '
      CUNIT(57) = 'SWT '
      CUNIT(58) = 'cfp '
      CUNIT(59) = 'PCGN'
      CUNIT(60) = '    '
      CUNIT(61) = 'FMP '
      CUNIT(62) = 'UPW '
      CUNIT(63) = 'NWT '

C     Get the name of the name file
*      CALL GETNAMFIL(FNAME)
      INUNIT = 0
      MAXUNIT= INUNIT

C     Open name file
      FNAME = 'modflow.mfn' !this is required
      OPEN (UNIT=INUNIT,FILE=FNAME,STATUS='OLD',ACTION=ACTION(1))
      NC=INDEX(FNAME,' ')
      WRITE(*,490)' Using NAME file: ',FNAME(1:NC)
      print *  !rtb
  490 FORMAT(A,A)

C     Get current date and time, assign to IBDT, and write to screen
      CALL DATE_AND_TIME(VALUES=IBDT)

C     Allocate and read (BAS,DIS files)
      IGRID=1
      NSOL=1
      CALL GWF2BAS7AR(INUNIT,VERSION,24,31,32,MAXUNIT,IGRID,12,
     &                HEADNG,26,MFVNAM)
      IF(IUNIT(50).GT.0 .AND. IUNIT(52).GT.0) THEN
        WRITE(IOUT,'(1X,/,1X,A)')
     &  'MNW1 and MNW2 cannot both be active in the same simulation'
        CALL USTOP(' ')
      END IF
      IF(IUNIT(50).GT.0 .AND. IUNIT(52).GT.0) THEN
        WRITE(IOUT,'(1X,/,1X,A)')
     &  'MNW1 and MNW2 cannot both be active in the same simulation'
        CALL USTOP(' ')
      END IF
      allocate(mf_SC2(ncol,nrow,9)) !rtb
      IF(IUNIT(1).GT.0) CALL GWF2BCF7AR(IUNIT(1),IGRID)
      IF(IUNIT(23).GT.0) CALL GWF2LPF7AR(IUNIT(23),IGRID)
      IF(IUNIT(37).GT.0) CALL GWF2HUF7AR(IUNIT(37),IUNIT(47),
     &                                   IUNIT(53),IGRID)

!     Allocate a variable for MODFLOW if UPW package isn't used, tcw
      IF(IUNIT(62).eq.0 ) THEN !tcw
        ALLOCATE(LAYTYPUPW(NLAY)) !tcw
        LAYTYPUPW = 0 !tcw
      END IF !tcw
      
!     Allocate a few variables for MODFLOW since it's subroutine (LMT7BAS7AR) isn't called, tcw
      ALLOCATE(ISSMT3D,IUMT3D,ILMTFMT) !tcw
      ISSMT3D = 0 !tcw
      IUMT3D = 0 !tcw
      ILMTFMT = 7 !tcw

C     Allocate arrays for Newton Solver
      IF(IUNIT(63).GT.0) CALL GWF2NWT1AR(IUNIT(63),MXITER, 
     &                                   IUNIT(22),IGRID)
      IF(IUNIT(62).GT.0) CALL GWF2UPW1AR(IUNIT(62), Igrid)
      IF(IUNIT(2).GT.0) CALL GWF2WEL7AR(IUNIT(2),IUNIT(63),IGRID)
      IF(IUNIT(3).GT.0) CALL GWF2DRN7AR(IUNIT(3),IGRID)
      IF(IUNIT(4).GT.0) CALL GWF2RIV7AR(IUNIT(4),IGRID)
      IF(IUNIT(5).GT.0) CALL GWF2EVT7AR(IUNIT(5),IGRID)
      IF(IUNIT(7).GT.0) CALL GWF2GHB7AR(IUNIT(7),IGRID)
      IF(IUNIT(8).GT.0) CALL GWF2RCH7AR(IUNIT(8),IGRID)
      IF(IUNIT(16).GT.0) CALL GWF2FHB7AR(IUNIT(16),IGRID)
      IF(IUNIT(17).GT.0) CALL GWF2RES7AR(IUNIT(17),IGRID)
      IF(IUNIT(18).GT.0) CALL GWF2STR7AR(IUNIT(18),IGRID)
      IF(IUNIT(19).GT.0) CALL GWF2IBS7AR(IUNIT(19),IUNIT(54),IGRID)
      IF(IUNIT(20).GT.0) CALL GWF2CHD7AR(IUNIT(20),IGRID)
      IF(IUNIT(21).GT.0) CALL GWF2HFB7AR(IUNIT(21),IGRID)
      IF(IUNIT(44).GT.0) CALL GWF2SFR7AR(IUNIT(44),IUNIT(1),IUNIT(23),
     &                           IUNIT(37),IUNIT(15),NSOL,IOUTS,
     &                           IUNIT(62),IUNIT(55),IGRID)
      IF(IUNIT(55).GT.0) CALL GWF2UZF1AR(IUNIT(55),IUNIT(1),
     &                                   IUNIT(23),IUNIT(37),
     &                                   IUNIT(63),IGRID)
      IF(IUNIT(22).GT.0 .OR. IUNIT(44).GT.0) CALL GWF2LAK7AR(
     &             IUNIT(22),IUNIT(44),IUNIT(15),IUNIT(55),NSOL,IGRID)
      IF(IUNIT(46).GT.0) CALL GWF2GAG7AR(IUNIT(46),IUNIT(44),
     &                                     IUNIT(22),IGRID)
      IF(IUNIT(39).GT.0) CALL GWF2ETS7AR(IUNIT(39),IGRID)
      IF(IUNIT(40).GT.0) CALL GWF2DRT7AR(IUNIT(40),IGRID)
      IF(IUNIT(9).GT.0) CALL SIP7AR(IUNIT(9),MXITER,IGRID)
      IF(IUNIT(10).GT.0) CALL DE47AR(IUNIT(10),MXITER,IGRID)
      IF(IUNIT(13).GT.0) CALL PCG7AR(IUNIT(13),MXITER,IGRID)
      IF(IUNIT(42).GT.0) CALL GMG7AR(IUNIT(42),MXITER,IGRID)
      IF(IUNIT(50).GT.0) CALL GWF2MNW27AR(IUNIT(50),IGRID)
      IF(IUNIT(51).GT.0) CALL GWF2MNW2I7AR(IUNIT(51),IUNIT(50),IGRID)
      IF(IUNIT(52).GT.0) CALL GWF2MNW17AR(IUNIT(52),IUNIT(9),
     &                     IUNIT(10),IUNIT(63),0,IUNIT(13),
     &                     0,IUNIT(42),FNAME,IGRID)
      IF(IUNIT(43).GT.0) CALL GWF2HYD7BAS7AR(IUNIT(43),IGRID)
      IF(IUNIT(43).GT.0 .AND. IUNIT(18).GT.0)
     &                   CALL GWF2HYD7STR7AR(IUNIT(43),IGRID)
      IF(IUNIT(43).GT.0 .AND. IUNIT(44).GT.0)
     &                   CALL GWF2HYD7SFR7AR(IUNIT(43),IGRID)

C     Observation allocate and read
      CALL OBS2BAS7AR(IUNIT(28),IGRID)
      IF(IUNIT(33).GT.0) CALL OBS2DRN7AR(IUNIT(33),IUNIT(3),IGRID)
      IF(IUNIT(34).GT.0) CALL OBS2RIV7AR(IUNIT(34),IUNIT(4),IGRID)
      IF(IUNIT(35).GT.0) CALL OBS2GHB7AR(IUNIT(35),IUNIT(7),IGRID)
      IF(IUNIT(38).GT.0) CALL OBS2CHD7AR(IUNIT(38),IGRID)

C     Modify conductance for HFB when using UPW.
      IF(IUNIT(62).GT.0 ) THEN
        IF(IUNIT(21).GT.0) CALL GWF2HFB7UPW(IGRID)
      END IF
      
C     Initialize the stress period
      kper = 1
      kkper = kper
      rt_kper = kper

C     Initialize daycount (for mf_read subroutine) and read_stress flag
      daycount = mf_interval
      daytotal = mf_interval
      read_stress = 1

      !Read in observation cells (for output at each time step)
      num_MF_obs = 0
      MF_obs = 0
      if(mf_obs_flag.eq.1) then
        open(30050,file='modflow.obs')
        read(30050,*)
        read(30050,*) num_MF_obs
        if(num_MF_obs.gt.0) then
          allocate(MF_obs(num_MF_obs,3))
          do i=1,num_MF_obs
            !1 = row, 2 = column, 3 = layer
            read(30050,*) (MF_obs(i,j),j=1,3)
          enddo
          open(30051,file='swatmf_out_MF_obs')
          write(30051,*) 'Head values for selected MODFLOW cells'
        endif
      endif

      return 
      end
