*     This file contains the modules for RT3D, in order to allocate variables dynamically

*     Modules contained in RT3D
*     GLOBAL
*     VST_ARRAYS !for VST package

*     Global module -------------------------------------------------------------------------------
      module rt_global
      
      !integer variables      
      integer:: LENX,LENIX,MXSTP,MXOBS,MXCOMP,MXPTS,MXFTIMES,
     &          MXCOL,MXROW,MXLAY,MXCROP,MXYR,MXCYCLE,MXPROC,
     &          INBTN,INADV,INDSP,INSSM,INRCT,INGCG,INHFF,
     &          OUTRES,FNAMES,
     &          IOUT,ICONC,ICBIN,ISSMASS,IOBS,IMAS,IPLOT,ISPC,
     &          INUMSPC,
     &          IXSECT,IXSBIN,ILCH,ILCHOBS,ILCHBIN,IMASSOUT,IMAPOUT,
     &          I1,I2,I3,
     &          NTRANS,NRC
      integer:: ISUMX,ISUMIX,ISUM,ISUM2,NCOMP,MCOMP,
     &          LCLAYC,LCDELR,LCDELC,LCDZ,LCPR,LCXBC,LCYBC,LCZBC,LCQX,
     &          LCQY,LCQZ,LCDH,LCIB,LCCOLD,LCCNEW,LCCADV,LCRETA,LCBUFF,
     &          MIXELM,MXPART,LCXP,LCYP,LCZP,LCCNPT,LCCHEK,
     &          LCAL,LCTRPT,LCTRPV,LCDM,LCDXX,LCDXY,
     &          LCDXZ,LCDYX,LCDYY,LCDYZ,LCDZX,LCDZY,LCDZZ,LCSSMC,
     &          LCIRCH,LCRECH,LCCRCH,LCIEVT,LCEVTR,LCCEVT,LCSS,
     &          ISOTHM,IREACT,isolver,NCRXNDATA,NVRXNDATA,
     &          LCRHOB,LCSP1,LCSP2,LCRC1,LCRC2,INTERP,
     &          ITRACK,NPL,NPH,NPMIN,NPMAX,NPLANE,NLSINK,NPSINK,
     &          NPRS,NOBS,NSS,KFLOW,NTSS,NPS,
     &          IFMTCN,IFMTNP,IFMTRF,IFMTDP,MXSTRN,
     &          NPER,NSTP,ISTAT,LCQSTO,LCHTOP,LCCWGT,LCSR,
     &          LCINDX,LCINDY,LCINDZ,ISS,IVER,NPROBS,NPRMAS,IRCTOP,
     &          RT_MXITER,IPRGCG,NADVFD,ITP,NODES,ICNVG,RT_ITER1,ITO,
     &          ISOLVE,LCA,LCQ,LCWK,LCCNCG,LCLRCH,LCRHS,LCIBND,LCCNEWN,
     &          IMPSOL,NCRS,ISPD,IGETSC,ICOMP,NPERFL,
     &          NIRRCELL,NWELL,NSCRNLYRS,NSAMPLE,PUMPSIM,PHASE(100),
     &          WRITEFLAG,WRITE_OUT(100)
      integer:: sum_tracker,rt_year
      integer*4::NUMARGS,RESULT,LEN

      !real variables      
      real::    PERCEL,HORIGN,XMAX,YMAX,ZMAX,CINACT,
     &          DCEPS,SRMULT,WD,DCHMOC,T1,T2,TIME1,
     &          TIME2,DT0,DTRACK,DTDISP,DT,THKMIN,
     &          DTSSM,DTRCT,DTRACK2,RFMIN,RT_ACCL,CCLOSE,
     &          TTSMULT,TTSMAX,SIMLEN,TIME,
     &          time_begin,time_end,time_total,
     &          DM_TOTAL,DM_SS,SOURCE_TOT,SINK_TOT,INITMASS,NEWMASS
      
      !arrays
      integer,allocatable,save :: LOCOBS(:,:),NCOUNT(:),NPINS(:)
      real,allocatable,save :: TIMPRS(:),TSLNGH(:),TMASIO(:,:,:),
     &                         RMASIO(:,:,:),TMASS(:,:),TMASIN(:),
     &                         TMASOUT(:),ERROR(:),ERROR2(:),
     &                         ERROR2_TS(:),SSCUMUL(:,:)
      
      !arrays for implicit solver
      integer,allocatable,save :: IX(:)
      real,allocatable,save :: X(:)

      !other variables      
      double precision,allocatable,save::RC(:),VRC(:,:,:,:)
      logical:: TRNOP(10),UNIDX,UNIDY,UNIDZ,SAVUCN,SAVCBM,CHKMAS,
     &          PRTOUT,UPDLHS,
     &          FAPP,FCNC
      character:: FLTYPE*4,FLNAME*80,FINDEX*30,TUNIT*4,LUNIT*4,MUNIT*4,
     &          FPRT*1,PATH*80,CONCFILE*80,CONCNAME*80,
     &          PREFIX*80,ARGV*80
      character,allocatable,save::SPECIES(:)*80

C     Arrays for general concentration calculations
        !Integer Arrays (GENERAL)      
        INTEGER, ALLOCATABLE, SAVE :: LAYCON(:),ICBUND(:,:,:,:)
        !Integer Arrays (ADV)      
        INTEGER, ALLOCATABLE, SAVE :: NPCHEK(:,:,:,:)
        INTEGER, ALLOCATABLE, SAVE :: INDEXX(:,:),INDEXY(:,:)
        INTEGER, ALLOCATABLE, SAVE :: INDEXZ(:,:)
        !Real Arrays (GENERAL)
        REAL, ALLOCATABLE, SAVE :: HTOP(:,:),DF(:,:,:)
        REAL, ALLOCATABLE, SAVE :: CELLTOP(:,:,:),DZ(:,:,:)
        REAL, ALLOCATABLE, SAVE :: XBC(:),YBC(:),ZBC(:,:,:)
        REAL, ALLOCATABLE, SAVE :: CNEW(:,:,:,:),COLD(:,:,:,:) 
        REAL, ALLOCATABLE, SAVE :: CWGT(:,:,:,:),CADV(:,:,:,:)
        REAL, ALLOCATABLE, SAVE :: CSR(:,:,:,:)      
        REAL, ALLOCATABLE, SAVE :: RETA(:,:,:,:),PRSITY(:,:,:)
        !Real Arrays (ADV)
        REAL, ALLOCATABLE, SAVE :: XP(:,:),YP(:,:),ZP(:,:),CNPT(:,:,:)
        !Real Arrays (DSP)
        REAL, ALLOCATABLE, SAVE :: ALPHAL(:,:,:),TRPT(:),TRPV(:),
     &                             DMCOEF(:)
        REAL, ALLOCATABLE, SAVE :: DXX(:,:,:),DXY(:,:,:),DXZ(:,:,:)
        REAL, ALLOCATABLE, SAVE :: DYX(:,:,:),DYY(:,:,:),DYZ(:,:,:) 
        REAL, ALLOCATABLE, SAVE :: DZX(:,:,:),DZY(:,:,:),DZZ(:,:,:)
        !Real Arrays (SSM)
        REAL, ALLOCATABLE, SAVE :: CEVT(:,:,:),CCNC(:,:,:,:)
        !Real Arrays (RCT)
        REAL, ALLOCATABLE, SAVE :: RHOB(:,:,:),SP1(:,:,:,:),SP2(:,:,:,:)
        !Real Arrays (shale, reaeration processes)
        REAL, ALLOCATABLE, SAVE :: MATSE(:)
        REAL, ALLOCATABLE, SAVE :: MATTYPE(:,:,:,:)
        INTEGER, ALLOCATABLE, SAVE :: SESOURCE(:)   

C       Source/Sink concentration decay boundary condition
        REAL, ALLOCATABLE, SAVE,DIMENSION(:,:,:,:) :: ccnh,cghb

C       Leaching arrays
        REAL, ALLOCATABLE, SAVE :: TLEACH(:,:,:,:),LEACH(:,:,:,:)

C       Species arrays
        LOGICAL, ALLOCATABLE, SAVE :: NTR(:)
        INTEGER, ALLOCATABLE, SAVE :: NRCSPECIES(:),NVCSPECIES(:)
        INTEGER, ALLOCATABLE, SAVE :: NTRSPECIES(:)
        INTEGER, ALLOCATABLE, SAVE :: RCINDEX(:,:),VCINDEX(:,:),SPI(:,:)
        REAL, ALLOCATABLE, SAVE :: RETASED(:,:,:,:),RETAOM(:,:,:,:)    

      endmodule rt_global


C     Arrays for VST package ----------------------------------------------------------------------
      MODULE VST_ARRAYS
        REAL, ALLOCATABLE, SAVE :: MCOLD(:,:,:)
        REAL, ALLOCATABLE, SAVE :: SATOLD(:,:,:),SATNEW(:,:,:)
        REAL, ALLOCATABLE, SAVE :: WATVOL(:,:,:),QZVAR(:,:,:)
        REAL, ALLOCATABLE, SAVE :: WT(:,:)
        REAL, ALLOCATABLE, SAVE :: APPL(:,:)
        REAL, ALLOCATABLE, SAVE :: CAPPL(:,:,:)
        REAL, ALLOCATABLE, SAVE :: CUZET(:,:,:,:), CGWET(:,:,:,:)
        REAL, ALLOCATABLE, SAVE :: MNEW(:,:,:,:),MOLD(:,:,:,:),
     &                             MTEMP(:,:,:,:)
      ENDMODULE VST_ARRAYS