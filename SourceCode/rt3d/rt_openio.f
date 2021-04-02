C**************************  OPENIO ************************************
C
C  open the input and output files
C
C***********************************************************************

      subroutine rt_openio

      use rt_global

      implicit none
      
      INTEGER      I,NSPECIES
      CHARACTER*30 FN1,FN2,FN3,FN4,FN5,FN6,FN7,FN8,FN9,FN10
      CHARACTER*30 FN20
      CHARACTER*30 COMP

*     Set the integer values for the input/output files
      !input files
      FNAMES=8000
      INBTN=8001
      INADV=8002
      INDSP=8003
      INSSM=8004
      INRCT=8005
      INGCG=8006
      INHFF=8007
      OUTRES=8010
      !output files
      IOUT=8020
      ISPC=8023
      ICONC=8200
      ICBIN=8250
      ISSMASS=8300
      IOBS=8400
      IMAS=8450
      IPLOT=8500
      IXSECT=8550
      IXSBIN=8600
      ILCH=9000
      ILCHBIN=9100
      ILCHOBS=9200
      IMASSOUT=9300

*     Open up the file containing the input and output file names      
      OPEN(FNAMES,FILE='rt3d_filenames',STATUS='OLD')

C	Read the input file names and the main output file, and open them up
C     BTN = Basic Transport: record 1
C     ADV = Advection: record 2
C     DSP = Dispersion: record 3
C     SSM = Source/Sink Mixing: record 4
C     RCT = Chemical Reaction: record 5
C     GCG = General Conjugate Gradient: record 6 
C     HHF = Flow simulation linker file (binary: unformatted file): record 7
C     OUT = General output file: record 20
     
*     Read the names of the files from the RT3DAG.fnames file      
      READ(FNAMES,*) FN1
      READ(FNAMES,*) FN2
      READ(FNAMES,*) FN3
      READ(FNAMES,*) FN4
      READ(FNAMES,*) FN5
      READ(FNAMES,*) FN6
c      READ(FNAMES,*) FN7 !linker file is not needed!
      READ(FNAMES,*) FN10

      OPEN(INBTN, FILE=FN1, STATUS='OLD')
      OPEN(INADV, FILE=FN2, STATUS='OLD')
      OPEN(INDSP, FILE=FN3, STATUS='OLD')
      OPEN(INSSM, FILE=FN4, STATUS='OLD')
      OPEN(INRCT, FILE=FN5, STATUS='OLD')
      OPEN(INGCG, FILE=FN6, STATUS='OLD')
c      OPEN(INHFF, FILE=FN7, STATUS='OLD')
      OPEN(OUTRES, FILE=FN10, STATUS='UNKNOWN')
      
*     Read (from the BTN file) the packages that are used for this simulation
      READ(INBTN,*)
      DO I=1,10
        TRNOP(I)=.FALSE.
      ENDDO
      READ(INBTN,*) (TRNOP(I),I=1,10)    

*     Open up file for CPU timing information
      OPEN(5000,FILE='CPU')

      OPEN(IOUT,FILE='swatmf_out_RT_out',STATUS='UNKNOWN')
	  
      return
      end
