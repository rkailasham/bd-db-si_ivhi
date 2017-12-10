C     STEADY STATE VISCOSITY AND FIRST 
C     NORMAL STRESS DIFFERENCE COEFFICIENT
C     FOR DUMBBELLS WITH FENE + INTERNAL VISCOSITY
c     + EXCLUDED VOLUME INTERACTIONS
C     SECOND ORDER SEMI-IMPLICIT PREDICTOR-CORRECTOR
C     ALGORITHM
C     30 MAY 2001
C     10TH MARCH 2017 : MODIFIED SO THAT
C     INITIAL CONFIGS ARE SAMPLED FROM A 
C     GAUSSIAN DISTRIBUTION
C     CONNECTOR VECTOR CALCULATED USING DEFINITION OF "MOLECULAR
C     "EXTENSION" GIVEN IN HUR AND SHAQFEH'S PAPER[J.RHEOL 2000]


C     MODIFIED ON 29TH APRIL 2017 : HF = AMPL = 1./(1-QL/B)

C     MODIFIED ON 2ND MAY 2017 : INCLUDE EV TERMS

c     MODIFIED ON 8TH JULY 2017 : CHANGING DEFINITION OF Q2
C
c     MODIFIED ON 7TH AUG 2017 : CHANGE DEFINITION OF HEAUX
C     TO INCLUDE AMPL2
C
c     25-SEP-2017 : FOR FENE DUMBBELLS, SAMPLING FROM A 
C     GAUSSIAN.
C     
c     ALSO ENSURING THAT EQUILIBRATION AND PRODUCTION STEPS
c     ARE DONE USING DIFFERENT TIME-STEP WIDTHS
c
C     HAVE ALSO IMPLEMENTED STRESS-TENSOR EXPRESSIONS
C
c     26-SEP-2017 : REALIZED THAT THE IDEA OF DIFFERENT TIMESTEPS
c     HAVE TO BE INCORPORATED IN DTH, DTQ, C1P AND C2P AS WELL
c
c     16-OCT-2017 : IMPLEMENTED REGULARIZED OSEEN BURGER
C
c     10-NOV-2017 : PROPERLY IMPLEMENTING DOUBLE PRECISION BY  
c     PUTTING D0 AT THE END OF ALL DECIMAL NUMBERS.
c     ALSO USING HIGHER NUMBER OF DIGITS FOR FOURPI
c     ALSO USING HIGHER NUMBER OF DIGITS FOR PI
c
c     NOTE : CHANGE SO THAT EQBN AND PRODUCTION ARE DONE USING THE 
c     SAME T-STEP WIDTHS
c
c     13-NOV-2017 : CHANGE Q SO THAT IT IS STORED AS VECTOR : Q(1),Q(2),Q(3)
c     RATHER THAN 3 SEPARATE NUMBERS Q1,Q2,Q3 
c
c     NOTE : WHEN RUNNING EQUILIBRIUM SIMULATIONS, DIVERT PROGRAM
c     FLOW FROM ENTERING TEXTRA. WHEN SR=0.0, MATERIAL FUNCTIONS BECOME UNDEFINED.
C     THIS CREATES HAVOC WITH TEXTRA
C
c     15-NOV-2017 : MODIFYING CODE SO THAT IT ACCEPTS INITIAL CONDITIONS FROM A
C     DATABASE OF EQUILIBRATED FENE DUMBBELL CONFIGURATIONS
C
c     NOTE : FPATH CONTAINS THE FULL PATH TO eqbconfigs.dat WHICH CONTAINS THE 
c     EQUILIBRATED CONFIGURATIONS. CHANGE THE PATH NAME BY CHANGING FPATH
c
C     16-NOV-2017 : TRYING TO IMPLEMENT RANDOM ACCESS OF eqbconfigs.dat.
c     DEFAULT IS SEQUENTIAL ACCESS
c
c     17-NOV-2017 : IMPLEMENTED RANDOM ACCESS BY LOADING DATABSE ON TO AN 
c     ARRAY. PROGRAM EXECUTION OVERHEAD OF 30s.
c     CHANGE FPATH ACCORDING TO PLATFORM
c
c     18-NOV-2017 : CHANGING PARAMETERS BEING WRITTEN TO TEXTRA OUTPUT
c
c     10-DEC-2017 : Writing configurations at the end of simulation run;
c     so that I can continue where I left off
c



      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NDATM=50,NOUT=10) 
      PARAMETER (NBINS=50)
      PARAMETER (NDB=10000000)
      REAL*8 FENFAC,TEMPB
      REAL*8 MPSI1,MPSI2
      REAL*8 AQ2,MQ2,VQ2
      REAL*8 QMAG,QLEN,TEMP12,TEMP23
      REAL*8 METE,METD,MET
      REAL*8 META
      REAL*8 AVQ2(NOUT),MEQ2,ERQ2(NOUT)
      REAL*8 AVTE(NOUT),ERTE(NOUT)
      REAL*8 AVTD(NOUT),ERTD(NOUT)
      REAL*8 AVT(NOUT),ERT(NOUT)
      REAL*8 AVETA(NOUT),ERETA(NOUT)
      REAL*8 AVPSI1(NOUT),ERPSI1(NOUT)
      REAL*8 AVPSI2(NOUT),ERPSI2(NOUT)
      REAL*8 DTARR(10) 
      REAL*8 XARR(NDATM),YARR(NDATM),SIGARR(NDATM) 
      REAL*8 TEMP1(NDATM, NOUT),TEMP2(NDATM, NOUT),TEMP3(NDATM, NOUT)
      REAL*8 OUTIME(NOUT)
      REAL*8 Q(3)
      REAL*8 DB(NDB,3)
      INTEGER NTIARR(10) 
      INTEGER TIME (8)
      CHARACTER (LEN=12) CLK(3)
      CHARACTER (LEN=512) FPATH
      CHARACTER (LEN=512) CPATH
      PARAMETER(FPATH="/home/rkailasham/Sdrive/Desktop/Kailash_Files/Com
     &bined_Code_Validation/equilibrated_database/eqbconfigs.dat")
      PARAMETER(CPATH="/home/rkailasham/Sdrive/Desktop/Kailash_Files/fin
     &configs.dat")
      COMMON /STEPBL/ THI,B,ZMU,RMU2,BAUXQ,BAUXR,DTH,DTQ,SRDT,
     &SRDTH,C1P,C2P,E,AMPL2,BETA,AB,A2,A4,AUX1,AUX2,
     &AUX3,AUX4,AUX5,S,QALPH,GEE1,GEE2,GEE3,GEE4 
      COMMON /EXTRP/ NOPT,NDUOPT,XOPT,YOPT,SIGOPT,ALIOPT,VLIOPT 
      LOGICAL THERE
      LOGICAL CTHERE
c      INTEGER :: THI
C     Excluded-volume and FENE parameters 
c     Z = Solvent Quality, RMU = EV parameter, B= FENE parameter
c     SR = Shear Rate, E = IV PARAMETER,
C     H0 = HYDRODYNAMIC INTERACTION PARAMETER
C     THI = EXPRESSION TO 
C     BE USED OFR THE HYDRODYNAMIC INTERACTION TENSOR.
C     THI = 1 FOR REGULARIZED OSEEN-BURGER
C     THI = 2 FOR ROTNE-PRAGER-YAMAKAWA 
      open (unit=1, file='inp.dat') 
      READ (1,*) Z, RMU, B, SR,E,H0,THI,INPAR 
      open (unit=2, file='tstep.dat') 
      open (unit=3, file='eta.dat',STATUS='UNKNOWN') 
      open (unit=4, file='psi1.dat',STATUS='UNKNOWN')
      open (unit=5, file='psi2.dat',STATUS='UNKNOWN')
      open (unit=7, file='q2.dat', STATUS='UNKNOWN')
      open (unit=8, file='tau.dat',STATUS='UNKNOWN')
      open (unit=9, file='tau_e.dat',STATUS='UNKNOWN')
      open (unit=10, file='tau_d.dat',STATUS='UNKNOWN')
      INQUIRE (FILE=FPATH,EXIST=THERE)
      INQUIRE (FILE=FPATH,EXIST=CTHERE)
      open (unit=112,file='finconfigs.dat',STATUS='UNKNOWN')
C      open (unit=113, file='eqbconfigs150DT.dat',STATUS='UNKNOWN')
c      OPEN (UNIT=89, file='bin_data.dat',STATUS='UNKNOWN')

C     RES_TIME DECIDES OUTIME BASED ON IF RUNS RESUME FROM PREVIOUS CONFIGS OR
C     EQB DATABASE 
C      RES_TIME=0.D0
      TMAX=10.D0
      READ (2,*) NTIWID, NTRAJ
      DO 15 I = 1, NTIWID
          READ (2,*) NTIARR(I)
          DTARR(I)=TMAX/NTIARR(I)
15    CONTINUE

C 
      ZMU=Z/RMU**(5.D0) 
      RMU2=(2.D0)*RMU*RMU
      SRTEMP=SR
      NTHI=NINT(THI) 
C 
      SELECT CASE (NTHI)
          CASE (1)
              WRITE(*,*) "THI : ",THI
              WRITE(*,*) "TYPE OF HI : ROBT"
              IF(H0.EQ.0) THEN
                  WRITE(*,*) "H0=0; SETTING THI=3"
                  THI=3.D0
              ENDIF
          CASE (2)
              WRITE(*,*) "THI : ",THI
              WRITE(*,*) "TYPE OF HI : RPY"
              IF(H0.EQ.0) THEN
                  WRITE(*,*) "H0=0; SETTING THI=3"
                  THI=3.D0
              ENDIF
          CASE DEFAULT
              WRITE(*,*) "HI OPTION NOT VALID"
              WRITE(*,*) "SETTING H0=0.0"
              H0=0.D0
      END SELECT

      CALL CPU_TIME(STARTTIME)

      SELECT CASE (INPAR)
          CASE (1)
              WRITE(*,*) "INPAR : ",INPAR
              IF(THERE)THEN
                  OPEN(UNIT=114,file=FPATH)
                  WRITE(*,*) "LOADING FROM EQUILIBRATED DATABASE.."
                  DO 12 I=1,NDB
                      READ(114,8) K,DB(I,1),DB(I,2),DB(I,3)
12                CONTINUE
                  CLOSE(UNIT=114)
              ELSE
                  STOP "eqbconfigs.dat file not found. Execution 
     &terminated"
              ENDIF  
              ISEED=20171113              
              WRITE(*,*) "LOADED DATABASE"
          CASE (2)
              WRITE(*,*) "INPAR : ",INPAR
              IF(CTHERE)THEN
                  OPEN(UNIT=115,file=CPATH)
                  WRITE(*,*) "READING FROM FINAL CONFIGS OF PREVIOUS 
     &RUN.."
                  DO 13 I=1,(3*NTRAJ)
                      READ(115,9) TIMEI,DB(I,1),DB(I,2),DB(I,3)
13                CONTINUE
                  READ(115,*) ISEED
                  CLOSE(UNIT=115)
              ELSE
                  STOP "finconfigs.dat file not found. Execution 
     &terminated"
              ENDIF
          CASE DEFAULT
              WRITE(*,*) "READ_FROM_INPUT OPTION NOT VALID"
              STOP "USE INPAR=1 FOR EQB DBASE, INPAR=2 FOR RESUMING 
     &FROM PREVIOUS RUN"
      END SELECT








      PI=3.1415926535897931D0
      AB=SQRT(PI)*H0
      A2=AB*AB
      A4=A2*A2
45    FORMAT(F20.16)
      WRITE(*,45) E
      WRITE(*,45) H0


   
C     Loop for different time step widths 
      CALL SRAND(ISEED)



8     FORMAT(I10,4X,F20.16,4X,F20.16,4X,F20.16)




c      CALL CPU_TIME(STARTTIME)
      DO 1000 IDT=1,NTIWID
C        Auxiliary parameters 
c         OPEN(UNIT=114,file=FPATH)
         WRITE(*,*) "DOING TIMESTEP WIDTH : ",DTARR(IDT)
         DELTAT=DTARR(IDT) 
         TEMPDT=DELTAT
         NTIME=NTIARR(IDT)/NOUT
         DTH=0.5D0*DELTAT 
         DTQ=0.25D0*DELTAT 
         SQDT=SQRT(DELTAT) 
         C1P=14.14855378D0*SQDT 
         C2P=1.21569221D0*SQDT 
         



C     Initializing averages and errors... 
C     AV... ARE AVERAGES OF CORRESPONDING MATL. FNS. OVER
C     ALL THE BLOCKS; ER... ARE STANDARD ERRORS FOR THE CORR-
C     -ESPONDING MATL. FNS.
         DO 25 I=1,NOUT
             AVQ2(I)=0.D0
             ERQ2(I)=0.D0
             AVTE(I)=0.D0
             ERTE(I)=0.D0
             AVTD(I)=0.D0
             ERTD(I)=0.D0
             AVT(I)=0.D0
             ERT(I)=0.D0
             AVETA(I)=0.D0
             ERETA(I)=0.D0
             AVPSI1(I)=0.D0
             ERPSI1(I)=0.D0
             AVPSI2(I)=0.D0
             ERPSI2(I)=0.D0
25       CONTINUE









C        A FRESH SEED IS GIVEN FOR EACH TIME-STEP WIDTH
C        TO ENSURE NON-OCCURENCE OF PERIOD EXHAUSTION

c         CALL DATE_AND_TIME(CLK(1),CLK(2),CLK(3),TIME)
c         ISEED=TIME(8)*100000+TIME(7)*1000+TIME(6)*10+TIME(5)
C        write (*,*) ISEED 
C         ISEED=ISEED+1         
         CALL RANILS(ISEED)

         DO 100 ITRAJ=1,NTRAJ 
        
             SR=0.D0
             SRDT=SR*DELTAT
             SRDTH=0.5D0*SRDT

C            Initial conditions taken from memory  
c             READ(114,8) J,Q(1),Q(2),Q(3)

C              
             PICK=RAND()
             NSEED=NSEED+1
             CHANGE=NDB*PICK
             IF(INPAR.EQ.1)THEN
                 NCHOOSE=NINT(CHANGE)
             ELSE
                 NCHOOSE=((IDT-1)*NTRAJ)+(ITRAJ)
             ENDIF 
             Q(1)=DB(NCHOOSE,1)
             Q(2)=DB(NCHOOSE,2)
             Q(3)=DB(NCHOOSE,3)



             IF(MODULO(ITRAJ,1000).EQ.0)THEN
             WRITE(*,*) "STATUS : EQB.TIME-STEP WIDTH : ",DELTAT,
     &"TRAJ # ",ITRAJ
             ENDIF
c             WRITE(*,*) "EQUILIBRATION AT SR = ",SR 
C            Relaxation of initial conditions
             NEQB=1.D0/DELTAT
             DO 50 ITIME=1,NEQB
                CALL SEMIMP(Q)
50           CONTINUE
C 
             SR=SRTEMP
             SRDT=SR*DELTAT
             SRDTH=0.5D0*SRDT

c             WRITE(*,*) "PRODUCTION AT SR = ",SR 
             IWAIT=0
             IOUT=0
C
C           Time integration: semi-implicit predictor-corrector scheme 
             DO 10 ITIME=1,NTIARR(IDT) 
                 CALL SEMIMP(Q) 
                 IWAIT=IWAIT+1
                 IF (IWAIT.EQ.NTIME) THEN
                     IWAIT=0
                     IOUT=IOUT+1
                     QMAG =Q(1)*Q(1) + Q(2)*Q(2) + Q(3)*Q(3)
                     TEMP12=(Q(1)*Q(1) - Q(2)*Q(2))
                     TEMP23=(Q(2)*Q(2) - Q(3)*Q(3))
                     AVQ2(IOUT)=AVQ2(IOUT)+(QMAG)
                     ERQ2(IOUT)=ERQ2(IOUT)+(QMAG*QMAG)
                     METD=(2.D0)*E*GEE3*SR*(Q(1)**(2.D0))*
     &(Q(2)**(2.D0))/QMAG
                     METE=((GEE3*Q(1)*Q(2))/(1.D0-(QMAG/B)))+
     &(E*GEE4*Q(1)*Q(2)/QMAG)
                     MET=METE+METD
                     META=MET/SR

                     MPSI1=((GEE3*TEMP12)/(1.D0-(QMAG/B)))+
     &(2.D0*E*GEE3*SR*Q(1)*Q(2)*TEMP12/QMAG)+(E*GEE4*TEMP12/QMAG)
                     MPSI1=(MPSI1)/(SR*SR)

                     MPSI2=((GEE3*TEMP23)/(1.D0-(QMAG/B)))+
     &(2.D0*E*GEE3*SR*Q(1)*Q(2)*TEMP23/QMAG)+(E*GEE4*TEMP23/QMAG)
                     MPSI2=(MPSI2)/(SR*SR)

                     AVTE(IOUT)=AVTE(IOUT)+METE
                     ERTE(IOUT)=ERTE(IOUT)+(METE*METE)
                     AVTD(IOUT)=AVTD(IOUT)+METD
                     ERTD(IOUT)=ERTD(IOUT)+(METD*METD)
                     AVT(IOUT)=AVT(IOUT)+MET
                     ERT(IOUT)=ERT(IOUT)+(MET*MET)
                     AVETA(IOUT)=AVETA(IOUT)+META
                     ERETA(IOUT)=ERETA(IOUT) + (META*META)
                     AVPSI1(IOUT)=AVPSI1(IOUT)+MPSI1
                     ERPSI1(IOUT)=ERPSI1(IOUT)+(MPSI1*MPSI1)
                     AVPSI2(IOUT)=AVPSI2(IOUT)+MPSI2
                     ERPSI2(IOUT)=ERPSI2(IOUT)+(MPSI2*MPSI2)
                 ENDIF
10           CONTINUE 

             WRITE(112,9) DELTAT,Q(1),Q(2),Q(3)
9            FORMAT(F11.8,4X,F20.16,4X,F20.16,4X,F20.16)
100      CONTINUE 
C 
C        Averages, statistical errors 

         DO 35 I=1,NOUT

             AVQ2(I)=AVQ2(I)/NTRAJ
             ERQ2(I)=ERQ2(I)/NTRAJ
             ERQ2(I)=(ERQ2(I)-AVQ2(I)*AVQ2(I))/(NTRAJ-1)
             ERQ2(I)=SQRT(ERQ2(I))

             AVPSI1(I)=AVPSI1(I)/NTRAJ
             ERPSI1(I)=ERPSI1(I)/NTRAJ
             ERPSI1(I)=(ERPSI1(I)-AVPSI1(I)*AVPSI1(I))/(NTRAJ-1)
             ERPSI1(I)=SQRT(ERPSI1(I))
         
             AVETA(I)=AVETA(I)/NTRAJ
             ERETA(I)=ERETA(I)/NTRAJ
             ERETA(I)=(ERETA(I)-AVETA(I)*AVETA(I))/(NTRAJ-1)
             ERETA(I)=SQRT(ERETA(I))

             AVPSI2(I)=AVPSI2(I)/NTRAJ
             ERPSI2(I)=ERPSI2(I)/NTRAJ
             ERPSI2(I)=(ERPSI2(I)-AVPSI2(I)*AVPSI2(I))/(NTRAJ-1)
             ERPSI2(I)=SQRT(ERPSI2(I))

             AVT(I)=AVT(I)/NTRAJ
             ERT(I)=ERT(I)/NTRAJ
             ERT(I)=(ERT(I)-AVT(I)*AVT(I))/(NTRAJ-1)
             ERT(I)=SQRT(ERT(I))
    
             AVTD(I)=AVTD(I)/NTRAJ
             ERTD(I)=ERTD(I)/NTRAJ
             ERTD(I)=(ERTD(I)-AVTD(I)*AVTD(I))/(NTRAJ-1)
             ERTD(I)=SQRT(ERTD(I))
         
             AVTE(I)=AVTE(I)/NTRAJ
             ERTE(I)=ERTE(I)/NTRAJ
             ERTE(I)=(ERTE(I)-AVTE(I)*AVTE(I))/(NTRAJ-1)
             ERTE(I)=SQRT(ERTE(I))


             OUTIME1=NTIME*I*DELTAT
c        Output of results 
             WRITE(3,4) DELTAT,OUTIME1,AVETA(I),ERETA(I)
             WRITE(4,4) DELTAT,OUTIME1,AVPSI1(I),ERPSI1(I)
             WRITE(5,4) DELTAT,OUTIME1,AVPSI2(I),ERPSI2(I)
             WRITE(7,4)  DELTAT,OUTIME1,AVQ2(I),ERQ2(I)
c             WRITE(7,*)  DELTAT,OUTIME,AVQ2(I),ERQ2(I),SIGSQ(I) 
             WRITE(8,4)  DELTAT,OUTIME1,AVT(I),ERT(I)
c             WRITE(*,1)  DELTAT,OUTIME,AVT(I),ERT(I)
             WRITE(9,4)  DELTAT,OUTIME1,AVTE(I),ERTE(I)
             WRITE(10,4) DELTAT,OUTIME1,AVTD(I),ERTD(I)
35       CONTINUE
c1        FORMAT(F11.8,4X,F6.2,2X,F16.5,2X,F16.5)
4        FORMAT(F11.8,4X,F10.5,4X,F24.12,4X,F24.12) 
c         CLOSE(UNIT=114)
1000  CONTINUE 
C 

      CALL CPU_TIME(ENDTIME)


      WRITE(*,*) "SHEAR RATE TO BE WRITTEN : ",SR
23    FORMAT(I12)
      WRITE(112,23) ISEED

      BACKSPACE (UNIT=1)
      WRITE (1,3) Z,RMU,B,SR,E,H0,THI,INPAR,ENDTIME-STARTTIME
3     FORMAT(F5.1,2X,F4.1,4X,F8.1,6X,F8.2,6X,F5.2,
     &2X,F5.2,2X,F5.2,4X,I1,4X,F10.1)
      CLOSE (UNIT=1)
      close (unit=2) 
      CLOSE (unit=3)
      CLOSE (UNIT=4)
      CLOSE (UNIT=5)
      CLOSE (UNIT=7)
      CLOSE (UNIT=8)
      CLOSE (UNIT=9)
      CLOSE (UNIT=10)
      CLOSE (UNIT=89)
      CLOSE (UNIT=112)
C      CLOSE (UNIT=114)
c1100  STOP

C 
      open (unit=15, file='q2.dat')
      open (unit=16, file='psi1.dat')
      open (unit=17, file='eta.dat')
      open (unit=18, file='psi2.dat')
      open (unit=19, file='tau.dat')
      open (unit=20, file='tau_d.dat')
      open (unit=21, file='tau_e.dat')

      open (unit=2, file='tstep.dat')
      open (unit=1, file='inp.dat')

      open (unit=40, file='q2x.dat',status='UNKNOWN')
      open (unit=41, file='psi1x.dat',status='UNKNOWN')
      open (unit=42, file='etax.dat',status='UNKNOWN')
      open (unit=43, file='psi2x.dat',status='UNKNOWN')
      open (unit=44, file='taux.dat',status='UNKNOWN')
      open (unit=45, file='tau_dx.dat',status='UNKNOWN')
      open (unit=46, file='tau_ex.dat',status='UNKNOWN')

      READ (1,*) Z,RMU,B,SR,E,H0,THI 
      READ (2,*) NTIWID, NTRAJ

1     FORMAT(F11.8,4X,F8.4,4X,F16.12,4X,F16.12)
2     FORMAT(2X,F8.4,3X,F5.2,3X,F5.2,4X,F8.4,4X,F16.5,4X,F16.5,
     &4X,I3,3X,I3,A4,I3,3X,F10.5,A4,F10.5)



C     EXTRAPOLATION TO ZERO STEP SIZE 
      DELTAT=0.D0
      IFLAG=0
C 
 
c     Q2 VALUES

      DO 31 J=1,NTIWID
        DO 31 I=1,NOUT
c         write(*,*) I,J
         READ(15,*) TEMP1(J,I),OUTIME(I),TEMP2(J,I),TEMP3(J,I)
c         WRITE(*,*) TEMP1(J,I),OUTIME(I),TEMP2(J,I),TEMP3(J,I)

31    CONTINUE

      DO 33 I=1,NOUT
        DO 34 J=1,NTIWID
          XARR(J)=TEMP1(J,I)
          YARR(J)=TEMP2(J,I)
          SIGARR(J)=TEMP3(J,I)
c          write(*,*) XARR(J),YARR(J),SIGARR(J)
34    CONTINUE
C     IF NOPT=-1, TEXTRA HAS NOT BEEN ABLE TO EXTRAPOLATE.
C     IN SUCH A CASE, THE ERRLEV PARAMETER IS REDUCED 10 TIMES
C     AND THE EXTRAPOLATION IS REPEATED. IF, EVEN AFTER 6 TIMES
C     THE RESULT IS THE SAME, ALL VALUES ARE REPORTED AS 0.0
C
      NOPT = -1
      ERRLEV = 0.25D0
      DO 32 WHILE ((NOPT.EQ.-1).AND.(ERRLEV.GT.0.0000025D0))
        CALL TEXTRA(XARR,YARR,SIGARR,NTIWID,NDATM,IFLAG,ERRLEV)
        ERRLEV = ERRLEV/10.D0
32    CONTINUE
      IF(NOPT.EQ.-1) THEN
         YOPT=0.D0
         SIGOPT=0.D0
         NOPT=0
         NDUOPT=0
         XOPT=0
         ALIOPT=0.D0
         VLIOPT=0.D0
      ENDIF
      IF(IFLAG.EQ.0) THEN
         ALIOPT=0.D0
         VLIOPT=0.D0
      ENDIF
      WRITE (40,2) SR, E , H0, OUTIME(I), YOPT, SIGOPT, NOPT, NDUOPT,
     &' of ', NTIWID, ALIOPT,' +- ',VLIOPT
33    CONTINUE

c1100  STOP 



C    FIRST NORMAL STRESS DIFFERENCE
      DO 41 J=1,NTIWID
        DO 41 I=1,NOUT
c         write(*,*) I,J
         READ(16,*) TEMP1(J,I),OUTIME(I),TEMP2(J,I),TEMP3(J,I)
c         WRITE(*,*) TEMP1(J,I),OUTIME(I),TEMP2(J,I),TEMP3(J,I)

41    CONTINUE

      DO 43 I=1,NOUT
        DO 44 J=1,NTIWID
          XARR(J)=TEMP1(J,I)
          YARR(J)=TEMP2(J,I)
          SIGARR(J)=TEMP3(J,I)
44    CONTINUE
C     IF NOPT=-1, TEXTRA HAS NOT BEEN ABLE TO EXTRAPOLATE.
C     IN SUCH A CASE, THE ERRLEV PARAMETER IS REDUCED 10 TIMES
C     AND THE EXTRAPOLATION IS REPEATED. IF, EVEN AFTER 6 TIMES
C     THE RESULT IS THE SAME, ALL VALUES ARE REPORTED AS 0.0
C
      NOPT = -1
      ERRLEV = 0.25D0
      DO 42 WHILE ((NOPT.EQ.-1).AND.(ERRLEV.GT.0.0000025D0))
        CALL TEXTRA(XARR,YARR,SIGARR,NTIWID,NDATM,IFLAG,ERRLEV)
        ERRLEV = ERRLEV/10.D0
42    CONTINUE
      IF(NOPT.EQ.-1) THEN
         YOPT=0.D0
         SIGOPT=0.D0
         NOPT=0
         NDUOPT=0
         XOPT=0
         ALIOPT=0.D0
         VLIOPT=0.D0
      ENDIF
      IF(IFLAG.EQ.0) THEN
         ALIOPT=0.D0
         VLIOPT=0.D0
      ENDIF
      WRITE (41,2) SR, E, H0, OUTIME(I), YOPT, SIGOPT, NOPT, NDUOPT,
     &' of ', NTIWID, ALIOPT,' +- ',VLIOPT
43    CONTINUE


C     VISCOSITY

      DO 51 J=1,NTIWID
        DO 51 I=1,NOUT
c         write(*,*) I,J
         READ(17,*) TEMP1(J,I),OUTIME(I),TEMP2(J,I),TEMP3(J,I)
c         WRITE(*,*) TEMP1(J,I),OUTIME(I),TEMP2(J,I),TEMP3(J,I)

51    CONTINUE

      DO 53 I=1,NOUT
        DO 54 J=1,NTIWID
          XARR(J)=TEMP1(J,I)
          YARR(J)=TEMP2(J,I)
          SIGARR(J)=TEMP3(J,I)
54    CONTINUE
C     IF NOPT=-1, TEXTRA HAS NOT BEEN ABLE TO EXTRAPOLATE.
C     IN SUCH A CASE, THE ERRLEV PARAMETER IS REDUCED 10 TIMES
C     AND THE EXTRAPOLATION IS REPEATED. IF, EVEN AFTER 6 TIMES
C     THE RESULT IS THE SAME, ALL VALUES ARE REPORTED AS 0.0
C
      NOPT = -1
      ERRLEV = 0.25D0
      DO 52 WHILE ((NOPT.EQ.-1).AND.(ERRLEV.GT.0.0000025))
        CALL TEXTRA(XARR,YARR,SIGARR,NTIWID,NDATM,IFLAG,ERRLEV)
        ERRLEV = ERRLEV/10.D0
52    CONTINUE
      IF(NOPT.EQ.-1) THEN
         YOPT=0.D0
         SIGOPT=0.D0
         NOPT=0
         NDUOPT=0
         XOPT=0
         ALIOPT=0.D0
         VLIOPT=0.D0
      ENDIF
      IF(IFLAG.EQ.0) THEN
         ALIOPT=0.D0
         VLIOPT=0.D0
      ENDIF
      WRITE (42,2) SR, E, H0, OUTIME(I), YOPT, SIGOPT, NOPT, NDUOPT,
     &' of ', NTIWID, ALIOPT,' +- ',VLIOPT
53    CONTINUE


C     SECOND NORMAL STRESS DIFFERENCE

      DO 61 J=1,NTIWID
        DO 61 I=1,NOUT
c         write(*,*) I,J
         READ(18,*) TEMP1(J,I),OUTIME(I),TEMP2(J,I),TEMP3(J,I)
c         WRITE(*,*) TEMP1(J,I),OUTIME(I),TEMP2(J,I),TEMP3(J,I)

61    CONTINUE

      DO 63 I=1,NOUT
        DO 64 J=1,NTIWID
          XARR(J)=TEMP1(J,I)
          YARR(J)=TEMP2(J,I)
          SIGARR(J)=TEMP3(J,I)
64    CONTINUE
C     IF NOPT=-1, TEXTRA HAS NOT BEEN ABLE TO EXTRAPOLATE.
C     IN SUCH A CASE, THE ERRLEV PARAMETER IS REDUCED 10 TIMES
C     AND THE EXTRAPOLATION IS REPEATED. IF, EVEN AFTER 6 TIMES
C     THE RESULT IS THE SAME, ALL VALUES ARE REPORTED AS 0.0
C
      NOPT = -1
      ERRLEV = 0.25D0
      DO 62 WHILE ((NOPT.EQ.-1).AND.(ERRLEV.GT.0.0000025D0))
        CALL TEXTRA(XARR,YARR,SIGARR,NTIWID,NDATM,IFLAG,ERRLEV)
        ERRLEV = ERRLEV/10.D0
62    CONTINUE
      IF(NOPT.EQ.-1) THEN
         YOPT=0.D0
         SIGOPT=0.D0
         NOPT=0
         NDUOPT=0
         XOPT=0
         ALIOPT=0.D0
         VLIOPT=0.D0
      ENDIF
      IF(IFLAG.EQ.0) THEN
         ALIOPT=0.D0
         VLIOPT=0.D0
      ENDIF
      WRITE (43,2) SR, E, H0, OUTIME(I), YOPT, SIGOPT, NOPT, NDUOPT,
     &' of ', NTIWID, ALIOPT,' +- ',VLIOPT
63    CONTINUE

C     TAU (XY COMPONENT OF STRESS TENSOR)

      DO 71 J=1,NTIWID
        DO 71 I=1,NOUT
c         write(*,*) I,J
         READ(19,*) TEMP1(J,I),OUTIME(I),TEMP2(J,I),TEMP3(J,I)
c         WRITE(*,*) TEMP1(J,I),OUTIME(I),TEMP2(J,I),TEMP3(J,I)

71    CONTINUE

      DO 73 I=1,NOUT
        DO 74 J=1,NTIWID
          XARR(J)=TEMP1(J,I)
          YARR(J)=TEMP2(J,I)
          SIGARR(J)=TEMP3(J,I)
74    CONTINUE
C     IF NOPT=-1, TEXTRA HAS NOT BEEN ABLE TO EXTRAPOLATE.
C     IN SUCH A CASE, THE ERRLEV PARAMETER IS REDUCED 10 TIMES
C     AND THE EXTRAPOLATION IS REPEATED. IF, EVEN AFTER 6 TIMES
C     THE RESULT IS THE SAME, ALL VALUES ARE REPORTED AS 0.0
C
      NOPT = -1
      ERRLEV = 0.25D0
      DO 72 WHILE ((NOPT.EQ.-1).AND.(ERRLEV.GT.0.0000025D0))
        CALL TEXTRA(XARR,YARR,SIGARR,NTIWID,NDATM,IFLAG,ERRLEV)
        ERRLEV = ERRLEV/10.D0
72    CONTINUE
      IF(NOPT.EQ.-1) THEN
         YOPT=0.D0
         SIGOPT=0.D0
         NOPT=0
         NDUOPT=0
         XOPT=0
         ALIOPT=0.D0
         VLIOPT=0.D0
      ENDIF
      IF(IFLAG.EQ.0) THEN
         ALIOPT=0.D0
         VLIOPT=0.D0
      ENDIF
      WRITE (44,2) SR, E, H0, OUTIME(I), YOPT, SIGOPT, NOPT, NDUOPT,
     &' of ', NTIWID, ALIOPT,' +- ',VLIOPT
73    CONTINUE


C     TAU_E (ELASTIC COMPONENT OF TAU)

      DO 91 J=1,NTIWID
        DO 91 I=1,NOUT
c         write(*,*) I,J
         READ(21,*) TEMP1(J,I),OUTIME(I),TEMP2(J,I),TEMP3(J,I)
c         WRITE(*,*) TEMP1(J,I),OUTIME(I),TEMP2(J,I),TEMP3(J,I)

91    CONTINUE

      DO 93 I=1,NOUT
        DO 94 J=1,NTIWID
          XARR(J)=TEMP1(J,I)
          YARR(J)=TEMP2(J,I)
          SIGARR(J)=TEMP3(J,I)
94    CONTINUE
C     IF NOPT=-1, TEXTRA HAS NOT BEEN ABLE TO EXTRAPOLATE.
C     IN SUCH A CASE, THE ERRLEV PARAMETER IS REDUCED 10 TIMES
C     AND THE EXTRAPOLATION IS REPEATED. IF, EVEN AFTER 6 TIMES
C     THE RESULT IS THE SAME, ALL VALUES ARE REPORTED AS 0.0
C
      NOPT = -1
      ERRLEV = 0.25D0
      DO 92 WHILE ((NOPT.EQ.-1).AND.(ERRLEV.GT.0.0000025D0))
        CALL TEXTRA(XARR,YARR,SIGARR,NTIWID,NDATM,IFLAG,ERRLEV)
        ERRLEV = ERRLEV/10.D0
92    CONTINUE
      IF(NOPT.EQ.-1) THEN
         YOPT=0.D0
         SIGOPT=0.D0
         NOPT=0
         NDUOPT=0
         XOPT=0
         ALIOPT=0.D0
         VLIOPT=0.D0
      ENDIF
      IF(IFLAG.EQ.0) THEN
         ALIOPT=0.D0
         VLIOPT=0.D0
      ENDIF
      WRITE (46,2) SR, E, H0, OUTIME(I), YOPT, SIGOPT, NOPT, NDUOPT,
     &' of ', NTIWID, ALIOPT,' +- ',VLIOPT
93    CONTINUE

C    TAU_D (DISSIPATIVE COMPONENT OF TAU)

      DO 81 J=1,NTIWID
        DO 81 I=1,NOUT
c         write(*,*) I,J
         READ(20,*) TEMP1(J,I),OUTIME(I),TEMP2(J,I),TEMP3(J,I)
c         WRITE(*,*) TEMP1(J,I),OUTIME(I),TEMP2(J,I),TEMP3(J,I)

81    CONTINUE

      DO 83 I=1,NOUT
        DO 84 J=1,NTIWID
          XARR(J)=TEMP1(J,I)
          YARR(J)=TEMP2(J,I)
          SIGARR(J)=TEMP3(J,I)
84    CONTINUE
C     IF NOPT=-1, TEXTRA HAS NOT BEEN ABLE TO EXTRAPOLATE.
C     IN SUCH A CASE, THE ERRLEV PARAMETER IS REDUCED 10 TIMES
C     AND THE EXTRAPOLATION IS REPEATED. IF, EVEN AFTER 6 TIMES
C     THE RESULT IS THE SAME, ALL VALUES ARE REPORTED AS 0.0
C
      NOPT = -1
      ERRLEV = 0.25D0
      DO 82 WHILE ((NOPT.EQ.-1).AND.(ERRLEV.GT.0.0000025D0))
        CALL TEXTRA(XARR,YARR,SIGARR,NTIWID,NDATM,IFLAG,ERRLEV)
        ERRLEV = ERRLEV/10.D0
82    CONTINUE
      IF(NOPT.EQ.-1) THEN
         YOPT=0.D0
         SIGOPT=0.D0
         NOPT=0
         NDUOPT=0
         XOPT=0
         ALIOPT=0.D0
         VLIOPT=0.D0
      ENDIF
      IF(IFLAG.EQ.0) THEN
         ALIOPT=0.D0
         VLIOPT=0.D0
      ENDIF
      WRITE (45,2) SR, E, H0, OUTIME(I), YOPT, SIGOPT, NOPT, NDUOPT,
     &' of ', NTIWID, ALIOPT,' +- ',VLIOPT
83    CONTINUE































































      close (unit=15)
      close (unit=16)
      close (unit=17)
      close (unit=18)
      close (unit=19)
      close (unit=20)
      close (unit=21)

      close (unit=40)
      close (unit=41)
      close (unit=42)
      close (unit=43)
      close (unit=44)
      close (unit=45)
      close (unit=46)

      close (unit=1)
      close (unit=2)




1100  STOP 


      END 
C 





      SUBROUTINE HISET(Y)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /STEPBL/ THI,B,ZMU,RMU2,BAUXQ,BAUXR,DTH,DTQ,SRDT,
     &SRDTH,C1P,C2P,E,AMPL2,BETA,AB,A2,A4,AUX1,AUX2,
     &AUX3,AUX4,AUX5,S,QALPH,GEE1,GEE2,GEE3,GEE4 
       
      DATA C23/0.6666666666666667D0/
      DATA C43/1.3333333333333333D0/
      DATA C143/4.6666666666666670D0/
      DATA C83/2.6666666666666665D0/   
    
      REAL*8 COMP,COMPD,COMPI,Y1,ALPHA
      REAL*8 AUXD1,AUXD2,AUXD3,RES
      Y1=SQRT(Y)
      ALPHA=0.75D0*AB
      NTHI=NINT(THI)

c     Depending on THI, variables for the 
c     hydrodynamic tensor will be defined
C
C     I am considering HI tensor
C     to be of the form : 
C     \zeta Omega = (\alpha)/(Q)[A\delta+BQQ/Q^2]
C     so A = AUX2, B=AUX3 for both cases
C     In Regularized Oseen Burger case,
C     AUX4=M,AUX5=N
C
      SELECT CASE(NTHI)
C     FOR REGULARIZED OSEEN BURGERS          
          CASE(1)
              AUX1=Y+C43*A2
              AUX4=Y**(3.D0)+C143*A2*Y*Y+8.D0*A4*Y
              AUX5=Y**(3.D0)+2.D0*A2*Y*Y-C83*A4*Y
              AUX2=AUX4/(AUX1**(3.D0))
              AUX3=AUX5/(AUX1**(3.D00))
              AUXD1=((C43*A2)+(7.D0*Y))/(AUX1)
              AUXD2=(12.D0*(Y**(3.D0))+10.D0*C83*A2*Y*Y+4.D0*C83*A4*Y)
     &/(AUX1**(3.D0))
              AUXD3=(AUX2+AUX3)*AUXD1
              S=0.5D0*(AUXD3-AUXD2)
c     FOR RPY CASE      
          CASE(2)
              COMP=Y1/(2.D0*AB)
              COMPD=2.D0*COMP
              COMPI=1.D0/(COMPD) 
              IF(COMP.GE.1)THEN
                  AUX2=1.D0+(C23*COMPI*COMPI)
                  AUX3=1.D0-(2.D0*COMPI*COMPI) 
              ELSE
                  AUX2=(C43*COMPD)-(0.375D0*COMPD*COMPD)
                  AUX3=(0.125D0*COMPD*COMPD)
              ENDIF  
              S=AUX3
          CASE DEFAULT
              AUX1=Y
              AUX2=1.D0 
              AUX3=1.D0 
              S=1.D0
      END SELECT
      RES=S-AUX3
C     RES=S-B=A-r
      BETA=1.D0-(((ALPHA)/(Y1))*(AUX2+AUX3)) 
      AMPL2=(BETA)/((E*BETA)+1.D0)
      QALPH=(Y1)-(ALPHA*AUX2)
      GEE1=AMPL2*(((ALPHA*AUX3)/(BETA*QALPH))+E)
      GEE2=2.D0*(((AMPL2*AMPL2*ALPHA*AUX3)/(BETA*BETA*Y))
     &-(QALPH/Y)*GEE1 
     &-(E*(AMPL2*AMPL2*ALPHA*(RES)/(BETA*BETA*Y)))
     &-(E*ALPHA*(RES)*AMPL2/Y))
      GEE3=AMPL2/BETA
      GEE4=((2.D0*AMPL2*AMPL2*ALPHA*S)/(BETA*BETA*Y1))+(3.D0*AMPL2)
 
      RETURN
      END


      SUBROUTINE SEMIMP(Q)
C     Time step: semi-implicit predictor-corrector scheme for FENE+IV model
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (FOURPI=12.5663706143591725D0) 
      COMMON /STEPBL/ THI,B,ZMU,RMU2,BAUXQ,BAUXR,DTH,DTQ,SRDT,
     &SRDTH,C1P,C2P,E,AMPL2,BETA,AB,A2,A4,AUX1,AUX2,
     &AUX3,AUX4,AUX5,S,QALPH,GEE1,GEE2,GEE3,GEE4 

      DATA C23/0.6666666666666667D0/
      DATA C43/1.3333333333333333D0/
      DATA C143/4.6666666666666670D0/
      DATA C83/2.6666666666666665D0/      

      REAL*8 Q(3)
      REAL*8 B_D(3,3)
      REAL*8 S1,S2,S3,IVH,KH,IVAUX,KAUX
      NTHI=NINT(THI)
c      WRITE(*,*) "REACHED SI"
C     B_D refers to the square root of the diffusion tensor
c     To be constructed at the beginning of the timestep
c     based on initial values of Q
C     S refers to dot product of diffusion tensor with 
c     the weiner numbers, W
 
C     Auxiliary parameters 
      QL=Q(1)*Q(1)+Q(2)*Q(2)+Q(3)*Q(3) 
      SQQL=SQRT(QL) 
      CALL HISET(QL)

C      WRITE(*,*) "AUX1 : ",AUX1
C      WRITE(*,*) "AUX2 : ",AUX2
C      WRITE(*,*) "AUX3 : ",AUX3
C      WRITE(*,*) "AUX4 : ",AUX4
C      WRITE(*,*) "AUX5 : ",AUX5


      HF=(1.D0)/(1.D0-(QL/B)) 
      GT=-(GEE2*DTH/SQQL)
      KH=(E*AMPL2*Q(1)*Q(2)*SRDT)/QL 
      HE=-ZMU*EXP(-QL/RMU2)
      T0=DTH*AMPL2*(HF+HE) + GT + KH 
      T1=1.D0-T0
C     Construction of suitable random numbers 
      W1=RANULS()-0.5D0 
      W2=RANULS()-0.5D0 
      W3=RANULS()-0.5D0 
      W1=W1*(C1P*W1*W1+C2P) 
      W2=W2*(C1P*W2*W2+C2P) 
      W3=W3*(C1P*W3*W3+C2P) 
C 

C     Construction of B_D : square root of Diffusion Tensor
C     See Eqn (31) in Prabhakar's 2002 paper in J.Rheology

      GEE = (QALPH/SQQL)
      GEETIL = -(QALPH/SQQL)*GEE1
c      GEETIL=-(E*AMPL2*BETA)-((AUX3*0.75*AB)/SQQL)
c      BOTH DEFINITIONS OF GEETIL ARE ESSENTIALLY THE SAME
c      JUST MAKING IT MATCH THE EXPRESSION IN MY DERIVATION

      AUXF = SQRT(GEE)
      AUXS = (SQRT(GEE+GEETIL)-SQRT(GEE))/QL            

      B_D(1,1) = AUXF + (AUXS*Q(1)*Q(1))
      B_D(1,2) = AUXS*Q(1)*Q(2)
      B_D(1,3) = AUXS*Q(1)*Q(3)

      B_D(2,1) = AUXS*Q(2)*Q(1)
      B_D(2,2) = AUXF + (AUXS*Q(2)*Q(2))
      B_D(2,3) = AUXS*Q(2)*Q(3)

      B_D(3,1) = AUXS*Q(3)*Q(1)
      B_D(3,2) = AUXS*Q(3)*Q(2)
      B_D(3,3) = AUXF + (AUXS*Q(3)*Q(3))

c     S = (B_D).W

      S1 = B_D(1,1)*W1 + B_D(1,2)*W2 + B_D(1,3)*W3
      S2 = B_D(2,1)*W1 + B_D(2,2)*W2 + B_D(2,3)*W3
      S3 = B_D(3,1)*W1 + B_D(3,2)*W2 + B_D(3,3)*W3       

C     Predictor step 

      QAUX1=T1*Q(1)+(SRDT*Q(2))+S1 
      QAUX2=T1*Q(2)+S2 
      QAUX3=T1*Q(3)+S3
      QAUXL=QAUX1*QAUX1+QAUX2*QAUX2+QAUX3*QAUX3


c     All the auxilary functions have to be
c     calculated again. Bit of a pain.
c     Better to do this through function calls.

      SQAUXQL=SQRT(QAUXL) 
      CALL HISET(QAUXL)

      BAUXQ=3.D0*B*(1.D0+DTQ*AMPL2)
      BAUXR=9.D0*B*(1.D0-0.5D0*DTQ*AMPL2)


      GTAUX=-(GEE2*DTQ/SQAUXQL)
      KAUX=(E*AMPL2*SRDTH*QAUX1*QAUX2)/QAUXL
      HEAUX=-DTQ*AMPL2*ZMU*EXP(-QAUXL/RMU2)
      T4=GTAUX+KAUX+HEAUX 
      T3=1.D0-T4


C     Corrector step 
      Q(1)=T3*QAUX1+(SRDTH*(QAUX2-Q(2)))+0.5D0*T0*Q(1)
      Q(2)=T3*QAUX2+0.5D0*T0*Q(2)
      Q(3)=T3*QAUX3+0.5D0*T0*Q(3)
C     Exact solution of the implicit equation for the 
C     length of the connector vector 
      QL=Q(1)*Q(1)+Q(2)*Q(2)+Q(3)*Q(3) 
      SQQL=SQRT(QL) 
      AUXQ=(BAUXQ+QL)/9.D0 
      AUXR=SQQL*(BAUXR-QL)/27.D0 
      SQAUXQ=SQRT(AUXQ) 
      AUX=ACOS(AUXR/(AUXQ*SQAUXQ)) 
      XL=-2.D0*SQAUXQ*COS((AUX+FOURPI)/3.D0)+SQQL/3.D0 
C     Rescaling to obtain the proper length 
      RED=XL/SQQL 
      Q(1)=RED*Q(1) 
      Q(2)=RED*Q(2) 
      Q(3)=RED*Q(3) 
      RETURN 
      END  
C 
      SUBROUTINE RANILS(ISEED) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     Choice of ISEED: 0 <= ISEED <= 2000000000 (2E+9); 
C     ISEED can, for example, be formed from current time and date: 
C     2 digits each for seconds, minutes, hours, day, month 
C                    or minutes, hours, day, month, year 
      PARAMETER (IN=2147483563,IK=40014,IQ=53668,IR=12211,NTAB=32) 
      INTEGER IV(NTAB) 
      COMMON /RANBLS/ IDUM,IDUM2,IY,IV 
C     Initial seeds for two random number generators 
      IDUM=ISEED+123456789 
      IDUM2=IDUM 
C     Load the shuffle table (after 8 warm-ups) 
      DO 10 J=NTAB+8,1,-1 
         K=IDUM/IQ 
         IDUM=IK*(IDUM-K*IQ)-K*IR 
         IF(IDUM.LT.0) IDUM=IDUM+IN 
         IF(J.LE.NTAB) IV(J)=IDUM 
10    CONTINUE 
      IY=IV(1) 
      RETURN 
      END 
C 
      FUNCTION RANULS() 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      PARAMETER (IN1=2147483563,IK1=40014,IQ1=53668,IR1=12211, 
     &           IN2=2147483399,IK2=40692,IQ2=52774,IR2=3791, 
     &           NTAB=32,AN=1.D0/IN1,INM1=IN1-1,NDIV=1+INM1/NTAB) 
      INTEGER IV(NTAB) 
      COMMON /RANBLS/ IDUM,IDUM2,IY,IV 
C     Linear congruential generator 1 
      K=IDUM/IQ1 
      IDUM=IK1*(IDUM-K*IQ1)-K*IR1 
      IF(IDUM.LT.0) IDUM=IDUM+IN1 
C     Linear congruential generator 2 
      K=IDUM2/IQ2 
      IDUM2=IK2*(IDUM2-K*IQ2)-K*IR2 
      IF(IDUM2.LT.0) IDUM2=IDUM2+IN2 
C     Shuffling and subtracting 
      J=1+IY/NDIV 
      IY=IV(J)-IDUM2 
      IV(J)=IDUM 
      IF(IY.LT.1) IY=IY+INM1 
      RANULS=AN*IY 
      RETURN 
      END 

C 
      FUNCTION RANGLS()
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      SAVE IFLAG,GAUSS2
      DATA IFLAG/0/
      IF(IFLAG.EQ.0) THEN
10       CONTINUE
C        Pair of uniform random numbers in [-1,1]x[-1,1] 
         X1=2.D0*RANULS()-1.D0
         X2=2.D0*RANULS()-1.D0
C        If not in the unit circle, try again 
         XSQ=X1*X1+X2*X2
         IF(XSQ.GE.1.D0) GOTO 10
C        Pair of Gaussian random numbers; return one and 
C        save the other for next time 
         AUX=SQRT(-2.D0*LOG(XSQ)/XSQ)
         RANGLS=X1*AUX
         GAUSS2=X2*AUX
         IFLAG=1
      ELSE
         RANGLS=GAUSS2
         IFLAG=0
      ENDIF
      RETURN
      END
C 


C 
      SUBROUTINE TEXTRA(XARR,YARR,SIGARR,NDAT,NDATM,IFLAG,ERRLEV) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NDATP=50) 
      REAL*8 XARR(NDATM),YARR(NDATM),SIGARR(NDATM)   
      REAL*8 A(NDATP),COVAR(NDATP,NDATP) 
      INTEGER LFLAG(NDATP) 
      COMMON /EXTRP/ NOPT,NDUOPT,XOPT,YOPT,SIGOPT,ALIOPT,VLIOPT 
      EXTERNAL fpoly 
C 
C     Sorting the data such that XARR(1).LE.XARR(2).LE.XARR(3) ... 
C     (by straight insertion) 
      DO 10 J=2,NDAT 
         X=XARR(J) 
         Y=YARR(J) 
C         IF((Y+1.D0).EQ.Y) THEN
C             WRITE(*,*) "ERROR"
C         ENDIF 
         SIG=SIGARR(J) 
         DO 20 I=J-1,1,-1 
            IF(XARR(I).LE.X) GOTO 30 
            XARR(I+1)=XARR(I) 
            YARR(I+1)=YARR(I) 
            SIGARR(I+1)=SIGARR(I) 
20       CONTINUE 
         I=0 
30       XARR(I+1)=X 
         YARR(I+1)=Y 
         SIGARR(I+1)=SIG 
10    CONTINUE 
C 
      NOPT=-1 
      SIGOPT=6.022D23 
C     Fitting polynomials of various degrees N.LE.NMAX 
      NMAX=NDAT-2 
      IF(IFLAG.EQ.0) NMAX=NMAX+1 
      DO 1000 N=0,NMAX 
         NDATMI=N+2 
         IF(IFLAG.EQ.0.AND.N.GE.1) NDATMI=NDATMI-1 
C        Discarding data with large XARR 
         DO 500 NDATU=NDAT,NDATMI,-1 
            NDF=NDATU-NDATMI+1 
C           Least squares fit 
            DO 40 I=1,N+1 
               A(I)=0.D0 
               LFLAG(I)=1 
40          CONTINUE 
            DO 50 I=N+2,NDATP 
               A(I)=0.D0 
               LFLAG(I)=0 
50          CONTINUE 
            IF(N.GT.0) LFLAG(2)=IFLAG 
            CALL lfit(XARR,YARR,SIGARR,NDATU,A,LFLAG,NDATP, 
     &                                     COVAR,NDATP,TEST,fpoly) 
C           Chi-squared test; smaller statistical error bars? 
            IF(gammq(0.5D0*NDF,0.5D0*TEST).GT.ERRLEV) THEN 
               IF(SQRT(COVAR(1,1)).LT.SIGOPT) THEN 
                  YOPT=A(1) 
                  SIGOPT=SQRT(COVAR(1,1)) 
                  NOPT=N 
                  NDUOPT=NDATU 
                  XOPT=XARR(NDATU) 
                  ALIOPT=A(2) 
                  VLIOPT=SQRT(COVAR(2,2)) 
               ENDIF 
            ENDIF 
500      CONTINUE 
1000  CONTINUE 
C 
      RETURN 
      END   
C  
      SUBROUTINE lfit(x,y,sig,ndat,a,ia,ma,covar,npc,chisq,funcs) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER ma,ia(ma),npc,ndat,MMAX 
      REAL*8 chisq,a(ma),covar(npc,npc),sig(ndat),x(ndat),y(ndat) 
      EXTERNAL funcs 
      PARAMETER (MMAX=50) 
CU    USES covsrt,gaussj 
      INTEGER i,j,k,l,m,mfit 
      REAL*8 sig2i,sum,wt,ym,afunc(MMAX),beta(MMAX) 
      mfit=0 
      do 11 j=1,ma 
        if(ia(j).ne.0) mfit=mfit+1 
11    continue 
      if(mfit.eq.0) pause 'lfit: no parameters to be fitted' 
      do 13 j=1,mfit 
        do 12 k=1,mfit 
          covar(j,k)=0.D0 
12      continue 
        beta(j)=0.D0 
13    continue 
      do 17 i=1,ndat 
        call funcs(x(i),afunc,ma) 
        ym=y(i) 
        if(mfit.lt.ma) then 
          do 14 j=1,ma 
            if(ia(j).eq.0) ym=ym-a(j)*afunc(j) 
14        continue 
        endif 
        sig2i=1.D0/sig(i)**(2.D0) 
        j=0 
        do 16 l=1,ma 
          if (ia(l).ne.0) then 
            j=j+1 
            wt=afunc(l)*sig2i 
            k=0 
            do 15 m=1,l 
              if (ia(m).ne.0) then 
                k=k+1 
                covar(j,k)=covar(j,k)+wt*afunc(m) 
              endif 
15          continue 
            beta(j)=beta(j)+ym*wt 
          endif 
16      continue 
17    continue 
      do 19 j=2,mfit 
        do 18 k=1,j-1 
          covar(k,j)=covar(j,k) 
18      continue 
19    continue 
      call gaussj(covar,mfit,npc,beta,1,1) 
      j=0 
      do 21 l=1,ma 
        if(ia(l).ne.0) then 
          j=j+1 
          a(l)=beta(j) 
        endif 
21    continue 
      chisq=0.D0 
      do 23 i=1,ndat 
        call funcs(x(i),afunc,ma) 
        sum=0. 
        do 22 j=1,ma 
          sum=sum+a(j)*afunc(j) 
22      continue 
        chisq=chisq+((y(i)-sum)/sig(i))**2 
23    continue 
      call covsrt(covar,npc,ma,ia,mfit) 
      return 
      END 
C   
      SUBROUTINE fpoly(x,p,np) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER np 
      REAL*8 x,p(np) 
      INTEGER j 
      p(1)=1.D0 
      do 11 j=2,np 
        p(j)=p(j-1)*x 
11    continue 
      return 
      END 
C   
      SUBROUTINE covsrt(covar,npc,ma,ia,mfit) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER ma,mfit,npc,ia(ma) 
      REAL*8 covar(npc,npc) 
      INTEGER i,j,k 
      REAL*8 swap 
      do 12 i=mfit+1,ma 
        do 11 j=1,i 
          covar(i,j)=0.D0 
          covar(j,i)=0.D0 
11      continue 
12    continue 
      k=mfit 
      do 15 j=ma,1,-1 
        if(ia(j).ne.0)then 
          do 13 i=1,ma 
            swap=covar(i,k) 
            covar(i,k)=covar(i,j) 
            covar(i,j)=swap 
13        continue 
          do 14 i=1,ma 
            swap=covar(k,i) 
            covar(k,i)=covar(j,i) 
            covar(j,i)=swap 
14        continue 
          k=k-1 
        endif 
15    continue 
      return 
      END 
C   
      SUBROUTINE gaussj(a,n,np,b,m,mp) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER m,mp,n,np,NMAX 
      REAL*8 a(np,np),b(np,mp) 
      PARAMETER (NMAX=50) 
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX) 
      REAL*8 big,dum,pivinv 
      do 11 j=1,n 
        ipiv(j)=0 
11    continue 
      do 22 i=1,n 
        big=0.D0 
        do 13 j=1,n 
          if(ipiv(j).ne.1)then 
            do 12 k=1,n 
              if (ipiv(k).eq.0) then 
                if (abs(a(j,k)).ge.big)then 
                  big=abs(a(j,k)) 
                  irow=j 
                  icol=k 
                endif 
              else if (ipiv(k).gt.1) then 
                pause 'singular matrix in gaussj' 
              endif 
12          continue 
          endif 
13      continue 
        ipiv(icol)=ipiv(icol)+1 
        if (irow.ne.icol) then 
          do 14 l=1,n 
            dum=a(irow,l) 
            a(irow,l)=a(icol,l) 
            a(icol,l)=dum 
14        continue 
          do 15 l=1,m 
            dum=b(irow,l) 
            b(irow,l)=b(icol,l) 
            b(icol,l)=dum 
15        continue 
        endif 
        indxr(i)=irow 
        indxc(i)=icol 
        if (a(icol,icol).eq.0.D0) pause 'singular matrix in gaussj' 
        pivinv=1.D0/a(icol,icol) 
        a(icol,icol)=1.D0 
        do 16 l=1,n 
          a(icol,l)=a(icol,l)*pivinv 
16      continue 
        do 17 l=1,m 
          b(icol,l)=b(icol,l)*pivinv 
17      continue 
        do 21 ll=1,n 
          if(ll.ne.icol)then 
            dum=a(ll,icol) 
            a(ll,icol)=0.D0 
            do 18 l=1,n 
              a(ll,l)=a(ll,l)-a(icol,l)*dum 
18          continue 
            do 19 l=1,m 
              b(ll,l)=b(ll,l)-b(icol,l)*dum 
19          continue 
          endif 
21      continue 
22    continue 
      do 24 l=n,1,-1 
        if(indxr(l).ne.indxc(l))then 
          do 23 k=1,n 
            dum=a(k,indxr(l)) 
            a(k,indxr(l))=a(k,indxc(l)) 
            a(k,indxc(l))=dum 
23        continue 
        endif 
24    continue 
      return 
      END 
C   
      FUNCTION gammq(a,x) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      REAL*8 a,gammq,x 
CU    USES gcf,gser 
      REAL*8 gammcf,gamser,gln 
      if(x.lt.0.D0.or.a.le.0.D0)pause 'bad arguments in gammq' 
      if(x.lt.a+1.D0)then 
        call gser(gamser,a,x,gln) 
        gammq=1.D0-gamser 
      else 
        call gcf(gammcf,a,x,gln) 
        gammq=gammcf 
      endif 
      return 
      END 
C  
      SUBROUTINE gcf(gammcf,a,x,gln) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER ITMAX 
      REAL*8 a,gammcf,gln,x,EPS,FPMIN 
      PARAMETER (ITMAX=100,EPS=3.D-07,FPMIN=1.D-30) 
CU    USES gammln 
      INTEGER i 
      REAL*8 an,b,c,d,del,h,gammln 
      gln=gammln(a) 
      b=x+1.D0-a 
      c=1.D0/FPMIN 
      d=1.D0/b 
      h=d 
      do 11 i=1,ITMAX 
        an=-i*(i-a) 
        b=b+2.D0 
        d=an*d+b 
        if(abs(d).lt.FPMIN)d=FPMIN 
        c=b+an/c 
        if(abs(c).lt.FPMIN)c=FPMIN 
        d=1.D0/d 
        del=d*c 
        h=h*del 
        if(abs(del-1.D0).lt.EPS)goto 1 
11    continue 
      pause 'a too large, ITMAX too small in gcf' 
1     gammcf=exp(-x+a*log(x)-gln)*h 
      return 
      END 
C   
      SUBROUTINE gser(gamser,a,x,gln) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER ITMAX 
      REAL*8 a,gamser,gln,x,EPS 
      PARAMETER (ITMAX=100,EPS=3.D-07) 
CU    USES gammln 
      INTEGER n 
      REAL*8 ap,del,sum,gammln 
      gln=gammln(a) 
      if(x.le.0.D0)then 
        if(x.lt.0.D0)pause 'x < 0 in gser' 
        gamser=0.D0 
        return 
      endif 
      ap=a 
      sum=1.D0/a 
      del=sum 
      do 11 n=1,ITMAX 
        ap=ap+1.D0 
        del=del*x/ap 
        sum=sum+del 
        if(abs(del).lt.abs(sum)*EPS)goto 1 
11    continue 
      pause 'a too large, ITMAX too small in gser' 
1     gamser=sum*exp(-x+a*log(x)-gln) 
      return 
      END 
C   
      FUNCTION gammln(xx) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      REAL*8 gammln,xx 
      INTEGER j 
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6) 
      SAVE cof,stp 
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0, 
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2, 
     *-.5395239384953d-5,2.5066282746310005d0/ 
      x=xx 
      y=x 
      tmp=x+5.5d0 
      tmp=(x+0.5d0)*log(tmp)-tmp 
      ser=1.000000000190015d0 
      do 11 j=1,6 
        y=y+1.d0 
        ser=ser+cof(j)/y 
11    continue 
      gammln=tmp+log(stp*ser/x) 
      return 
      END 
C  

      FUNCTION betafn(x,y)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      REAL*8 betafnans,x,y
C     USES gammln
c     Returns  the  value  of  the  beta  function B(x,y).
      REAL*8 gammln
      betafnans=exp(gammln(x)+gammln(y)-gammln(x+y))
      return
      END


      FUNCTION EXPPHI(R)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      REAL*8 R
      COMMON /STEPBL/ THI,B,ZMU,RMU2,BAUXQ,BAUXR,DTH,DTQ,SRDT,
     &SRDTH,C1P,C2P,E,AMPL2,BETA,AB,A2,A4,AUX1,AUX2,
     &AUX3,AUX4,AUX5,S,QALPH,GEE1,GEE2,GEE3,GEE4

      VAR=1.D0-((R*R)**(B/2.D0))
      EXPPHI=VAR
      RETURN
      END

      FUNCTION FENGEN()
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(LN=21)
      PARAMETER(UPBOUND=2.D0)
      COMMON /STEPBL/ THI,B,ZMU,RMU2,BAUXQ,BAUXR,DTH,DTQ,SRDT,
     &SRDTH,C1P,C2P,E,AMPL2,BETA,AB,A2,A4,AUX1,AUX2,
     &AUX3,AUX4,AUX5,S,QALPH,GEE1,GEE2,GEE3,GEE4

      INTEGER DONE 
      REAL*8 L0S
      REAL*8 NORMEQPR
      REAL*8 R,CHECK
      REAL*8 X(LN),W(LN),SCRATCH(LN),XFACT,RG,FB,M

      DONE=0      
      NSEED=2233445
      MSEED=9988776
      RV=RAND(NSEED)
      RVCHK=RAND(MSEED)
      L0S=SQRT(B)
      X(01) = -0.993752170620389;   W(01) = 0.016017228257774
      X(02) = -0.967226838566306;   W(02) = 0.03695378977085309
      X(03) = -0.920099334150401;   W(03) = 0.05713442542685689
      X(04) = -0.853363364583318;   W(04) = 0.0761001136283793
      X(05) = -0.768439963475678;   W(05) = 0.09344442345603385
      X(06) = -0.667138804197413;   W(06) = 0.1087972991671478
      X(07) = -0.55161883588722;   W(07) = 0.1218314160537286
      X(08) = -0.424342120207439;   W(08) = 0.1322689386333376
      X(09) = -0.288021316802401;   W(09) = 0.1398873947910733
      X(10) = -0.145561854160895;   W(10) = 0.14452440398997
      X(11) = -2.4782829604619e-16;   W(11) = 0.1460811336496907
      X(12) = 0.145561854160895;   W(12) = 0.1445244039899697
      X(13) = 0.288021316802401;   W(13) = 0.1398873947910732
      X(14) = 0.424342120207439;   W(14) = 0.1322689386333371
      X(15) = 0.55161883588722;   W(15) = 0.1218314160537288
      X(16) = 0.667138804197413;   W(16) = 0.1087972991671492
      X(17) = 0.768439963475678;   W(17) = 0.09344442345603389
      X(18) = 0.853363364583317;   W(18) = 0.07610011362837897
      X(19) = 0.920099334150401;   W(19) = 0.05713442542685797
      X(20) = 0.967226838566306;   W(20) = 0.03695378977085323
      X(21) = 0.993752170620389;   W(21) = 0.01601722825777395

      M=18

      IF(M.GT.(B/2.D0)) THEN
          M=B/2.D0
      ENDIF


10    IF(DONE.NE.1) THEN
          NSEED=NSEED+1
          MSEED=MSEED+1
          RV=RAND(NSEED)
          RVCHK=RAND(MSEED)
          XFACT=SQRT(2.D0*M/B)
          NORM_B=0.D0
          FB=0.D0
          DO 11 N=1,LN
              RG = ((X(N)+1.D0)/2.D0)*XFACT
              NORM_B=NORM_B + W(N)*EXPPHI(RG)*RG*RG
              FB=FB + W(N)*EXPPHI(RG)*(RG**4.D0)
11        CONTINUE
          NORMEQPR=(RV*RV*EXPPHI(RV))/NORM_B
          CHECK=NORMEQPR/UPBOUND
          IF(RVCHK.LE.CHECK) THEN
              DONE=1
          ENDIF
      GOTO 10
      ENDIF
      FENGEN=NORMEQPR
      RETURN
      END 



 
