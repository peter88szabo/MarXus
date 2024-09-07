      PROGRAM SACM
C**********************************************************************
C Author: M. OLZMANN
C This program calculates microcanonical (E,J-resolved) and canonical
C rate constants for dissociation and recombination reactions.
C
C Remark: Centrifugal distortion of the reactant molecule is included
C in the case of a prolate symm. top. The treatment of linear reactant
C molecules is not possible in this version.
C Version: (02/11/95)
C**********************************************************************
      DIMENSION RHOJ(0:100005),AIRR(10),DR(100),D1(100),D2(100),
     *          ENULL(0:2000),SUM(0:100005),SUMA(0:100005),FAMEC(4),
     *          FAMUC(5)
      INTEGER SIGMAR,SIGMA1,SIGMA2,SUNEND
      REAL NYRK,MX,MY,MZ,LENGTH,NYR(100),NY1(100),NY2(100),KEJ,KB,KT,
     *     KDIPST,KREPST
      CHARACTER*70 TITLE
      DOUBLE PRECISION DGAMMA,ENVIBR,RHOJ,RHO,SUM,WPST,SUMM
      COMMON/A/ RHO(0:100005)
      COMMON/RHOI/ NROTI1,NROTI2,AIR1(10),AIR2(10)
      COMMON/BCEN/BE,DMORSE,A1,A2
      COMMON/SPOLL/NFG,EPSE(110),EPSP(110),XE(110),XP(110),
     *             EPSTS(110),XTS(110)
      COMMON/WCOUN/WPST(0:100005)
      COMMON/SCOUN/SUMM(0:100005)
      DATA FAMEC,FAMUC /1.,2.,2.,2.,1.2732,2.5465,1.6211,3.2423,6.4846/
C
C  Read in data
C
      OPEN(7,FILE='sacm.inp')
      OPEN(8,FILE='sacm.out')
      OPEN(9,FILE='sacmke.out')
      OPEN(10,FILE='sacmkt.out')
	OPEN(11,FILE='sacmrhoe.out')
      OPEN(12,FILE='sacmwe.out')
  100 FORMAT(1H )
      READ(5,*)NUMBER
      IC=0
 2001 READ(5,101) TITLE
  101 FORMAT(A70)
      WRITE(6,*)'SACM Program is running'
      READ(5,*)KKAN,KZR1,KZR2,KZP,KZE,KZBETA,KANH
      IF(KZR1 .NE. 1) GOTO 2003
      WRITE(6,*)'Linear reactants not possible in this program!'
      WRITE(8,*)'Linear reactants not possible in this program!'
      WRITE(9,*)'Linear reactants not possible in this program!'
      WRITE(10,*)'Linear reactants not possible in this program!'
      GOTO 2002
 2003 READ(5,*)NATOM,NROTIR,NVIBR,NROTI1,NVIB1,NROTI2,NVIB2
      NFG=3*NATOM-6
      IF(KZE)2,1,2
    1 READ(5,*)DMORSE
      GOTO3
    2 READ(5,*)H0K
    3 IF(KZBETA)5,4,5
    4 READ(5,*)ALBE,BETA,BETAEF
      GOTO6
    5 READ(5,*)ALBE,NYRK,BETAEF
    6 READ(5,*)MX,MY,MZ,QE,Q2E,THETAE
      READ(5,*)SIGMAR,SIGMAQ,SIGMA1,SIGMA2,SUNEND
      READ(5,*)EMAX,LENGTH,STEP
      READ(5,*)JOTA,JOTE,JSTEP
      READ(5,*)AR,BR,CR,B1,B2
      IF(NROTIR)70,69,70
   70 READ(5,*)(AIRR(I),I=1,NROTIR),V0R,NMR
      GOTO 691
   69 V0R=0.
  691 IF(NROTI1)72,71,72
   72 READ(5,*)(AIR1(I),I=1,NROTI1),V01,NM1
      GOTO 721
   71 V01=0.
  721 IF(NROTI2)74,73,74
   74 READ(5,*)(AIR2(I),I=1,NROTI2),V02,NM2
      GOTO 741
   73 V02=0.
  741 READ(5,*)(NYR(I),I=1,NVIBR)
      READ(5,*)(NY1(I),I=1,NVIB1)
      IF(NVIB2)76,75,76
   76 READ(5,*)(NY2(I),I=1,NVIB2)
   75 IF(KANH)86,85,86
   86 READ(5,*)(DR(I),I=1,NVIBR)
      READ(5,*)(D1(I),I=1,NVIB1)
      IF(NVIB2)87,85,87
   87 READ(5,*)(D2(I),I=1,NVIB2)
   85 READ(5,*)(EPSE(I),XE(I),I=1,NFG)
      READ(5,*)(EPSP(I),XP(I),I=1,NFG)
      READ(5,*)QELA,QELB,QELC
      READ(5,*)TEMPA,TEMPE,TEMPST
      EZP=0.
      EZN=0.

      DO 57 I=1,NVIBR
   57 EZN=EZN+.5*NYR(I)

      DO 58 I=1,NVIB1
   58 EZP=EZP+.5*NY1(I)

      DO 59 I=1,NVIB2
   59 EZP=EZP+.5*NY2(I)

      SEPSE=0.
      SEPSEQ=0.
      DO 67 I=1,NFG-1
      SEPSE=SEPSE+EPSE(I)
   67 SEPSEQ=SEPSEQ+EPSE(I)*EPSE(I)

      ESS=NVIBR-1.
      BETWR=(ESS-1.)*(ESS+NROTIR/2.)/ESS*SEPSEQ/(SEPSE*SEPSE)
      IF(KZE)56,92,56
   56 DMORSE=H0K+EZN-EZP
   92 H0K=DMORSE-EZN+EZP
      BE=(BR+CR)/2.
      GOTO(60,60,61),KZR2
   60 REDMAS=MX*MY/(MX+MY)

      GOTO62

   61 REDMAS=(MX+MY)*MZ/(MX+MY+MZ)
   62 IF(KZBETA)63,64,63
   63 BETA=.1218*SQRT(REDMAS/DMORSE)*NYRK
   64 A1=2./(BETA*QE)
      A2=1./(BETA*QE*BETA*QE)
      GOTO(66,66,65),KZR2
   65 ANEN=MX*(MY+MZ)*Q2E*Q2E-2.*MX*MZ*Q2E*QE*COS(THETAE)+MZ*(MX+MY)
     **QE*QE
      A1=A1*(MZ*(MX+MY)*QE*QE-MX*MZ*Q2E*QE*COS(THETAE))/ANEN
      A2=A2*MZ*(MX+MY)*QE*QE/ANEN
C
C Print of the input data
C
   66 WRITE(8,100)
      WRITE(8,101)TITLE
      WRITE(9,101)TITLE
      WRITE(9,100)
      WRITE(8,100)
      WRITE(8,100)
      GOTO (7,8,9,10),KZR1
    7 WRITE(8,102)
      GOTO11
    8 WRITE(8,103)
      GOTO11
    9 WRITE(8,104)
      GOTO11
   10 WRITE(8,105)
  102 FORMAT(10X,16HREACTANT: LINEAR)
  103 FORMAT(10X,23HREACTANT: SPHERICAL TOP)
  104 FORMAT(10X,26HREACTANT: OBLATE SYMM. TOP)
  105 FORMAT(10X,27HREACTANT: PROLATE SYMM. TOP)
   11 IF(KZR2-2)12,12,13
   12 WRITE(8,106)
      GOTO14
   13 WRITE(8,107)
  106 FORMAT(20X,15H(QUASIDIATOMIC),/)
  107 FORMAT(20X,16H(QUASITRIATOMIC),/)
   14 GOTO (15,16,17,18,19),KZP
   15 WRITE(8,108)
      GOTO20
   16 WRITE(8,109)
      GOTO20
   17 WRITE(8,110)
      GOTO20
   18 WRITE(8,111)
      GOTO20
   19 WRITE(8,112)
  108 FORMAT(10X,32HFRAGMENT1/FRAGMENT2: LINEAR/ATOM,//)
  109 FORMAT(10X,39HFRAGMENT1/FRAGMENT2: SPHERICAL TOP/ATOM,//)
  110 FORMAT(10X,34HFRAGMENT1/FRAGMENT2: LINEAR/LINEAR,//)
  111 FORMAT(10X,41HFRAGMENT1/FRAGMENT2: SPHERICAL TOP/LINEAR,//)
  112 FORMAT(10X,48HFRAGMENT1/FRAGMENT2: SPHERICAL TOP/SPHERICAL TOP,//)
   20 WRITE(8,113)
  113 FORMAT(10X,9HREACTANT:)
      WRITE(8,114)NVIBR
      WRITE(8,115)NROTIR
  114 FORMAT(15X,20H        VIBRATIONS  ,I3,/)
  115 FORMAT(15X,21HINTERNAL ROTATIONS   ,I2,/)
      WRITE(8,116)
  116 FORMAT(10X,11HFRAGMENT 1:)
      WRITE(8,114)NVIB1
      WRITE(8,115)NROTI1
      WRITE(8,117)
  117 FORMAT(10X,11HFRAGMENT 2:)
      WRITE(8,114)NVIB2
      WRITE(8,115)NROTI2
      WRITE(8,100)
      WRITE(8,100)
      IF(KZE)22,21,22
   21 WRITE(8,118)DMORSE
  118 FORMAT(10X,19HD(MORSE) IN 1/CM:  ,F8.1,/)
      GOTO23
   22 WRITE(8,119)H0K
  119 FORMAT(10X,37HENTHALPY OF REACTION AT 0K IN 1/CM:  ,F8.1,/)
   23 WRITE(8,120)BETA
  120 FORMAT(10X,37HMORSE-PARAMETER BETA IN 1/(1E-10M):  ,F5.3,/)
      IF(KZBETA)25,26,25
   25 WRITE(8,121)NYRK
  121 FORMAT(10X,34HNY(REACTION COORDINATE) IN 1/CM:  ,F6.1,/)
   26 WRITE(8,152)ALBE
  152 FORMAT(10X,13HALPHA/BETA:  ,F5.3,/)
      WRITE(8,171)BETAEF
  171 FORMAT(10X,12HBETA-EFF.:  ,F5.3,/)
      WRITE(8,150)EZN
      WRITE(8,151)EZP
  150 FORMAT(10X,48HVIB. ZERO-POINT ENERGY OF THE REACTANT IN 1/CM: ,
     *F8.1,/)
  151 FORMAT(10X,49HVIB. ZERO-POINT ENERGY OF THE FRAGMENTS IN 1/CM: ,
     *F8.1,/)
      WRITE(8,122)MX,MY,MZ
  122 FORMAT(10X,16HMASSES IN G/MOL:,//,10X,4HMX= ,F5.1,5X,
     *4HMY= ,F5.1,5X,4HMZ= ,F5.1,/)
      IF(KZR2-2)27,27,28
   27 WRITE(8,123)QE
  123 FORMAT(10X,14HQE IN 1E-10M: ,F5.3,/)
      GOTO29
   28 WRITE(8,124)QE,Q2E,THETAE
  124 FORMAT(10X,14HQE IN 1E-10M: ,F5.3,//,10X,15HQ2E IN 1E-10M: ,
     *F5.3,//,10X,19HTHETAE IN DEGREES: ,F5.1,/)
   29 WRITE(8,125)SIGMAR,SIGMA1,SIGMA2
  125 FORMAT(10X,17HSYMMETRY NUMBERS:,//,10X,10HREACTANT: ,I2,5X,
     *12HFRAGMENT 1: ,I2,5X,12HFRAGMENT 2: ,I2,/)
      IF(SUNEND-1)30,30,31
   30 WRITE(8,126)
      GOTO 32
   31 WRITE(8,127)
  126 FORMAT(10X,25H(NON-IDENTICAL FRAGMENTS),//)
  127 FORMAT(10X,21H(IDENTICAL FRAGMENTS),//) 
   32 WRITE(8,128)
  128 FORMAT(10X,50HROTATIONAL CONSTANTS AND FREQUENCIES, BOTH IN 1/CM,
     *//)
      WRITE(8,129)
  129 FORMAT(30X,8HREACTANT,/)
      WRITE(8,130)
  130 FORMAT(10X,19HEXTERNAL ROTATIONS:,/)
      WRITE(8,131)AR,BR,CR
  131 FORMAT(10X,3HA= ,F6.3,5X,3HB= ,F6.3,5X,3HC= ,F6.3,/)
      IF(NROTIR)77,78,77
   77 WRITE(8,132)
  132 FORMAT(10X,19HINTERNAL ROTATIONS:,/)
      WRITE(8,133)(AIRR(I),I=1,NROTIR)
      IF(V0R)771,78,771
  771 WRITE(8,100)
      WRITE(8,1331)V0R,NMR
 1331 FORMAT(11X,22H(HINDRANCE POTENTIAL: ,F6.1,2H, ,I1,6H-FOLD))      
  133 FORMAT((10X,F6.3,4(4X,F6.3)))
   78 WRITE(8,100)
      WRITE(8,134)
  134 FORMAT(10X,12HFREQUENCIES:,/)
      WRITE(8,135)(NYR(I),I=1,NVIBR)
  135 FORMAT((10X,F6.1,5(3X,F6.1)))
      WRITE(8,100)
      IF(KANH)88,89,88
   88 WRITE(8,147)
  147 FORMAT(10X,'INDIVIDUAL DISSOCIATION ENERGIES OF THE OSCILLATORS:
     *',/)
      WRITE(8,148)(DR(I),I=1,NVIBR)
  148 FORMAT((9X,F8.1,5(1X,F8.1)))   
   89 WRITE(8,100)
      WRITE(8,100)
      WRITE(8,136)
  136 FORMAT(30X,10HFRAGMENT 1,/)
      WRITE(8,130)
      WRITE(8,137)B1
  137 FORMAT(10X,3HB= ,F6.3,/)
      IF(NROTI1)79,80,79
   79 WRITE(8,132)
      WRITE(8,133)(AIR1(I),I=1,NROTI1)
      IF(V01)791,80,791
  791 WRITE(8,100)
      WRITE(8,1331)V01,NM1
   80 WRITE(8,100)
      WRITE(8,134)
      WRITE(8,135)(NY1(I),I=1,NVIB1)
      WRITE(8,100)
      IF(KANH)93,94,93
   93 WRITE(8,147)
      WRITE(8,148)(D1(I),I=1,NVIB1)   
   94 WRITE(8,100)
      WRITE(8,100)
      IF(KZP.EQ.1.OR.KZP.EQ.2) GOTO 84
      WRITE(8,138)
  138 FORMAT(30X,10HFRAGMENT 2,/)
      WRITE(8,130)
      WRITE(8,137)B2
      IF(NROTI2)81,82,81
   81 WRITE(8,132)
      WRITE(8,133)(AIR2(I),I=1,NROTI2)
      IF(V02)811,82,811
  811 WRITE(8,100)
      WRITE(8,1331)V02,NM2
   82 WRITE(8,100)
      IF(NVIB2)83,84,83
   83 WRITE(8,134)
      WRITE(8,135)(NY2(I),I=1,NVIB2)
      WRITE(8,100)
      IF(KANH)95,84,95
   95 WRITE(8,147)
      WRITE(8,148)(D2(I),I=1,NVIB2)   
   84 WRITE(8,100)
      WRITE(8,100)
      WRITE(8,100)
C=============================================================================      

      EPSEAV=0.
      KAV=0
      SSTAR=1.

      DO 184 I=1,NFG-1
      SSTAR=SSTAR+XE(I)
      IF(XP(I).GT.0.51) GOTO 184
      IF(XE(I).LT.0.51) GOTO 184
      KAV=KAV+1
      EPSEAV=EPSEAV+EPSE(I)
  184 CONTINUE

      ALBEEF=ALBE*BETA/BETAEF
      DMOREF=BETA/BETAEF*BETA/BETAEF*DMORSE
      EPSEAV=EPSEAV/KAV
      C3=4.-1.29*LOG(EPSEAV/DMOREF)
C============================================================
C interpolates between reactant and product eigenvalues

      CALL SPOL(ALBEEF,DMOREF)
C============================================================

      WRITE(8,153)
  153 FORMAT(10X,34HCORRELATIONS (VIB: 1.0, ROT: 0.5):,//
     *14X,54HREACTANT               CHANNEL                 PRODUCT,/)
      WRITE(8,154)(EPSE(I),XE(I),EPSTS(I),XTS(I),EPSP(I),XP(I),I=1,NFG)
  154 FORMAT(10X,F7.2,3X,1H(,F3.1,1H),2X,5H---->,2X,F7.2,3X,
     *       1H(,F4.2,1H),2X,5H---->,2X,F7.2,3X,1H(,F3.1,1H))
      WRITE(8,100)
      WRITE(8,100)
      IF(KKAN.EQ.1) GOTO 99
      WRITE(8,139)
  139 FORMAT(30X,17HCOUNTING-ROUTINES,//)
      WRITE(8,140)EMAX,LENGTH
  140 FORMAT(10X,24HMAXIMUM ENERGY IN 1/CM: ,F8.1,//,10X,
     *18HSTEPSIZE IN 1/CM: ,F6.1,/)
      WRITE(8,145)
  145 FORMAT(10X,9H    E    ,7X,6HRHO(E),13X,4HW(E),13X,4HK(E),//
     *10X,9H [1/CM]  ,8X,6H[CM]  ,12X,4H    ,12X,7H[1/SEC],//)
C
C Planck's constant/Boltzmann's constant
C
   99 HCM=3.3356E-11
      KB=.69503
      SIGMAP=SIGMA1*SIGMA2*SUNEND
      SIGMAW=SIGMAP+(SIGMAQ-SIGMAP)*EXP(-C3*ALBEEF)
      DEZ=EZN-EZP-NYR(NVIBR-2)/2.
      IF(KKAN.EQ.1) GOTO 1500
      FACTOR=SIGMAR/(SIGMAW*HCM) 
      ISTEP=INT(STEP/LENGTH+.5)
C
C Calculation of the number of open channels
C
      N=INT(EMAX/LENGTH+.5)
      CALL SCOUNT(NFG,EPSTS,XTS,EMAX,LENGTH)
      DO 1527 I=0,N
 1527 SUMA(I)=SUMM(I)
      CALL SCOUNT(NFG,EPSP,XP,EMAX,LENGTH)
C
C Initiation of the J-loop
C
      DO 1000 JOT=JOTA,JOTE,JSTEP
      WRITE(6,155)JOT
  155 FORMAT(1X,4HJ = ,I4)
C
C Estimation of E0(J)
C
      CALL ETERAT(DMORSE,EZP,DEZ,ALBE,BE,A1,A2,JOT,ENULLJ)
      IF(ENULLJ)90,90,91
   90 WRITE(8,146)JOT
  146 FORMAT(10X,3HJ= ,I4,3X,35HTHIS STATE IS ROTATIONALLY INSTABLE,/)
      GOTO2000
C
C Calculation of the density of states
C
   91 GOTO (36,37,38,39),KZR1
C
C Reactant = linear
C
C The following code for the linear reactant is not correct.
C Do not make access to this part of the program.
C
   36 IF(JOT.NE.0) GOTO40
      NVIBR=NVIBR-2
      CALL DCOUNT(NVIBR,EMAX,LENGTH,NYR,NROTIR,AIRR)
      N=INT(EMAX/LENGTH+.5)
      DO 1521 I=0,N
      EI=I*LENGTH
 1521 RHO(I)=RHO(I)*FANHRO(KANH,EI,NVIBR,NYR,DR)      
      NVIBR=NVIBR+2
      DO 41 I=0,N
   41 RHOJ(I)=RHO(I)
      GOTO 54
   40 IF(JOT.EQ.JOTA.OR.JOT.EQ.JOTA+JSTEP) GOTO43
      GOTO44
   43 CALL DCOUNT(NVIBR,EMAX,LENGTH,NYR,NROTIR,AIRR)
      N=INT(EMAX/LENGTH+.5)
      DO 1522 I=0,N
      EI=I*LENGTH
 1522 RHO(I)=RHO(I)*FANHRO(KANH,EI,NVIBR,NYR,DR)      
      ENVIBR=DBLE(NVIBR)
   44 F1A=REAL(DGAMMA(ENVIBR)/DGAMMA(ENVIBR-2))*(1+JOT)*(1+JOT/2)
      BEQ=BCENT(JOT)
      BRJOT=INT(BEQ*JOT*(JOT+1)/LENGTH+.5)*LENGTH
      MINIMI=INT(BRJOT/LENGTH)
      DO 42 I=MINIMI,N
      EKORR=I*LENGTH-BRJOT
      F1=F1A*(NYR(NVIBR)/EKORR)*(NYR(NVIBR)/EKORR)
      F1=F1/(F1+1.)
   42 RHOJ(I)=RHO(INT(EKORR/LENGTH))*F1
      GOTO 54
C
C Reactant = spherical top
C
   37 IF(JOT-JOTA)45,45,46
   45 CALL DCOUNT(NVIBR,EMAX,LENGTH,NYR,NROTIR,AIRR)
      N=INT(EMAX/LENGTH+.5)
      DO 1523 I=0,N
      EI=I*LENGTH
 1523 RHO(I)=RHO(I)*FANHRO(KANH,EI,NVIBR,NYR,DR)
     *             *RRHIND(NVIBR,V0R,EI)
   46 BEQ=BR
      BRJOT=INT(BEQ*JOT*(JOT+1)/LENGTH+.5)*LENGTH
      MINIMI=INT(BRJOT/LENGTH)
      IRHOJ=2*JOT+1
      DO 47 I=MINIMI,N
      EKORR=I*LENGTH-BRJOT
   47 RHOJ(I)=IRHOJ*RHO(INT(EKORR/LENGTH))
      GOTO 54
C
C Reactant = oblate symmetrical top
C
   38 IF(JOT-JOTA)48,48,49
   48 CALL DCOUNT(NVIBR,EMAX,LENGTH,NYR,NROTIR,AIRR)
      N=INT(EMAX/LENGTH+.5)
      DO 1524 I=0,N
      EI=I*LENGTH
 1524 RHO(I)=RHO(I)*FANHRO(KANH,EI,NVIBR,NYR,DR)
     *             *RRHIND(NVIBR,V0R,EI)
   49 BEQ=BR
      BECE=BEQ-CR
      BRJOT=INT(BEQ*JOT*(JOT+1)/LENGTH+.5)*LENGTH
      MINIMI=INT((BEQ*JOT*(JOT+1)-BECE*JOT*JOT)/LENGTH+.5)
      DO 187 I=0,N
      IF(JOT.EQ.0) THEN
         RHOJ(I)=RHO(I)
         GOTO 187
      ELSE
      IF(I.LT.MINIMI) THEN
         RHOJ(I)=0.
      ELSE
      IKORR=INT((I*LENGTH-BEQ*JOT*(JOT+1)+BECE*JOT*JOT)/LENGTH+.5)
      RHOJ(I)=2.*RHO(IKORR)
      IF(I*LENGTH-BRJOT) 189,185,185
  185 KKMIN=1
      KMIN=1
      GOTO 186
  189 KKMIN=2
      KMIN=INT(SQRT((BEQ*JOT*(JOT+1)-I*LENGTH)/BECE)+1.)
  186 DO 50 K=JOT-1,KMIN,-1
      IKORRB=INT((I*LENGTH-BEQ*JOT*(JOT+1)+BECE*K*K)/LENGTH+.5)
   50 RHOJ(I)=RHOJ(I)+2.*RHO(IKORRB)
      IF(KKMIN-1) 188,188,187
  188 RHOJ(I)=RHOJ(I)+RHO(INT((I*LENGTH-BRJOT)/LENGTH))
      ENDIF
      ENDIF
  187 CONTINUE
      GOTO 54
C
C Reactant = prolate symmetrical top
C
   39 IF(JOT-JOTA)51,51,52
   51 CALL DCOUNT(NVIBR,EMAX,LENGTH,NYR,NROTIR,AIRR)
      N=INT(EMAX/LENGTH+.5)
      DO 1525 I=0,N
      EI=I*LENGTH
 1525 RHO(I)=RHO(I)*FANHRO(KANH,EI,NVIBR,NYR,DR)
     *             *RRHIND(NVIBR,V0R,EI)
   52 BEQ=BCENT(JOT)
      ABE=AR-BEQ
      BRJOT=INT(BEQ*JOT*(JOT+1)/LENGTH+.5)*LENGTH
      MINIMI=INT(BRJOT/LENGTH)
      DO 53 I=MINIMI,N
      EKORR=I*LENGTH-BRJOT
      RHOJ(I)=RHO(INT(EKORR/LENGTH))
      KX=INT(SQRT(EKORR/ABE))
      KMAX=MIN(JOT,KX)
      DO 53 K=1,KMAX
      EKORRA=EKORR-INT(ABE*K*K/LENGTH+.5)*LENGTH
   53 RHOJ(I)=RHOJ(I)+2.*RHO(INT(EKORRA/LENGTH))
   54 CONTINUE
C
C Calculation of the angular momentum correction factor Faminf
C
      CALL WCOUNT(KZP,0,EMAX,LENGTH,JOT,B1,B2,NVIB1,NVIB2,NY1,
     *            NY2,NROTI1,NROTI2)
      DO 1529 I=1,N
 1529 WPST(I)=WPST(I)/SUMM(I)
C Now, Faminf is contained in WPST
      V0F=V01+V02
C Either V01 or V02 or both must be zero
      NVIBG=NVIB1+NVIB2
      DO 55 I=1,N
      EI=I*LENGTH
      EIF=EI+AWR(2.*EI/SEPSE,BETWR)*SEPSE/2.
      FAMES=FAME(KZR1,EIF,JOT,SSTAR,EPSE(NFG),SQRT(BR*CR))
      WPST(I)=WPST(I)+(FAMES-WPST(I))*EXP(-C3*ALBEEF)
   55 SUM(I)=SUMA(I)*FANHWW(KANH,EI,NVIB1,NVIB2,NY1,NY2,D1,D2)
     *              *RWHIND(NVIBG,V0F,EI)
     *              *WPST(I)
C
C Print of the results
C
      WRITE(8,143)JOT,ENULLJ-EZN,BEQ
  143 FORMAT(10X,3HJ= ,I4,5X,15HE0(J) IN 1/CM: ,F7.1,5X,
     *15HBCENT IN 1/CM: ,F5.3,/)
      WRITE(9,100)
      WRITE(9,180)JOT,ENULLJ-EZN
  180 FORMAT(10X,3HJ= ,I4,5X,15HE0(J) IN 1/CM: ,F7.1,/)
      ISTART=INT((ENULLJ-EZN)/LENGTH+.5)
      EI=0.
      WRITE(8,144)ENULLJ-EZN,RHOJ(ISTART)/SIGMAR,SUMA(0),FACTOR*SIGMAW
     */RHOJ(ISTART)
      WRITE(9,181)ENULLJ-EZN,FACTOR*SIGMAW/RHOJ(ISTART)
      WRITE(12,'(I6)')INT((N-ISTART-ISTEP)/ISTEP)+2
      WRITE(12,181)ENULLJ-EZN,SUMA(0)
  181 FORMAT(10X,F8.1,5X,1PE12.6)
      DO 68 I=ISTEP,N-ISTART,ISTEP
      EI=I*LENGTH
      KEJ=FACTOR*SUM(I)/RHOJ(I+ISTART)
      WRITE(8,144)EI+ENULLJ-EZN,RHOJ(I+ISTART)/SIGMAR,SUM(I)/SIGMAW,KEJ
  144 FORMAT(10X,F8.1,3(5X,1PE12.6))
      WRITE(9,181)EI+ENULLJ-EZN,KEJ
      WRITE(12,181)EI+ENULLJ-EZN,SUM(I)/SIGMAW
   68 CONTINUE
      WRITE(11,'(I6)')INT((N-MINIMI)/ISTEP)+1
      DO 699 I=MINIMI,N,ISTEP
	EI=I*LENGTH
      WRITE(11,181)EI,RHOJ(I)/SIGMAR
  699 CONTINUE
      WRITE(8,100)
 1000 CONTINUE
C
C Calculation of the high-pressure rate constants
C
 1500 WRITE(8,100)
      WRITE(8,100)
      WRITE(8,156)
  156 FORMAT(24X,28HHIGH-PRESSURE RATE CONSTANTS,/)
      WRITE(10,182)
  182 FORMAT(10X,1HT,9X,5HKDISS,7X,4HKREC,8X,4HKEQU,/)
      JOT=-1
 1501 JOT=JOT+1
      CALL ETERAT(DMORSE,EZP,DEZ,ALBE,BE,A1,A2,JOT,ENULL(JOT))
      IF(ENULL(JOT))1503,1503,1502
 1502 IF((ENULL(JOT)-ENULL(0))/(KB*TEMPE).GT.87.) GOTO 1503
      TERM=(2*JOT+1)*EXP(-(ENULL(JOT)-ENULL(0))/(KB*TEMPE))
      IF(TERM.GE.1E-2) GOTO 1501
 1503 JOTMAX=JOT-1
      DE0Z=ENULL(0)-DMORSE-EZP
      FAMC=FAMUC(KZP)+(FAMEC(KZR1)-FAMUC(KZP))*EXP(-C3*ALBEEF)
      DO 1600 TEMP=TEMPA,TEMPE,TEMPST
      KT=KB*TEMP
      QCENT=0.
      DO 1504 JOT=0,JOTMAX
      IF((ENULL(JOT)-ENULL(0))/KT.GT.87.) GOTO 1504
      QCENT=QCENT+(2*JOT+1)*EXP(-(ENULL(JOT)-ENULL(0))/KT)
 1504 CONTINUE
      QVIBA=1.
      DO 1505 I=1,NVIBR
 1505 QVIBA=QVIBA/(1.-EXP(-NYR(I)/KT))
      QVIBB=1.
      DO 1506 I=1,NVIB1
 1506 QVIBB=QVIBB/(1.-EXP(-NY1(I)/KT))
      QVIBC=1.
      DO 1507 I=1,NVIB2
 1507 QVIBC=QVIBC/(1.-EXP(-NY2(I)/KT))
      QIRA=1.478**NROTIR*TEMP**(NROTIR/2.)
      DO 1508 I=1,NROTIR
 1508 QIRA=QIRA/SQRT(AIRR(I))
      QIRA=QIRA*FQHIND(KT,V0R,QIRA,NMR)
      QIRB=1.478**NROTI1*TEMP**(NROTI1/2.)
      DO 1509 I=1,NROTI1
 1509 QIRB=QIRB/SQRT(AIR1(I))
      QIRB=QIRB*FQHIND(KT,V01,QIRB,NM1)
      QIRC=1.478**NROTI2*TEMP**(NROTI2/2.)
      DO 1510 I=1,NROTI2
 1510 QIRC=QIRC/SQRT(AIR2(I))
      QIRC=QIRC*FQHIND(KT,V02,QIRC,NM2)
      GOTO (1511,1512,1513,1513),KZR1
 1511 QERA=.695*TEMP/BR
      GOTO 1514
 1512 QERA=1.027*TEMP/BR*SQRT(TEMP/BR)
      GOTO 1514
 1513 QERA=1.027*SQRT(TEMP*TEMP*TEMP/(AR*BR*CR))
 1514 GOTO (1515,1516,1517,1518,1519),KZP
 1515 QERB=.695*TEMP/B1
      QERC=1.
      GOTO 1520
 1516 QERB=1.027*TEMP/B1*SQRT(TEMP/B1)
      QERC=1.
      GOTO 1520
 1517 QERB=.695*TEMP/B1
      QERC=.695*TEMP/B2
      GOTO 1520
 1518 QERB=1.027*TEMP/B1*SQRT(TEMP/B1)
      QERC=.695*TEMP/B2
      GOTO 1520
 1519 QERB=1.027*TEMP/B1*SQRT(TEMP/B1)
      QERC=1.027*TEMP/B2*SQRT(TEMP/B2)
 1520 QVRA=QVIBA*QIRA*QERA/SIGMAR
      QVRB=QVIBB*QIRB*QERB/SIGMA1
      QVRC=QVIBC*QIRC*QERC/SIGMA2
      QSTARP=1.
      DO 1530 I=1,NFG
 1530 QSTARP=QSTARP*QSTAR(EPSTS(I),XTS(I),KT)
      IF(NROTI1+NROTI2)1526,1528,1526
 1526 QSTARP=QSTARP*(NROTI1+NROTI2)*2.
C The partition function of internal rotors is underestimated
C in QSTAR by a factor of 2.
C
 1528 QPROD=QCENT*FAMC/SIGMAW*QSTARP*EXP(-DE0Z/KT)
      KDIPST=KT/HCM*QPROD/QVRA*EXP(-H0K/KT)
      KREPST=KT/HCM*3.0835E-21/(REDMAS*KT*SQRT(REDMAS*KT))
     **QELA/(QELB*QELC)*QPROD/(QVRB*QVRC)
      WRITE(8,157)TEMP
  157 FORMAT(10X,18HTEMPERATURE IN K: ,F6.1,/)
      WRITE(8,158)
  158 FORMAT(10X,47HPARTITION FUNCTIONS (WITHOUT SYMMETRY NUMBERS):,/)
      WRITE(8,159)
  159 FORMAT(33X,38HREACTANT      FRAGMENT1      FRAGMENT2,/)
      WRITE(8,160)QERA,QERB,QERC
  160 FORMAT(10X,16HQROT (EXTERNAL):,3(5X,1PE10.4))
      WRITE(8,161)QIRA,QIRB,QIRC
  161 FORMAT(10X,16HQROT (INTERNAL):,3(5X,1PE10.4))
      WRITE(8,162)QVIBA,QVIBB,QVIBC
  162 FORMAT(10X,16HQVIB:           ,3(5X,1PE10.4))
      WRITE(8,163)QELA,QELB,QELC
  163 FORMAT(10X,16HQEL:            ,3(5X,1PE10.4))
      WRITE(8,100)
      WRITE(8,164)QCENT
  164 FORMAT(10X,7HQCENT: ,1PE10.4,/)
      WRITE(8,165)
  165 FORMAT(10X,13HDISSOCIATION:,/)
      WRITE(8,166)KDIPST
  166 FORMAT(15X,12HK IN 1/SEC: ,1PE10.4,/)
      WRITE(8,167)
  167 FORMAT(10X,14HRECOMBINATION:,/)
      WRITE(8,168)KREPST
  168 FORMAT(15X,16HK IN CM**3/SEC: ,1PE10.4,/)
      WRITE(8,169)
  169 FORMAT(10X,21HEQUILIBRIUM CONSTANT:,/)
      WRITE(8,170)KDIPST/KREPST
  170 FORMAT(15X,18HK(EQ) IN 1/CM**3: ,1PE10.4,/)
C#######      WRITE(10,183)TEMP,KDIPST,KREPST,KDIPST/KREPST
      WRITE(10,183)TEMP,KDIPST,KREPST*6.02E23,KDIPST/KREPST
  183 FORMAT(8X,F6.1,2X,1PE10.4,2X,1PE10.4,2X,1PE10.4)  
 1600 CONTINUE
 2000 IC=IC+1
      IF(IC.EQ.NUMBER) GOTO 2002
      GOTO 2001
 2002 CONTINUE
      WRITE(6,*)'Program finished'
      END


c=================================================================
      SUBROUTINE DCOUNT(K,AEMAX,LENGTH,NY,MROT,B)
c=================================================================
c This subroutine calculates the density of states for an ensemble of 
c k harmonic vibrations and mrot free one-dimensional rotations. It
c constructs the input data required by the subroutine COUNT.
c The densities are contained in the array mrho.
c (In pure vibrational systems densities are smoothed.)
c********************************************************************
c k: number of vibrational degrees of freedom
c aemax: energy up to which the density is counted
c length: step size of counting
c ny: array of the frequencies
c mrot: number of rotational degrees of freedom
c b: array of the rotational constants
c********************************************************************
      DIMENSION B(10)
      COMMON/A/ MRHO(0:100005)
      COMMON/C/ MT(0:100005)
      COMMON/D/ MR(200)
      DOUBLE PRECISION DGAMMA,MRHO,MT
      REAL NY(100),LENGTH
      N=INT(AEMAX/LENGTH+.5)
      IF(MROT)1,1,2
    1 MT(0)=1.
      DO 5 I=1,N
    5 MT(I)=0.
      DO 3 I=1,K
    3 MR(I)=INT(NY(I)/LENGTH+.5)
      CALL COUNT(K,N)
c
c Smoothing of the densities in pure vibrational systems
c
      I=0
   50 J=0
   20 I=I+1
      IF(MT(I).GT.0..OR.I.EQ.N) GOTO 30
      J=J+1
      GOTO 20
   30 DO 40 IK=I-1,I-(J+1),-1
   40 MT(IK)=MT(I-(J+1))/(J+1)
      IF(I.LT.N) GOTO 50
      MT(N)=MT(N-1)
c
c End of smoothing
c
      DO 60 I=0,N
   60 MRHO(I)=MT(I)/LENGTH
      RETURN
    2 GD=MROT/2.0
      GDZ=GD+1.0
      PRODBI=1.
      DO 9 I=1,MROT
    9 PRODBI=PRODBI*SQRT(B(I))  
      CON=1.772454**MROT/(DGAMMA(DBLE(GDZ))*PRODBI)
      MT(0)=1.
      DO 6 I=1,N
      RIND=1.0*I*LENGTH
      RIN=RIND-LENGTH
    6 MT(I)=CON*(RIND**GD-RIN**GD)
      DO 8 I=1,K
    8 MR(I)=INT(NY(I)/LENGTH+.5)
      CALL COUNT(K,N)
      DO 7 I=0,N
    7 MRHO(I)=MT(I)/LENGTH
      RETURN
      END
C======================================================================      

C======================================================================      
      SUBROUTINE SCOUNT(NF,EPS,X,EMAX,LENGTH)
C**********************************************************************
C This subroutine calculates the number of open channels by convolution
C of Wdisapp(E,J) with the conserved degrees of freedom via the Beyer-
C Swinehart algorithm.
C NF: Number of the degrees of freedom,
C EPS(I): Energy quanta of the degrees of freedom,
C X(I): "Exponents" of the degrees of freedom.
C**********************************************************************
      REAL LENGTH
      DOUBLE PRECISION MT,MW
      DIMENSION EPS(110),X(110)
      COMMON/SCOUN/ MW(0:100005)
      COMMON/C/ MT(0:100005)
      COMMON/D/ MR(200)
      COMMON/SPOLL/NFG,EPSE(110),EPSP(110),XE(110),XP(110),
     *             EPSTS(110),XTS(110)
      N=INT(EMAX/LENGTH+.5)
      MT(0)=1.
      DO 1 I=1,N
      EI=1.0*I*LENGTH
      EIM=EI-LENGTH
    1 MT(I)=WLOOSE(NF,EPS,X,EI)-WLOOSE(NF,EPS,X,EIM)
      K=0
      DO 2 I=1,NF
      IF(XP(I).GT.0.51) GOTO 3
      GOTO 2
    3 K=K+1
      MR(K)=INT(EPS(I)/LENGTH+.5)
    2 CONTINUE
      CALL COUNT(K,N)
      DO 5 I=1,N
    5 MT(I)=MT(I)+MT(I-1)
      DO 4 I=0,N
    4 MW(I)=MT(I)
      RETURN
      END
C======================================================================      
C
C
C======================================================================      
      SUBROUTINE COUNT(K,M)
C*********************************************************************
C BEYER-SWINEHART-count
C This subroutine counts the number of complexions for an ensemble of
C k harmonic oscillators up to an energy of m*deltaE; k, m, and r as
C well as the starting vector t (and of course deltaE) have to be 
C specified in the calling routine.
C
C ref.: S.E.STEIN, B.S.RABINOVITCH: J.Chem.Phys. 58, 2438 (1973)
C*********************************************************************
      COMMON/C/ T(0:100005)
      COMMON/D/ R(200)
      INTEGER R
      DOUBLE PRECISION T
      DO 1 I=1,K
      DO 1 J=R(I),M
    1 T(J)=T(J)+T(J-R(I))
      RETURN
      END
C======================================================================      
 
 
C======================================================================      
      FUNCTION WLOOSE(NF,EPS,X,E)
C********************************************************************
C This function calculates the integrated density of states of the
C disappearing oscillators and conserved internal rotors according to
C eq. (3.3) in JCP 86, 6171 (1987).
C********************************************************************
      COMMON/SPOLL/NFG,EPSE(110),EPSP(110),XE(110),XP(110),
     *             EPSTS(110),XTS(110)
      DOUBLE PRECISION DGAMMA
      DIMENSION EPS(110),X(110)
      GH=1.772454
      WLOOSE=1.
      R=0.
      SUMGAM=1.
      PROGAM=1.
      DO 1 I=1,NF
      IF(XP(I).GT.0.51) GOTO 1
      IF(XE(I).LT.0.51.AND.XP(I).LT.0.51) THEN
          R=R+1.
          SUMGAM=SUMGAM+.5
          WLOOSE=WLOOSE*SQRT(E/EPS(I)+.5)
      ELSE
          PROGAM=PROGAM*DGAMMA(DBLE(1.+X(I)))
          SUMGAM=SUMGAM+X(I)
          WLOOSE=WLOOSE*(E/EPS(I)+.5)**X(I)
      ENDIF
    1 CONTINUE
      WLOOSE=.5*GH**R*PROGAM/DGAMMA(DBLE(SUMGAM))*WLOOSE
C The factor .5 corrects for the fact that the K-rotor is formally
C taken in the subroutine as free internal rotor but has to be taken
C as a disappearing degree of freedom.
C
      RETURN
      END
C
C
      FUNCTION FAME(KZR1,E,JOT,SSTAR,AE,BE)
C**********************************************************************
C This function calculates angular momentum coupling factors Fame(E,J).
C Ref.: JCP 79, 6017 (1983), eqs. A17, C1, C11, prolate symmetric top
C expression (C10) modified after Ber. Bunsenges. Phys. Chem. 98, 1563
C (1994).
C**********************************************************************
      DOUBLE PRECISION DGAMMA,QGAM
      QGAM=DGAMMA(DBLE(SSTAR))*DGAMMA(DBLE(1.5))/DGAMMA(DBLE(SSTAR+.5))
      ARG=AE/E
      ARGP=AE/(AE-BE)
      GOTO(3,4,4,5),KZR1
    3 SSTAR=SSTAR+1.
C Linear molecule, SSTAR one to small, because EPSE(NFG) is not the
C K-rotor, see main program.
      F1=DGAMMA(DBLE(SSTAR+1.))/DGAMMA(DBLE(SSTAR-1.))
     *   *(1+JOT)*(1.+JOT/2.)*ARG*ARG
      FAME=F1/(F1+1.)
      RETURN
    4 FAME=(2.*JOT+1.)*SQRT(ARG)/QGAM
      RETURN
    5 IF(JOT.EQ.0) THEN
          FAME=SQRT(ARG)/(QGAM)
          GOTO 6
      END IF
      X=(AE-BE)*(JOT+1)*(JOT+1)/E
      FAME=SQRT(ARGP)
     *     *TANH(1./(2.*QGAM)*(1.+(2.+SSTAR)/10.*X)*SQRT(X)*2.)
    6 RETURN
      END
C======================================================================      
 
       
C======================================================================      
      SUBROUTINE ETERAT(DMORSE,EZP,DEZ,ALBE,BE,A1,A2,J,ENULLJ)
C********************************************************************
C This subroutine estimates by an iteration procedure the J-dependent
C threshold energy E0(J) with a MORSE type radial part. If enullj=-1
C then the molecule is in a rotationally repulsive state for the J 
C under consideration.
C********************************************************************
C dmorse: dissociation energy of the MORSE-potential
C ezp: zero-point energy of the reaction products
C dez:ezn-ezp-erc/2
C albe: alpha/beta
C be: mean of the two smallest rotational constants of the reactant
C a1, a2: centrifugal parameters (see ref.)
C j: rotational quantum number
C enullj: E0(J)
C
C ref.: J.TROE, J.Chem.Phys. 75, 226 (1981)
C*******************************************************************
C znew: starting parameter for the iteration, should be as large as
C       possible, restriction: the compiler must be able to treat
C       exp(-znew)
C
      B=DEZ/DMORSE
      ZNEW=85.
      IF(J)1004,1005,1004
 1005 ZOLD=ZNEW
      ZNEW=LOG(2.*(1.-EXP(-ZOLD))/(B*ALBE*EXP(-ALBE*ZOLD)))
      IF(ZNEW.GT.ZOLD) THEN
      ENULLJ=DMORSE+EZP
      RETURN
      ELSE
      TOL=ABS(ZNEW-ZOLD)
      IF(TOL.LE..001) GOTO 1002
      GOTO 1005
      ENDIF
 1004 A=BE*J*(J+1)/DMORSE
 1001 ZOLD=ZNEW
      ZNEW=-LOG(A*(A1+2.*A2*ZOLD)/(2.*(1.-EXP(-ZOLD))*(1.+A1*ZOLD
     *+A2*ZOLD*ZOLD)*(1.+A1*ZOLD+A2*ZOLD*ZOLD))+
     *B*ALBE*EXP(-ALBE*ZOLD)/(2.*(1.-EXP(-ZOLD))))
      IF(ZNEW.LT.0.) GOTO 1003
      TOL=ABS(ZNEW-ZOLD)
      IF(TOL.LE..001) GOTO 1002
      GOTO 1001
 1002 ENULLJ=(1.-EXP(-ZNEW))**2+BE*J*(J+1)/(DMORSE*(1.+A1*ZNEW+
     *A2*ZNEW*ZNEW))+B*EXP(-ALBE*ZNEW)
      ENULLJ=ENULLJ*DMORSE+EZP
      RETURN
 1003 ENULLJ=-1.
      RETURN
      END

C======================================================================      
      FUNCTION FANHRO(KANH,E,NVIBR,NYR,DR)
C********************************************************************      
C This function calculates anharmonicity factors for the densities of     
C states of the reactant molecule.
C********************************************************************      
      REAL NYR(100),DR(100)

      GOTO (2,3), KANH+1

    3 NENNER=2*NVIBR-3
      FANHRO=1.
      DO 1 I=1,NVIBR
    1 FANHRO=FANHRO*(1.+(E+NYR(I)/2.)/DR(I)/NENNER)

      RETURN
    2 FANHRO=1.
      RETURN
      END
C======================================================================      


C======================================================================      
      FUNCTION FANHWW(KANH,EI,NVIB1,NVIB2,NY1,NY2,D1,D2)
C********************************************************************
C This function calculates anharmonicity factors for the sums of
C states of the product molecules.
C********************************************************************
      REAL NY1(100),NY2(100),D1(100),D2(100)
      GOTO(2,3),KANH+1
    3 NENNER=2*(NVIB1+NVIB2)-1
      FANHWW=1.
      DO 1 I=1,NVIB1
    1 FANHWW=FANHWW*(1.+(EI+NY1(I)/2.)/D1(I)/NENNER)
      DO 4 I=1,NVIB2
    4 FANHWW=FANHWW*(1.+(EI+NY2(I)/2.)/D2(I)/NENNER)
      RETURN
    2 FANHWW=1.
      RETURN
      END
C======================================================================      
C     
C     
C======================================================================      
      FUNCTION BCENT(JOT)
C**********************************************************************
C This function calculates for a given J the minimum of the lowest
C channel potential, estimating in this way the J-dependence of the 
C equilibrium value for Bcent. The minimum is determined by a Newton-
C iteration of the first derivative. 
C *********************************************************************
      COMMON/BCEN/BE,DMORSE,A1,A2
      COMMON/CENT/ ECENTD
      IF(JOT)1002,1001,1002
 1001 BCENT=BE
      RETURN
 1002 ECENTD=BE*JOT*(JOT+1)/DMORSE
      ZNEW=1E-10
 1003 ZOLD=ZNEW
      ZNEW=ZOLD-F(ZOLD)/FDERIV(ZOLD)
      TOL=ABS(ZNEW-ZOLD)
      IF(TOL.LE.1E-4) GOTO 1005
      GOTO 1003
 1005 BCENT=BE/(1.+A1*ZNEW+A2*ZNEW*ZNEW)
      RETURN
      END
C======================================================================      
 
 
C======================================================================      
      FUNCTION F(Z)
C********************************************************************
C First derivative of the lowest channel potential.
C********************************************************************
      COMMON/BCEN/BE,DMORSE,A1,A2
      COMMON/CENT/ ECENTD
      BRACKT=1.+A1*Z+A2*Z*Z
      F=2.*(1.-EXP(-Z))*EXP(-Z)-ECENTD*(A1+2.*A2*Z)/(BRACKT*BRACKT)
      RETURN
      END
C======================================================================      
 
 
C======================================================================      
      FUNCTION FDERIV(Z)
C********************************************************************
C Second derivative of the lowest channel potential.
C********************************************************************
      COMMON/BCEN/BE,DMORSE,A1,A2
      COMMON/CENT/ ECENTD
      BRACKT=1.+A1*Z+A2*Z*Z
      FACTOR=((A1+2.*A2*Z)**2-A2*BRACKT)/(BRACKT*BRACKT*BRACKT)
      FDERIV=4.*EXP(-2.*Z)-2.*EXP(-Z)+2.*ECENTD*FACTOR
      RETURN
      END
C======================================================================      
 
 
C======================================================================      
      FUNCTION RRHIND(S,V0,E)
C*********************************************************************
C This function calculates the correction factor rhohind(E)/rhofree(E)
C for the densities of states of an ensemble of S oscillators and
C one hindered rotor with an hindrance potential V0.
C*********************************************************************
      INTEGER S
      DOUBLE PRECISION PI,VE,A,GS,B(0:100),RHIND,DGAMMA
      SAVE IC,PI,GS,A,B
      DATA IC,PI /0,1.772453850905516/
      IF(V0.GT.0..AND.E.GT.0.) GOTO 5
      RRHIND=1.
      RETURN
    5 IC=IC+1
      IF(IC.GT.1) GOTO 2
      GS=DGAMMA(DBLE(S+.5))
      A=GS/(PI*DGAMMA(DBLE(S+1.)))
      DO 1 NY=0,S-2
    1 B(NY)=(-1)**NY*GS*(NY+2.5)/(DGAMMA(DBLE(NY+1.))
     *    *DGAMMA(DBLE(S-1.-NY))*PI*(NY+2.)*(NY+1.5))
    2 VE=DBLE(V0)/DBLE(E)
      IF(VE.LT.1.) GOTO 3
      RHIND=A/DSQRT(VE)
      RRHIND=REAL(RHIND)
      RETURN
    3 RHIND=1.
      DO 4 I=0,S-2
      RHIND=RHIND-B(I)*VE**(I+1.5)
    4 RRHIND=REAL(RHIND)
      RETURN
      END
C======================================================================      
       
       
C======================================================================      
      FUNCTION RWHIND(S,V0,E)
C*********************************************************************
C This function calculates the correction factor Whind(E)/Wfree(E)
C for the sums of states of an ensemble of S oscillators and
C one hindered rotor with an hindrance potential V0.
C*********************************************************************
      INTEGER S
      DOUBLE PRECISION PI,VE,A,GS,B(0:100),WHIND,DGAMMA
      SAVE IC,PI,GS,A,B
      DATA IC,PI /0,1.772453850905516/
      IF(V0.GT.0..AND.E.GT.0.) GOTO 5
      RWHIND=1.
      RETURN
    5 IC=IC+1
      IF(IC.GT.1) GOTO 2
      GS=DGAMMA(DBLE(S+1.5))
      A=GS/(PI*DGAMMA(DBLE(S+2.)))
      DO 1 NY=0,S-1
    1 B(NY)=(-1)**NY*GS*(NY+2.5)/(DGAMMA(DBLE(NY+1.))
     *    *DGAMMA(DBLE(S-NY))*PI*(NY+2.)*(NY+1.5))
    2 VE=DBLE(V0)/DBLE(E)
      IF(VE.LT.1.) GOTO 3
      WHIND=A/DSQRT(VE)
      RWHIND=REAL(WHIND)
      RETURN
    3 WHIND=1.
      DO 4 I=0,S-1
      WHIND=WHIND-B(I)*VE**(I+1.5)
    4 RWHIND=REAL(WHIND)
      RETURN
      END
C======================================================================      
       
       
C======================================================================      
      FUNCTION DGAMMA(X)
C********************************************************************
C Calculation of the Gamma-function with double precision.
C******************************************************************** 
      DOUBLE PRECISION PI,DGAMMA,Z,GAMMLN,X
      DATA PI /3.141592654/
      IF (X-1.) 1,2,2
    1 Z=1.-X
      DGAMMA=PI*Z/(DEXP(GAMMLN(1.+Z))*SIN(PI*Z))
      RETURN
    2 DGAMMA=DEXP(GAMMLN(X))
      RETURN
      END
C======================================================================      

C======================================================================      
      FUNCTION GAMMLN(XX)
C********************************************************************
C Calculation of ln(gamma(x)), following W. H. Press et al.,
C 'Numerical Recipes', Cambridge Univ. Press 1989
C********************************************************************
      DOUBLE PRECISION COF(6),STP,HALF,ONE,FPF,X,TMP,SER,XX,
     *GAMMLN
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *-1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
      X=X+ONE
   11 SER=SER+COF(J)/X
      GAMMLN=TMP+LOG(STP*SER)
      RETURN
      END
C======================================================================      
 
 
C======================================================================      
      FUNCTION FQHIND(KT,V0,QFREE,N)
C********************************************************************
C This routine calculates the correction factor between the partition
C function of a hindered and a free internal rotor 
C********************************************************************
      REAL KT,KTV0
      IF(V0)1,2,1
    1 KTV0=KT/V0
      HNYKT=N*SQRT(3.1416/KTV0)/QFREE
      FQHIND=EXP(-1.2*KTV0)/(QFREE*(1.-EXP(-HNYKT)))
     *    +(1.-EXP(-KTV0))**1.2
      RETURN
    2 FQHIND=1.
      RETURN
      END
C======================================================================      

C======================================================================      
      SUBROUTINE SPOL(ALBE,DMORSE)
C*******************************************************************
C This subroutine interpolates between reactant and product 
C eigenvalues.
C Reactant: EPSE,XE
C Channel: EPSTS,XTS
C Product: EPSP,XP
C [Ref.: J. TROE, JPC 79, 6017(1983)]
C*******************************************************************
      COMMON/SPOLL/NFG,EPSE(110),EPSP(110),XE(110),XP(110),
     *             EPSTS(110),XTS(110)
      DO 1 I=1,NFG
      IF(XE(I)-XP(I))3,2,3
    2 C3=4.-1.29*LOG(EPSE(I)/DMORSE)
      EPSTS(I)=EPSP(I)+(EPSE(I)-EPSP(I))*EXP(-C3*ALBE)
      XTS(I)=XE(I)
      GOTO 1
    3 C2=1.24+55.*EPSP(I)/EPSE(I)
      C3=4.-1.29*LOG(EPSE(I)/DMORSE)
      C4=2.8-5.19*LOG(EPSE(I)/DMORSE)
      EN=2.25+.005*EPSE(I)/EPSP(I)
      Y=C2*ALBE
      CALBE=C3*ALBE+C4*ALBE*ALBE*ALBE*ALBE
      YN=EXP(-Y-Y**EN)
      EPSTS(I)=EPSP(I)+(EPSE(I)-EPSP(I))*EXP(-CALBE)
      XTS(I)=XE(I)*YN+XP(I)*(1.-YN)**EN
    1 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE WCOUNT(KZP,NCE,EMAX,LENGTH,J,B1,B2,NVIB1,NVIB2,NY1,
     *NY2,NROTI1,NROTI2)
C**********************************************************************
C This subroutine calculates the overall sum of states by convolution
C of Wl(E,J) with the vibrational degrees of freedom of the fragments
C via the Beyer-Swinehart algorithm followed by the convolution with
C eventually existing free internal rotors of the fragments.
C**********************************************************************
      REAL LENGTH,NY1(100),NY2(100),NY(200)
      DOUBLE PRECISION MT,MW
      COMMON/WCOUN/ MW(0:100005)
      COMMON/C/ MT(0:100005)
      COMMON/D/ MR(200)
      COMMON/ROUT/ RHOIR(0:100005)
      IF(KZP.EQ.1) GOTO 1
      IF(KZP.EQ.2) GOTO 2
      IF(KZP.EQ.3) GOTO 3
      IF(KZP.EQ.4) GOTO 4
      IF(KZP.EQ.5) GOTO 5
    1 CALL WLLA(EMAX,LENGTH,J,B1)
      GOTO 6
    2 CALL WLSA(EMAX,LENGTH,J,B1)
      GOTO 6
    3 CALL WLLL(NCE,EMAX,LENGTH,J,B1,B2)
      GOTO 6
    4 CALL WLSL(NCE,EMAX,LENGTH,J,B1,B2)
      GOTO 6
    5 CALL WLSS(NCE,EMAX,LENGTH,J,B1,B2)
    6 NVIBG=NVIB1+NVIB2
      DO 7 I=1,NVIB1
    7 NY(I)=NY1(I)
      DO 8 I=NVIB1+1,NVIBG
    8 NY(I)=NY2(I-NVIB1)
      DO 9 I=1,NVIBG
    9 MR(I)=INT(NY(I)/LENGTH+.5)
      N=INT(EMAX/LENGTH+.5)
      CALL COUNT(NVIBG,N)
      DO 10 I=1,N
   10 MT(I)=MT(I)+MT(I-1)
      IF(NROTI1.EQ.0.AND.NROTI2.EQ.0) GOTO13
      CALL RHOIRO(N,LENGTH)
      MW(0)=1.
      DO 12 I=1,N
      MW(I)=0.
      DO 11 L=0,I
   11 MW(I)=MW(I)+MT(L)*RHOIR(I-L)
   12 MW(I)=MW(I)*LENGTH
      RETURN
   13 DO 14 I=0,N
   14 MW(I)=MT(I)
      RETURN
      END
C
C
      SUBROUTINE WLLA(AEMAX,LENGTH,J,B)
C**********************************************************************
C This subroutine calculates Wl(E,J) for systems linear/atom and
C creates the starting vector for the Beyer-Swinehart count.
C**********************************************************************
      DIMENSION MT(0:100005)
      COMMON/C/ MT
      REAL LENGTH
      DOUBLE PRECISION MT
      N=INT(AEMAX/LENGTH+.5)
      MT(0)=1.
      DO 1 I=1,N
      E=REAL(I*LENGTH)
      EM=E-LENGTH
    1 MT(I)=WLLIAT(E,J,B)-WLLIAT(EM,J,B)
      RETURN
      END
C
C
      SUBROUTINE WLSA(AEMAX,LENGTH,J,B)
C**********************************************************************
C This subroutine calculates Wl(E,J) for systems spherical top/atom and
C creates the starting vector for the Beyer-Swinehart count.
C**********************************************************************
      DIMENSION MT(0:100005)
      COMMON/C/ MT
      REAL LENGTH
      DOUBLE PRECISION MT
      N=INT(AEMAX/LENGTH+.5)
      MT(0)=1.
      DO 1 I=1,N
      E=REAL(I*LENGTH)
      EM=E-LENGTH
    1 MT(I)=WLSTAT(E,J,B)-WLSTAT(EM,J,B)
      RETURN
      END
C
C
      SUBROUTINE WLLL(NCE,AEMAX,LENGTH,JG,B1,B2)
C**********************************************************************
C This subroutine estimates Wl(E,J) for systems linear/linear and
C creates the starting vector for the Beyer-Swinehart count.
C**********************************************************************
      DIMENSION MT(0:100005)
      COMMON/C/ MT
      REAL LENGTH
      DOUBLE PRECISION MT
      N=INT(AEMAX/LENGTH+.5)
      MT(0)=1.
      IF(NCE)6,6,7
    6 NSCH=0
      DO 1 I=1,N
      E=REAL(I*LENGTH)
      EM=E-LENGTH
      IF(NSCH)2,2,3
    2 EXACT=WLLILI(E,JG,B1,B2)
      ORM=WLHILL(E,B1,B2)
      ARG=(2*JG+1)*WLLOLL(E,B1,B2)/ORM
      CLASS=ORM*POLAT(ARG)
      ABW=ABS(EXACT-CLASS)/EXACT*100.
      IF(ABW.LE.2.) GOTO4
      GOTO5
    4 NSCH=1
    5 MT(I)=EXACT-WLLILI(EM,JG,B1,B2)
      GOTO1
    3 ORM=WLHILL(E,B1,B2)
      ARG=(2*JG+1)*WLLOLL(E,B1,B2)/ORM
      CLASS=ORM*POLAT(ARG)
      ORM=WLHILL(EM,B1,B2)
      ARG=(2*JG+1)*WLLOLL(EM,B1,B2)/ORM
      CLASSN=ORM*POLAT(ARG)
      MT(I)=CLASS-CLASSN
    1 CONTINUE
      RETURN
    7 DO 8 I=1,N
      E=REAL(I*LENGTH)
      EM=E-LENGTH  
    8 MT(I)=WLLILI(E,JG,B1,B2)-WLLILI(EM,JG,B1,B2)
      RETURN
      END
C
C
      SUBROUTINE WLSL(NCE,AEMAX,LENGTH,JG,BST,BLI)
C**********************************************************************
C This subroutine estimates Wl(E,J) for systems spherical top/linear
C and creates the starting vector for the Beyer-Swinehart count.
C**********************************************************************
      DIMENSION MT(0:100005)
      COMMON/C/ MT
      REAL LENGTH
      DOUBLE PRECISION MT
      N=INT(AEMAX/LENGTH+.5)
      MT(0)=1.
      IF(NCE)6,6,7
    6 NSCH=0
      DO 1 I=1,N
      E=REAL(I*LENGTH)
      EM=E-LENGTH
      IF(NSCH)2,2,3
    2 EXACT=WLSTLI(E,JG,BST,BLI)
      ORM=WLHISL(E,BST,BLI)
      ARG=(2*JG+1)*WLLOSL(E,BST,BLI)/ORM
      CLASS=ORM*POLAT(ARG)
      ABW=ABS(EXACT-CLASS)/EXACT*100.
      IF(ABW.LE.2.) GOTO4
      GOTO5
    4 NSCH=1
    5 MT(I)=EXACT-WLSTLI(EM,JG,BST,BLI)
      GOTO1
    3 ORM=WLHISL(E,BST,BLI)
      ARG=(2*JG+1)*WLLOSL(E,BST,BLI)/ORM
      CLASS=ORM*POLAT(ARG)
      ORM=WLHISL(EM,BST,BLI)
      ARG=(2*JG+1)*WLLOSL(EM,BST,BLI)/ORM
      CLASSN=ORM*POLAT(ARG)
      MT(I)=CLASS-CLASSN
    1 CONTINUE
      RETURN
    7 DO 8 I=1,N
      E=REAL(I*LENGTH)
      EM=E-LENGTH  
    8 MT(I)=WLSTLI(E,JG,BST,BLI)-WLSTLI(EM,JG,BST,BLI)
      RETURN
      END
C
C
      SUBROUTINE WLSS(NCE,AEMAX,LENGTH,JG,B1,B2)
C**********************************************************************
C This subroutine estimates Wl(E,J) for systems spherical top/spherical
C top and creates the starting vector for the Beyer-Swinehart count.
C**********************************************************************
      DIMENSION MT(0:100005)
      COMMON/C/ MT
      REAL LENGTH
      DOUBLE PRECISION MT
      N=INT(AEMAX/LENGTH+.5)
      MT(0)=1.
      IF(NCE)6,6,7
    6 NSCH=0
      DO 1 I=1,N
      E=REAL(I*LENGTH)
      EM=E-LENGTH
      IF(NSCH)2,2,3
    2 EXACT=WLSTST(E,JG,B1,B2)
      ORM=WLHISS(E,B1,B2)
      ARG=(2*JG+1)*WLLOSS(E,B1,B2)/ORM
      CLASS=ORM*POLAT(ARG)
      ABW=ABS(EXACT-CLASS)/EXACT*100.
      IF(ABW.LE.2.) GOTO4
      GOTO5
    4 NSCH=1
    5 MT(I)=EXACT-WLSTST(EM,JG,B1,B2)
      GOTO1
    3 ORM=WLHISS(E,B1,B2)
      ARG=(2*JG+1)*WLLOSS(E,B1,B2)/ORM
      CLASS=ORM*POLAT(ARG)
      ORM=WLHISS(EM,B1,B2)
      ARG=(2*JG+1)*WLLOSS(EM,B1,B2)/ORM
      CLASSN=ORM*POLAT(ARG)
      MT(I)=CLASS-CLASSN
    1 CONTINUE
      RETURN
    7 DO 8 I=1,N
      E=REAL(I*LENGTH)
      EM=E-LENGTH  
    8 MT(I)=WLSTST(E,JG,B1,B2)-WLSTST(EM,JG,B1,B2)
      RETURN
      END
C
C
      FUNCTION WLLIAT (E,J,B)
C*********************************************************************
C This function calculates exact values of Wl(E,J) for systems linear/
C atom.
C*********************************************************************
      JMAX=INT(SQRT(.25+E/B)-.5)
      IF(J.GT.JMAX) GOTO 1
      WLLIAT=REAL((2*J+1)*(JMAX+1)-J*(J+1))
      RETURN
    1 WLLIAT=REAL((JMAX+1)*(JMAX+1))
      RETURN
      END
C
C
      FUNCTION WLSTAT(E,J,B)
C********************************************************************
C This function calculates exact values of Wl(E,J) for systems
C spherical top/atom.
C********************************************************************
      JMAX=INT(SQRT(.25+E/B)-.5)
      IF(J.GT.JMAX) GOTO 1
      WLSTAT=REAL((2*J+1)*(JMAX+1)*(JMAX+1)-(2*J+1)*(J+1)*J/3)
      RETURN
    1 WLSTAT=REAL((JMAX+1)*(2*JMAX+3)*(2*JMAX+1)/3)
      RETURN
      END
C
C
      FUNCTION WLLILI(E,JG,B1,B2)
C********************************************************************
C This function counts exact values of Wl(E,J) for systems linear/
C linear.
C********************************************************************
      WLLILI=0.
      J1MAX=INT(SQRT(.25+E/B1)-.5)
      DO 1 J1=0,J1MAX
      J2MAX=INT(SQRT(.25+(E-B1*J1*(J1+1))/B2)-.5)
      DO 1 J2=0,J2MAX
      DO 1 J=IABS(J1-J2),J1+J2
      DO 1 L=IABS(JG-J),JG+J
    1 WLLILI=WLLILI+1.
      RETURN
      END
C
C
      FUNCTION WLLOLL(E,B1,B2)
C********************************************************************
C This function calculates classically  Wl(E,J=0) for systems linear/
C linear.
C********************************************************************
      WLLOLL=2.*E*SQRT(E)*(SQRT(B1)+SQRT(B2)-SQRT(B1+B2))/
     *(3.*B1*B2)
      RETURN
      END
C
C
      FUNCTION WLHILL(E,B1,B2)
C*********************************************************************
C This function calculates classically Wl(E,highJ) for systems linear/
C linear.
C*********************************************************************
      WLHILL=E*E/(2.*B1*B2)
      RETURN
      END
C
C
      FUNCTION WLSTLI(E,JG,BST,BLI)
C********************************************************************
C This function counts exact values of Wl(E,J) for systems spherical
C top/linear.
C********************************************************************
      WLSTLI=0.
      J1MAX=INT(SQRT(.25+E/BST)-.5)
      DO 1 J1=0,J1MAX
      J2MAX=INT(SQRT(.25+(E-BST*J1*(J1+1))/BLI)-.5)
      SWL=REAL(2*J1+1)
      DO 1 J2=0,J2MAX
      DO 1 J=IABS(J1-J2),J1+J2
      DO 1 L=IABS(JG-J),JG+J
    1 WLSTLI=WLSTLI+SWL
      RETURN
      END
C
C
      FUNCTION WLLOSL(E,BST,BLI)
C**********************************************************************
C This function calculates classically Wl(E,J=0) for systems
C spherical top/linear.
C**********************************************************************
      WLLOSL=E*E*ASIN(SQRT(BST/(BST+BLI)))/(2.*BST*SQRT(BST)*SQRT(BLI))
      RETURN
      END
C
C
      FUNCTION WLHISL(E,BST,BLI)
C**********************************************************************
C This function calculates classically Wl(E,highJ) for systems
C spherical top/linear.
C**********************************************************************
      WLHISL=8.*E*E*SQRT(E)/(15.*BST*SQRT(BST)*BLI)
      RETURN
      END
C
C
      FUNCTION WLSTST(E,JG,BST1,BST2)
C**********************************************************************
C This function counts exact values of Wl(E,J) for systems spherical
C top/spherical top.
C**********************************************************************
      WLSTST=0.
      J1MAX=INT(SQRT(.25+E/BST1)-.5)
      DO 1 J1=0,J1MAX
      J2MAX=INT(SQRT(.25+(E-BST1*J1*(J1+1))/BST2)-.5)
      SWL1=REAL(2*J1+1)
      DO 1 J2=0,J2MAX
      SWL=REAL(2*J2+1)*SWL1
      DO 1 J=IABS(J1-J2),J1+J2
      DO 1 L=IABS(JG-J),JG+J
    1 WLSTST=WLSTST+SWL
      RETURN
      END
C
C
      FUNCTION WLLOSS(E,BST1,BST2)
C**********************************************************************
C This function calculates classically Wl(E,J=0) for systems spherical
C top/spherical top.
C**********************************************************************
      WLLOSS=8.*E*E*SQRT(E)/(15.*BST1*BST2*SQRT(BST1+BST2))
      RETURN
      END
C
C
      FUNCTION WLHISS(E,BST1,BST2)
C**********************************************************************
C This function calculates classically Wl(E,highJ) for systems
C spherical top/spherical top.
C**********************************************************************
      WLHISS=3.1416*E*E*E/(6.*BST1*BST2*SQRT(BST1*BST2))
      RETURN
      END
C
C
      FUNCTION POLAT(X)
C**********************************************************************
C This function interpolates Wl(E,J) between Wl(E,J=0) and Wl(E,highJ)
C in a doubly reduced representation.
C**********************************************************************
      POLAT=TANH(X)*(.08*SQRT(X)*EXP(-1.1*(X-1.)*(X-1.))+1.)
      RETURN
      END
      SUBROUTINE RHOIRO(N,LENGTH)
C***********************************************************************
C This function calculates the classical density of states of the
C internal rotations of the fragments.
C***********************************************************************
      REAL LENGTH
      DOUBLE PRECISION DGAMMA
      COMMON/RHOI/ NROTI1,NROTI2,AIR1(10),AIR2(10)
      COMMON/ROUT/ RHOIR(0:100005)
      R=REAL(NROTI1+NROTI2)
      RHOIR(0)=0.
      SQRB=1.
      DO 1 I=1,NROTI1
    1 SQRB=SQRB*SQRT(AIR1(I))
      DO 2 I=1,NROTI2
    2 SQRB=SQRB*SQRT(AIR2(I))
      GD=R/2.-1.
      DO 3 I=1,N
      E=I*LENGTH
    3 RHOIR(I)=1.7725**R/(DGAMMA(DBLE(R/2.))*SQRB)*E**GD
      RETURN
      END
C
C
      FUNCTION QSTAR(EPS,X,KT)
C*********************************************************************
C This function calculates (pseudo) partition functions for the
C "transition state", see JCP 75, 226 (1981), eqs. (6.4) and (6.6).
C*********************************************************************
      REAL KT
      DOUBLE PRECISION DGAMMA
      QSTAR=REAL(DGAMMA(DBLE(1.+X)))*(1.-EXP(-EPS/KT))**(-X)
      RETURN
      END
C
C
      FUNCTION AWR(E,BETA)
C*********************************************************************
C This function calculates the Whitten-Rabinovitch a(E).
C*********************************************************************
      IF (E.GT.1.) THEN
                W=10**(-1.0506*E**.25)
          ELSE
                W=1/(5.*E+2.73*SQRT(E)+3.51)
      ENDIF
      AWR=1.-BETA*W
      RETURN
      END

