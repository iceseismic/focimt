//---------------------------------------------------------------------------
#include <trilib/fortranmath.h>
#include <trilib/georoutines.h>
#include "usmtcore.h"
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
using namespace Taquart::UsmtCore;
using namespace Taquart;

namespace Taquart {
  namespace UsmtCore {
    int NDAE[10] = { 0, 36, 36, 32, 32, 24, 24, 16, 8, 4 };
    double U[MAXCHANNEL + 1];
    double AZM[MAXCHANNEL + 1];
    double TKF[MAXCHANNEL + 1];
    double GA[MAXCHANNEL + 1][3 + 1];
    double A[MAXCHANNEL + 1][6 + 1];
    double FIJ[3 + 1][3 + 1][MAXCHANNEL + 1];
    double RM[6 + 1][3 + 1];
    double COV[6 + 1][6 + 1][3 + 1];
    int RO[MAXCHANNEL + 1];
    int VEL[MAXCHANNEL + 1];
    int R[MAXCHANNEL + 1];
    double UTH[MAXCHANNEL + 1];
    int N = 0;
    double TROZ = 0.0;
    double QSD = 0.0;
    double QF = 0.0;
    bool FSTCLL = true;
    int ICOND = 0;
    Taquart::FaultSolution Solution[4];
    int ISTA = 1;
    int * ThreadProgress;
  } // namespace UsmtCore
} // namespace Foci

//---------------------------------------------------------------------------
void TransferSolution(Taquart::SolutionType AType,
    std::list<Taquart::FaultSolution> &ASolution) {
  ASolution.push_back(Taquart::UsmtCore::Solution[int(AType)]);
}

//---------------------------------------------------------------------------
void TransferSolution(Taquart::SolutionType AType,
    Taquart::FaultSolution &ASolution) {
  ASolution = Taquart::UsmtCore::Solution[int(AType)];
}

//---------------------------------------------------------------------------
void USMTCore(Taquart::NormType ANormType, int QualityType,
    Taquart::SMTInputData &InputData, int * const AThreadProgress) {
  int IEXP = 0;
  ThreadProgress = AThreadProgress;
  PROGRESS(0, 350);
  RDINP(InputData);
  ANGGA();
  JEZ();
  switch (ANormType) {
    case Taquart::ntL1:
      MOM2(false, QualityType);
      SIZEMM(IEXP);
      MOM1(IEXP, QualityType);
      break;
    case Taquart::ntL2:
      MOM2(true, QualityType);
      break;
  }
  PROGRESS(360, 350);
}

//---------------------------------------------------------------------------
void Taquart::UsmtCore::MOM1(int &IEXP, int QualityType) {
  //      SUBROUTINE MOM1(IEXP)
  //      CHARACTER PS(80),TITLE*40
  //      REAL U(80),ARR(80),AZM(80),TKF(80),LLA(3),HA(2)
  //      INTEGER RO(80),VEL(80),R(80)
  //      COMMON/MDATA/ PS,U,ARR,AZM,TKF,RO,VEL,R,TITLE,N,TROZ
  //      COMMON/GAGAGA/ GA(80,3)
  //      COMMON/RESTEX/ RESULT
  //      CHARACTER*56 RESULT(34)
  //      COMMON/MOMNT/ RM(6,3)
  //      DOUBLE PRECISION DHELP1,DHELP2,COV
  //      COMMON/PDATA/ A(80,6)
  //      DIMENSION EQM(3),AA(80,6),B(6),H(5),IW(80),PA(3),BB(6)
  //      REAL RM0(3),RMT(3),PCLVD(2),PDBCP(2)
  double AA[MAXCHANNEL + 1][6 + 1];
  unsigned int IW[MAXCHANNEL + 1];
  double PA[3 + 1];
  Zero(PA, 4);
  double EQM[3 + 1];
  Zero(EQM, 4);
  double BB[6 + 1];
  Zero(BB, 7);
  double RM0[3 + 1];
  Zero(RM0, 4);
  double RMT[3 + 1];
  Zero(RMT, 4);
  double RMERR[3 + 1];
  Zero(RMERR, 4);
  double PEXPL[3 + 1];
  Zero(PEXPL, 4);
  double PCLVD[3 + 1];
  Zero(PCLVD, 4);
  double PDBCP[3 + 1];
  Zero(PDBCP, 4);
  double B[7];
  Zero(B, 7);
  double H[5 + 1];
  Zero(H, 6);
  double ALF = 0.0;
  double HELP = 0.0;
  double PEXPLO = 0.0;
  double RMAG = 0.0;
  double MAGN[4];

  //      PI=4.*ATAN(1.)
  const double PI = 4.0 * atan(1.0);

  //      DO 3 I=1,N
  for (int i = 1; i <= N; i++) {
    //      IW(I)=0
    //      IF((PS(I).EQ.'S').OR.(PS(I).EQ.'s')) IW(I)=1
    //      IF((PS(I).EQ.'H').OR.(PS(I).EQ.'h')) IW(I)=2
    //      ALF=FLOAT(VEL(I))
    //      HELP=4.*PI*FLOAT(RO(I))*ALF*ALF*ALF*FLOAT(R(I))*TROZ
    IW[i] = 0;
    //if(PS[i] == 'S' || PS[i] == 's') IW[i] = 1;
    //if(PS[i] == 'H' || PS[i] == 'h') IW[i] = 2;
    ALF = VEL[i];
    /* DONE  5 -c2.4.14 : Problem z kalibracj� (USMTCORE) */
    HELP = 4.0 * PI * double(RO[i]) * ALF * ALF * ALF * double(R[i]);

    //      IF(IW(I).NE.0) GO TO 5
    if (IW[i] == 0) {
      //C     For P:
      //      A(I,1)=GA(I,1)*GA(I,1)/HELP
      //      A(I,2)=2.*GA(I,1)*GA(I,2)/HELP
      //      A(I,3)=2.*GA(I,1)*GA(I,3)/HELP
      //      A(I,4)=GA(I,2)*GA(I,2)/HELP
      //      A(I,5)=2.*GA(I,2)*GA(I,3)/HELP
      //      A(I,6)=GA(I,3)*GA(I,3)/HELP
      //      GO TO 3
      A[i][1] = GA[i][1] * GA[i][1] / HELP;
      A[i][2] = 2.0 * GA[i][1] * GA[i][2] / HELP;
      A[i][3] = 2.0 * GA[i][1] * GA[i][3] / HELP;
      A[i][4] = GA[i][2] * GA[i][2] / HELP;
      A[i][5] = 2.0 * GA[i][2] * GA[i][3] / HELP;
      A[i][6] = GA[i][3] * GA[i][3] / HELP;
    }

    /* TODO -o3.1.19 : Code for SV and SH is switched off by default. */
    /*
     else
     {
     //C     For SV:
     //    5 SVANG=REAL(ARR(I))*PI/180.
     //      PA(3)=-SIN(SVANG)
     //      SKAL=-COS(SVANG)/SQRT(GA(I,1)*GA(I,1)+GA(I,2)*GA(I,2))
     //      PA(1)=GA(I,1)*SKAL
     //      PA(2)=GA(I,2)*SKAL
     double SVANG = double(ARR[i]) * PI / 180.0;
     PA[3] = -sin(SVANG);
     const double SKAL = -cos(SVANG) / sqrt(GA[i][1] * GA[i][1] + GA[i][2] * GA[i][2]);
     PA[1] = GA[i][1] * SKAL;
     PA[2] = GA[i][2] * SKAL;

     //      IF(IW(I).EQ.2) GO TO 6
     if(IW[i] != 2)
     {
     //      A(I,1)=GA(I,1)*PA(1)/HELP
     //      A(I,2)=(GA(I,1)*PA(2)+GA(I,2)*PA(1))/HELP
     //      A(I,3)=(GA(I,1)*PA(3)+GA(I,3)*PA(1))/HELP
     //      A(I,4)=GA(I,2)*PA(2)/HELP
     //      A(I,5)=(GA(I,2)*PA(3)+GA(I,3)*PA(2))/HELP
     //      A(I,6)=GA(I,3)*PA(3)/HELP
     //      GO TO 3
     A[i][1]=GA[i][1]*PA[1]/HELP;
     A[i][2]=(GA[i][1]*PA[2]+GA[i][2]*PA[1])/HELP;
     A[i][3]=(GA[i][1]*PA[3]+GA[i][3]*PA[1])/HELP;
     A[i][4]=GA[i][2]*PA[2]/HELP;
     A[i][5]=(GA[i][2]*PA[3]+GA[i][3]*PA[2])/HELP;
     A[i][6]=GA[i][3]*PA[3]/HELP;
     }
     else
     {
     //C     For SH:
     //    6 LLA(3)=-COS(SVANG)
     //      SKAL=SIN(SVANG)/SQRT(GA(I,1)*GA(I,1)+GA(I,2)*GA(I,2))
     //      LLA(1)=GA(I,1)*SKAL
     //      LLA(2)=GA(I,2)*SKAL
     //      HA(1)=LLA(2)*PA(3)-LLA(3)*PA(2)
     //      HA(2)=-LLA(1)*PA(3)+LLA(3)*PA(1)
     //      A(I,1)=GA(I,1)*HA(1)/HELP+LLA(3)*GA(I,1)*PA(1)/HELP
     //      A(I,2)=(GA(I,1)*HA(2)+GA(I,2)*HA(1))/HELP
     //     $+LLA(3)*(GA(I,1)*PA(2)+GA(I,2)*PA(1))/HELP
     //      A(I,3)=GA(I,3)*HA(1)/HELP+LLA(3)*(GA(I,1)*PA(3)+GA(I,3)*PA(1))
     //     $/HELP
     //      A(I,4)=GA(I,2)*HA(2)/HELP+LLA(3)*GA(I,2)*PA(2)/HELP
     //      A(I,5)=GA(I,3)*HA(2)/HELP+LLA(3)*(GA(I,2)*PA(3)+GA(I,3)*PA(2))
     //     $/HELP
     //      A(I,6)=+LLA(3)*GA(I,3)*PA(3)/HELP
     double LLA[4];
     double HA[4];

     LLA[3]=-cos(SVANG);
     double SKAL = sin(SVANG)/sqrt(GA[i][1]*GA[i][1]+GA[i][2]*GA[i][2]);
     LLA[1]=GA[i][1]*SKAL;
     LLA[2]=GA[i][2]*SKAL;
     HA[1]=LLA[2]*PA[3]-LLA[3]*PA[2];
     HA[2]=-LLA[1]*PA[3]+LLA[3]*PA[1];
     A[i][1]=GA[i][1]*HA[1]/HELP+LLA[3]*GA[i][1]*PA[1]/HELP;
     A[i][2]=(GA[i][1]*HA[2]+GA[i][2]*HA[1])/HELP+LLA[3]*(GA[i][1]*PA[2]+GA[i][2]*PA[1])/HELP;
     A[i][3]=GA[i][3]*HA[1]/HELP+LLA[3]*(GA[i][1]*PA[3]+GA[i][3]*PA[1])/HELP;
     A[i][4]=GA[i][2]*HA[2]/HELP+LLA[3]*GA[i][2]*PA[2]/HELP;
     A[i][5]=GA[i][3]*HA[2]/HELP+LLA[3]*(GA[i][2]*PA[3]+GA[i][3]*PA[2])/HELP;
     A[i][6]=+LLA[3]*GA[i][3]*PA[3]/HELP;
     }
     }
     */
  }
  //    3 CONTINUE
  //      CALL GSOL(B,iexp)
  //C     Find and sort eigenvalues:
  //      CALL EIG3(B,0,EQM)
  GSOL(B, IEXP);
  EIG3(B, 0, EQM);

  //      EQQ1=EQM(1)
  //      EQQ2=EQM(2)
  //      EQQ3=EQM(3)
  double EQQ1 = EQM[1];
  double EQQ2 = EQM[2];
  double EQQ3 = EQM[3];

  //      do 1004 i=1,6
  // 1004 rm(i,1)=b(i)
  for (int i = 1; i <= 6; i++)
    RM[i][1] = B[i];

  //      RMM=(EQM(1)+EQM(2)+EQM(3))/3.
  //      RMY=EQM(1)
  //      IF(EQM(2).LT.RMY) RMY=EQM(2)
  //      IF(EQM(3).LT.RMY) RMY=EQM(3)
  //      RMX=EQM(3)
  //double RMM = (EQM[1]+EQM[2]+EQM[3]) / 3.0;
  double RMY = EQM[1];
  if (EQM[2] < RMY) RMY = EQM[2];
  if (EQM[3] < RMY) RMY = EQM[3];
  double RMX = EQM[3];

  //      IF(EQM(1).GT.RMX) RMX=EQM(1)
  //      IF(EQM(2).GT.RMX) RMX=EQM(2)
  //      IF((EQM(1).NE.RMX).AND.(EQM(1).NE.RMY)) RMZ=EQM(1)
  //      IF((EQM(2).NE.RMX).AND.(EQM(2).NE.RMY)) RMZ=EQM(2)
  //      IF((EQM(3).NE.RMX).AND.(EQM(3).NE.RMY)) RMZ=EQM(3)
  if (EQM[1] > RMX) RMX = EQM[1];
  if (EQM[2] > RMX) RMX = EQM[2];
  //double RMZ = 0.0;
  //if(EQM[1] != RMX && EQM[1] != RMY) RMZ = EQM[1];
  //if(EQM[2] != RMX && EQM[2] != RMY) RMZ = EQM[2];
  //if(EQM[3] != RMX && EQM[3] != RMY) RMZ = EQM[3];

  //C     Finds scalar seismic moment:
  //      DO 9 I=1,6
  //    9 BB(I)=RM(I,1)
  for (int i = 1; i <= 6; i++)
    BB[i] = RM[i][1];

  //      CALL EIG3(BB,0,EQM)
  //      EQM(1)=ABS(EQM(1))*1.E-12
  //      EQM(2)=ABS(EQM(2))*1.E-12
  //      EQM(3)=ABS(EQM(3))*1.E-12
  //      HELP1=ABS(EQM(1)-EQM(2))
  //      HELP2=ABS(EQM(1)-EQM(3))
  //      HELP3=ABS(EQM(2)-EQM(3))
  EIG3(BB, 0, EQM);
  EQM[1] = fabs(EQM[1]) * 1.0E-12;
  EQM[2] = fabs(EQM[2]) * 1.0E-12;
  EQM[3] = fabs(EQM[3]) * 1.0E-12;
  double HELP1 = fabs(EQM[1] - EQM[2]);
  double HELP2 = fabs(EQM[1] - EQM[3]);
  double HELP3 = fabs(EQM[2] - EQM[3]);

  //      EU=SQRT(.5*(EQM(1)*EQM(1)+EQM(2)*EQM(2)+EQM(3)*EQM(3)))
  //      EC=AMAX1(HELP1,HELP2)
  //      EC=AMAX1(EC,HELP3)
  //      RMT(1)=AMAX1(EC,EU)*1.E+12
  //      RM0(1)=AMIN1(EC,EU)*1.E+12
  double EU = sqrt(0.5 * (EQM[1] * EQM[1] + EQM[2] * EQM[2] + EQM[3] * EQM[3]));
  double EC = Taquart::amax1(HELP1, HELP2);
  EC = amax1(EC, HELP3);
  RMT[1] = amax1(EC, EU) * 1.0e+12;
  RM0[1] = amin1(EC, EU) * 1.0e+12;

  //C     Error calculation:
  //      DO 401 I=1,N
  for (int i = 1; i <= N; i++) {
    //      EPS=U(I)
    double EPS = U[i];

    //      DO 406 J=1,6
    //  406 EPS=EPS-A(I,J)*RM(J,1)
    for (int j = 1; j <= 6; j++)
      EPS = EPS - A[i][j] * RM[j][1];

    //      SAI=0.
    double SAI = 0.0;

    //      DO 403 J=1,6
    //  403 SAI=SAI+ABS(A(I,J))
    for (int j = 1; j <= 6; j++)
      SAI = SAI + fabs(A[i][j]);

    //      DO 404 J=1,6
    //  404 AA(I,J)=EPS/SAI*SIGN(1.,A(I,J))
    for (int j = 1; j <= 6; j++)
      AA[i][j] = EPS / SAI * sign(1.0, A[i][j]);

    //      DO 401 J=1,6
    //      AA(I,J)=AA(I,J)+RM(J,1)
    for (int j = 1; j <= 6; j++)
      AA[i][j] = AA[i][j] + RM[j][1];

  }
  //  401 CONTINUE

  //      COV=DBLE(0.) 
  double COV = 0.0;
  //      DO 4384 I=1,6
  //      DO 4384 J=1,6
  //      DO 405 K=1,N
  for (int i = 1; i <= 6; i++) {
    for (int j = 1; j <= 6; j++) {
      for (int k = 1; k <= N; k++) {
        //      DHELP1=DBLE(AA(K,I)-RM(I,1))*1.E-12
        //      DHELP2=DBLE(AA(K,J)-RM(J,1))*1.E-12
        //  405 COV=COV+DHELP1*DHELP2
        double DHELP1 = (AA[k][i] - RM[i][1]) * 1.0E-12;
        double DHELP2 = (AA[k][j] - RM[j][1]) * 1.0E-12;
        COV = COV + DHELP1 * DHELP2;
      }
    }
  }
  // 4384 CONTINUE

  //      COV=COV/36.
  //      DHELP1=DABS(COV)/DBLE(FLOAT(N-1))
  //      SIG=SNGL(DSQRT(DHELP1))*1.E+12
  COV = COV / 36.0;
  double DHELP1 = fabs(COV) / double(N - 1);
  double SIG = sqrt(DHELP1) * 1.0e+12;
  RMERR[1] = SIG;

  //C     Output:
  //      WRITE(75,4999)
  // 4999 FORMAT(1X,/,72(1H-))
  //      WRITE(RESULT(1),5011) TITLE
  // 5011 FORMAT('L1 tensor for : ',A40)
  //      WRITE(75,5000) RESULT(1)
  // 5000 FORMAT(1H ,A56)
  //      WRITE(RESULT(2),5012)
  // 5012 FORMAT('  Full solution :',39X)
  //      WRITE(75,5000) RESULT(2)
  //      WRITE(RESULT(3),5013) RM(1,1),RM(2,1),RM(3,1)
  //      WRITE(75,5000) RESULT(3)
  //      WRITE(RESULT(4),5013) RM(2,1),RM(4,1),RM(5,1)
  //      WRITE(75,5000) RESULT(4)
  //      WRITE(RESULT(5),5013) RM(3,1),RM(5,1),RM(6,1)
  //      WRITE(75,5000) RESULT(5)
  // 5013 FORMAT(11X,E10.3,2X,E10.3,2X,E10.3,11X)
  //      WRITE(RESULT(6),5014) TROZ,RM0(1),RMT(1),SIG
  // 5014 FORMAT(2X,'T0=',F8.4,2X,'M0=',E9.3,2X,'MT=',E9.3,2X,'ERR=',E9.3)
  //      WRITE(75,5000) RESULT(6)
#ifdef USMTCORE_DEBUG
  std::cout << "L1 Tensor:" << std::endl;
  std::cout << RM[1][1] << RM[2][1] << RM[3][1] << std::endl;
  std::cout << RM[2][1] << RM[4][1] << RM[5][1] << std::endl;
  std::cout << RM[3][1] << RM[5][1] << RM[6][1] << std::endl;
  std::cout << "T0 = " << TROZ << " M0 = " << RM0[1] << " MT = " << RMT[1] << " ERR = " << SIG;
#endif

  //      CALL EIGGEN(EQQ1,EQQ2,EQQ3,PEXPL,PCLVD(1),PDBCP(1))
  //      RMAG=.6667*ALOG10(RM0(1))-6.0
  //      IF(RMAG.GE.100.) RMAG=99.9
  //      IF(RMAG.LE.-10.) RMAG=-9.9
  EIGGEN(EQQ1, EQQ2, EQQ3, PEXPLO, PCLVD[1], PDBCP[1]);
  RMAG = 0.6667 * alog10(RM0[1]) - 6.0;
  if (RMAG >= 100.0) RMAG = 99.9;
  if (RMAG <= -10.0) RMAG = -9.9;
  MAGN[1] = RMAG;

  //      WRITE(RESULT(7),5015) PEXPL,PCLVD(1),PDBCP(1),RMAG
  // 5015 FORMAT(2X,'Expl.=',F6.1,'%   CLVD.=',F6.1,'%   DBCP.=',F6.1,
  //     $'%   M=',F4.1)
  //      WRITE(75,5000) RESULT(7)
#ifdef USMTCORE_DEBUG
  std::cout << " EXPL = " << PEXPL << " PCLVD = " << PCLVD[1];
  std::cout << " PDBCP = " << PDBCP[1] << " Mw = " << RMAG << std::endl;
#endif
  PEXPL[1] = PEXPLO;

  //C     Postep wyswietlany przez rutyny GSOL?
  //C     ### TRACE NULL TENSOR:

  //      CALL GSOL5(H,iexp)
  //      DO 55 I=1,5
  //   55 RM(I,2)=H(I)
  //      RM(6,2)=-RM(1,2)-RM(4,2)
  GSOL5(H, IEXP);
  for (int i = 1; i <= 5; i++)
    RM[i][2] = H[i];
  RM[6][2] = -RM[1][2] - RM[4][2];

  //C     Find and sort eigenvalues:
  //      DO 1003 I=1,6
  // 1003 B(I)=RM(I,2)
  //      CALL EIG3(B,0,EQM)
  for (int i = 1; i <= 6; i++)
    B[i] = RM[i][2];
  EIG3(B, 0, EQM);

  //      EQQ1=EQM(1)
  //      EQQ2=EQM(2)
  //      EQQ3=EQM(3)
  //      RMM=(EQM(1)+EQM(2)+EQM(3))/3.
  //      RMY=EQM(1)
  //      IF(EQM(2).LT.RMY) RMY=EQM(2)
  //      IF(EQM(3).LT.RMY) RMY=EQM(3)
  //      RMX=EQM(3)
  EQQ1 = EQM[1];
  EQQ2 = EQM[2];
  EQQ3 = EQM[3];
  //RMM = (EQM[1] + EQM[2] + EQM[3]) / 3.0;
  RMY = EQM[1];
  if (EQM[2] < RMY) RMY = EQM[2];
  if (EQM[3] < RMY) RMY = EQM[3];
  RMX = EQM[3];

  //      IF(EQM(1).GT.RMX) RMX=EQM(1)
  //      IF(EQM(2).GT.RMX) RMX=EQM(2)
  //      IF((EQM(1).NE.RMX).AND.(EQM(1).NE.RMY)) RMZ=EQM(1)
  //      IF((EQM(2).NE.RMX).AND.(EQM(2).NE.RMY)) RMZ=EQM(2)
  //      IF((EQM(3).NE.RMX).AND.(EQM(3).NE.RMY)) RMZ=EQM(3)
  if (EQM[1] > RMX) RMX = EQM[1];
  if (EQM[2] > RMX) RMX = EQM[2];
  //if(EQM[1] != RMX && EQM[1] != RMY) RMZ = EQM[1];
  //if(EQM[2] != RMX && EQM[2] != RMY) RMZ = EQM[2];
  //if(EQM[3] != RMX && EQM[3] != RMY) RMZ = EQM[3];

  //C     Finds scalar seismic moment:
  //      DO 2009 I=1,6
  // 2009 BB(I)=RM(I,2)
  for (int i = 1; i <= 6; i++)
    BB[i] = RM[i][2];

  //      CALL EIG3(BB,0,EQM)
  //      EQM(1)=ABS(EQM(1))*1.E-12
  //      EQM(2)=ABS(EQM(2))*1.E-12
  //      EQM(3)=ABS(EQM(3))*1.E-12
  //      HELP1=ABS(EQM(1)-EQM(2))
  //      HELP2=ABS(EQM(1)-EQM(3))
  //      HELP3=ABS(EQM(2)-EQM(3))
  //      EU=SQRT(.5*(EQM(1)*EQM(1)+EQM(2)*EQM(2)+EQM(3)*EQM(3)))
  //      EC=AMAX1(HELP1,HELP2)
  //      EC=AMAX1(EC,HELP3)
  //      RMT(2)=AMAX1(EC,EU)*1.E+12
  //      RM0(2)=AMIN1(EC,EU)*1.E+12

  EIG3(BB, 0, EQM);
  EQM[1] = fabs(EQM[1]) * 1.0E-12;
  EQM[2] = fabs(EQM[2]) * 1.0E-12;
  EQM[3] = fabs(EQM[3]) * 1.0E-12;
  HELP1 = fabs(EQM[1] - EQM[2]);
  HELP2 = fabs(EQM[1] - EQM[3]);
  HELP3 = fabs(EQM[2] - EQM[3]);
  EU = sqrt(0.5 * (EQM[1] * EQM[1] + EQM[2] * EQM[2] + EQM[3] * EQM[3]));
  EC = amax1(HELP1, HELP2);
  EC = amax1(EC, HELP3);
  RMT[2] = amax1(EC, EU) * 1.0E+12;
  RM0[2] = amin1(EC, EU) * 1.0E+12;

  //C     Error calculation:
  //      DO 2401 I=1,N
  for (int i = 1; i <= N; i++) {
    //      EPS=U(I)
    double EPS = U[i];
    //      DO 2406 J=1,6
    // 2406 EPS=EPS-A(I,J)*RM(J,2)
    for (int j = 1; j <= 6; j++)
      EPS = EPS - A[i][j] * RM[j][2];

    //      SAI=0.
    double SAI = 0.0;

    //      DO 2403 J=1,6
    // 2403 SAI=SAI+ABS(A(I,J))
    for (int j = 1; j <= 6; j++)
      SAI = SAI + fabs(A[i][j]);

    //      DO 2404 J=1,6
    // 2404 AA(I,J)=EPS/SAI*SIGN(1.,A(I,J))
    for (int j = 1; j <= 6; j++)
      AA[i][j] = EPS / SAI * sign(1.0, A[i][j]);

    //      DO 2401 J=1,6
    //      AA(I,J)=AA(I,J)+RM(J,2)
    for (int j = 1; j <= 6; j++)
      AA[i][j] = AA[i][j] + RM[j][2];
  }
  // 2401 CONTINUE

  //      COV=DBLE(0.)
  //      DO 4482 I=1,6
  //      DO 4482 J=1,6
  //      DO 2405 K=1,N
  //      DHELP1=DBLE(AA(K,I)-RM(I,2))*1.E-12
  //      DHELP2=DBLE(AA(K,J)-RM(J,2))*1.E-12
  // 2405 COV=COV+DHELP1*DHELP2
  // 4482 CONTINUE

  COV = 0.0;
  for (int i = 1; i <= 6; i++) {
    for (int j = 1; j <= 6; j++) {
      for (int k = 1; k <= N; k++) {
        double DHELP1 = (AA[k][i] - RM[i][2]) * 1.0E-12;
        double DHELP2 = (AA[k][j] - RM[j][2]) * 1.0E-12;
        COV = COV + DHELP1 * DHELP2;
      }
    }
  }

  //      COV=COV/36.
  //      DHELP1=DABS(COV)/DBLE(FLOAT(N-1))
  //      SIG=SNGL(DSQRT(DHELP1))*1.E+12
  COV = COV / 36.0;
  DHELP1 = fabs(COV) / double(N - 1);
  SIG = sqrt(DHELP1) * 1.0e+12;
  RMERR[2] = SIG;

  //C     Output:
  //      WRITE(75,4998)
  // 4998 FORMAT(1H )
  //      WRITE(RESULT(13),5016)
  // 5016 FORMAT('  Trace null solution :',33X)
  //      WRITE(75,5000) RESULT(13)
  //      WRITE(RESULT(14),5013) RM(1,2),RM(2,2),RM(3,2)
  //      WRITE(75,5000) RESULT(14)
  //      WRITE(RESULT(15),5013) RM(2,2),RM(4,2),RM(5,2)
  //      WRITE(75,5000) RESULT(15)
  //      WRITE(RESULT(16),5013) RM(3,2),RM(5,2),RM(6,2)
  //      WRITE(75,5000) RESULT(16)
  //      WRITE(RESULT(17),5014) TROZ,RM0(2),RMT(2),SIG
  //      WRITE(75,5000) RESULT(17)
#ifdef USMTCORE_DEBUG
  std::cout << RM[1][2] << RM[2][2] << RM[3][2] << std::endl;
  std::cout << RM[2][2] << RM[4][2] << RM[5][2] << std::endl;
  std::cout << RM[3][2] << RM[5][2] << RM[6][2] << std::endl;
  std::cout << "T0 = " << TROZ << " M0 = " << RM0[2] << " MT = " << RMT[2] << " ERR = " << SIG;
#endif

  //      CALL EIGGEN(EQQ1,EQQ2,EQQ3,DUM,PCLVD(2),PDBCP(2))
  //      RMAG=.6667*ALOG10(RM0(2))-6.0
  //      IF(RMAG.GE.100.) RMAG=99.9
  //      IF(RMAG.LE.-10.) RMAG=-9.9

  double DUM = 0.0;
  EIGGEN(EQQ1, EQQ2, EQQ3, DUM, PCLVD[2], PDBCP[2]);
  RMAG = 0.6667 * alog10(RM0[2]) - 6.0;
  if (RMAG >= 100.0) RMAG = 99.9;
  if (RMAG <= -10.0) RMAG = -9.9;
  MAGN[2] = RMAG;

  //      WRITE(RESULT(18),5017) PCLVD(2),PDBCP(2),RMAG
  // 5017 FORMAT(18X,'CLVD.=',F6.1,'%   DBCP.=',F6.1,'%   M=',F4.1)
  //      WRITE(75,5000) RESULT(18)
#ifdef USMTCORE_DEBUG
  std::cout << " PCLVD = " << PCLVD[2];
  std::cout << " PDBCP = " << PDBCP[2] << " Mw = " << RMAG << std::endl;
#endif
  PEXPL[2] = 0.0;

  //C     ### DOUBLE COUPLE
  //      CALL GSOLA(H,IEXP)
  //      DO 3055 I=1,5
  // 3055 RM(I,3)=H(I)
  //      RM(6,3)=-RM(1,3)-RM(4,3)
  GSOLA(H, IEXP);
  for (int i = 1; i <= 5; i++)
    RM[i][3] = H[i];
  RM[6][3] = -RM[1][3] - RM[4][3];

  //C     Finds scalar seismic moment:
  //      DO 3009 I=1,6
  // 3009 BB(I)=RM(I,3)
  for (int i = 1; i <= 6; i++)
    BB[i] = RM[i][3];

  //      CALL EIG3(BB,0,EQM)
  //      EQQ1=EQM(1)
  //      EQQ2=EQM(2)
  //      EQQ3=EQM(3)
  //      EQM(1)=ABS(EQM(1))*1.E-12
  //      EQM(2)=ABS(EQM(2))*1.E-12
  //      EQM(3)=ABS(EQM(3))*1.E-12
  //      HELP1=ABS(EQM(1)-EQM(2))
  //      HELP2=ABS(EQM(1)-EQM(3))
  //      HELP3=ABS(EQM(2)-EQM(3))
  //      EC=AMAX1(HELP1,HELP2)
  //      EC=AMAX1(EC,HELP3)
  //      RMT(3)=EC*1.E+12
  //      RM0(3)=RMT(3)

  EIG3(BB, 0, EQM);
  EQQ1 = EQM[1];
  EQQ2 = EQM[2];
  EQQ3 = EQM[3];
  EQM[1] = fabs(EQM[1]) * 1.0E-12;
  EQM[2] = fabs(EQM[2]) * 1.0E-12;
  EQM[3] = fabs(EQM[3]) * 1.0E-12;
  HELP1 = fabs(EQM[1] - EQM[2]);
  HELP2 = fabs(EQM[1] - EQM[3]);
  HELP3 = fabs(EQM[2] - EQM[3]);
  EC = amax1(HELP1, HELP2);
  EC = amax1(EC, HELP3);
  RMT[3] = EC * 1.0E+12;
  RM0[3] = RMT[3];

  //C     Error calculation:
  //      DO 3401 I=1,N
  for (int i = 1; i <= N; i++) {
    //      EPS=U(I)
    double EPS = U[i];
    //      DO 3406 J=1,6
    // 3406 EPS=EPS-A(I,J)*RM(J,3)
    for (int j = 1; j <= 6; j++)
      EPS = EPS - A[i][j] * RM[j][3];

    //      SAI=0.
    //      DO 3403 J=1,6
    // 3403 SAI=SAI+ABS(A(I,J))
    double SAI = 0.0;
    for (int j = 1; j <= 6; j++)
      SAI = SAI + fabs(A[i][j]);

    //      DO 3404 J=1,6
    // 3404 AA(I,J)=EPS/SAI*SIGN(1.,A(I,J))
    for (int j = 1; j <= 6; j++)
      AA[i][j] = EPS / SAI * sign(1.0, A[i][j]);

    //      DO 3401 J=1,6
    //      AA(I,J)=AA(I,J)+RM(J,3)
    for (int j = 1; j <= 6; j++)
      AA[i][j] = AA[i][j] + RM[j][3];
  }
  // 3401 CONTINUE

  //      COV=DBLE(0.)
  //      DO 4492 I=1,6
  //      DO 4492 J=1,6
  //      DO 3405 K=1,N
  //      DHELP1=DBLE(AA(K,I)-RM(I,3))*1.E-12
  //      DHELP2=DBLE(AA(K,J)-RM(J,3))*1.E-12
  // 3405 COV=COV+DHELP1*DHELP2
  // 4492 CONTINUE
  COV = 0.0;
  for (int i = 1; i <= 6; i++) {
    for (int j = 1; j <= 6; j++) {
      for (int k = 1; k <= N; k++) {
        double DHELP1 = (AA[k][i] - RM[i][3]) * 1.0E-12;
        double DHELP2 = (AA[k][j] - RM[j][3]) * 1.0E-12;
        COV = COV + DHELP1 * DHELP2;
      }
    }
  }

  //      COV=COV/36.
  //      DHELP1=DABS(COV)/DBLE(FLOAT(N-1))
  //      SIG=SNGL(DSQRT(DHELP1))*1.E+12
  COV = COV / 36.0;
  DHELP1 = fabs(COV) / double(N - 1);
  SIG = sqrt(DHELP1) * 1.0e+12;
  RMERR[3] = SIG;

  //C     Output:
  //      WRITE(75,4998)
  //      WRITE(RESULT(24),5018)
  // 5018 FORMAT('  Double couple solution :',30X)
  //      WRITE(75,5000) RESULT(24)
  //      WRITE(RESULT(25),5013) RM(1,3),RM(2,3),RM(3,3)
  //      WRITE(75,5000) RESULT(25)
  //      WRITE(RESULT(26),5013) RM(2,3),RM(4,3),RM(5,3)
  //      WRITE(75,5000) RESULT(26)
  //      WRITE(RESULT(27),5013) RM(3,3),RM(5,3),RM(6,3)
  //      WRITE(75,5000) RESULT(27)
  //      WRITE(RESULT(28),5014) TROZ,RM0(3),RMT(3),SIG
  //      WRITE(75,5000) RESULT(28)
#ifdef USMTCORE_DEBUG
  std::cout << RM[1][3] << RM[2][3] << RM[3][3] << std::endl;
  std::cout << RM[2][3] << RM[4][3] << RM[5][3] << std::endl;
  std::cout << RM[3][3] << RM[5][3] << RM[6][3] << std::endl;
  std::cout << "T0 = " << TROZ << " M0 = " << RM0[3] << " MT = " << RMT[3] << " ERR = " << SIG;
#endif

  //      RMAG=.6667*ALOG10(RM0(3))-6.0
  //      IF(RMAG.GE.100.) RMAG=99.9
  //      IF(RMAG.LE.-10.) RMAG=-9.9
  //      WRITE(RESULT(29),5019) RMAG
  RMAG = 0.6667 * alog10(RM0[3]) - 6.0;
  if (RMAG >= 100.0) RMAG = 99.9;
  if (RMAG <= -10.0) RMAG = -9.9;
  MAGN[3] = RMAG;

  // 5019 FORMAT(34X,'DBCP.= 100.0%   M=',F4.1)
#ifdef USMTCORE_DEBUG
  std::cout << " PDBCP = 100% Mw = " << RMAG << std::endl;
#endif
  PEXPL[3] = 0.0;
  PCLVD[3] = 0.0;
  PDBCP[3] = 100.0;

  //      WRITE(75,5000) RESULT(29)
  //      WRITE(75,4997)
  // 4997 FORMAT(72(1H-))
  //      ICOND=1
  //      CALL XTRINF(ICOND)
  //      NORM=1
  //      CALL WROUT(NORM)

  //---- Transfer data to Taquart::FaultSolution structures

  for (int i = 1; i <= 3; i++) {
    int z = 0;
    int v[] = { 1, 2, 3, 2, 4, 5, 3, 5, 6 };
    for (int m = 1; m <= 3; m++)
      for (int n = 1; n <= 3; n++) {
        Solution[i].M[m][n] = RM[v[z]][i];
        z++;
      }
    Solution[i].T0 = TROZ;
    Solution[i].M0 = RM0[i];
    Solution[i].MT = RMT[i];
    Solution[i].ERR = RMERR[i];
    Solution[i].EXPL = PEXPL[i];
    Solution[i].CLVD = PCLVD[i];
    Solution[i].DBCP = PDBCP[i];
    Solution[i].MAGN = MAGN[i];
    for (int m = 1; m <= 6; m++)
      for (int n = 1; n <= 6; n++)
        Solution[i].Covariance[m][n] = 0.0;
  }

  ICOND = 1;

  XTRINF(ICOND, 1, RM0, RMERR);

  //      RETURN
  //      END
}

//-----------------------------------------------------------------------------
void Taquart::UsmtCore::EIG3(double RM[], int ISTER, double E[]) {
  //      SUBROUTINE EIG3(RM,ISTER,E)
  //C     Routine finds eigenvalues of matrix 3*3 by solving characteristic
  //C     equation.
  //      DIMENSION RM(6),E(3)
  //      DOUBLE PRECISION A2,A1,A0,VAL,ONE,TWO,ZERO,A,B,X,Q(6)
  //      ZERO=DBLE(0.)
  //      ONE=DBLE(1.)
  //      TWO=DBLE(2.)
  double ZERO = 0.0, TWO = 2.0;
  double Q[6 + 1];
  Zero(Q, 6);
  double A2 = 0.0, A1 = 0.0, A0 = 0.0, VAL = 0.0, A = 0.0, B = 0.0, X = 0.0;
  double F = 0.0, G = 0.0, OX = 0.0, FF = 0.0, FG = 0.0, QX = 0.0, HELP = 0.0;
  double FX = 0.0;

  //      DO 1 I=1,6
  //    1 Q(I)=DBLE(RM(I))
  for (int i = 1; i <= 6; i++)
    Q[i] = RM[i];

  //      A2=-(Q(1)+Q(4)+Q(6))
  //      A1=Q(1)*Q(6)+Q(4)*Q(6)+Q(1)*Q(4)-Q(2)*Q(2)-Q(3)*Q(3)-Q(5)*Q(5)
  //      A0=-TWO*Q(2)*Q(3)*Q(5)-Q(1)*Q(4)*Q(6)+Q(3)*Q(3)*Q(4)+Q(2)*Q(2)*Q(6)+Q(5)*Q(5)*Q(1)
  //      A=-1.D+30
  //      B=1.D+30
  //      X=ZERO

  A2 = -(Q[1] + Q[4] + Q[6]);
  A1 = Q[1] * Q[6] + Q[4] * Q[6] + Q[1] * Q[4] - Q[2] * Q[2] - Q[3] * Q[3]
      - Q[5] * Q[5];
  A0 = -TWO * Q[2] * Q[3] * Q[5] - Q[1] * Q[4] * Q[6] + Q[3] * Q[3] * Q[4]
      + Q[2] * Q[2] * Q[6] + Q[5] * Q[5] * Q[1];
  A = -1.0e+30;
  B = 1.0e+30;
  X = ZERO;

  //      DO 2 KROK=1,200
  for (int KROK = 1; KROK <= 200; KROK++) {
    //      VAL=X**3+A2*X**2+A1*X+A0
    VAL = pow(X, 3.0) + A2 * pow(X, 2) + A1 * X + A0;
    //      IF(VAL.EQ.ZERO) GO TO 3
    //      IF(VAL.LT.ZERO) A=X
    //      IF(VAL.GT.ZERO) B=X
    //      X=(A+B)/TWO
    if (VAL == ZERO) goto p3;
    if (VAL < ZERO) A = X;
    if (VAL > ZERO) B = X;
    X = (A + B) / TWO;
  }
  //    2 CONTINUE

  //    3 A0=A1+A2*X+X*X
  //      A1=A2+X
  //      VAL=A1*A1-TWO*TWO*A0

  p3: A0 = A1 + A2 * X + X * X;
  A1 = A2 + X;
  VAL = A1 * A1 - TWO * TWO * A0;

  //      IF(VAL.GE.ZERO) GO TO 5
  if (VAL < ZERO) VAL = ZERO;

  //    5 VAL=DSQRT(VAL)
  //      A=(-A1-VAL)/TWO
  //      B=(-A1+VAL)/TWO
  //      E(1)=SNGL(X)
  //      E(2)=SNGL(A)
  //      E(3)=SNGL(B)
  VAL = sqrt(VAL);
  A = (-A1 - VAL) / TWO;
  B = (-A1 + VAL) / TWO;
  E[1] = X;
  E[2] = A;
  E[3] = B;

  //      IF(ISTER.EQ.0) RETURN
  if (ISTER == 0) return;

  //      DO 6 I=1,3
  for (int i = 1; i <= 3; i++) {
    //      iter=0
    int iter = 0;

    //      IF(E(I).EQ.0.) THEN
    //       F=-1.E-06
    //       G=1.E-06
    //      ELSE
    //       F=E(I)*.999
    //       G=E(I)*1.001
    //      ENDIF
    F = 0.0;
    G = 0.0;
    if (E[i] == 0.0) {
      F = -1.0e-06;
      G = +1.0e-06;
    }
    else {
      F = E[i] * 0.999;
      G = E[i] * 1.001;
    }

    //      OX=F
    OX = F;

    p7:
    //    7 FF=DETR(RM,F)
    //      FG=DETR(RM,G)
    //      QX=(F+G)/2.
    //      iter=iter+1
    FF = DETR(RM, F);
    FG = DETR(RM, G);
    QX = (F + G) / 2.0;
    iter = iter + 1;

    //      if(iter.gt.1000) go to 6
    if (iter > 1000) continue;

    //      IF(ABS((QX-OX)/E(I)).LT.1.E-06) GO TO 8

    /* DONE 5 -c2.4.8 : Correction for the /div0 error. */
    double ee = (E[i] == 0.0) ? 1.0e-6 : E[i];
    //--------------------------------------
    if (fabs((QX - OX) / ee) < 1.0e-6) goto p8;

    //      HELP=FF*FG
    HELP = FF * FG;

    //      IF(HELP.GT.0.) GO TO 9
    if (HELP > 0.0) goto p9;

    //      FX=DETR(RM,QX)
    //      HELP=FF*FX
    FX = DETR(RM, QX);
    HELP = FF * FX;

    //      IF(HELP.LT.0.) THEN
    //       G=QX
    //      ELSE
    //       F=QX
    //      ENDIF
    if (HELP < 0.0) {
      G = QX;
    }
    else {
      F = QX;
    }

    //      OX=QX
    //      GO TO 7
    OX = QX;
    goto p7;

    //    9 F=.999*F
    //      G=1.001*G
    //      GO TO 7
    p9: F = 0.999 * F;
    G = 1.001 * G;
    goto p7;

    //    8 E(I)=QX
    p8: E[i] = QX;
  }
  //    6 CONTINUE
  //      RETURN
  //      END
}

/*
 disp('skladowe tensora liczone wedlug Jost and Herrmann, 1989');
 W=eig(M);
 Wa=eig(M)-trace(M)/3;
 V1=(Wa(find(abs(Wa)==max(abs(Wa)))));
 V3=Wa(find(abs(Wa)==min(abs(Wa))));
 FF=-V3/V1;
 iso1=sum(eig(M))/3;
 clvd1=V1*FF;
 dc1=V1*(1-2*FF);
 dc1=abs(dc1);
 ISO=iso1/(abs(iso1)+abs(clvd1)+abs(dc1))
 CLVD=clvd1/(abs(iso1)+abs(clvd1)+abs(dc1))
 DC=dc1/(abs(iso1)+abs(clvd1)+abs(dc1))
 void Taquart::UsmtCore::EIGGEN_NEW(double &e1, double &e2, double &e3, double &ALFA, double &BETA,
 {
 }
 */
//-----------------------------------------------------------------------------
void Taquart::UsmtCore::EIGGEN_NEW(double e1, double e2, double e3, double &iso,
    double &clvd, double &dbcp) {
  const double trace = (e1 + e2 + e3) / 3.0;
  const double e1a = fabs(e1 - trace);
  const double e2a = fabs(e2 - trace);
  const double e3a = fabs(e3 - trace);
  double v1 = 0.0, v3 = 0.0;

  // Find the largest absolute deviatoric eigenvalue.
  if (e1a >= e2a && e1a >= e3a) {
    v1 = e1 - trace;
  }
  if (e2a >= e1a && e2a >= e3a) {
    v1 = e2 - trace;
  }
  if (e3a >= e1a && e3a >= e2a) {
    v1 = e3 - trace;
  }

  // Fint the smallest...
  if (e1a <= e2a && e1a <= e3a) {
    v3 = e1 - trace;
  }
  if (e2a <= e1a && e2a <= e3a) {
    v3 = e2 - trace;
  }
  if (e3a <= e1a && e3a <= e2a) {
    v3 = e3 - trace;
  }

  const double ff = -v3 / v1;
  iso = trace;
  clvd = 2.0 * v1 * ff; /* TODO 5 -c3.4.2 : Correction of the error, now it is according to Jost & Hermann */
  dbcp = fabs(v1 * (1.0 - 2.0 * ff));
  const double s = fabs(iso) + fabs(clvd) + fabs(dbcp);

  iso = iso / s * 100.0;
  clvd = clvd / s * 100.0;
  dbcp = dbcp / s * 100.0;
}

//-----------------------------------------------------------------------------
void Taquart::UsmtCore::EIGGEN(double &E1, double &E2, double &E3, double &ALFA,
    double &BETA, double &GAMA) {
  EIGGEN_NEW(E1, E2, E3, ALFA, BETA, GAMA);
  /*
   int IA = 1;
   double EA = E1;
   double EB = 0.0;
   int IB = 0;
   double EC = 0.0;
   int IC = 0;

   if(EA >= E2)
   {
   EA = E2;
   IA = 2;
   }
   if(EA >= E3)
   {
   EA = E3;
   IA = 3;
   }
   if(IA == 1)
   {
   EB = E2;
   IB = 2;
   if(EB > E3) goto p3;
   EB = E3;
   IB = 3;
   goto p3;
   }

   if(IA != 3)
   {
   EB = E1;
   IB = 1;
   if(EB > E3) goto p3;
   EB = E3;
   IB = 3;
   goto p3;
   }
   else
   {
   EB = E1;
   IB = 1;
   if(EB > E2) goto p3;
   EB = E2;
   IB = 2;
   }

   p3: IC = 6 - IA - IB;
   if(IC == 1) EC = E1;
   if(IC == 2) EC = E2;
   if(IC == 3) EC = E3;

   ALFA = (E1 + E2 + E3) / 3.0;
   BETA = ALFA - EC;
   GAMA = EC - EA;
   ALFA = ALFA * 1.0e-12;
   BETA = BETA * 1.0e-12;
   GAMA = GAMA * 1.0e-12;

   const double skal = fabs(ALFA) + fabs(BETA) + fabs(GAMA);
   ALFA = ALFA / skal * 100.0;
   BETA = BETA / skal * 100.0;
   GAMA = GAMA / skal * 100.0;
   */
}

//-----------------------------------------------------------------------------
void Taquart::UsmtCore::XTRINF(int &ICOND, int LNORM, double Moment0[],
    double MomentErr[]) {
  //      SUBROUTINE XTRINF(ICOND)
  //      CHARACTER PS(80),TITLE*40
  //      REAL U(80),ARR(80),AZM(80),TKF(80),ETA(3),RKAPPA(3),VALKAP(3),
  //     $B(6),EQM(3),V(3,3,3),RNU(3),PLUNGE(3,3),TREND(3,3),X(3),Y(3),
  //     $STRIKE(3,2),DIP(3,2)
  //      INTEGER RO(80),VEL(80),R(80)
  //      INTEGER*2 QF
  //      COMMON/MDATA/ PS,U,ARR,AZM,TKF,RO,VEL,R,TITLE,N,TROZ
  //      COMMON/MOMNT/ RM(6,3)
  //      COMMON/JEZQUA/ QSD,QF
  //      COMMON/RESTEX/ RESULT
  //      CHARACTER*56 RESULT(34)
  //      DOUBLE PRECISION COV(6,6,3)
  //      COMMON/VCOVAR/ COV
  //      CHARACTER*13 INFO(3),INFFTD(3)
  //      DATA VALKAP/1.,.707,.5/
  //      DATA INFO/'normal  fault',' strike slip ',
  //     $'reverse fault'/

  //      SQR2=SQRT(2.)
  //      PI=4.*ATAN(1.)
  double SQR2 = sqrt(2.0);
  double PI = 4.0 * atan(1.0);
  double B[6];
  Zero(B, 6);
  double V[4][4][4];
  Zero(&V[0][0][0], 64);
  double EQM[4];
  Zero(EQM, 4);
  double X[4], Y[4];
  Zero(X, 4);
  Zero(Y, 4);
  double STRIKE[4][3];
  Zero(&STRIKE[0][0], 12);
  double DIP[4][3];
  Zero(&DIP[0][0], 12);
  double RAKE[4][3];
  Zero(&RAKE[0][0], 12);
  double PLUNGE[4][4];
  Zero(&PLUNGE[0][0], 16);
  double TREND[4][4];
  Zero(&TREND[0][0], 16);
  double AMP[4][4];
  Zero(&AMP[0][0], 16);

  String INFO[4] = { "", "Normal fault", "Strike fault", "Reverse fault" };
  String INFFTD[4];

  double ETA[4];
  double HELP = 0, HELP1 = 0.0, HELP2 = 0.0, HELP3 = 0.0;
  int IMIN = 0;
  int IMAX = 0;
  int j = 0;

  //C     Quality factor calculation
  //      READ(RESULT(1),1) LNORM
  //    1 FORMAT(1X,I1)

  //      IF(LNORM.EQ.2) GO TO 2
  if (LNORM != 2) {
    //      DO 4 I=1,3
    //      READ(RESULT(11*I-5),3) HELP1,HELP2
    //    3 FORMAT(18X,E9.3,20X,E9.3)
    //    4 ETA(I)=1.-HELP2/HELP1
    //      GO TO 5
    for (int i = 1; i <= 3; i++)
      ETA[i] = 1.0 - MomentErr[i] / Moment0[i];
  }
  else {
    //    2 DO 6 I=1,3
    for (int i = 1; i <= 3; i++) {
      //      READ(RESULT(11*I-5),7) HELP1
      //    7 FORMAT(18X,E9.3)
      //      HELP2=SNGL(DSQRT(COV(1,1,I)))
      HELP1 = Moment0[i];
      HELP3 = 0.0;
      HELP2 = sqrt(COV[1][1][i]);
      //      DO 8 J=2,6
      for (int j = 2; j <= 6; j++) {
        //      HELP3=SNGL(DSQRT(COV(J,J,I)))
        //    8 HELP2=AMAX1(HELP2,HELP3)
        HELP3 = sqrt(COV[j][j][i]);
        HELP2 = amax1(HELP2, HELP3);
      }
      ETA[i] = 1.0 - HELP2 / HELP1;
    }
    //    6 ETA(I)=1.-HELP2/HELP1
  }
  //    5 RNU(1)=FLOAT(N-7)/FLOAT(N)
  //      RNU(2)=FLOAT(N-6)/FLOAT(N)
  //      RNU(3)=FLOAT(N-5)/FLOAT(N)
  double RNU[4];
  RNU[1] = double(N - 7) / double(N);
  RNU[2] = double(N - 6) / double(N);
  RNU[3] = double(N - 5) / double(N);

  //      IF((ICOND.LE.0).OR.(ICOND.GT.3)) ICOND=1
  //      RKAPPA(1)=1.
  //      RKAPPA(2)=1.
  //      RKAPPA(3)=VALKAP(ICOND)
  if (ICOND <= 0 || ICOND > 3) ICOND = 1;
  const double VALKAP[4] = { 0.0, 1.0, 0.707, 0.5 };
  double RKAPPA[4] = { 0.0, 1.0, 1.0, VALKAP[ICOND] };

  //      DO 10 I=1,3
  for (int i = 1; i <= 3; i++) {
    //      IF(ETA(I).LT.0.) ETA(I)=.01
    //      ETA(I)=ETA(I)*RNU(I)*RKAPPA(I)*100.
    //      IF(QF.NE.0) ETA(I)=ETA(I)*QSD
    if (ETA[i] < 0.0) ETA[i] = 0.01;
    ETA[i] = ETA[i] * RNU[i] * RKAPPA[i] * 100.0;
    if (QF != 0.0) ETA[i] = ETA[i] * QSD;
  }
  //   10 CONTINUE

  //C     Calculate eigenvectors for three solutions
  //      DO 11 I=1,3
  //      DO 12 J=1,6
  //   12 B(J)=RM(J,I)
  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 6; j++)
      B[j] = RM[j][i];

    //      CALL EIG3(B,0,EQM)
    EIG3(B, 0, EQM);

    //      IMIN=1
    //      IMAX=3
    IMIN = 1;
    IMAX = 3;

    //      IF(EQM(1).GT.EQM(2)) IMIN=2
    //      IF(EQM(IMIN).GT.EQM(3)) IMIN=3
    //      IF(EQM(3).LT.EQM(2)) IMAX=2
    //      IF(EQM(IMAX).LT.EQM(1)) IMAX=1
    if (EQM[1] > EQM[2]) IMIN = 2;
    if (EQM[IMIN] > EQM[3]) IMIN = 3;
    if (EQM[3] < EQM[2]) IMAX = 2;
    if (EQM[IMAX] < EQM[1]) IMAX = 1;

    //      INUL=6-IMIN-IMAX
    //      HELP1=EQM(IMIN)
    //      HELP2=EQM(INUL)
    //      HELP3=EQM(IMAX)
    //      EQM(1)=HELP1
    //      EQM(2)=HELP2
    //      EQM(3)=HELP3
    int INUL = 6 - IMIN - IMAX;
    HELP1 = EQM[IMIN];
    HELP2 = EQM[INUL];
    HELP3 = EQM[IMAX];
    EQM[1] = HELP1;
    EQM[2] = HELP2;
    EQM[3] = HELP3;

    //      DO 13 J=1,3
    //      B(1)=B(1)-EQM(J)
    //      B(4)=B(4)-EQM(J)
    //      B(6)=B(6)-EQM(J)
    //      HELP1=B(1)*B(4)-B(2)*B(2)
    for (int j = 1; j <= 3; j++) {
      B[1] = B[1] - EQM[j];
      B[4] = B[4] - EQM[j];
      B[6] = B[6] - EQM[j];
      HELP1 = B[1] * B[4] - B[2] * B[2];

      //      IF(ABS(HELP1).LT..001) GO TO 14
      if (fabs(HELP1) >= 0.001) {
        //      V(3,J,I)=1.
        //      V(1,J,I)=(-B(3)*B(4)+B(2)*B(5))/HELP1
        //      V(2,J,I)=(-B(1)*B(5)+B(2)*B(3))/HELP1
        //      GO TO 15
        V[3][j][i] = 1.0;
        V[1][j][i] = (-B[3] * B[4] + B[2] * B[5]) / HELP1;
        V[2][j][i] = (-B[1] * B[5] + B[2] * B[3]) / HELP1;
      }
      else {
        //   14 HELP1=B(1)*B(5)-B(2)*B(3)
        HELP1 = B[1] * B[5] - B[2] * B[3];
        //      IF(ABS(HELP1).LT..001) GO TO 16
        if (fabs(HELP1) >= 0.001) {
          //      V(2,J,I)=1.
          //      V(1,J,I)=(-B(2)*B(5)+B(4)*B(3))/HELP1
          //      V(3,J,I)=(-B(1)*B(4)+B(2)*B(2))/HELP1
          //      GO TO 15
          V[2][j][i] = 1.0;
          V[1][j][i] = (-B[2] * B[5] + B[4] * B[3]) / HELP1;
          V[3][j][i] = (-B[1] * B[4] + B[2] * B[2]) / HELP1;
        }
        else {
          //   16 HELP1=B(2)*B(5)-B(4)*B(3)
          HELP1 = B[2] * B[5] - B[4] * B[3];
          //      IF(ABS(HELP1).LT..001) GO TO 17
          if (fabs(HELP1) >= 0.001) {
            //      V(1,J,I)=1.
            //      V(2,J,I)=(-B(1)*B(5)+B(2)*B(3))/HELP1
            //      V(3,J,I)=(-B(2)*B(2)+B(4)*B(1))/HELP1
            V[1][j][i] = 1.0;
            V[2][j][i] = (-B[1] * B[5] + B[2] * B[3]) / HELP1;
            V[3][j][i] = (-B[2] * B[2] + B[4] * B[1]) / HELP1;
            //      GO TO 15
          }
          else {
            //   17 V(3,J,I)=1.
            //      V(1,J,I)=0.
            //      V(2,J,I)=0.
            V[3][j][i] = 1.0;
            V[1][j][i] = 0.0;
            V[2][j][i] = 0.0;
          }
        }
      }

      //   15 HELP1=SQRT(V(1,J,I)*V(1,J,I)+V(2,J,I)*V(2,J,I)+V(3,J,I)*V(3,J,I))
      //      IF(V(3,J,I).LT.0.) HELP1=-1.*HELP1
      HELP1 = sqrt(
          V[1][j][i] * V[1][j][i] + V[2][j][i] * V[2][j][i]
              + V[3][j][i] * V[3][j][i]);
      if (V[3][j][i] < 0.0) HELP1 = -1.0 * HELP1;

      //      DO 18 K=1,3
      //   18 V(K,J,I)=V(K,J,I)/HELP1
      for (int k = 1; k <= 3; k++)
        V[k][j][i] = V[k][j][i] / HELP1;

      //      B(1)=B(1)+EQM(J)
      //      B(4)=B(4)+EQM(J)
      //   13 B(6)=B(6)+EQM(J)
      //   11 CONTINUE
      B[1] = B[1] + EQM[j];
      B[4] = B[4] + EQM[j];
      B[6] = B[6] + EQM[j];
    }
  }

  //C     Now we have EIGENVECTORS V(K,J,I), J pertaining to P,B,T,
  //C     I pertaining to F,T,D, K pertaining to 1,2,3
  //      DO 19 I=1,3
  //      DO 19 J=1,3
  //      HELP=AMIN1(V(3,J,I),1.)
  //      PLUNGE(J,I)=ASIN(AMAX1(HELP,-1.))
  //      TREND(J,I)=ATAN2(V(2,J,I),V(1,J,I))
  //   19 CONTINUE
  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      HELP = amin1(V[3][j][i], 1.0);
      PLUNGE[j][i] = asin(amax1(HELP, -1.0));
      //TREND[j][i] = atan2(V[2][j][i],V[1][j][i]);

      /* TODO 5 -c3.1.19 : Correction of ATAN2 error. */
      TREND[j][i] = datan2(V[2][j][i], V[1][j][i]) * M_PI / 180.0;

      AMP[j][i] = sqrt(
          V[1][j][i] * V[1][j][i] + V[2][j][i] * V[2][j][i]
              + V[3][j][i] * V[3][j][i]);
    }
  }

  //C     Plunge & Trend, first ccordinate is PBT, second is FTD
  //      DO 20 I=1,3
  HELP1 = 0.0;
  j = 0;
  double FaultType[4] = { 0.0, 0.0, 0.0, 0.0 };
  for (int i = 1; i <= 3; i++) {
    //      J=1
    //      HELP1=PLUNGE(1,I)
    j = 1;
    HELP1 = PLUNGE[1][i];
    FaultType[i] = -1.0;
    //      IF(PLUNGE(2,I).LE.HELP1) GO TO 21
    //      J=2
    //      HELP1=PLUNGE(2,I)
    //   21 IF(PLUNGE(3,I).GT.HELP1) J=3
    //   20 INFFTD(I)=INFO(J)
    if (PLUNGE[2][i] > HELP1) {
      j = 2;
      HELP1 = PLUNGE[2][i];
    }
    if (PLUNGE[3][i] > HELP1) j = 3;
    INFFTD[i] = INFO[j];

    if (PLUNGE[3][i] > PLUNGE[1][i]) {
      FaultType[i] = 1.0;
    }
  }

  //C     INFFTD contains information on fault type.
  //      DO 22 I=1,3
  for (int i = 1; i <= 3; i++) {
    //      DO 23 K=1,3
    for (int k = 1; k <= 3; k++) {
      //      X(K)=(V(K,3,I)-V(K,1,I))/SQR2
      X[k] = (V[k][3][i] - V[k][1][i]) / SQR2;
      Y[k] = (V[k][3][i] + V[k][1][i]) / SQR2;
    }
    //   23 Y(K)=(V(K,3,I)+V(K,1,I))/SQR2

    //      IF(X(3).GE.0.) GO TO 24
    if (X[3] < 0.0) {
      //      DO 25 K=1,3
      //   25 X(K)=-X(K)
      for (int k = 1; k <= 3; k++)
        X[k] = -X[k];
    }
    //   24 IF(Y(3).GE.0.) GO TO 26
    else if (Y[3] < 0.0) {
      //      DO 27 K=1,3
      //   27 Y(K)=-Y(K)
      for (int k = 1; k <= 3; k++)
        Y[k] = -Y[k];
    }
    //   26 STRIKE(I,1)=0.
    //      STRIKE(I,2)=0.
    STRIKE[i][1] = 0.0;
    STRIKE[i][2] = 0.0;
    //      X(3)=AMIN1(X(3),1.)
    //      Y(3)=AMIN1(Y(3),1.)
    X[3] = amin1(X[3], 1.0);
    Y[3] = amin1(Y[3], 1.0);
    //      DIP(I,1)=ACOS(X(3))
    //      DIP(I,2)=ACOS(Y(3))
    DIP[i][1] = acos(X[3]);
    DIP[i][2] = acos(Y[3]);
    //      IF(DIP(I,1).GT..01) STRIKE(I,1)=ATAN2(X(1),-X(2))
    //      IF(DIP(I,2).GT..01) STRIKE(I,2)=ATAN2(Y(1),-Y(2))
    //      IF(DIP(I,1).GE.DIP(I,2)) GO TO 22

    /* TODO 5 -c3.1.19 : Correction of ATAN2 error */
    /*
     if(DIP[i][1] > 0.01) STRIKE[i][1] = atan2(X[1],-X[2]);
     if(DIP[i][2] > 0.01) STRIKE[i][2] = atan2(Y[1],-Y[2]);
     */
    if (DIP[i][1] > 0.01) STRIKE[i][1] = datan2(X[1], -X[2]) * M_PI / 180.0;
    if (DIP[i][2] > 0.01) STRIKE[i][2] = datan2(Y[1], -Y[2]) * M_PI / 180.0;

    if (DIP[i][1] >= DIP[i][2]) continue;

    //      HELP1=DIP(I,1)
    //      DIP(I,1)=DIP(I,2)
    //      DIP(I,2)=HELP1
    //      HELP1=STRIKE(I,1)
    //      STRIKE(I,1)=STRIKE(I,2)
    //      STRIKE(I,2)=HELP1
    HELP1 = DIP[i][1];
    DIP[i][1] = DIP[i][2];
    DIP[i][2] = HELP1;
    HELP1 = STRIKE[i][1];
    STRIKE[i][1] = STRIKE[i][2];
    STRIKE[i][2] = HELP1;
  }
  //   22 CONTINUE

  //---- RAKE angles calculation.
  /*
   for(int i=1; i<=3; i++)
   {
   for(int k=1; k<=2; k++)
   {
   double Value = 0.0;
   double Divisor = 0.0;

   Divisor = sin(2.0 * DIP[i][k]) * Moment0[i];
   if(Divisor == 0.0 && RM[6][i] > 0.0)
   Value = 1.0;
   else if(Divisor == 0.0 && RM[6][i] < 0.0)
   Value = -1.0;
   else Value = RM[6][i] / Divisor;
   if(Value < -1.0) Value = -1.0;
   if(Value > 1.0) Value = 1.0;
   RAKE[i][k] = asin(Value);
   }
   }
   !!! SEE LATER >>>>>
   */

  //      DO 28 I=1,3
  for (int i = 1; i <= 3; i++) {
    //      DO 29 K=1,2
    for (int k = 1; k <= 2; k++) {
      //      STRIKE(I,K)=STRIKE(I,K)*180./PI
      //      IF(STRIKE(I,K).LT.0.) STRIKE(I,K)=360.+STRIKE(I,K)
      //   29 DIP(I,K)=DIP(I,K)*180./PI
      STRIKE[i][k] = STRIKE[i][k] * 180.0 / PI;
      if (STRIKE[i][k] < 0.0) STRIKE[i][k] = 360.0 + STRIKE[i][k];
      DIP[i][k] = DIP[i][k] * 180.0 / PI;
      //RAKE[i][k] = RAKE[i][k]*180.0/PI; // Latter calculations >>>>>

      /* DONE  5 -c2.4.13 : Usuniete w 2.4.13 */
      //if(RAKE[i][k] < -180.0) RAKE[i][k] = RAKE[i][k] + 360.0;
      //if(RAKE[i][k] > 180.0) RAKE[i][k] = RAKE[i][k] - 360.0;
      //----
    }

    /* DONE 5 -c2.4.13 :
     Oblicz warto�ci RAKE przy pomocy metody z programu PSMECA. Korekta
     krytycznego b��du a� do wersji 2.4.13 dotycz�ca b��dnej warto�ci rake
     (nie mia�a wp�ywu na rozwi�zanie graficzne tensora, tylko na warto��
     rake w plikach wynikowych).*/

    /*
     double Strike1 = 78.9; //STRIKE[i][1];
     double Strike2 = 255.6; //STRIKE[i][2];
     double Dip1 = 48.8; //DIP[i][1];
     double Dip2 = 41.2; //DIP[i][2];
     FaultType[i] = 1.0;
     */

    double Strike1 = STRIKE[i][1];
    double Strike2 = STRIKE[i][2];
    double Dip1 = DIP[i][1];
    double Dip2 = DIP[i][2];
    double Rake1 = Taquart::computed_rake2(Strike2, Dip2, Strike1, Dip1,
        FaultType[i]);
    double Rake2 = Taquart::computed_rake2(Strike1, Dip1, Strike2, Dip2,
        FaultType[i]);
    RAKE[i][1] = Rake1;
    RAKE[i][2] = Rake2;
    //------

    //      DO 28 K=1,3
    for (int k = 1; k <= 3; k++) {
      //      TREND(I,K)=TREND(I,K)*180./PI
      //      IF(TREND(I,K).LT.0.) TREND(I,K)=360.+TREND(I,K)
      //   28 PLUNGE(I,K)=PLUNGE(I,K)*180./PI
      TREND[i][k] = TREND[i][k] * 180.0 / PI;
      if (TREND[i][k] < 0.0) TREND[i][k] = 360.0 + TREND[i][k];
      PLUNGE[i][k] = PLUNGE[i][k] * 180.0 / PI;
    }
  }

  //      WRITE(RESULT(8),31) STRIKE(1,1),DIP(1,1),
  //     $STRIKE(1,2),DIP(1,2)
  //   31 FORMAT(12X,3HfA=,F7.2,2X,3HdA=,F6.2,2X,3HfB=,F7.2,2X,3HdB=,F6.2)

  //      WRITE(RESULT(9),32) TREND(1,1),PLUNGE(1,1)
  //   32 FORMAT(12X,14HP-axis: trend=,F7.2,2X,7Hplunge=,F6.2,8X)
  //      WRITE(RESULT(10),33) TREND(3,1),PLUNGE(3,1)
  //   33 FORMAT(12X,14HT-axis: trend=,F7.2,2X,7Hplunge=,F6.2,8X)
  //      WRITE(RESULT(11),34) TREND(2,1),PLUNGE(2,1)
  //   34 FORMAT(12X,14HB-axis: trend=,F7.2,2X,7Hplunge=,F6.2,8X)
  //      WRITE(RESULT(12),37) ETA(1),INFFTD(1)
#ifdef USMTCORE_DEBUG
  std::cout << "Full solution:" << std::endl;
  std::cout << "--------------" << std::endl;
  std::cout << "P-axis: trend = " << TREND[1][1] << " plunge = " << PLUNGE[1][1] << std::endl;
  std::cout << "T-axis: trend = " << TREND[3][1] << " plunge = " << PLUNGE[3][1] << std::endl;
  std::cout << "B-axis: trend = " << TREND[2][1] << " plunge = " << PLUNGE[2][1] << "  " << INFFTD[1].c_str() << std::endl;
  std::cout << "           fA = " << STRIKE[1][1] << " dA = " << DIP[1][1];
  std::cout << "           fB = " << STRIKE[1][2] << " dB = " << DIP[1][2] << std::endl;
  std::cout << "Quality factor = " << ETA[1] << std::endl;
  std::cout << std::endl;

  std::cout << "Trace-null solution:" << std::endl;
  std::cout << "--------------------" << std::endl;
  std::cout << "P-axis: trend = " << TREND[1][2] << " plunge = " << PLUNGE[1][2] << std::endl;
  std::cout << "T-axis: trend = " << TREND[3][2] << " plunge = " << PLUNGE[3][2] << std::endl;
  std::cout << "B-axis: trend = " << TREND[2][2] << " plunge = " << PLUNGE[2][2] << "  " << INFFTD[2].c_str() << std::endl;
  std::cout << "           fA = " << STRIKE[2][1] << " dA = " << DIP[2][1];
  std::cout << "           fB = " << STRIKE[2][2] << " dB = " << DIP[2][2] << std::endl;
  std::cout << "Quality factor = " << ETA[2] << std::endl;
  std::cout << std::endl;

  std::cout << "Double-couple solution:" << std::endl;
  std::cout << "-----------------------" << std::endl;
  std::cout << "P-axis: trend = " << TREND[1][3] << " plunge = " << PLUNGE[1][3] << std::endl;
  std::cout << "T-axis: trend = " << TREND[3][3] << " plunge = " << PLUNGE[3][3] << std::endl;
  std::cout << "B-axis: trend = " << TREND[2][3] << " plunge = " << PLUNGE[2][3] << "  " << INFFTD[3].c_str() << std::endl;
  std::cout << "           fA = " << STRIKE[3][1] << " dA = " << DIP[3][1];
  std::cout << "           fB = " << STRIKE[3][2] << " dB = " << DIP[3][2] << std::endl;
  std::cout << "Quality factor = " << ETA[3] << std::endl;
  std::cout << std::endl;
#endif

  for (int i = 1; i <= 3; i++) {
    Solution[i].FIA = STRIKE[i][1];
    Solution[i].DLA = DIP[i][1];
    Solution[i].RAKEA = RAKE[i][1];
    Solution[i].FIB = STRIKE[i][2];
    Solution[i].DLB = DIP[i][2];
    Solution[i].RAKEB = RAKE[i][2];
    Solution[i].PXTR = TREND[1][i];
    Solution[i].PXPL = PLUNGE[1][i];
    Solution[i].PXAM = AMP[1][i];
    Solution[i].TXTR = TREND[3][i];
    Solution[i].TXPL = PLUNGE[3][i];
    Solution[i].TXAM = AMP[3][i];
    Solution[i].BXTR = TREND[2][i];
    Solution[i].BXPL = PLUNGE[2][i];
    Solution[i].BXAM = AMP[2][i];
    Solution[i].QI = ETA[i];
    Solution[i].Type = INFFTD[i];
  }

  //   37 FORMAT(12X,15HQuality index =,F6.1,2X,A13,8X)
  //      WRITE(RESULT(19),31) STRIKE(2,1),DIP(2,1),
  //     $STRIKE(2,2),DIP(2,2)
  //      WRITE(RESULT(20),32) TREND(1,2),PLUNGE(1,2)
  //      WRITE(RESULT(21),33) TREND(3,2),PLUNGE(3,2)
  //      WRITE(RESULT(22),34) TREND(2,2),PLUNGE(2,2)
  //      WRITE(RESULT(23),37) ETA(2),INFFTD(2)
  //      WRITE(RESULT(30),31) STRIKE(3,1),DIP(3,1),
  //     $STRIKE(3,2),DIP(3,2)
  //      WRITE(RESULT(31),32) TREND(1,3),PLUNGE(1,3)
  //      WRITE(RESULT(32),33) TREND(3,3),PLUNGE(3,3)
  //      WRITE(RESULT(33),34) TREND(2,3),PLUNGE(2,3)
  //      WRITE(RESULT(34),37) ETA(3),INFFTD(3)
  //      RETURN
  //      END
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void Taquart::UsmtCore::MOM2(bool REALLY, int QualityType) {
  //      SUBROUTINE MOM2(REALLY)
  //      CHARACTER PS(80),TITLE*40
  //      REAL U(80),ARR(80),AZM(80),TKF(80),A(80,6),ATA(6,6),
  //     $ATAINV(6,6),B(6),EQM(3),AA(80,6),Z1(9,9),Z2(9,9),H(80,5),
  //     $IW(80),PA(3),BB(6),RM0(3),RMT(3),PCLVD(2),PDBCP(2),LLA(3),HA(2)
  //      INTEGER RO(80),VEL(80),R(80)
  //      COMMON/MDATA/ PS,U,ARR,AZM,TKF,RO,VEL,R,TITLE,N,TROZ
  //      COMMON/GAGAGA/ GA(80,3)
  //      COMMON/MOMNT/ RM(6,3)
  //      COMMON/RESTEX/ RESULT
  //      CHARACTER*56 RESULT(34)
  //      DOUBLE PRECISION COV(6,6,3),SAI22
  //      COMMON/VCOVAR/ COV
  //      LOGICAL REALLY

  //      PI=4.*ATAN(1.)
  double PI = 4.0 * atan(1.0);
  int IW[MAXCHANNEL + 1];
  double PA[3 + 1];
  Zero(&PA[0], 4);
  double ATA[6 + 1][6 + 1];
  Zero(&ATA[0][0], 49);
  double ATAINV[6 + 1][6 + 1];
  Zero(&ATAINV[0][0], 49);
  double Z1[9 + 1][9 + 1];
  Zero(&Z1[0][0], 100);
  double Z2[9 + 1][9 + 1];
  Zero(&Z2[0][0], 100);
  double BB[6 + 1];
  Zero(&BB[0], 7);
  double B[6 + 1];
  Zero(&B[0], 7);
  double EQM[3 + 1];
  Zero(EQM, 4);
  double EQQ1 = 0.0, EQQ2 = 0.0, EQQ3 = 0.0;
  double HELP1 = 0.0, HELP2 = 0.0, HELP3 = 0.0;
  double EU = 0.0, EC = 0.0;
  double RM0[3 + 1];
  Zero(RM0, 4);
  double RMT[3 + 1];
  Zero(RMT, 4);
  double RMERR[3 + 1];
  Zero(RMERR, 4);
  double EPS = 0.0;
  double SAI22 = 0.0;
  double SIG = 0.0;
  double AA[MAXCHANNEL + 1][6 + 1];
  Zero(&AA[0][0], (MAXCHANNEL + 1) * 7);
  double DUM = 0.0;
  double H[MAXCHANNEL + 1][5 + 1];
  Zero(&H[0][0], (MAXCHANNEL + 1) * 6);
  double RMX = 0.0, RMY = 0.0, RMZ = 0.0;
  double RMAG = 0.0;
  double PEXPL[4], PCLVD[3 + 1], PDBCP[3 + 1];
  double ALF = 0.0;
  ;
  double MAGN[4];
  double HELP = 0.0;
  //double SVANG = 0.0;
  //double SKAL = 0.0;
  double LLA[4];
  Zero(LLA, 4);
  double HA[4];
  Zero(HA, 4);

  //      DO 3 I=1,N
  for (int i = 1; i <= N; i++) {
    //      IW(I)=0
    //      IF((PS(I).EQ.'S').OR.(PS(I).EQ.'s')) IW(I)=1
    //      IF((PS(I).EQ.'H').OR.(PS(I).EQ.'h')) IW(I)=2
    IW[i] = 0;
    //if(PS[i] == 'S' || PS[i] == 's') IW[i] = 1;
    //if(PS[i] == 'H' || PS[i] == 'h') IW[i] = 2;

    //      ALF=FLOAT(VEL(I))
    //      HELP=4.*PI*FLOAT(RO(I))*ALF*ALF*ALF*FLOAT(R(I))*TROZ*1.E-12
    ALF = VEL[i];

    /* DONE  5 -c2.4.14 : Problem z kalibracj� (USMTCORE) */
    HELP = 4.0 * PI * double(RO[i]) * ALF * ALF * ALF * double(R[i]) * 1.0e-12;

    //      IF(IW(I).NE.0) GO TO 5
    if (IW[i] == 0) {
      //C     For P:
      //      A(I,1)=GA(I,1)*GA(I,1)/HELP
      //      A(I,2)=2.*GA(I,1)*GA(I,2)/HELP
      //      A(I,3)=2.*GA(I,1)*GA(I,3)/HELP
      //      A(I,4)=GA(I,2)*GA(I,2)/HELP
      //      A(I,5)=2.*GA(I,2)*GA(I,3)/HELP
      //      A(I,6)=GA(I,3)*GA(I,3)/HELP
      //      GO TO 3
      A[i][1] = GA[i][1] * GA[i][1] / HELP;
      A[i][2] = 2.0 * GA[i][1] * GA[i][2] / HELP;
      A[i][3] = 2.0 * GA[i][1] * GA[i][3] / HELP;
      A[i][4] = GA[i][2] * GA[i][2] / HELP;
      A[i][5] = 2. * GA[i][2] * GA[i][3] / HELP;
      A[i][6] = GA[i][3] * GA[i][3] / HELP;
    }
    /* TODO -o3.1.19 : Code for SV and SH is switched off by default. */
    /*
     else
     {
     //C     For SV:
     //    5 SVANG=REAL(ARR(I))*PI/180.
     //      PA(3)=-SIN(SVANG)
     //      SKAL=-COS(SVANG)/SQRT(GA(I,1)*GA(I,1)+GA(I,2)*GA(I,2))
     //      PA(1)=GA(I,1)*SKAL
     //      PA(2)=GA(I,2)*SKAL
     SVANG = double(ARR[i]) * PI / 180.0;
     PA[3] = -sin(SVANG);
     SKAL = -cos(SVANG) / sqrt(GA[i][1] * GA[i][1] + GA[i][2] * GA[i][2]);
     PA[1] = GA[i][1] * SKAL;
     PA[2] = GA[i][2] * SKAL;

     //      IF(IW(I).EQ.2) GO TO 6
     if(IW[i] != 2)
     {
     //      A(I,1)=GA(I,1)*PA(1)/HELP
     //      A(I,2)=(GA(I,1)*PA(2)+GA(I,2)*PA(1))/HELP
     //      A(I,3)=(GA(I,1)*PA(3)+GA(I,3)*PA(1))/HELP
     //      A(I,4)=GA(I,2)*PA(2)/HELP
     //      A(I,5)=(GA(I,2)*PA(3)+GA(I,3)*PA(2))/HELP
     //      A(I,6)=GA(I,3)*PA(3)/HELP
     //      GO TO 3
     A[i][1]=GA[i][1]*PA[1]/HELP;
     A[i][2]=(GA[i][1]*PA[2]+GA[i][2]*PA[1])/HELP;
     A[i][3]=(GA[i][1]*PA[3]+GA[i][3]*PA[1])/HELP;
     A[i][4]=GA[i][2]*PA[2]/HELP;
     A[i][5]=(GA[i][2]*PA[3]+GA[i][3]*PA[2])/HELP;
     A[i][6]=GA[i][3]*PA[3]/HELP;
     }
     else
     {
     //C     For SH:
     //    6 LLA(3)=-COS(SVANG)
     //      SKAL=SIN(SVANG)/SQRT(GA(I,1)*GA(I,1)+GA(I,2)*GA(I,2))
     //      LLA(1)=GA(I,1)*SKAL
     //      LLA(2)=GA(I,2)*SKAL
     //      HA(1)=LLA(2)*PA(3)-LLA(3)*PA(2)
     //      HA(2)=-LLA(1)*PA(3)+LLA(3)*PA(1)
     //      A(I,1)=GA(I,1)*HA(1)/HELP+LLA(3)*GA(I,1)*PA(1)/HELP
     //      A(I,2)=(GA(I,1)*HA(2)+GA(I,2)*HA(1))/HELP
     //     $+LLA(3)*(GA(I,1)*PA(2)+GA(I,2)*PA(1))/HELP
     //      A(I,3)=GA(I,3)*HA(1)/HELP+LLA(3)*(GA(I,1)*PA(3)+GA(I,3)*PA(1))
     //     $/HELP
     //      A(I,4)=GA(I,2)*HA(2)/HELP+LLA(3)*GA(I,2)*PA(2)/HELP
     //      A(I,5)=GA(I,3)*HA(2)/HELP+LLA(3)*(GA(I,2)*PA(3)+GA(I,3)*PA(2))
     //     $/HELP
     //      A(I,6)=+LLA(3)*GA(I,3)*PA(3)/HELP

     LLA[3]=-cos(SVANG);
     SKAL = sin(SVANG)/sqrt(GA[i][1]*GA[i][1]+GA[i][2]*GA[i][2]);
     LLA[1]=GA[i][1]*SKAL;
     LLA[2]=GA[i][2]*SKAL;
     HA[1]=LLA[2]*PA[3]-LLA[3]*PA[2];
     HA[2]=-LLA[1]*PA[3]+LLA[3]*PA[1];
     A[i][1]=GA[i][1]*HA[1]/HELP+LLA[3]*GA[i][1]*PA[1]/HELP;
     A[i][2]=(GA[i][1]*HA[2]+GA[i][2]*HA[1])/HELP+LLA[3]*(GA[i][1]*PA[2]+GA[i][2]*PA[1])/HELP;
     A[i][3]=GA[i][3]*HA[1]/HELP+LLA[3]*(GA[i][1]*PA[3]+GA[i][3]*PA[1])/HELP;
     A[i][4]=GA[i][2]*HA[2]/HELP+LLA[3]*GA[i][2]*PA[2]/HELP;
     A[i][5]=GA[i][3]*HA[2]/HELP+LLA[3]*(GA[i][2]*PA[3]+GA[i][3]*PA[2])/HELP;
     A[i][6]=+LLA[3]*GA[i][3]*PA[3]/HELP;
     }
     }
     */
  }
  //    3 CONTINUE

  //      IF(.NOT.REALLY) GO TO 1011
  if (REALLY) {
    //C     #################### FULL TENSOR PART:
    //C     Solve equation AM=U for M:
    //      DO 101 I=1,6
    //      DO 101 J=1,6
    //      ATA(I,J)=0.
    //      DO 101 K=1,N
    //  101 ATA(I,J)=ATA(I,J)+A(K,J)*A(K,I)
    for (int i = 1; i <= 6; i++) {
      for (int j = 1; j <= 6; j++) {
        ATA[i][j] = 0.0;
        for (int k = 1; k <= N; k++)
          ATA[i][j] = ATA[i][j] + A[k][j] * A[k][i];
      }
    }

    //      DO 1002 I=1,6
    //      DO 1002 J=1,6
    // 1002 Z1(I,J)=ATA(I,J)
    for (int i = 1; i <= 6; i++) {
      for (int j = 1; j <= 6; j++) {
        Z1[i][j] = ATA[i][j];
      }
    }

    //      CALL INVMAT(Z1,Z2,6)
    INVMAT(Z1, Z2, 6);

    //      DO 1001 I=1,6
    //      DO 1001 J=1,6
    // 1001 ATAINV(I,J)=Z2(I,J)
    for (int i = 1; i <= 6; i++) {
      for (int j = 1; j <= 6; j++) {
        ATAINV[i][j] = Z2[i][j];
      }
    }

    //      DO 102 I=1,6
    //      B(I)=0.
    //      DO 102 J=1,N
    //  102 B(I)=B(I)+A(J,I)*U(J)*1.E+12
    for (int i = 1; i <= 6; i++) {
      B[i] = 0.0;
      for (int j = 1; j <= N; j++) {
        B[i] = B[i] + A[j][i] * U[j] * 1.0e+12;
      }
    }

    //      DO 103 I=1,6
    //      RM(I,1)=0.
    //      DO 103 J=1,6
    //  103 RM(I,1)=RM(I,1)+ATAINV(I,J)*B(J)
    for (int i = 1; i <= 6; i++) {
      RM[i][1] = 0.0;
      for (int j = 1; j <= 6; j++) {
        RM[i][1] = RM[i][1] + ATAINV[i][j] * B[j];
      }
    }

    //C     Finds scalar seismic moment:
    //      DO 9 I=1,6
    //    9 BB(I)=RM(I,1)
    for (int i = 1; i <= 6; i++)
      BB[i] = RM[i][1];

    //      CALL EIG3(BB,0,EQM)
    EIG3(BB, 0, EQM);

    //      EQQ1=EQM(1)
    //      EQQ2=EQM(2)
    //      EQQ3=EQM(3)
    //      EQM(1)=ABS(EQM(1))*1.E-12
    //      EQM(2)=ABS(EQM(2))*1.E-12
    //      EQM(3)=ABS(EQM(3))*1.E-12
    //      HELP1=ABS(EQM(1)-EQM(2))
    //      HELP2=ABS(EQM(1)-EQM(3))
    //      HELP3=ABS(EQM(2)-EQM(3))
    EQQ1 = EQM[1];
    EQQ2 = EQM[2];
    EQQ3 = EQM[3];
    EQM[1] = fabs(EQM[1]) * 1.0E-12;
    EQM[2] = fabs(EQM[2]) * 1.0E-12;
    EQM[3] = fabs(EQM[3]) * 1.0E-12;
    HELP1 = fabs(EQM[1] - EQM[2]);
    HELP2 = fabs(EQM[1] - EQM[3]);
    HELP3 = fabs(EQM[2] - EQM[3]);

    //      EU=SQRT(.5*(EQM(1)*EQM(1)+EQM(2)*EQM(2)+EQM(3)*EQM(3)))
    //      EC=AMAX1(HELP1,HELP2)
    //      EC=AMAX1(EC,HELP3)
    //      RMT(1)=AMAX1(EC,EU)*1.E+12
    //      RM0(1)=AMIN1(EC,EU)*1.E+12
    EU = sqrt(0.5 * (EQM[1] * EQM[1] + EQM[2] * EQM[2] + EQM[3] * EQM[3]));
    EC = amax1(HELP1, HELP2);
    EC = amax1(EC, HELP3);
    RMT[1] = amax1(EC, EU) * 1.0E+12;
    RM0[1] = amin1(EC, EU) * 1.0E+12;

    //C     Covariance calculation:
    //      DO 401 I=1,N
    for (int i = 1; i <= N; i++) {
      //      EPS=U(I)
      //      DO 406 J=1,6
      //  406 EPS=EPS-A(I,J)*RM(J,1)*1.E-12
      EPS = U[i];

      // Calculate theoretical displacement.
      UTH[i] = 0.0f;
      for (int j = 1; j <= 6; j++)
        UTH[i] += A[i][j] * RM[j][1] * 1.0e-12;

      for (int j = 1; j <= 6; j++)
        EPS = EPS - A[i][j] * RM[j][1] * 1.0e-12;

      //      SAI22=DBLE(0.)
      //      DO 403 J=1,6
      //  403 SAI22=SAI22+A(I,J)*A(I,J)
      SAI22 = 0.0;
      for (int j = 1; j <= 6; j++)
        SAI22 = SAI22 + A[i][j] * A[i][j];

      //      SAI2=SNGL(SAI22)
      double SAI2 = double(SAI22);

      //      DO 404 J=1,6
      //  404 AA(I,J)=EPS*(A(I,J)/SAI2)*1.E+12
      for (int j = 1; j <= 6; j++)
        AA[i][j] = EPS * (A[i][j] / SAI2) * 1.0e+12;

      //      DO 401 J=1,6
      //      AA(I,J)=AA(I,J)+RM(J,1)
      for (int j = 1; j <= 6; j++)
        AA[i][j] = AA[i][j] + RM[j][1];
    }
    //  401 CONTINUE

    //      DO 402 I=1,6
    //      DO 402 J=1,6
    //      COV(I,J,1)=DBLE(0.)
    //      DO 405 K=1,N
    //  405 COV(I,J,1)=COV(I,J,1)+(AA(K,I)-RM(I,1))*(AA(K,J)-RM(J,1))
    //  402 COV(I,J,1)=COV(I,J,1)/DBLE(FLOAT(N-6)**2)
#ifdef USMTCORE_DEBUG
    std::cout << std::endl;
#endif
    for (int i = 1; i <= 6; i++) {
      for (int j = 1; j <= 6; j++) {
        COV[i][j][1] = 0.0;
        for (int k = 1; k <= N; k++)
          COV[i][j][1] = COV[i][j][1]
              + (AA[k][i] - RM[i][1]) * (AA[k][j] - RM[j][1]);
        COV[i][j][1] = COV[i][j][1] / double((N - 6) * (N - 6));
#ifdef USMTCORE_DEBUG
        std::cout << FormatFloat("0.000e+00",COV[i][j][1]).c_str() << " ";
#endif
      }
#ifdef USMTCORE_DEBUG
      std::cout << std::endl;
#endif
    }
#ifdef USMTCORE_DEBUG
    std::cout << std::endl;
#endif

    //      SIG=0.
    //      DO 201 I=1,6
    //      SIGH=SNGL(DSQRT(COV(I,I,1)))
    //      IF(SIG.LT.SIGH) SIG=SIGH
    //  201 CONTINUE
    SIG = 0.0;
    for (int i = 1; i <= 6; i++) {
      double SIGH = sqrt(COV[i][i][1]);
      if (SIG < SIGH) SIG = SIGH;
    }
    RMERR[1] = SIG;

    //C     Output:
    //      WRITE(75,4999)
    // 4999 FORMAT(72(1H=))
    //      WRITE(RESULT(1),5001) TITLE
    // 5001 FORMAT('L2 tensor for : ',A40)
    //      WRITE(75,5000) RESULT(1)
    // 5000 FORMAT(1X,A56)
    //      WRITE(RESULT(2),5002)
    // 5002 FORMAT('  Full solution :',39X)
    //      WRITE(75,5000) RESULT(2)
    //      WRITE(RESULT(3),5003) RM(1,1),RM(2,1),RM(3,1)
    //      WRITE(75,5000) RESULT(3)
    //      WRITE(RESULT(4),5003) RM(2,1),RM(4,1),RM(5,1)
    //      WRITE(75,5000) RESULT(4)
    //      WRITE(RESULT(5),5003) RM(3,1),RM(5,1),RM(6,1)
    //      WRITE(75,5000) RESULT(5)
    // 5003 FORMAT(11X,E10.3,2X,E10.3,2X,E10.3,11X)
    //      WRITE(RESULT(6),5004) TROZ,RM0(1),RMT(1),SIG
    // 5004 FORMAT(2X,'T0=',F8.4,2X,'M0=',E9.3,2X,'MT=',E9.3,2X,'ERR=',E9.3)
    //      WRITE(75,5000) RESULT(6)
#ifdef USMTCORE_DEBUG
    std::cout << "L2 Tensor:" << std::endl << std::endl;
    std::cout << "Full solution:" << std::endl;
    std::cout << RM[1][1] << '\t' << RM[2][1] << '\t' << RM[3][1] << std::endl;
    std::cout << RM[2][1] << '\t' << RM[4][1] << '\t' << RM[5][1] << std::endl;
    std::cout << RM[3][1] << '\t' << RM[5][1] << '\t' << RM[6][1] << std::endl;
    std::cout << "T0 = " << TROZ << " M0 = " << RM0[1] << " MT = " << RMT[1] << " ERR = " << SIG;
    std::cout << std::endl;
#endif

    //      CALL EIGGEN(EQQ1,EQQ2,EQQ3,PEXPL,PCLVD(1),PDBCP(1))
    //      RMAG=.6667*ALOG10(RM0(1))-6.0
    //      IF(RMAG.GE.100.) RMAG=99.9
    //      IF(RMAG.LE.-10.) RMAG=-9.9
    //      WRITE(RESULT(7),5005) PEXPL,PCLVD(1),PDBCP(1),RMAG
    // 5005 FORMAT(2X,'Expl.=',F6.1,'%   CLVD.=',F6.1,'%   DBCP.=',F6.1,
    //     $'%   M=',F4.1)
    //      WRITE(75,5000) RESULT(7)
    EIGGEN(EQQ1, EQQ2, EQQ3, PEXPL[1], PCLVD[1], PDBCP[1]);
    //   EIGGEN_NEW(EQQ1,EQQ2,EQQ3,PEXPL[1],PCLVD[1],PDBCP[1]);
    RMAG = 0.6667 * alog10(RM0[1]) - 6.0;
    if (RMAG >= 100.0) RMAG = 99.9;
    if (RMAG < -10.0) RMAG = -9.9;
    MAGN[1] = RMAG;

#ifdef USMTCORE_DEBUG
    std::cout << " EXPL = " << PEXPL << " PCLVD = " << PCLVD[1];
    std::cout << " PDBCP = " << PDBCP[1] << " Mw = " << RMAG << std::endl;
    std::cout << std::endl;
#endif
  } // if REALLY

  // 1011 DO 55 I=1,N
  //      H(I,1)=A(I,1)-A(I,6)
  //      H(I,2)=A(I,2)
  //      H(I,3)=A(I,3)
  //      H(I,4)=A(I,4)-A(I,6)
  //   55 H(I,5)=A(I,5)
  for (int i = 1; i <= N; i++) {
    H[i][1] = A[i][1] - A[i][6];
    H[i][2] = A[i][2];
    H[i][3] = A[i][3];
    H[i][4] = A[i][4] - A[i][6];
    H[i][5] = A[i][5];
  }

  //C     Solve equation HM=U for M:
  //      DO 2101 I=1,5
  //      DO 2101 J=1,5
  //      ATA(I,J)=0.
  //      DO 2101 K=1,N
  // 2101 ATA(I,J)=ATA(I,J)+H(K,J)*H(K,I)
  for (int i = 1; i <= 5; i++) {
    for (int j = 1; j <= 5; j++) {
      ATA[i][j] = 0.0;
      for (int k = 1; k <= N; k++)
        ATA[i][j] = ATA[i][j] + H[k][j] * H[k][i];
    }
  }

  //      DO 1005 I=1,5
  //      DO 1005 J=1,5
  // 1005 Z1(I,J)=ATA(I,J)
  for (int i = 1; i <= 5; i++) {
    for (int j = 1; j <= 5; j++) {
      Z1[i][j] = ATA[i][j];
    }
  }

  //      CALL INVMAT(Z1,Z2,5)
  INVMAT(Z1, Z2, 5);

  //      DO 1004 I=1,5
  //      DO 1004 J=1,5
  // 1004 ATAINV(I,J)=Z2(I,J)
  for (int i = 1; i <= 5; i++) {
    for (int j = 1; j <= 5; j++) {
      ATAINV[i][j] = Z2[i][j];
    }
  }

  //      DO 2102 I=1,5
  //      B(I)=0.
  //      DO 2102 J=1,N
  // 2102 B(I)=B(I)+H(J,I)*U(J)*1.E+12
  for (int i = 1; i <= 5; i++) {
    B[i] = 0.0;
    for (int j = 1; j <= N; j++) {
      B[i] = B[i] + H[j][i] * U[j] * 1.0e+12;
    }
  }

  //      DO 2103 I=1,5
  //      RM(I,2)=0.
  //      DO 2103 J=1,5
  // 2103 RM(I,2)=RM(I,2)+ATAINV(I,J)*B(J)
  for (int i = 1; i <= 5; i++) {
    RM[i][2] = 0.0;
    for (int j = 1; j <= 5; j++) {
      RM[i][2] = RM[i][2] + ATAINV[i][j] * B[j];
    }
  }

  //      RM(6,2)=-RM(1,2)-RM(4,2)
  RM[6][2] = -RM[1][2] - RM[4][2];

  //C     Finds scalar seismic moment:
  //      DO 2009 I=1,6
  // 2009 BB(I)=RM(I,2)
  for (int i = 1; i <= 6; i++)
    BB[i] = RM[i][2];

  //      CALL EIG3(BB,0,EQM)
  EIG3(BB, 0, EQM);

  //      EQQ1=EQM(1)
  //      EQQ2=EQM(2)
  //      EQQ3=EQM(3)
  //      EQM(1)=ABS(EQM(1))*1.E-12
  //      EQM(2)=ABS(EQM(2))*1.E-12
  //      EQM(3)=ABS(EQM(3))*1.E-12
  //      HELP1=ABS(EQM(1)-EQM(2))
  //      HELP2=ABS(EQM(1)-EQM(3))
  //      HELP3=ABS(EQM(2)-EQM(3))
  EQQ1 = EQM[1];
  EQQ2 = EQM[2];
  EQQ3 = EQM[3];
  EQM[1] = fabs(EQM[1]) * 1.0E-12;
  EQM[2] = fabs(EQM[2]) * 1.0E-12;
  EQM[3] = fabs(EQM[3]) * 1.0E-12;
  HELP1 = fabs(EQM[1] - EQM[2]);
  HELP2 = fabs(EQM[1] - EQM[3]);
  HELP3 = fabs(EQM[2] - EQM[3]);

  //      EU=SQRT(.5*(EQM(1)*EQM(1)+EQM(2)*EQM(2)+EQM(3)*EQM(3)))
  //      EC=AMAX1(HELP1,HELP2)
  //      EC=AMAX1(EC,HELP3)
  //      RMT(2)=AMAX1(EC,EU)*1.E+12
  //      RM0(2)=AMIN1(EC,EU)*1.E+12
  EU = sqrt(0.5 * (EQM[1] * EQM[1] + EQM[2] * EQM[2] + EQM[3] * EQM[3]));
  EC = amax1(HELP1, HELP2);
  EC = amax1(EC, HELP3);
  RMT[2] = amax1(EC, EU) * 1.0E+12;
  RM0[2] = amin1(EC, EU) * 1.0E+12;

  //C     Error calculation:
  //      DO 2401 I=1,N
  for (int i = 1; i <= N; i++) {
    //      EPS=U(I)
    //      DO 2406 J=1,6
    // 2406 EPS=EPS-A(I,J)*RM(J,2)*1.E-12
    EPS = U[i];
    for (int j = 1; j <= 6; j++)
      EPS = EPS - A[i][j] * RM[j][2] * 1.0e-12;

    //      SAI22=DBLE(0.)
    //      DO 2403 J=1,6
    // 2403 SAI22=SAI22+a(I,J)*A(I,J)
    SAI22 = 0.0;
    for (int j = 1; j <= 6; j++)
      SAI22 = SAI22 + A[i][j] * A[i][j];

    //      SAI2=SNGL(SAI22)
    double SAI2 = double(SAI22);

    //      DO 2404 J=1,6
    // 2404 AA(I,J)=EPS*(A(I,J)/SAI2)*1.E+12
    for (int j = 1; j <= 6; j++)
      AA[i][j] = EPS * (A[i][j] / SAI2) * 1.0e+12;

    //      DO 2401 J=1,6
    //      AA(I,J)=AA(I,J)+RM(J,2)
    for (int j = 1; j <= 6; j++)
      AA[i][j] = AA[i][j] + RM[j][2];
  }
  // 2401 CONTINUE

  //      DO 2402 I=1,6
  //      DO 2402 J=1,6
  //      COV(I,J,2)=DBLE(0.)
  //      DO 2405 K=1,N
  // 2405 COV(I,J,2)=COV(I,J,2)+(AA(K,I)-RM(I,2))*(AA(K,J)-RM(J,2))
  // 2402 COV(I,J,2)=COV(I,J,2)/DBLE(FLOAT(N-6)**2)
#ifdef USMTCORE_DEBUG
  std::cout << std::endl;
#endif
  for (int i = 1; i <= 6; i++) {
    for (int j = 1; j <= 6; j++) {
      COV[i][j][2] = 0.0;
      for (int k = 1; k <= N; k++)
        COV[i][j][2] = COV[i][j][2]
            + (AA[k][i] - RM[i][2]) * (AA[k][j] - RM[j][2]);
      COV[i][j][2] = COV[i][j][2] / double((N - 6) * (N - 6));
#ifdef USMTCORE_DEBUG
      std::cout << FormatFloat("0.000e+00",COV[i][j][2]).c_str() << " ";
#endif
    }
#ifdef USMTCORE_DEBUG
    std::cout << std::endl;
#endif
  }
#ifdef USMTCORE_DEBUG
  std::cout << std::endl;
#endif
  //      IF(.NOT.REALLY) RETURN
  if (!REALLY) return;

  //      SIG=0.
  //      DO 202 I=1,6
  //      SIGH=SNGL(DSQRT(COV(I,I,2)))
  //      IF(SIG.LT.SIGH) SIG=SIGH
  //  202 CONTINUE
  SIG = 0.0;
  for (int i = 1; i <= 6; i++) {
    double SIGH = sqrt(COV[i][i][2]);
    if (SIG < SIGH) SIG = SIGH;
  }
  RMERR[2] = SIG;

  //C     in L1 calculation only the RM(?,2) is important
  //C     Output:
  //      WRITE(RESULT(13),5006)
  // 5006 FORMAT('  Trace null solution :',33X)
  //      WRITE(75,5000) RESULT(13)
  //      WRITE(RESULT(14),5003) RM(1,2),RM(2,2),RM(3,2)
  //      WRITE(75,5000) RESULT(14)
  //      WRITE(RESULT(15),5003) RM(2,2),RM(4,2),RM(5,2)
  //      WRITE(75,5000) RESULT(15)
  //      WRITE(RESULT(16),5003) RM(3,2),RM(5,2),RM(6,2)
  //      WRITE(75,5000) RESULT(16)
  //      WRITE(RESULT(17),5004) TROZ,RM0(2),RMT(2),SIG
  //      WRITE(75,5000) RESULT(17)
#ifdef USMTCORE_DEBUG
  std::cout << "Trace-null solution:" << std::endl;
  std::cout << RM[1][2] << '\t' << RM[2][2] << '\t' << RM[3][2] << std::endl;
  std::cout << RM[2][2] << '\t' << RM[4][2] << '\t' << RM[5][2] << std::endl;
  std::cout << RM[3][2] << '\t' << RM[5][2] << '\t' << RM[6][2] << std::endl;
  std::cout << "T0 = " << TROZ << " M0 = " << RM0[2] << " MT = " << RMT[2] << " ERR = " << SIG;
  std::cout << std::endl;
#endif

  //      CALL EIGGEN(EQQ1,EQQ2,EQQ3,DUM,PCLVD(2),PDBCP(2))
  //      RMAG=.6667*ALOG10(RM0(2))-6.0
  //      IF(RMAG.GE.100.) RMAG=99.9
  //      IF(RMAG.LE.-10.) RMAG=-9.9
  //      WRITE(RESULT(18),5007) PCLVD(2),PDBCP(2),RMAG
  // 5007 FORMAT(18X,'CLVD.=',F6.1,'%   DBCP.=',F6.1,'%   M=',F4.1)
  //      WRITE(75,5000) RESULT(18)
  EIGGEN(EQQ1, EQQ2, EQQ3, DUM, PCLVD[2], PDBCP[2]);
  RMAG = 0.6667 * alog10(RM0[2]) - 6.0;
  if (RMAG >= 100.0) RMAG = 99.9;
  if (RMAG < -10.0) RMAG = -9.9;
  MAGN[2] = RMAG;

#ifdef USMTCORE_DEBUG
  std::cout << " PCLVD = " << PCLVD[2];
  std::cout << " PDBCP = " << PDBCP[2] << " Mw = " << RMAG << std::endl;
  std::cout << std::endl;
#endif
  PEXPL[2] = 0.0;

  //==== DOUBLE COUPLE SOLUTION ==============================================

  //C     #################### TRACE 0, DET 0 TENSOR PART:
  //      ISTER=1
  //      CALL FIJGEN
  //ISTER = 1;
  FIJGEN();

  //      DO 1003 I=1,6
  // 1003 B(I)=RM(I,2)
  for (int i = 1; i <= 6; i++)
    B[i] = RM[i][2];

  //      CALL EIG3(B,ISTER,EQM)
  EIG3(B, 1, EQM);

  //      RMM=(EQM(1)+EQM(2)+EQM(3))/3.
  //      RMX=EQM(1)
  //      IF(ABS(EQM(2)).LT.ABS(RMx)) RMX=EQM(2)
  //      IF(ABS(EQM(3)).LT.ABS(RMx)) RMX=EQM(3)
  //RMM = (EQM[1] + EQM[2] + EQM[3]) / 3.0;
  RMX = EQM[1];
  if (fabs(EQM[2]) < fabs(RMX)) RMX = EQM[2];
  if (fabs(EQM[3]) < fabs(RMX)) RMX = EQM[3];

  //      RMY=EQM(3)
  //      IF(ABS(EQM(1)).GT.ABS(RMy)) RMY=EQM(1)
  //      IF(ABS(EQM(2)).GT.ABS(RMy)) RMY=EQM(2)
  RMY = EQM[3];
  if (fabs(EQM[1]) > fabs(RMY)) RMY = EQM[1];
  if (fabs(EQM[2]) > fabs(RMY)) RMY = EQM[2];

  //      IF((EQM(1).NE.RMX).AND.(EQM(1).NE.RMY)) RMZ=EQM(1)
  //      IF((EQM(2).NE.RMX).AND.(EQM(2).NE.RMY)) RMZ=EQM(2)
  //      IF((EQM(3).NE.RMX).AND.(EQM(3).NE.RMY)) RMZ=EQM(3)
  if (EQM[1] != RMX && EQM[1] != RMY) RMZ = EQM[1];
  if (EQM[2] != RMX && EQM[2] != RMY) RMZ = EQM[2];
  if (EQM[3] != RMX && EQM[3] != RMY) RMZ = EQM[3];

  //      CALL BETTER(RMY,RMZ,RM0(3),RMT(3),ICOND)
  BETTER(RMY, RMZ, RM0[3], RMT[3], ICOND);

  //      DO 2003 I=1,6
  // 2003 B(I)=RM(I,3)
  for (int i = 1; i <= 6; i++)
    B[i] = RM[i][3];

  //      CALL EIG3(B,ISTER,EQM)
  EIG3(B, 1, EQM);

  //      EQQ1=EQM(1)
  //      EQQ2=EQM(2)
  //      EQQ3=EQM(3)
  EQQ1 = EQM[1];
  EQQ2 = EQM[2];
  EQQ3 = EQM[3];

  //C     Error calculation:
  //      DO 3401 I=1,N
  for (int i = 1; i <= N; i++) {
    //      EPS=U(I)
    //      DO 3406 J=1,6
    // 3406 EPS=EPS-A(I,J)*RM(J,3)*1.E-12
    EPS = U[i];
    for (int j = 1; j <= 6; j++)
      EPS = EPS - A[i][j] * RM[j][3] * 1.0e-12;

    //      SAI22=DBLE(0.)
    //      DO 3403 J=1,6
    // 3403 SAI22=SAI22+a(I,J)*A(I,J)
    SAI22 = 0.0;
    for (int j = 1; j <= 6; j++)
      SAI22 = SAI22 + A[i][j] * A[i][j];

    //      SAI2=SNGL(SAI22)
    double SAI2 = double(SAI22);

    //      DO 3404 J=1,6
    // 3404 AA(I,J)=EPS*(A(I,J)/SAI2)*1.E+12
    for (int j = 1; j <= 6; j++)
      AA[i][j] = EPS * (A[i][j] / SAI2) * 1.0e+12;

    //      DO 3401 J=1,6
    //      AA(I,J)=AA(I,J)+RM(J,3)
    for (int j = 1; j <= 6; j++)
      AA[i][j] = AA[i][j] + RM[j][3];
  }
  // 3401 CONTINUE

  //      DO 3402 I=1,6
  //      DO 3402 J=1,6
  //      COV(I,J,3)=DBLE(0.)
  //      DO 3405 K=1,N
  // 3405 COV(I,J,3)=COV(I,J,3)+(AA(K,I)-RM(I,3))*(AA(K,J)-RM(J,3))
  // 3402 COV(I,J,3)=COV(I,J,3)/DBLE(FLOAT(N-6)**2)
#ifdef USMTCORE_DEBUG
  std::cout << std::endl;
#endif
  for (int i = 1; i <= 6; i++) {
    for (int j = 1; j <= 6; j++) {
      COV[i][j][3] = 0.0;
      for (int k = 1; k <= N; k++)
        COV[i][j][3] = COV[i][j][3]
            + (AA[k][i] - RM[i][3]) * (AA[k][j] - RM[j][3]);
      COV[i][j][3] = COV[i][j][3] / double((N - 6) * (N - 6));
#ifdef USMTCORE_DEBUG
      std::cout << FormatFloat("0.000e+00",COV[i][j][3]).c_str() << " ";
#endif
    }
#ifdef USMTCORE_DEBUG
    std::cout << std::endl;
#endif
  }
#ifdef USMTCORE_DEBUG
  std::cout << std::endl;
#endif
  //      SIG=0.
  //      DO 203 I=1,6
  //      SIGH=SNGL(DSQRT(COV(I,I,3)))
  //      IF(SIG.LT.SIGH) SIG=SIGH
  //  203 CONTINUE
  SIG = 0.0;
  for (int i = 1; i <= 6; i++) {
    double SIGH = sqrt(COV[i][i][3]);
    if (SIG < SIGH) SIG = SIGH;
  }
  RMERR[3] = SIG;

  //C     Output:
  //      WRITE(RESULT(24),5008)
  // 5008 FORMAT('  Double couple solution :',30X)
  //      WRITE(75,5000) RESULT(24)
  //      WRITE(RESULT(25),5003) RM(1,3),RM(2,3),RM(3,3)
  //      WRITE(75,5000) RESULT(25)
  //      WRITE(RESULT(26),5003) RM(2,3),RM(4,3),RM(5,3)
  //      WRITE(75,5000) RESULT(26)
  //      WRITE(RESULT(27),5003) RM(3,3),RM(5,3),RM(6,3)
  //      WRITE(75,5000) RESULT(27)
  //      WRITE(RESULT(28),5004) TROZ,RM0(3),RMT(3),SIG
  //      WRITE(75,5000) RESULT(28)
#ifdef USMTCORE_DEBUG
  std::cout << "Double-couple solution:" << std::endl;
  std::cout << RM[1][3] << '\t' << RM[2][3] << '\t' << RM[3][3] << std::endl;
  std::cout << RM[2][3] << '\t' << RM[4][3] << '\t' << RM[5][3] << std::endl;
  std::cout << RM[3][3] << '\t' << RM[5][3] << '\t' << RM[6][3] << std::endl;
  std::cout << "T0 = " << TROZ << " M0 = " << RM0[3] << " MT = " << RMT[3] << " ERR = " << SIG;
  std::cout << std::endl;
#endif

  //      RMAG=.6667*ALOG10(RM0(3))-6.0
  //      IF(RMAG.GE.100.) RMAG=99.9
  //      IF(RMAG.LE.-10.) RMAG=-9.9
  //      WRITE(RESULT(29),5009) RMAG
  // 5009 FORMAT(34X,'DBCP.= 100.0%   M=',F4.1)
  //      WRITE(75,5000) RESULT(29)
  RMAG = 0.6667 * alog10(RM0[3]) - 6.0;
  if (RMAG >= 100.0) RMAG = 99.9;
  if (RMAG < -10.0) RMAG = -9.9;
  MAGN[3] = RMAG;
  PDBCP[3] = 100.0;
  PCLVD[3] = 0.0;
  PEXPL[3] = 0.0;

#ifdef USMTCORE_DEBUG
  std::cout << " PDBCP = 100% Mw = " << RMAG << std::endl;
  std::cout << std::endl;
#endif

  //      CALL XTRINF(ICOND)
  //      NORM=2
  //      CALL WROUT(NORM)
  //      CALL WRCOV

  //---- Transfer data to Taquart::FaultSolution structures

  // Theoretical displacement (taken from full solution, partially incorrect).
  for (int i = 1; i <= N; i++) {
    Solution[1].U_n = N;
    Solution[2].U_n = N;
    Solution[3].U_n = N;
    Solution[1].U_th[i - 1] = UTH[i];
    Solution[2].U_th[i - 1] = UTH[i];
    Solution[3].U_th[i - 1] = UTH[i];
    Solution[1].U_measured[i - 1] = U[i];
    Solution[2].U_measured[i - 1] = U[i];
    Solution[3].U_measured[i - 1] = U[i];
    /*
     Solution[1].U_th.push_back(UTH[i]);
     Solution[2].U_th.push_back(UTH[i]);
     Solution[3].U_th.push_back(UTH[i]);
     Solution[1].U_measured.push_back(U[i]);
     Solution[2].U_measured.push_back(U[i]);
     Solution[3].U_measured.push_back(U[i]);
     */
  }

  double uerr = 0.0;
  double umax = -1.0e300;
  double umin = +1.0e300;
  double d = 0.0;
  for (int i = 1; i <= N; i++) {
    d = UTH[i] - U[i];
    uerr = uerr + d * d;
    if (umax < d) umax = d;
    if (umin > d) umin = d;
  }
  uerr = sqrt(uerr / N) / (umax - umin);
  Solution[1].UERR = uerr;
  Solution[2].UERR = uerr;
  Solution[3].UERR = uerr;

  for (int i = 1; i <= 3; i++) {
    int z = 0;
    int v[] = { 1, 2, 3, 2, 4, 5, 3, 5, 6 };
    for (int m = 1; m <= 3; m++)
      for (int n = 1; n <= 3; n++) {
        Solution[i].M[m][n] = RM[v[z]][i];
        z++;
      }
    Solution[i].T0 = TROZ;
    Solution[i].M0 = RM0[i];
    Solution[i].MT = RMT[i];
    Solution[i].ERR = RMERR[i];
    Solution[i].EXPL = PEXPL[i];
    Solution[i].CLVD = PCLVD[i];
    Solution[i].DBCP = PDBCP[i];
    Solution[i].MAGN = MAGN[i];

    for (int m = 1; m <= 6; m++)
      for (int n = 1; n <= 6; n++)
        Solution[i].Covariance[m][n] = COV[m][n][i];
  }

  ICOND = QualityType;
  XTRINF(ICOND, 2, RM0, RMERR);
  // WROUT(NORM);
  // WRCOV();
  //      RETURN
  //      END
}

//-----------------------------------------------------------------------------
void Taquart::UsmtCore::INVMAT(double A[][10], double B[][10], int NP) {
  //      SUBROUTINE INVMAT(A,B,NP)
  //C     Routine inverts matrix A of dimensions 9*9
  //      DIMENSION A(9,9),B(9,9),Y(9,9),ALU(9,9),INDX(9)
  //      EXTERNAL LUDCMP,LUBKSB

  double Y[10][10];
  Zero(&Y[0][0], 100);
  double ALU[10][10];
  Zero(&ALU[0][0], 100);
  int INDX[10];
  for (int i = 0; i < 10; i++)
    INDX[i] = 0;
  double D = 0;

  //      DO 12 I=1,NP
  //       DO 11 J=1,NP
  //        Y(I,J)=0.
  //   11   CONTINUE
  //       Y(I,I)=1.
  //   12 CONTINUE
  for (int i = 1; i <= NP; i++) {
    for (int j = 1; j <= NP; j++) {
      Y[i][j] = 0.0;
    }
    Y[i][i] = 1.0;
  }

  //      CALL LUDCMP(A,ALU,INDX,D,NP)
  LUDCMP(A, ALU, INDX, D, NP);

  //      DO 13 J=1,NP
  //       CALL LUBKSB(ALU,INDX,Y(1,J),B(1,J),NP)
  //   13  CONTINUE
  for (int j = 1; j <= NP; j++) {
    LUBKSB2(ALU, INDX, Y, B, NP, j);
  }

  //      RETURN
  //      END
}

//-----------------------------------------------------------------------------
void Taquart::UsmtCore::LUBKSB2(double A[][10], int INDX[], double C[][10],
    double B[][10], int &NP, int jj) {
  //      SUBROUTINE LUBKSB(A,INDX,C,B,NP)
  //C     Routine solves set of linear equations AX=B. Matrix A is input as its
  //C     LU decomposition. INDX is permutation vector. C is input as right-hand
  //C     side vector and solution is returned as B. A, and INDX are not modified.
  //      DIMENSION A(9,9),INDX(9),B(9),C(9)

  //      DO 1 I=1,NP
  //    1  B(I)=C(I)
  for (int i = 1; i <= NP; i++)
    B[i][jj] = C[i][jj];

  //      II=0
  int II = 0;

  //      DO 12 I=1,NP
  for (int i = 1; i <= NP; i++) {
    //       LL=INDX(I)
    //       SUM=B(LL)
    //       B(LL)=B(I)
    int LL = INDX[i];
    double SUM = B[LL][jj];
    B[LL][jj] = B[i][jj];

    //       IF(II.NE.0) THEN
    if (II != 0) {
      //        DO 11 J=II,I-1
      //         SUM=SUM-A(I,J)*B(J)
      //   11    CONTINUE
      for (int j = II; j <= i - 1; j++)
        SUM = SUM - A[i][j] * B[j][jj];
    }
    //        ELSE IF(SUM.NE.0.) THEN
    else if (SUM != 0.0) {
      //         II=I
      II = i;
    }
    //        ENDIF
    //       B(I)=SUM
    B[i][jj] = SUM;
  }
  //   12  CONTINUE

  //      DO 14 I=NP,1,-1
  for (int i = NP; i >= 1; i--) {
    //       SUM=B(I)
    double SUM = B[i][jj];
    //       IF(I.LT.NP) THEN
    if (i < NP) {
      //        DO 13 J=I+1,NP
      //         SUM=SUM-A(I,J)*B(J)
      //   13    CONTINUE
      for (int j = i + 1; j <= NP; j++)
        SUM = SUM - A[i][j] * B[j][jj];
    }
    //        ENDIF
    //       B(I)=SUM/A(I,I)
    B[i][jj] = SUM / A[i][i];
  }
  //   14  CONTINUE

  //      RETURN
  //      END
  //
  //
}

//-----------------------------------------------------------------------------
void Taquart::UsmtCore::LUDCMP(double B[][10], double A[][10], int INDX[],
    double &D, int &NP) {
  //      SUBROUTINE LUDCMP(B,A,INDX,D,NP)
  //C     Given a matrix A , this routine replaces it by LU decomposition of
  //C     a rowwise permutation of itself. B is input, copied into A. INDX is an
  //C     output vector which records row permutation. D is output as +-1 depending
  //C     on the odd or even number of row interchanges.
  //      DIMENSION A(9,9),B(9,9),INDX(9),VV(100)
  //      DATA TINY/1.E-20/
  double TINY = 1.0e-20;
  double VV[100];
  Zero(&VV[0], 100);
  double AAMAX = 0.0;
  int IMAX = 0;
  double SUM = 0.0;
  double DUM = 0.0;

  //      DO 1 I=1,NP
  //       DO 1 J=1,NP
  //    1   A(I,J)=B(I,J)
  for (int i = 1; i <= NP; i++)
    for (int j = 1; j <= NP; j++)
      A[i][j] = B[i][j];

  //      D=1.
  D = 1.0;

  //      DO 12 I=1,NP
  for (int i = 1; i <= NP; i++) {
    //       AAMAX=0.
    AAMAX = 0.0;
    //       DO 11 j=1,NP
    for (int j = 1; j <= NP; j++) {
      //        IF(ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
      if (fabs(A[i][j]) > AAMAX) AAMAX = fabs(A[i][j]);
    }
    //   11   CONTINUE
    //       IF(AAMAX.EQ.0.) AAMAX=TINY
    if (AAMAX == 0.0) AAMAX = TINY;
    //       VV(I)=1./AAMAX
    VV[i] = 1.0 / AAMAX;
  }
  //   12  CONTINUE

  //      DO 19 J=1,NP
  for (int j = 1; j <= NP; j++) {
    //       DO 14 I=1,J-1
    for (int i = 1; i <= j - 1; i++) {
      //        SUM=A(I,J)
      SUM = A[i][j];
      //        DO 13 K=1,I-1
      for (int k = 1; k <= i - 1; k++) {
        //         SUM=SUM-A(I,K)*A(K,J)
        SUM = SUM - A[i][k] * A[k][j];
      }
      //   13    CONTINUE
      //        A(I,J)=SUM
      A[i][j] = SUM;
    }
    //   14   CONTINUE
    //       AAMAX=0.
    AAMAX = 0.0;
    IMAX = 0;

    //       DO 16 I=J,NP
    for (int i = j; i <= NP; i++) {
      //        SUM=A(I,J)
      SUM = A[i][j];
      //        DO 15 K=1,J-1
      for (int k = 1; k <= j - 1; k++) {
        //         SUM=SUM-A(I,K)*A(K,J)
        SUM = SUM - A[i][k] * A[k][j];
      }
      //   15    CONTINUE
      //        A(I,J)=SUM
      //        DUM=VV(I)*ABS(SUM)
      A[i][j] = SUM;
      DUM = VV[i] * fabs(SUM);
      //        IF(DUM.GE.AAMAX) THEN
      //         IMAX=I
      //         AAMAX=DUM
      //         ENDIF
      if (DUM >= AAMAX) {
        IMAX = i;
        AAMAX = DUM;
      }
    }
    //   16   CONTINUE

    //       IF(J.NE.IMAX) THEN
    if (j != IMAX) {
      //        DO 17 K=1,NP
      for (int k = 1; k <= NP; k++) {
        //         DUM=A(IMAX,K)
        //         A(IMAX,K)=A(J,K)
        //         A(J,K)=DUM
        DUM = A[IMAX][k];
        A[IMAX][k] = A[j][k];
        A[j][k] = DUM;
      }
      //   17    CONTINUE
      //        D=-D
      D = -D;
      //        VV(IMAX)=VV(J)
      VV[IMAX] = VV[j];
    }
    //        ENDIF

    //       INDX(J)=IMAX
    //       IF(A(J,J).EQ.0.) A(J,J)=TINY
    INDX[j] = IMAX;
    if (A[j][j] == 0.0) A[j][j] = TINY;

    //       IF(J.NE.NP) THEN
    //        DUM=1./A(J,J)
    //        DO 18 I=J+1,NP
    //         A(I,J)=A(I,J)*DUM
    //   18    CONTINUE
    //        ENDIF
    if (j != NP) {
      DUM = 1.0 / A[j][j];
      for (int i = j + 1; i <= NP; i++) {
        A[i][j] = A[i][j] * DUM;
      }
    }
  }
  //   19  CONTINUE

  //      RETURN
  //      END
}

//-----------------------------------------------------------------------------
void Taquart::UsmtCore::FIJGEN(void) {
  //      SUBROUTINE FIJGEN
  //      CHARACTER PS(80),TITLE*40
  //      REAL U(80),ARR(80),AZM(80),TKF(80)
  //      INTEGER RO(80),VEL(80),R(80)
  //      COMMON/MDATA/ PS,U,ARR,AZM,TKF,RO,VEL,R,TITLE,N,TROZ
  //      COMMON/GAGAGA/ GA(80,3)
  //      COMMON/AUXIL/ FIJ(3,3,80)

  //      do 3003 i=1,n
  for (int i = 1; i <= N; i++) {
    //      THE=ACOS(GA(I,3))
    //      PHI=ATAN2(GA(I,1),GA(I,2))
    double THE = acos(GA[i][3]);
    /* TODO 5 -c3.1.19 : Correction of ATAN2 error */
    //double PHI = atan2(GA[i][1],GA[i][2);
    double PHI = datan2(GA[i][1], GA[i][2]) * M_PI / 180.0;

    //      FIJ(1,1,I)=(SIN(THE)*SIN(PHI))**2
    //      FIJ(1,2,I)=.5*SIN(THE)**2*SIN(2.*PHI)
    //      FIJ(1,3,I)=-.5*SIN(2.*THE)*SIN(PHI)
    //      FIJ(2,2,I)=(SIN(THE)*COS(PHI))**2
    //      FIJ(2,3,I)=-.5*SIN(2.*THE)*COS(PHI)
    //      FIJ(3,3,I)=COS(PHI)**2
    //      fij(2,1,i)=fij(1,2,i)
    //      fij(3,1,i)=fij(1,3,i)
    // 3003 fij(3,2,i)=fij(2,3,i)
    FIJ[1][1][i] = pow(sin(THE) * sin(PHI), 2.0);
    FIJ[1][2][i] = 0.5 * sin(THE) * sin(THE) * sin(2.0 * PHI);
    FIJ[1][3][i] = -0.5 * sin(2.0 * THE) * sin(PHI);
    FIJ[2][2][i] = pow(sin(THE) * cos(PHI), 2.0);
    FIJ[2][3][i] = -0.5 * sin(2.0 * THE) * cos(PHI);
    FIJ[3][3][i] = cos(PHI) * cos(PHI);
    FIJ[2][1][i] = FIJ[1][2][i];
    FIJ[3][1][i] = FIJ[1][3][i];
    FIJ[3][2][i] = FIJ[2][3][i];
  }

  //      RETURN
  //      END
}

//-----------------------------------------------------------------------------
void Taquart::UsmtCore::Zero(double * Address, int C) {
  for (int i = 0; i < C; i++)
    *(Address + i) = 0.0;
}

//-----------------------------------------------------------------------------
void Taquart::UsmtCore::BETTER(double &RMY, double &RMZ, double &RM0,
    double &RMT, int &ICOND) {
  //      SUBROUTINE BETTER(RMY,RMZ,RM0,RMT,ICOND)
  //      CHARACTER PS(80),TITLE*40
  //      REAL U(80),ARR(80),AZM(80),TKF(80),vn(3),ve(3),Z1(9,9),
  //     $Z2(9,9),DD(6),EQM(3)
  //      INTEGER RO(80),VEL(80),R(80)
  //      COMMON/MDATA/ PS,U,ARR,AZM,TKF,RO,VEL,R,TITLE,N,TROZ
  //      COMMON/GAGAGA/ GA(80,3)
  //      COMMON/MOMNT/ RM(6,3)
  //      COMMON/AUXIL/ FIJ(3,3,80)
  //      COMMON/PRIPAR/ IWUL,ILAS,ISTA,ISHD,IDCSHW,IAXS,IUND
  //      DOUBLE PRECISION U0,U1,C(6,80),BB(8,8),BBINV(8,8),DE(3),EPS,
  //     $DN(3),DU(80),CTDU(8),X(8),dhelp1,dhelp2,sume,sumn,AM0(80)
  double VN[3 + 1];
  Zero(VN, 4);
  double VE[3 + 1];
  Zero(VE, 4);
  double AM0[MAXCHANNEL + 1];
  Zero(AM0, MAXCHANNEL + 1);
  double DD[6 + 1];
  Zero(DD, 7);
  double EQM[3 + 1];
  Zero(EQM, 4);
  double EC = 0.0;
  double C[6 + 1][MAXCHANNEL + 1];
  Zero(&C[0][0], 7 * (MAXCHANNEL + 1));
  double BB[8 + 1][8 + 1];
  Zero(&BB[0][0], 81);
  double BBINV[8 + 1][8 + 1];
  Zero(&BBINV[0][0], 81);
  double Z1[9 + 1][9 + 1];
  Zero(&Z1[0][0], 100);
  double Z2[9 + 1][9 + 1];
  Zero(&Z2[0][0], 100);
  double DE[3 + 1];
  Zero(DE, 4);
  double DU[MAXCHANNEL + 1];
  Zero(DU, MAXCHANNEL + 1);
  double DN[3 + 1];
  Zero(DN, 4);
  double CTDU[8 + 1];
  Zero(CTDU, 9);
  double X[8 + 1];
  Zero(X, 9);
  double DHELP1 = 0.0;
  double DHELP2 = 0.0;

  //      DO 336 I=1,6
  //  336 RM(I,3)=RM(I,2)
  //      RM2=RM(2,3)
  //      RM3=RM(3,3)
  //      RM5=RM(5,3)
  //      RMMAX=0.
  //      IRMAX=0
  for (int i = 1; i <= 6; i++)
    RM[i][3] = RM[i][2];

  double RM2 = RM[2][3];
  double RM3 = RM[3][3];
  double RM5 = RM[5][3];
  double RMMAX = 0.0;
  int IRMAX = 0;

  //      DO 333 I=1,6
  //      IF(ABS(RM(I,3)).LT.ABS(RMMAX)) GO TO 333
  //      RMMAX=RM(I,3)
  //      IRMAX=I
  //  333 CONTINUE
  for (int i = 1; i <= 6; i++) {
    if (fabs(RM[i][3]) < fabs(RMMAX)) continue;
    RMMAX = RM[i][3];
    IRMAX = i;
  }

  //      HELP1=RM(1,3)-RMY
  //      HELP2=RM(4,3)-RMY
  //      HELP3=RM(6,3)-RMY
  //      CALL VEIG(HELP1,RM2,RM3,HELP2,RM5,HELP3,VN)
  //      HELP1=RM(1,3)-RMZ
  //      HELP2=RM(4,3)-RMZ
  //      HELP3=RM(6,3)-RMZ
  //      CALL VEIG(HELP1,RM2,RM3,HELP2,RM5,HELP3,VE)
  double HELP1 = RM[1][3] - RMY;
  double HELP2 = RM[4][3] - RMY;
  double HELP3 = RM[6][3] - RMY;
  VEIG(HELP1, RM2, RM3, HELP2, RM5, HELP3, VN);
  HELP1 = RM[1][3] - RMZ;
  HELP2 = RM[4][3] - RMZ;
  HELP3 = RM[6][3] - RMZ;
  VEIG(HELP1, RM2, RM3, HELP2, RM5, HELP3, VE);

  //      SE=0.
  //      SN=0.
  //      DO 3101 I=1,3
  //      SE=SE+VE(I)*VE(I)
  // 3101 SN=SN+VN(I)*VN(I)
  //      SE=SQRT(SE)
  //      SN=SQRT(SN)
  //      DO 3201 I=1,3
  //      VE(I)=VE(I)/SE
  // 3201 VN(I)=VN(I)/SN
  double SE = 0.0;
  double SN = 0.0;
  for (int i = 1; i <= 3; i++) {
    SE += VE[i] * VE[i];
    SN += VN[i] * VN[i];
  }
  SE = sqrt(SE);
  SN = sqrt(SN);
  for (int i = 1; i <= 3; i++) {
    VE[i] = VE[i] / SE;
    VN[i] = VN[i] / SN;
  }

  //      DO 3203 L=1,N
  //      U0=dble(0.)
  //      U1=dble(0.)
  //      DO 3202 I=1,3
  //      DO 3202 J=1,3
  // 3202 U0=U0+dble(FIJ(I,J,L)*(VN(I)*VE(J)+VN(J)*VE(I)))
  //      U1=U1+dble(U(L))
  // 3203 AM0(L)=U1/U0
  for (int l = 1; l <= N; l++) {
    double U0 = 0.0;
    double U1 = 0.0;
    for (int i = 1; i <= 3; i++) {
      for (int j = 1; j <= 3; j++) {
        U0 = U0 + double(FIJ[i][j][l] * (VN[i] * VE[j] + VN[j] * VE[i]));
      }
    }
    U1 = U1 + double(U[l]);
    AM0[l] = U1 / U0;
  }

  //C     Finds scalar seismic moment:
  //      DO 9 I=1,6
  //    9 DD(I)=RM(I,3)
  //      CALL EIG3(DD,0,EQM)
  for (int i = 1; i <= 6; i++)
    DD[i] = RM[i][3];
  EIG3(DD, 0, EQM);

  //      EQM(1)=ABS(EQM(1))*1.E-12
  //      EQM(2)=ABS(EQM(2))*1.E-12
  //      EQM(3)=ABS(EQM(3))*1.E-12
  //      HELP1=ABS(EQM(1)-EQM(2))
  //      HELP2=ABS(EQM(1)-EQM(3))
  //      HELP3=ABS(EQM(2)-EQM(3))
  //      EC=AMAX1(HELP1,HELP2)
  //      EC=AMAX1(EC,HELP3)
  //      RMT=EC*1.E+12
  //      RM0=RMT
  EQM[1] = fabs(EQM[1]) * 1.0e-12;
  EQM[2] = fabs(EQM[2]) * 1.0e-12;
  EQM[3] = fabs(EQM[3]) * 1.0e-12;
  HELP1 = fabs(EQM[1] - EQM[2]);
  HELP2 = fabs(EQM[1] - EQM[3]);
  HELP3 = fabs(EQM[2] - EQM[3]);
  EC = amax1(HELP1, HELP2);
  EC = amax1(EC, HELP3);
  RMT = EC * 1.0e+12;
  RM0 = RMT;

  //      do 3001 i=1,n
  for (int i = 1; i <= N; i++) {
    //      do 3103 j=1,3
    // 3103 c(j,i)=dble(0.)
    for (int j = 1; j <= 6; j++) //???? Poprawka b��du!
      C[j][i] = 0.0;

    //      do 3001 j=1,3
    for (int j = 1; j <= 3; j++) {
      //      do 3104 k=1,3
      // 3104 c(k,i)=c(k,i)+dble(ve(j)*(fij(k,j,i)+fij(j,k,i)))*AM0(I)
      for (int k = 1; k <= 3; k++)
        C[k][i] = C[k][i]
            + double(VE[j] * (FIJ[k][j][i] + FIJ[j][k][i])) * AM0[i];
      //      do 3105 k=1,3
      // 3105 c(k+3,i)=c(k+3,i)+dble(vn(j)*(fij(k,j,i)+fij(j,k,i)))*AM0(I)
      for (int k = 1; k <= 3; k++)
        C[k + 3][i] = C[k + 3][i]
            + double(VN[j] * (FIJ[k][j][i] + FIJ[j][k][i])) * AM0[i];
    }
  }
  // 3001 continue

  //C     calculate matrix bb and its inverse
  //      do 3008 i=1,6
  //      do 3008 j=1,6
  //      bb(i,j)=dble(0.)
  //      do 3008 k=1,n
  // 3008 bb(i,j)=bb(i,j)+c(i,k)*c(j,k)
  for (int i = 1; i <= 6; i++) {
    for (int j = 1; j <= 6; j++) {
      BB[i][j] = 0.0;
      for (int k = 1; k <= N; k++)
        BB[i][j] = BB[i][j] + C[i][k] * C[j][k];
    }
  }

  //      DO 3204 I=1,3
  //      BB(7,I)=dble(VE(I))
  //      bb(7,I+3)=dble(VN(I))
  //      bb(8,I)=dble(VN(I))
  // 3204 BB(8,I+3)=-dble(VE(I))
  for (int i = 1; i <= 3; i++) {
    BB[7][i] = VE[i];
    BB[7][i + 3] = VN[i];
    BB[8][i] = VN[i];
    BB[8][i + 3] = -VE[i];
  }

  //      bb(7,7)=dble(0.)
  //      bb(7,8)=dble(0.)
  //      bb(8,7)=dble(0.)
  //      bb(8,8)=dble(0.)
  BB[7][7] = 0.0;
  BB[7][8] = 0.0;
  BB[8][7] = 0.0;
  BB[8][8] = 0.0;

  //      do 3009 i=1,8
  //      bb(i,7)=bb(7,i)
  // 3009 bb(i,8)=bb(8,i)
  for (int i = 1; i <= 8; i++) {
    BB[i][7] = BB[7][i];
    BB[i][8] = BB[8][i];
  }

  //      DO 3371 I=1,8
  //      DO 3371 J=1,8
  // 3371 Z1(I,J)=BB(I,J)
  for (int i = 1; i <= 8; i++) {
    for (int j = 1; j <= 8; j++) {
      Z1[i][j] = BB[i][j];
    }
  }

  //      call invmat(Z1,Z2,8)
  INVMAT(Z1, Z2, 8);

  //      DO 3372 I=1,8
  //      DO 3372 J=1,8
  // 3372 BBINV(I,J)=Z2(I,J)
  for (int i = 1; i <= 8; i++) {
    for (int j = 1; j <= 8; j++) {
      BBINV[i][j] = Z2[i][j];
    }
  }

  //C     Initial de and dn are given by rotation of .01 radians:
  //      ICOND=1
  //      SF=SIN(.01)
  //      CF=COS(.01)
  ICOND = 1;
  double SF = sin(0.01);
  double CF = cos(0.01);

  //      DE(1)=DBLE(VE(1)*CF+VE(2)*SF)
  //      DE(2)=DBLE(-VE(1)*SF+VE(2)*CF)
  //      DE(3)=DBLE(VE(3))
  //      DN(1)=DBLE(VN(1)*CF+VN(2)*SF)
  //      DN(2)=DBLE(-VN(1)*SF+VN(2)*CF)
  //      DN(3)=DBLE(VN(3))
  DE[1] = VE[1] * CF + VE[2] * SF;
  DE[2] = -VE[1] * SF + VE[2] * CF;
  DE[3] = VE[3];
  DN[1] = VN[1] * CF + VN[2] * SF;
  DN[2] = -VN[1] * SF + VN[2] * CF;
  DN[3] = VN[3];

  //      DO 3004 i=1,3
  //      DE(i)=DE(i)-DBLE(VE(i))
  // 3004 DN(i)=DN(I)-DBLE(VN(I))
  for (int i = 1; i <= 3; i++) {
    DE[i] = DE[i] - VE[i];
    DN[i] = DN[i] - VN[i];
  }

  //      EPS=1.d-6
  //      iter=1
  double EPS = 1.0e-06;
  int ITER = 1;

  // 3012 do 3002 i=1,n
  p3012: for (int i = 1; i <= N; i++) {
    //      du(i)=dble(0.)
    DU[i] = 0.0;
    //      do 3005 j=1,3
    //      do 3005 k=1,3
    // 3005 du(i)=du(i)+dble(vn(k)*(fij(j,k,i)+fij(k,j,i)))*de(j)
    for (int j = 1; j <= 3; j++)
      for (int k = 1; k <= 3; k++)
        DU[i] = DU[i] + double(VN[k] * (FIJ[j][k][i] + FIJ[k][j][i])) * DE[j];
    //      do 3006 j=1,3
    //      do 3006 k=1,3
    // 3006 du(i)=du(i)+dble(ve(k)*(fij(j,k,i)+fij(k,j,i)))*dn(j)
    for (int j = 1; j <= 3; j++)
      for (int k = 1; k <= 3; k++)
        DU[i] = DU[i] + double(VE[k] * (FIJ[j][k][i] + FIJ[k][j][i])) * DN[j];
    //      du(i)=du(i)*am0(I)
    DU[i] = DU[i] * AM0[i];
  }
  // 3002 continue

  //      do 3007 i=1,6
  //      ctdu(i)=dble(0.)
  //      do 3007 j=1,n
  // 3007 ctdu(i)=ctdu(i)+c(i,j)*du(j)

  for (int i = 1; i <= 6; i++) {
    CTDU[i] = 0.0;
    for (int j = 1; j <= N; j++) {
      CTDU[i] = CTDU[i] + C[i][j] * DU[j];
    }
  }

  //      ctdu(7)=dble(0.)
  //      ctdu(8)=dble(0.)
  CTDU[7] = 0.0;
  CTDU[8] = 0.0;

  //      do 3010 i=1,8
  //      x(i)=dble(0.)
  //      do 3010 j=1,8
  // 3010 x(i)=x(i)+bbinv(i,j)*ctdu(j)
  for (int i = 1; i <= 8; i++) {
    X[i] = 0.0;
    for (int j = 1; j <= 8; j++)
      X[i] = X[i] + BBINV[i][j] * CTDU[j];
  }

  //      IF(ITER.EQ.1000) go TO 3418
  //      IF(ITER.GE.500) EPS=1.d-3
  if (ITER == 1000) goto p3418;
  if (ITER > 500) EPS = 1.0e-03;

  //      DO 3428 i=1,8
  //      if(dabs(x(i)).gt.1.d+99) go to 3417
  // 3428 continue
  for (int i = 1; i <= 8; i++) {
    if (fabs(X[i]) > 1.0e+99) goto p3417;
  }

  //      do 3013 i=1,3
  //      dhelp1=dabs(x(i+3)-de(i))
  //      dhelp2=dabs(de(i))*eps
  //      if((dhelp1.gt.dhelp2).and.(dabs(de(i)).gt.1.d-6)) go to 3014
  // 3013 continue
  for (int i = 1; i <= 3; i++) {
    DHELP1 = fabs(X[i + 3] - DE[i]);
    DHELP2 = fabs(DE[i]) * EPS;
    if (DHELP1 > DHELP2 && fabs(DE[i]) > 1.0e-06) goto p3014;
  }

  //      do 3015 i=1,3
  //      dhelp1=dabs(x(i)-dn(i))
  //      dhelp2=dabs(dn(i))*eps
  //      if((dhelp1.gt.dhelp2).and.(dabs(dn(i)).gt.1.d-6)) go to 3014
  // 3015 continue
  for (int i = 1; i <= 3; i++) {
    DHELP1 = fabs(X[i] - DN[i]);
    DHELP2 = fabs(DN[i]) * EPS;
    if (DHELP1 > DHELP2 && fabs(DN[i]) > 1.0e-06) goto p3014;
  }
  //      go to 3017
  goto p3017;

  // 3014 dhelp1=dble(0.)
  //      dhelp2=dble(0.)
  p3014: DHELP1 = 0.0;
  DHELP2 = 0.0;

  //      do 3011 i=1,3
  //      de(i)=x(i+3)
  //      dn(i)=x(i)
  //      dhelp1=dhelp1+de(i)*de(i)
  // 3011 dhelp2=dhelp2+dn(i)*dn(i)
  for (int i = 1; i <= 3; i++) {
    DE[i] = X[i + 3];
    DN[i] = X[i];
    DHELP1 = DHELP1 + DE[i] * DE[i];
    DHELP2 = DHELP2 + DN[i] * DN[i];
  }

  //      if((sngl(dhelp1).ge.1.).or.(sngl(dhelp2).ge.1.)) go to 3417
  if (DHELP1 >= 1.0 || DHELP2 >= 1.0) goto p3417;
  //      iter=iter+1
  //      go to 3012
  ITER++;
  goto p3012;

  // 3418 ICOND=2
  //      GO TO 3419
  p3418: ICOND = 2;
  goto p3419;

  // 3417 ICOND=3
  p3417: ICOND = 3;
  // 3419 IF(ISTA.EQ.0) GO TO 8887
  p3419: if (ISTA == 0) goto p8887;
  //      WRITE(75,400)
  //  400 FORMAT(' Warning: iteration does not converge - initial eigen',
  //     $'vectors used',/,'          Double couple solution may be art',
  //     $'ificial.')
  //      go to 8887
#ifdef USMTCORE_DEBUG
  std::cout << "Iteration does not converge - initial eigenvectors used. Double "
  "couple solution may be artificial.";
#endif
  goto p8887;

  // 3017 CALL ORT(VE,VN,DE,DN)
  p3017: ORT(VE, VN, DE, DN);

  //c     These ve and vn correspond to force coordinates XY;
  //c     they need to be converted to displacement axes PT.
  // 8887 do 8888 i=1,3
  //      de(i)=dble(ve(i)+vn(i))
  // 8888 dn(i)=dble(vn(i)-ve(i))
  p8887: for (int i = 1; i <= 3; i++) {
    DE[i] = VE[i] + VN[i];
    DN[i] = VN[i] - VE[i];
  }

  //      sume=dble(0.)
  //      sumn=dble(0.)
  double SUME = 0.0;
  double SUMN = 0.0;

  //      do 8889 i=1,3
  //      sume=sume+de(i)*de(i)
  // 8889 sumn=sumn+dn(i)*dn(i)
  for (int i = 1; i <= 3; i++) {
    SUME = SUME + DE[i] * DE[i];
    SUMN = SUMN + DN[i] * DN[i];
  }

  //      sume=dsqrt(sume)
  //      sumn=dsqrt(sumn)
  SUME = sqrt(SUME);
  SUMN = sqrt(SUMN);

  //      do 8890 i=1,3
  //      ve(i)=sngl(de(i)/sume)
  // 8890 vn(i)=sngl(dn(i)/sumn)
  for (int i = 1; i <= 3; i++) {
    VE[i] = DE[i] / SUME;
    VN[i] = DN[i] / SUMN;
  }

  //      r1=rm(1,3)
  //      r2=rm(4,3)
  //      r3=rm(6,3)
  double R1 = RM[1][3];
  double R2 = RM[4][3];
  double R3 = RM[6][3];

  //      rm(1,3)=2.*ve(1)*vn(1)
  //      rm(2,3)=ve(1)*vn(2)+ve(2)*vn(1)
  //      rm(3,3)=ve(1)*vn(3)+ve(3)*vn(1)
  //      rm(4,3)=2.*ve(2)*vn(2)
  //      rm(5,3)=ve(2)*vn(3)+ve(3)*vn(2)
  //      rm(6,3)=2.*ve(3)*vn(3)
  RM[1][3] = 2.0 * VE[1] * VN[1];
  RM[2][3] = VE[1] * VN[2] + VE[2] * VN[1];
  RM[3][3] = VE[1] * VN[3] + VE[3] * VN[1];
  RM[4][3] = 2.0 * VE[2] * VN[2];
  RM[5][3] = VE[2] * VN[3] + VE[3] * VN[2];
  RM[6][3] = 2.0 * VE[3] * VN[3];

  //      j=1
  int j = 1;
  //      if(abs(r2).gt.abs(r1)) then
  if (fabs(R2) > fabs(R1)) {
    //      j=4
    //      r1=r2
    //      endif
    j = 4;
    R1 = R2;
  }
  //      if(abs(r3).gt.abs(r1)) then
  if (fabs(R3) > fabs(R1)) {
    //      j=6
    //      r1=r3
    //      endif
    j = 6;
    R1 = R3;
  }
  //      help=r1*rm(j,3)
  //      if(help.ge.0.) go to 8899
  double HELP = R1 * RM[j][3];
  if (HELP < 0.0) {
    //      do 8898 i=1,6
    // 8898 rm(i,3)=-rm(i,3)
    for (int i = 1; i <= 6; i++)
      RM[i][3] = -RM[i][3];
  }
  // 8899 do 3020 i=1,6
  // 3020 rm(i,3)=rm(i,3)*rm0
  for (int i = 1; i <= 6; i++)
    RM[i][3] = RM[i][3] * RM0;

  //      help=sign(1.,rmmax)*rm(irmax,3)
  //      if(help.gt.0.) go to 334
  HELP = sign(1.0, RMMAX) * RM[IRMAX][3];
  if (HELP <= 0.0) {
    //      do 335 i=1,6
    //  335 rm(i,3)=-rm(i,3)
    for (int i = 1; i <= 6; i++)
      RM[i][3] = -RM[i][3];
  }
  //  334 return
  //      end
}

//-----------------------------------------------------------------------------
void Taquart::UsmtCore::VEIG(double &s1, double &s2, double &s3, double &s4,
    double &s5, double &s6, double v[]) {
  //      SUBROUTINE VEIG(S1,S2,S3,S4,S5,S6,V)
  //      dimension v(3)
  //      double precision help0(6),help1,help2
  double help0[6 + 1];
  double help1 = 0.0;
  double help2 = 0.0;

  //      v(1)=1.
  //      v(2)=0.
  //      v(3)=0.
  v[1] = 1.0;
  v[2] = 0.0;
  v[3] = 0.0;

  //      help0(1)=s4*s6-s5*s5
  //      help0(2)=s1*s6-s3*s3
  //      help0(3)=s1*s4-s2*s2
  //      help0(4)=s2*s5-s4*s3
  //      help0(5)=s2*s6-s3*s5
  //      help0(6)=s1*s5-s2*s3
  help0[1] = s4 * s6 - s5 * s5;
  help0[2] = s1 * s6 - s3 * s3;
  help0[3] = s1 * s4 - s2 * s2;
  help0[4] = s2 * s5 - s4 * s3;
  help0[5] = s2 * s6 - s3 * s5;
  help0[6] = s1 * s5 - s2 * s3;

  //      IMAX=1
  int IMAX = 1;

  //      DO 600 I=2,6
  //      IF(DABS(HELP0(I)).LE.DABS(HELP0(1))) GO TO 600
  //      HELP0(1)=HELP0(I)
  //      IMAX=I
  //  600 CONTINUE
  for (int i = 2; i <= 6; i++) {
    if (fabs(help0[i]) <= fabs(help0[1])) continue;
    help0[1] = help0[i];
    IMAX = i;
  }

  //      GO TO (601,602,603,604,605,606),IMAX
  //  601 help1=s5*s3-s2*s6
  //      help2=s5*s2-s4*s3
  //      v(2)=help1/help0(1)
  //      v(3)=help2/help0(1)
  //      GO TO 610
  //  602 v(2)=1.
  //      help1=s3*s5-s2*s6
  //      help2=s3*s2-s1*s5
  //      v(1)=help1/help0(1)
  //      v(3)=help2/help0(1)
  //      GO TO 610
  //  603 v(3)=1.
  //      help1=s2*s5-s3*s4
  //      help2=s2*s3-s1*s5
  //      v(1)=help1/help0(1)
  //      v(2)=help2/help0(1)
  //      GO TO 610
  //  604 v(3)=1.
  //      help1=s4*s6-s5*s5
  //      help2=s3*s5-s2*s6
  //      v(1)=help1/help0(1)
  //      v(2)=help2/help0(1)
  //      GO TO 610
  //  605 v(2)=1.
  //      help1=s5*s5-s4*s6
  //      help2=s3*s4-s2*s5
  //      v(1)=help1/help0(1)
  //      v(3)=help2/help0(1)
  //      GO TO 610
  //  606 v(2)=1.
  //      help1=s3*s4-s2*s5
  //      help2=s2*s2-s1*s4
  //      v(1)=help1/help0(1)
  //      v(3)=help2/help0(1)
  switch (IMAX) {
    case 1:
      help1 = s5 * s3 - s2 * s6;
      help2 = s5 * s2 - s4 * s3;
      v[2] = help1 / help0[1];
      v[3] = help2 / help0[1];
      break;
    case 2:
      v[2] = 1.0;
      help1 = s3 * s5 - s2 * s6;
      help2 = s3 * s2 - s1 * s5;
      v[1] = help1 / help0[1];
      v[3] = help2 / help0[1];
      break;
    case 3:
      v[3] = 1.0;
      help1 = s2 * s5 - s3 * s4;
      help2 = s2 * s3 - s1 * s5;
      v[1] = help1 / help0[1];
      v[2] = help2 / help0[1];
      break;
    case 4:
      v[3] = 1.0;
      help1 = s4 * s6 - s5 * s5;
      help2 = s3 * s5 - s2 * s6;
      v[1] = help1 / help0[1];
      v[2] = help2 / help0[1];
      break;
    case 5:
      v[2] = 1.0;
      help1 = s5 * s5 - s4 * s6;
      help2 = s3 * s4 - s2 * s5;
      v[1] = help1 / help0[1];
      v[3] = help2 / help0[1];
      break;
    case 6:
      v[2] = 1.0;
      help1 = s3 * s4 - s2 * s5;
      help2 = s2 * s2 - s1 * s4;
      v[1] = help1 / help0[1];
      v[3] = help2 / help0[1];
      break;
    default:
      throw;
  }
}

//-----------------------------------------------------------------------------
void Taquart::UsmtCore::ORT(double VE[], double VN[], double DE[],
    double DN[]) {
  //      SUBROUTINE ORT(VE,VN,DE,DN)
  //      DIMENSION VE(3),VN(3),WE(3),WN(3),B(3),C(3),G(3),V1(3),V2(3)
  //      DOUBLE PRECISION DE(3),DN(3)
  double WE[3 + 1], WN[3 + 1], C[3 + 1], B[3 + 1], V1[3 + 1], V2[3 + 1],
      G[3 + 1];
  Zero(WE, 4);
  Zero(WN, 4);
  Zero(C, 4);
  Zero(B, 4);
  Zero(V1, 4);
  Zero(V2, 4);
  Zero(G, 4);
  double HELP = 0.0;
  double HELP2 = 0.0;

  //      DO 1 I=1,3
  //      WE(I)=VE(I)+SNGL(DE(I))
  //    1 WN(I)=VN(I)+SNGL(DN(I))
  for (int i = 1; i <= 3; i++) {
    WE[i] = VE[i] + DE[i];
    WN[i] = VN[i] + DN[i];
  }

  //      HELP=SQRT(WE(1)*WE(1)+WE(2)*WE(2)+WE(3)*WE(3))
  HELP = sqrt(WE[1] * WE[1] + WE[2] * WE[2] + WE[3] * WE[3]);

  //      DO 2 I=1,3
  //    2 WE(I)=WE(I)/HELP
  for (int i = 1; i <= 3; i++)
    WE[i] = WE[i] / HELP;

  //      HELP=SQRT(WN(1)*WN(1)+WN(2)*WN(2)+WN(3)*WN(3))
  HELP = sqrt(WN[1] * WN[1] + WN[2] * WN[2] + WN[3] * WN[3]);

  //      DO 4 I=1,3
  //    4 WN(I)=WN(I)/HELP
  for (int i = 1; i <= 3; i++)
    WN[i] = WN[i] / HELP;

  //      HELP=ABS(WE(1)*WN(1)+WE(2)*WN(2)+WE(3)*WN(3))
  //      IF(HELP.GT..99) RETURN
  HELP = fabs(WE[1] * WN[1] + WE[2] * WN[2] + WE[3] * WN[3]);
  if (HELP > 0.99) return;

  //      DO 3 I=1,3
  //    3 C(I)=.5*(WN(I)+WE(I))
  for (int i = 1; i <= 3; i++)
    C[i] = 0.5 * (WN[i] + WE[i]);

  //      B(1)=WN(2)*WE(3)-WN(3)*WE(2)
  //      B(2)=-WN(1)*WE(3)+WN(3)*WE(1)
  //      B(3)=WN(1)*WE(2)-WN(2)*WE(1)
  //      HELP=SQRT(B(1)*B(1)+B(2)*B(2)+B(3)*B(3))
  B[1] = WN[2] * WE[3] - WN[3] * WE[2];
  B[2] = -WN[1] * WE[3] + WN[3] * WE[1];
  B[3] = WN[1] * WE[2] - WN[2] * WE[1];
  HELP = sqrt(B[1] * B[1] + B[2] * B[2] + B[3] * B[3]);

  //      DO 5 I=1,3
  //    5 B(I)=B(I)/HELP
  for (int i = 1; i <= 3; i++)
    B[i] = B[i] / HELP;

  //      G(1)=B(2)*C(3)-B(3)*C(2)
  //      G(2)=-B(1)*C(3)+B(3)*C(1)
  //      G(3)=B(1)*C(2)-B(2)*C(1)
  G[1] = B[2] * C[3] - B[3] * C[2];
  G[2] = -B[1] * C[3] + B[3] * C[1];
  G[3] = B[1] * C[2] - B[2] * C[1];

  //      DO 6 I=1,3
  //      V1(I)=C(I)+G(I)
  //    6 V2(I)=C(I)-G(I)
  for (int i = 1; i <= 3; i++) {
    V1[i] = C[i] + G[i];
    V2[i] = C[i] - G[i];
  }

  //      HELP=SQRT(V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3))
  HELP = sqrt(V1[1] * V1[1] + V1[2] * V1[2] + V1[3] * V1[3]);

  //      DO 10 I=1,3
  //   10 V1(I)=V1(I)/HELP
  for (int i = 1; i <= 3; i++)
    V1[i] = V1[i] / HELP;

  //      HELP=SQRT(V2(1)*V2(1)+V2(2)*V2(2)+V2(3)*V2(3))
  HELP = sqrt(V2[1] * V2[1] + V2[2] * V2[2] + V2[3] * V2[3]);
  //      DO 11 I=1,3
  //   11 V2(I)=V2(I)/HELP
  for (int i = 1; i <= 3; i++)
    V2[i] = V2[i] / HELP;

  //      HELP=V1(1)*WN(1)+V1(2)*WN(2)+V1(3)*WN(3)
  //      HELP2=V2(1)*WN(1)+V2(2)*WN(2)+V2(3)*WN(3)
  HELP = V1[1] * WN[1] + V1[2] * WN[2] + V1[3] * WN[3];
  HELP2 = V2[1] * WN[1] + V2[2] * WN[2] + V2[3] * WN[3];

  //      IF(HELP.LT.HELP2) GO TO 9
  if (HELP >= HELP2) {
    //      DO 8 I=1,3
    //      VN(I)=V1(I)
    //    8 VE(I)=V2(I)
    for (int i = 1; i <= 3; i++) {
      VN[i] = V1[i];
      VE[i] = V2[i];
    }
    //      RETURN
  }
  else {
    //    9 DO 7 I=1,3
    //      VN(I)=V2(I)
    //    7 VE(I)=V1(I)
    for (int i = 1; i <= 3; i++) {
      VN[i] = V2[i];
      VE[i] = V1[i];
    }
  }
  //      RETURN
  //      END
}

//-----------------------------------------------------------------------------
/*
 inline double Taquart::UsmtCore::alog(double Value)
 {
 return log(Value);
 }

 //-----------------------------------------------------------------------------
 inline double Taquart::UsmtCore::alog10(double Value)
 {
 return log10(Value);
 }

 //-----------------------------------------------------------------------------
 inline double Taquart::UsmtCore::amax1(double value1, double value2)
 {
 return value1 > value2 ? value1 : value2;
 }

 //-----------------------------------------------------------------------------
 inline double Taquart::UsmtCore::amin1(double value1, double value2)
 {
 return value1 > value2 ? value2 : value1;
 }

 //-----------------------------------------------------------------------------
 inline double Taquart::UsmtCore::sign(double value1, double value2)
 {
 if(value2 >= 0.0)
 return fabs(value1);
 else return -fabs(value1);
 }
 */

//-----------------------------------------------------------------------------
bool Taquart::UsmtCore::ANGGA(void) {
  //      DETOPI=4.*ATAN(1.)/180.
  //      IOK=.TRUE.
  //      IF(N.LT.8) GO TO 2
  const double DETOPI = 4.0 * atan(1.0) / 180.0;
  if (N >= MIN_ALLOWED_CHANNELS) {
    //      DO 1 I=1,N
    for (int i = 1; i <= N; i++) {
      //      HELP=TKF(I)
      //      IF(HELP.EQ.90.0) HELP=89.75
      //      GA(I,3)=COS(HELP*DETOPI)
      //      HELP=SQRT(1.-GA(I,3)*GA(I,3))
      //      GA(I,1)=COS(AZM(I)*DETOPI)*HELP
      //      GA(I,2)=SIN(AZM(I)*DETOPI)*HELP
      double HELP = TKF[i];
      if (HELP == 90.0) HELP = 89.75;
      GA[i][3] = cos(HELP * DETOPI);
      HELP = sqrt(1.0 - GA[i][3] * GA[i][3]);
      GA[i][1] = cos(AZM[i] * DETOPI) * HELP;
      GA[i][2] = sin(AZM[i] * DETOPI) * HELP;
    }
    return true;
  }
  else {
    //    2 IOK=.FALSE.
    //      WRITE(75,3) N
    //    3 FORMAT(' Error: Angle_calc failed: N = ',I3,' (must be > 7).')
    //      RETURN
    //      END
    return false;
  }
}

//-----------------------------------------------------------------------------
bool Taquart::UsmtCore::JEZ(void) {
  //      SUBROUTINE JEZ(IOK)
  //      CHARACTER YN
  //      INTEGER*2 QF
  //      LOGICAL IOK,FSTCLL,USEDAE(212)
  //      COMMON/MDATA/ N,TROZ
  //      COMMON/GAGAGA/ GA(80,3)
  //      COMMON/JEZQUA/ QSD,QF
  //      DIMENSION NDAE(9),DAE(212,3)
  //      DATA FSTCLL,NDAE/.TRUE.,2*36,2*32,2*24,16,8,4/

  bool USEDAE[213];
  double DAE[212 + 1][3 + 1];
  Zero(&DAE[0][0], 213 * 4);

  //      IOK=.TRUE.

  //      IF(.NOT.FSTCLL) GO TO 10
  if (FSTCLL) {
    //      FSTCLL=.FALSE.
    //      PI=4.*ATAN(1.)
    //      DETOPI=PI/180.
    //      K=0
    FSTCLL = false;
    const double PI = 4.0 * atan(1.0);
    const double DETOPI = PI / 180.0;
    int k = 0;

    //      DO 1 I=1,9
    for (int i = 1; i <= 9; i++) {
      //      ANG=(FLOAT(I-1)*10.+5.)*DETOPI
      double ANG = (double(i - 1) * 10.0 + 5.0) * DETOPI;
      //      DO 2 J=1,NDAE(I)
      for (int j = 1; j <= NDAE[i]; j++) {
        //      K=K+1
        //      DAE(K,3)=SIN(ANG)
        //      SKAL=COS(ANG)
        //      HELP=FLOAT(J)/FLOAT(NDAE(I))*2.*PI
        //      DAE(K,1)=COS(HELP)*SKAL
        //    2 DAE(K,2)=SIN(HELP)*SKAL
        k = k + 1;
        DAE[k][3] = sin(ANG);
        double SKAL = cos(ANG);
        double HELP = double(j) / double(NDAE[i]) * 2.0 * PI;
        DAE[k][1] = cos(HELP) * SKAL;
        DAE[k][2] = sin(HELP) * SKAL;
      }
    }
    //    1 CONTINUE
  }

  //   10 DO 3 I=1,212
  //    3 USEDAE(I)=.TRUE.
  for (int i = 1; i <= 213; i++)
    USEDAE[i] = true;

  //      NLIVE=212
  int NLIVE = 212;

  //      DO 5 I=1,N
  //      DO 5 J=1,212
  for (int i = 1; i <= N; i++) {
    for (int j = 1; j <= 212; j++) {
      //      IF(.NOT.USEDAE(J)) GO TO 5
      if (USEDAE[j] == false) continue;

      //      HELP=ABS(GA(I,1)*DAE(J,1)+GA(I,2)*DAE(J,2)+GA(I,3)*DAE(J,3))
      double HELP = fabs(
          GA[i][1] * DAE[j][1] + GA[i][2] * DAE[j][2] + GA[i][3] * DAE[j][3]);

      //      IF(HELP.LT..9659) GO TO 5
      if (HELP < 0.9659) continue;

      //      USEDAE(J)=.FALSE.
      //      NLIVE=NLIVE-1
      USEDAE[j] = false;
      NLIVE = NLIVE - 1;
    }
  }
  //    5 CONTINUE

  //      IF(NLIVE.LE.200) GO TO 4
  if (NLIVE <= 200) {
    //    4 SUM1=0.
    //      SUM2=0.
    //      SUM3=0.
    //      SUM4=0.
    //      SUM5=0.
    //      SUM6=0.
    double SUM1 = 0.0, SUM2 = 0.0, SUM3 = 0.0, SUM4 = 0.0, SUM5 = 0.0, SUM6 =
        0.0;

    //      DO 50 I=1,N
    for (int i = 1; i <= N; i++) {
      //      SUM1=SUM1+GA(I,1)*GA(I,1)
      //      SUM2=SUM2+GA(I,1)*GA(I,2)
      //      SUM3=SUM3+GA(I,1)*GA(I,3)
      //      SUM4=SUM4+GA(I,2)*GA(I,2)
      //      SUM5=SUM5+GA(I,2)*GA(I,3)
      //   50 SUM6=SUM6+GA(I,3)*GA(I,3)
      SUM1 = SUM1 + GA[i][1] * GA[i][1];
      SUM2 = SUM2 + GA[i][1] * GA[i][2];
      SUM3 = SUM3 + GA[i][1] * GA[i][3];
      SUM4 = SUM4 + GA[i][2] * GA[i][2];
      SUM5 = SUM5 + GA[i][2] * GA[i][3];
      SUM6 = SUM6 + GA[i][3] * GA[i][3];
    }
    //      QSD=SUM1*SUM4*SUM6+2*SUM2*SUM3*SUM5-SUM3*SUM3*SUM4-
    //     $SUM1*SUM5*SUM5-SUM2*SUM2*SUM6
    //      QSD=ABS(QSD)/FLOAT(N**3)
    //      IF(QF.EQ.2) QSD=1./(1.-ALOG(QSD))
    //      IF(QSD.LT..01) QSD=.01
    QSD = SUM1 * SUM4 * SUM6 + 2 * SUM2 * SUM3 * SUM5 - SUM3 * SUM3 * SUM4
        - SUM1 * SUM5 * SUM5 - SUM2 * SUM2 * SUM6;
    QSD = fabs(QSD) / double(N * N * N);
    if (QF == 2) QSD = 1.0 / (1.0 - alog(QSD));
    if (QSD < 0.01) QSD = 0.01;
    //      RETURN
    return true;
  }
  else {
    //      WRITE(75,100)
    //  100 FORMAT(' Error: Poor Station distribution.',$)
    //      IOK=.FALSE.
    //	RETURN
    //IOK = false;
    return false; // Poor station distribution!
  }
  //      END
}

//-----------------------------------------------------------------------------
void Taquart::UsmtCore::PROGRESS(double Progress, double Max) {
  //std::cout << '\r' << Progress / Max * 100 << std::endl;
#ifdef USMTCORE_DEBUG
  std::cout << "|";
  if(Progress == Max)
  std::cout << std::endl;
#endif

  if (ThreadProgress) *ThreadProgress = int(Progress);

}

//-----------------------------------------------------------------------------
void Taquart::UsmtCore::GSOL(double x[], int &iexp) {
  //      subroutine gsol(x,iexp)
  //      dimension x(6),ix(6)
  //      double precision xlo(6),xhi(6),xstep(6),six,size,xtry(6),
  //     $val,try
  //      DATA six,METH/6.D+00,1/
  int ix[6 + 1]; //dimension
  double xlo[6 + 1];
  Zero(xlo, 7);
  double xhi[6 + 1];
  Zero(xhi, 7);
  double xstep[6 + 1];
  Zero(xstep, 7);
  double six = 6.0e+00;
  double xtry[6 + 1];
  double val = 0.0;
  double tryy = 0.0;
  //int METH = 1;
  //double size = 0.0;

  //      IF((IEXP.LT.10).OR.(IEXP.GT.30)) IEXP=20
  if (iexp < 10 || iexp > 30) iexp = 20;

  //      iter=0
  //      jter=1
  //      val=1.d+30
  int iter = 0;
  //int jter = 1;
  val = 1.0e+30;

  //      do 1 i=1,6
  //      xlo(i)=DBLE(-1.*10.**IEXP)
  //      xhi(i)=DBLE(10.**IEXP)
  //    1 ix(i)=0
  for (int i = 1; i <= 6; i++) {
    xlo[i] = -1.0 * pow(10.0, iexp);
    xhi[i] = pow(10.0, iexp);
    ix[i] = 0;
  }

  //      DO 8 L=1,50
  for (int l = 1; l <= 50; l++) {
    PROGRESS(l, 350);
    //      do 2 i=1,6
    //    2 xstep(i)=(xhi(i)-xlo(i))/SIX
    for (int i = 1; i <= 6; i++)
      xstep[i] = (xhi[i] - xlo[i]) / six;

    //      size=xhi(1)-xlo(1)
    //      iter=iter+1
    //size = xhi[1] - xlo[1];
    iter = iter + 1;

    //      do 3 j1=1,7
    for (int j1 = 1; j1 <= 7; j1++) {
      //      xtry(1)=xlo(1)+DBLE(j1-1)*xstep(1)
      //      CALL POSTEP(METH,JTER,J1)
      xtry[1] = xlo[1] + double(j1 - 1) * xstep[1];
      //POSTEP(METH,jter,j1);

      //      do 3 j2=1,7
      for (int j2 = 1; j2 <= 7; j2++) {
        //      xtry(2)=xlo(2)+DBLE(j2-1)*xstep(2)
        xtry[2] = xlo[2] + double(j2 - 1) * xstep[2];

        //      do 3 j3=1,7
        for (int j3 = 1; j3 <= 7; j3++) {
          //      xtry(3)=xlo(3)+DBLE(j3-1)*xstep(3)
          xtry[3] = xlo[3] + double(j3 - 1) * xstep[3];

          //      do 3 j4=1,7
          for (int j4 = 1; j4 <= 7; j4++) {
            //      xtry(4)=xlo(4)+DBLE(j4-1)*xstep(4)
            xtry[4] = xlo[4] + double(j4 - 1) * xstep[4];

            //      do 3 j5=1,7
            for (int j5 = 1; j5 <= 7; j5++) {
              //      xtry(5)=xlo(5)+DBLE(j5-1)*xstep(5)
              xtry[5] = xlo[5] + double(j5 - 1) * xstep[5];

              //      do 3 j6=1,7
              for (int j6 = 1; j6 <= 7; j6++) {
                //      xtry(6)=xlo(6)+DBLE(j6-1)*xstep(6)
                xtry[6] = xlo[6] + double(j6 - 1) * xstep[6];
                //      call f1(xtry,try)
                f1(xtry, tryy);

                //      if(try.gt.val) go to 3
                if (tryy > val) continue;

                //      val=try
                //      ix(1)=j1
                //      ix(2)=j2
                //      ix(3)=j3
                //      ix(4)=j4
                //      ix(5)=j5
                //      ix(6)=j6
                val = tryy;
                ix[1] = j1;
                ix[2] = j2;
                ix[3] = j3;
                ix[4] = j4;
                ix[5] = j5;
                ix[6] = j6;

                //      do 12 i=1,6
                //   12 x(i)=SNGL(xtry(i))
                for (int i = 1; i <= 6; i++)
                  x[i] = xtry[i];
              }
            }
          }
        }
      }
    }
    //    3 CONTINUE

    //      DO 4 I=1,6
    //      xhi(i)=xlo(i)+DBLE(ix(i)+1)*xstep(i)
    //    4 xlo(i)=xlo(i)+DBLE(ix(i)-3)*xstep(i)
    for (int i = 1; i <= 6; i++) {
      xhi[i] = xlo[i] + (ix[i] + 1) * xstep[i];
      xlo[i] = xlo[i] + (ix[i] - 3) * xstep[i];
    }
  }
  //    8 CONTINUE
  //      return
  //      end
}

//-----------------------------------------------------------------------------
void Taquart::UsmtCore::f1(double X[], double &fff) {
  //      SUBROUTINE F1(X,FFF)
  //      CHARACTER PS(80),TITLE*40
  //      REAL U(80),ARR(80),AZM(80),TKF(80)
  //      INTEGER RO(80),VEL(80),R(80)
  //      COMMON/MDATA/ PS,U,ARR,AZM,TKF,RO,VEL,R,TITLE,N,TROZ
  //      COMMON/PDATA/ A(80,6)
  //      DOUBLE PRECISION SUM,X(6),FFF

  //      FFF=DBLE(0.)
  fff = 0.0;

  //      DO 1 I=1,N
  for (int i = 1; i <= N; i++) {
    //      SUM=DBLE(0.)
    double sum = 0.0;
    //      DO 2 J=1,6
    //    2 SUM=SUM+DBLE(A(I,J))*X(J)
    for (int j = 1; j <= 6; j++) {
      sum = sum + A[i][j] * X[j];
    }
    //      FFF=FFF+DABS(SUM-DBLE(U(I)))
    fff = fff + fabs(sum - U[i]);
  }
  //    1 CONTINUE

  //      IF(DABS(FFF).GT.1.D+30) FFF=1.D+30
  if (fabs(fff) > 1e+30) fff = 1e+30;
  //      RETURN
  //      END
}

//-----------------------------------------------------------------------------
void Taquart::UsmtCore::GSOL5(double x[], int &IEXP) {
  //      subroutine gsol5(x,IEXP)
  //      dimension x(5),ix(5)
  //      double precision xlo(5),xhi(5),xstep(5),six,size,xtry(5),VAL,TRY
  //      DATA SIX,METH/6.D+0,2/
  double xlo[5 + 1], xhi[5 + 1], xstep[5 + 1], six = 6.0, xtry[5], VAL = 0.0,
      TRY = 0.0;
  //int METH = 2;
  int ix[5 + 1];

  //      IF((IEXP.LT.10).OR.(IEXP.GT.30)) IEXP=20
  //      iter=0
  //      jter=1
  //      val=1.d+30
  if (IEXP < 10 || IEXP > 30) IEXP = 20;
  int iter = 0;
  //int jter = 1;
  VAL = 1.0e+30;

  //      do 1 i=1,5
  for (int i = 1; i <= 5; i++) {
    //      xlo(i)=DBLE(-1.*10.**IEXP)
    //      xhi(i)=DBLE(10.**IEXP)
    xlo[i] = -1.0 * pow(10.0, double(IEXP));
    xhi[i] = pow(10.0, double(IEXP));
    ix[i] = 0;
  }
  //    1 ix(i)=0

  //      DO 8 L=1,50
  for (int l = 1; l <= 50; l++) {
    PROGRESS(l + 50, 350);

    //      do 2 i=1,5
    //    2 xstep(i)=(xhi(i)-xlo(i))/six
    for (int i = 1; i <= 5; i++)
      xstep[i] = (xhi[i] - xlo[i]) / six;

    //      size=xhi(1)-xlo(1)
    //      iter=iter+1
    //      do 3 j1=1,7
    //size = xhi[1] - xlo[1];
    iter = iter + 1;
    for (int j1 = 1; j1 <= 7; j1++) {
      //      xtry(1)=xlo(1)+DBLE(j1-1)*xstep(1)
      //      CALL POSTEP(METH,JTER,J1)
      xtry[1] = xlo[1] + double(j1 - 1) * xstep[1];
      //POSTEP(METH,jter,j1);

      //      do 3 j2=1,7
      for (int j2 = 1; j2 <= 7; j2++) {
        //      xtry(2)=xlo(2)+DBLE(j2-1)*xstep(2)
        //      do 3 j3=1,7
        xtry[2] = xlo[2] + double(j2 - 1) * xstep[2];
        for (int j3 = 1; j3 <= 7; j3++) {
          //      xtry(3)=xlo(3)+DBLE(j3-1)*xstep(3)
          //      do 3 j4=1,7
          xtry[3] = xlo[3] + double(j3 - 1) * xstep[3];
          for (int j4 = 1; j4 <= 7; j4++) {
            //      xtry(4)=xlo(4)+DBLE(j4-1)*xstep(4)
            //      do 3 j5=1,7
            xtry[4] = xlo[4] + double(j4 - 1) * xstep[4];
            for (int j5 = 1; j5 <= 7; j5++) {
              //      xtry(5)=xlo(5)+DBLE(j5-1)*xstep(5)
              //      call f2(xtry,try)
              xtry[5] = xlo[5] + double(j5 - 1) * xstep[5];
              f2(xtry, TRY);

              //      if(try.gt.val) go to 3
              if (TRY > VAL) continue;

              //      val=try
              //      ix(1)=j1
              //      ix(2)=j2
              //      ix(3)=j3
              //      ix(4)=j4
              //      ix(5)=j5
              //      do 12 i=1,5
              //   12 x(i)=SNGL(Xtry(i))
              VAL = TRY;
              ix[1] = j1;
              ix[2] = j2;
              ix[3] = j3;
              ix[4] = j4;
              ix[5] = j5;
              for (int i = 1; i <= 5; i++) {
                x[i] = xtry[i];
              }
            }
          }
        }
      }
    }
    //    3 CONTINUE
    //      DO 4 I=1,5
    //      xhi(i)=xlo(i)+DBLE(ix(i)+1)*xstep(i)
    //    4 xlo(i)=xlo(i)+DBLE(ix(i)-3)*xstep(i)
    for (int i = 1; i <= 5; i++) {
      xhi[i] = xlo[i] + double(ix[i] + 1) * xstep[i];
      xlo[i] = xlo[i] + double(ix[i] - 3) * xstep[i];
    }
  }
  //    8 CONTINUE
  //      return
  //      end
}

//-----------------------------------------------------------------------------
void Taquart::UsmtCore::GSOLA(double x[], int &IEXP) {
  //      subroutine gsola(x,IEXP)
  //      double precision xlo(4),xhi(4),xstep(4),xtry(5),FOUR,SIZE,SIX,val,try,help,y(5),del,two,ZERO
  //      dimension x(5),ix(4),xmem(5,5),vmem(5)
  //      DATA TWO,FOUR,METH,SIX/2.D+0,4.D+0,3,6.D+0/

  double xlo[4 + 1];
  Zero(xlo, 5);
  double xhi[4 + 1];
  Zero(xhi, 5);
  double xstep[4 + 1];
  Zero(xstep, 5);
  double xtry[5 + 1];
  Zero(xtry, 6);
  double FOUR = 4.0, SIX = 6.0, val = 0.0, tryy = 0.0;
  double help, y[5 + 1], DEL = 0.0, TWO = 2.0;
  int ix[4 + 1];
  double xmem[5 + 1][5 + 1], vmem[5 + 1];
  Zero(&xmem[0][0], 36);
  Zero(vmem, 6);
  //int METH = 3;

  //      IF((IEXP.LT.10).OR.(IEXP.GT.30)) IEXP=20
  //      ZERO=dble(0.)
  //      iter=0
  //      jter=1
  if (IEXP < 10 || IEXP > 30) IEXP = 20;
  double ZERO = 0.0;
  int iter = 0;
  //int jter = 1;
  //      do 5 i=1,5
  //      vmem(i)=1.e+30
  //      do 5 j=1,5
  //    5 xmem(i,j)=0.
  for (int i = 1; i <= 5; i++) {
    vmem[i] = 1.0e+30;
    for (int j = 1; j <= 5; j++)
      xmem[i][j] = 0.0;
  }

  //C     Search for X1,X2,X3,X4; X5 is calculated:
  //      val=1.d+30
  //      do 1 i=1,4
  //      xlo(i)=DBLE(-1.*10.**IEXP)
  //      xhi(i)=DBLE(10.**IEXP)
  //    1 ix(i)=0
  val = 1e+30;
  for (int i = 1; i <= 4; i++) {
    xlo[i] = -1.0 * pow(10.0, IEXP);
    xhi[i] = pow(10.0, IEXP);
    ix[i] = 0;
  }

  //int j7 = 0;

  //      DO 8 L=1,50
  for (int l = 1; l <= 50; l++) {
    PROGRESS(l + 100, 350);
    //      do 2 i=1,4
    //    2 xstep(i)=(xhi(i)-xlo(i))/SIX
    for (int i = 1; i <= 4; i++)
      xstep[i] = (xhi[i] - xlo[i]) / SIX;

    //      size=xhi(1)-xlo(1)
    //      iter=iter+1
    //      xtry(5)=DBLE(0.)
    //      J7=7
    //SIZE = xhi[1] - xlo[1];
    iter++;
    xtry[5] = 0.0;
    //j7 = 7;

    //      CALL POSTEP(METH,JTER,J7)
    //POSTEP(METH,jter,j7);

    //      do 3 j1=1,7
    for (int j1 = 1; j1 <= 7; j1++) {
      //      xtry(1)=xlo(1)+DBLE(j1-1)*xstep(1)
      xtry[1] = xlo[1] + double(j1 - 1) * xstep[1];
      //      do 3 j2=1,7
      for (int j2 = 1; j2 <= 7; j2++) {
        //      xtry(2)=xlo(2)+DBLE(j2-1)*xstep(2)
        xtry[2] = xlo[2] + double(j2 - 1) * xstep[2];
        //      do 3 j3=1,7
        for (int j3 = 1; j3 <= 7; j3++) {
          //      xtry(3)=xlo(3)+DBLE(j3-1)*xstep(3)
          xtry[3] = xlo[3] + double(j3 - 1) * xstep[3];
          //      do 3 j4=1,7
          for (int j4 = 1; j4 <= 7; j4++) {
            //      xtry(4)=xlo(4)+DBLE(j4-1)*xstep(4)
            xtry[4] = xlo[4] + double(j4 - 1) * xstep[4];

            //      if(dabs(xtry(1)).lt.1.d-6) go to 23
            if (fabs(xtry[1]) < 1.0e-6) goto p23;

            //      do 21 i=1,5
            //   21 y(i)=xtry(i)*1.d-10
            for (int i = 1; i <= 5; i++)
              y[i] = xtry[i] * 1.0e-10;

            //      DEL=(TWO*y(2)*y(3))**2+FOUR*Y(1)*(Y(4)*(-Y(1)**2-Y(1)*Y(4)-
            //     $Y(3)**2+Y(2)**2)+Y(2)**2*Y(1))
            DEL = pow(TWO * y[2] * y[3], 2.0)
                + FOUR * y[1]
                    * (y[4]
                        * (-pow(y[1], 2.0) - y[1] * y[4] - pow(y[3], 2.0)
                            + pow(y[2], 2.0)) + pow(y[2], 2.0) * y[1]);

            //      IF(DEL.LT.ZERO) GO TO 3
            if (DEL < ZERO) continue;

            //      DEL=DSQRT(DEL)
            //      Y(5)=(TWO*Y(2)*Y(3)+DEL)/TWO/Y(1)
            //      XTRY(5)=Y(5)*1.D+10
            //      call f2(xtry,try)
            DEL = sqrt(DEL);
            y[5] = (TWO * y[2] * y[3] + DEL) / TWO / y[1];
            xtry[5] = y[5] * 1.0e+10;
            f2(xtry, tryy);

            //      IF(TRY.GT.VAL) GO TO 22
            //      CALL RENUM(TRY,VAL,IX,J1,J2,J3,J4)
            if (tryy > val) goto p22;
            RENUM(tryy, val, ix, j1, j2, j3, j4);

            //      DO 12 i=1,5
            //   12 x(i)=SNGL(XTRY(I))
            for (int i = 1; i <= 5; i++)
              x[i] = double(xtry[i]);

            //   22 Y(5)=(TWO*Y(2)*Y(3)-DEL)/TWO/Y(1)
            //      XTRY(5)=Y(5)*1.D+10
            //      call f2(xtry,try)
            p22: y[5] = (TWO * y[2] * y[3] - DEL) / TWO / y[1];
            xtry[5] = y[5] * 1.0e+10;
            f2(xtry, tryy);

            //      IF(TRY.GT.VAL) GO TO 3
            //      CALL RENUM(TRY,VAL,IX,J1,J2,J3,J4)
            if (tryy > val) continue;
            RENUM(tryy, val, ix, j1, j2, j3, j4);

            //      DO 25 i=1,5
            //   25 x(i)=SNGL(xtry(i))
            for (int i = 1; i <= 5; i++)
              x[i] = double(xtry[i]);
            //      go to 3
            continue;

            //   23 if((DABS(XTRY(2)).lt.1.d-6).OR.(DABS(XTRY(3)).LT.1.d-6)) GO TO 3
            p23: if (fabs(xtry[2]) < 1.0e-6 || fabs(xtry[3]) < 1.0e-6) continue;

            //      DO 27 i=1,5
            //   27 y(i)=xtry(i)*1.d-10
            for (int i = 1; i <= 5; i++)
              y[i] = xtry[i] * 1.0e-10;

            //      DEL=Y(4)*(-Y(1)**2-Y(1)*Y(4)-Y(3)**2+Y(2)**2)+Y(2)**2*Y(1)

            DEL = y[4]
                * (-pow(y[1], 2.0) - y[1] * y[4] - pow(y[3], 2.0)
                    + pow(y[2], 2.0)) + pow(y[2], 2.0) * y[1];

            //      Y(5)=-DEL/TWO/Y(3)/Y(2)
            //      XTRY(5)=Y(5)*1.d+10
            //      call f2(xtry,try)
            y[5] = -DEL / TWO / y[3] / y[2];
            xtry[5] = y[5] * 1.0e+10;
            f2(xtry, tryy);

            //      if(try.gt.val) go to 3
            //      CALL RENUM(TRY,VAL,IX,J1,J2,J3,J4)
            //      do 28 i=1,5
            //   28 x(i)=SNGL(xtry(i))
            if (tryy > val) continue;
            RENUM(tryy, val, ix, j1, j2, j3, j4);
            for (int i = 1; i <= 5; i++)
              x[i] = double(xtry[i]);
          }
        }
      }
    }
    //    3 continue

    //      if(val.eq.1.d+30) go to 30
    if (val == 1e+30) goto p30;
    //      DO 4 I=1,4
    //      xhi(i)=xlo(i)+DBLE(ix(i)+1)*xstep(i)
    //    4 xlo(i)=xlo(i)+DBLE(ix(i)-3)*xstep(i)
    for (int i = 1; i <= 4; i++) {
      xhi[i] = xlo[i] + double(ix[i] + 1) * xstep[i];
      xlo[i] = xlo[i] + double(ix[i] - 3) * xstep[i];
    }
  }
  //    8 CONTINUE

  //      DO 29 I=1,5
  //   29 XMEM(1,I)=X(I)
  //      VMEM(1)=SNGL(VAL)
  for (int i = 1; i <= 5; i++)
    xmem[1][i] = x[i];
  vmem[1] = val;

  //C     Search for X1,X2,X3,X5; X4 is calculated:
  //   30 val=1.d+30
  //      do 101 i=1,4
  //      xlo(i)=DBLE(-1.*10.**IEXP)
  //      xhi(i)=DBLE(10.**IEXP)
  //  101 ix(i)=0
  p30: val = 1.0e+30;
  for (int i = 1; i <= 4; i++) {
    xlo[i] = -1.0 * pow(10, IEXP);
    xhi[i] = pow(10, IEXP);
    ix[i] = 0;
  }

  //      DO 108 l=1,50
  for (int l = 1; l <= 50; l++) {
    PROGRESS(l + 150, 350);
    //      do 102 i=1,4
    //  102 xstep(i)=(xhi(i)-xlo(i))/SIX
    for (int i = 1; i <= 4; i++)
      xstep[i] = (xhi[i] - xlo[i]) / SIX;

    //      size=xhi(1)-xlo(1)
    //      iter=iter+1
    //      xtry(4)=DBLE(0.)
    //      CALL POSTEP(METH,JTER,J7)
    //SIZE = xhi[1] - xlo[1];
    iter++;
    xtry[4] = 0.0;
    //POSTEP(METH,jter,j7);

    //      do 103 j1=1,7
    for (int j1 = 1; j1 <= 7; j1++) {
      //      xtry(1)=xlo(1)+DBLE(j1-1)*xstep(1)
      xtry[1] = xlo[1] + double(j1 - 1) * xstep[1];
      //      do 103 j2=1,7
      for (int j2 = 1; j2 <= 7; j2++) {
        //      xtry(2)=xlo(2)+DBLE(j2-1)*xstep(2)
        xtry[2] = xlo[2] + double(j2 - 1) * xstep[2];
        //      do 103 j3=1,7
        for (int j3 = 1; j3 <= 7; j3++) {
          //      xtry(3)=xlo(3)+DBLE(j3-1)*xstep(3)
          xtry[3] = xlo[3] + double(j3 - 1) * xstep[3];
          //      do 103 j4=1,7
          for (int j4 = 1; j4 <= 7; j4++) {
            //      xtry(5)=xlo(4)+DBLE(j4-1)*xstep(4)  
            //      if(dabs(xtry(1)).lt.1.d-6) go to 123
            xtry[5] = xlo[4] + double(j4 - 1) * xstep[4];
            if (fabs(xtry[1]) < 1.0e-06) goto p123;

            //      do 121 i=1,5
            //  121 y(i)=xtry(i)*1.d-10
            for (int i = 1; i <= 5; i++)
              y[i] = xtry[i] * 1.0e-10;

            //      DEL=(-y(1)**2-y(3)**2+y(2)**2)**2+FOUR*Y(1)*(TWO*Y(2)*Y(3)*
            //     $Y(5)-Y(1)*Y(5)**2+Y(1)*Y(2)**2)
            DEL = pow(-pow(y[1], 2.0) - pow(y[3], 2.0) + pow(y[2], 2.0), 2.0)
                + FOUR * y[1]
                    * (TWO * y[2] * y[3] * y[5] - y[1] * pow(y[5], 2.0)
                        + y[1] * pow(y[2], 2.0));

            //      IF(DEL.LT.ZERO) GO TO 103
            if (DEL < ZERO) continue;

            //      DEL=DSQRT(DEL)
            //      Y(4)=(-Y(1)**2-Y(3)**2+Y(2)**2+DEL)/TWO/Y(1)
            //      XTRY(4)=Y(4)*1.D+10
            //      call f2(xtry,try)
            DEL = sqrt(DEL);
            y[4] = (-pow(y[1], 2.0) - pow(y[3], 2.0) + pow(y[2], 2.0) + DEL)
                / TWO / y[1];
            xtry[4] = y[4] * 1.0e+10;
            f2(xtry, tryy);

            //      IF(TRY.GT.VAL) GO TO 122
            //      CALL RENUM(TRY,VAL,IX,J1,J2,J3,J4)
            //      DO 112 i=1,5
            //  112 x(i)=SNGL(XTRY(I))
            if (tryy <= val) {
              RENUM(tryy, val, ix, j1, j2, j3, j4);
              for (int i = 1; i <= 5; i++)
                x[i] = xtry[i];
            }

            //  122 Y(4)=(-Y(1)**2-Y(3)**2+Y(2)**2-DEL)/TWO/Y(1)
            y[4] = (-pow(y[1], 2.0) - pow(y[3], 2.0) + pow(y[2], 2.0) - DEL)
                / TWO / y[1];

            //      XTRY(4)=Y(4)*1.D+10
            //      call f2(xtry,try)
            //      IF(TRY.GT.VAL) GO TO 103
            //      CALL RENUM(TRY,VAL,IX,J1,J2,J3,J4)
            xtry[4] = y[4] * 1.0e+10;
            f2(xtry, tryy);
            if (tryy > val) continue;
            RENUM(tryy, val, ix, j1, j2, j3, j4);

            //      DO 125 i=1,5
            //  125 x(i)=SNGL(xtry(i))
            //      go to 103
            for (int i = 1; i <= 5; i++)
              x[i] = xtry[i];

            continue;

            //  123 DO 127 i=1,5
            //  127 y(i)=xtry(i)*1.d-10
            p123: for (int i = 1; i <= 5; i++)
              y[i] = xtry[i] * 1.0e-10;

            //      HELP=-y(1)**2-y(3)**2+y(2)**2
            //      IF(DABS(HELP).LT.1.d-20) GO TO 103
            //      DEL=TWO*y(2)*y(3)*y(5)-y(5)**2*y(1)+y(2)**2*y(1)
            help = -pow(y[1], 2.0) - pow(y[3], 2.0) + pow(y[2], 2.0);
            if (fabs(help) < 1.0e-20) continue;
            DEL = TWO * y[2] * y[3] * y[5] - y[5] * y[5] * y[1]
                + y[2] * y[2] * y[1];

            //      Y(4)=-DEL/HELP
            //      XTRY(4)=Y(4)*1.d+10
            //      call f2(xtry,try)
            y[4] = -DEL / help;
            xtry[4] = y[4] * 1.0e+10;
            f2(xtry, tryy);

            //      if(try.gt.val) go to 103
            //      CALL RENUM(TRY,VAL,IX,J1,J2,J3,J4)
            //      do 128 i=1,5
            //  128 x(i)=SNGL(xtry(i))
            if (tryy <= val) {
              RENUM(tryy, val, ix, j1, j2, j3, j4);
              for (int i = 1; i <= 5; i++)
                x[i] = xtry[i];
            }
          }
        }
      }
    }
    //  103 continue
    //      if(val.eq.1.d+30) go to 130
    //      DO 104 I=1,4
    //      xhi(i)=xlo(i)+DBLE(ix(i)+1)*xstep(i)
    //  104 xlo(i)=xlo(i)+DBLE(ix(i)-3)*xstep(i)
    if (val == 1.0e+30) goto p130;
    for (int i = 1; i <= 4; i++) {
      xhi[i] = xlo[i] + double(ix[i] + 1) * xstep[i];
      xlo[i] = xlo[i] + double(ix[i] - 3) * xstep[i];
    }
  }
  //  108 CONTINUE
  //      DO 129 I=1,5
  //  129 XMEM(2,I)=X(I)
  //      VMEM(2)=SNGL(VAL)
  for (int i = 1; i <= 5; i++)
    xmem[2][i] = x[i];
  vmem[2] = val;

  //C     Search for X1,X2,X4,X5; X3 is calculated:
  //  130 val=1.d+30
  p130: val = 1.0e+30;

  //      do 201 i=1,4
  //      xlo(i)=DBLE(-1.*10.**IEXP)
  //      xhi(i)=DBLE(10.**IEXP)
  //  201 ix(i)=0
  for (int i = 1; i <= 4; i++) {
    xlo[i] = -1.0 * pow(10, IEXP);
    xhi[i] = pow(10, IEXP);
    ix[i] = 0;
  }

  //      DO 208 l=1,50
  for (int l = 1; l <= 50; l++) {
    PROGRESS(l + 200, 350);
    //      do 202 i=1,4
    //  202 xstep(i)=(xhi(i)-xlo(i))/SIX
    for (int i = 1; i <= 4; i++)
      xstep[i] = (xhi[i] - xlo[i]) / SIX;

    //      size=xhi(1)-xlo(1)
    //      iter=iter+1
    //      xtry(3)=ZERO
    //      CALL POSTEP(METH,JTER,J7)
    //SIZE = xhi[1] - xlo[1];
    iter++;
    xtry[3] = 0.0;
    //POSTEP(METH,jter,j7);

    //      do 203 j1=1,7
    for (int j1 = 1; j1 <= 7; j1++) {
      //      xtry(1)=xlo(1)+DBLE(j1-1)*xstep(1)
      xtry[1] = xlo[1] + double(j1 - 1) * xstep[1];
      //      do 203 j2=1,7
      for (int j2 = 1; j2 <= 7; j2++) {
        //      xtry(2)=xlo(2)+DBLE(j2-1)*xstep(2)
        xtry[2] = xlo[2] + double(j2 - 1) * xstep[2];
        //      do 203 j3=1,7
        for (int j3 = 1; j3 <= 7; j3++) {
          //      xtry(4)=xlo(3)+DBLE(j3-1)*xstep(3)
          xtry[4] = xlo[3] + double(j3 - 1) * xstep[3];
          //      do 203 j4=1,7
          for (int j4 = 1; j4 <= 7; j4++) {
            //      xtry(5)=xlo(4)+DBLE(j4-1)*xstep(4)
            xtry[5] = xlo[4] + double(j4 - 1) * xstep[4];

            //      if(abs(xtry(4)).lt.1.d-6) go to 223
            if (fabs(xtry[4]) < 1.0e-06) goto p223;

            //      do 221 i=1,5
            //  221 y(i)=xtry(i)*1.d-10
            for (int i = 1; i <= 5; i++)
              y[i] = xtry[i] * 1.0e-10;

            //      DEL=(TWO*Y(2)*Y(5))**2+FOUR*Y(4)*(-y(1)**2*y(4)-Y(1)*Y(4)**2-
            //     $Y(5)**2*Y(1)+Y(2)**2*Y(1)+Y(2)**2*Y(4))
            DEL = pow(TWO * y[2] * y[5], 2.0)
                + FOUR * y[4]
                    * (-y[1] * y[1] * y[4] - y[1] * y[4] * y[4]
                        - y[5] * y[5] * y[1] + y[2] * y[2] * y[1]
                        + y[2] * y[2] * y[4]);

            //      IF(DEL.LT.ZERO) GO TO 203
            if (DEL < ZERO) continue;

            //      DEL=dSQRT(DEL)
            //      Y(3)=(TWO*Y(2)*Y(5)+DEL)/TWO/Y(4)
            //      XTRY(3)=Y(3)*1.D+10
            //      call f2(xtry,try)
            DEL = sqrt(DEL);
            y[3] = (TWO * y[2] * y[5] + DEL) / TWO / y[4];
            xtry[3] = y[3] * 1.0e+10;
            f2(xtry, tryy);

            //      IF(TRY.GT.VAL) GO TO 222
            //      CALL RENUM(TRY,VAL,IX,J1,J2,J3,J4)
            //      DO 212 i=1,5
            //  212 x(i)=SNGL(XTRY(I))
            if (tryy <= val) {
              RENUM(tryy, val, ix, j1, j2, j3, j4);
              for (int i = 1; i <= 5; i++)
                x[i] = xtry[i];
            }

            //  222 Y(3)=(TWO*Y(2)*Y(5)-DEL)/TWO/Y(4)
            y[3] = (TWO * y[2] * y[5] - DEL) / TWO / y[4];

            //      XTRY(3)=Y(3)*1.d+10
            //      call f2(xtry,try)
            //      IF(TRY.GT.VAL) GO TO 203
            //      CALL RENUM(TRY,VAL,IX,J1,J2,J3,J4)
            xtry[3] = y[3] * 1.0e+10;
            f2(xtry, tryy);
            if (tryy > val) continue;
            RENUM(tryy, val, ix, j1, j2, j3, j4);

            //      DO 225 i=1,5
            //  225 x(i)=sngl(xtry(i))
            //      go to 203
            for (int i = 1; i <= 5; i++)
              x[i] = xtry[i];
            continue;

            //  223 IF((dABS(XTRY(2)).LT.1.d-6).OR.(dABS(XTRY(5)).LT.1.d-6))
            //     $ GO TO 203
            p223: if (fabs(xtry[2]) < 1.0e-6 || fabs(xtry[5]) < 1.0e-06)
              continue;

            //      DO 227 i=1,5
            //  227 y(i)=xtry(i)*1.d-10
            for (int i = 1; i <= 5; i++)
              y[i] = xtry[i] * 1.0e-10;

            //      DEL=-y(1)**2*y(4)-Y(1)*y(4)**2-Y(5)**2*Y(1)+Y(2)**2*(Y(1)+Y(4))
            DEL = -y[1] * y[1] * y[4] - y[1] * y[4] * y[4] - y[5] * y[5] * y[1]
                + y[2] * y[2] * (y[1] + y[4]);

            //      Y(3)=-DEL/TWO/Y(2)/Y(5)
            //      XTRY(3)=Y(3)*1.d+10
            //      call f2(xtry,try)
            y[3] = -DEL / TWO / y[2] / y[5];
            xtry[3] = y[3] * 1.0e+10;
            f2(xtry, tryy);

            //      if(try.gt.val) go to 203
            //      CALL RENUM(TRY,VAL,IX,J1,J2,J3,J4)
            //      do 228 i=1,5
            //  228 x(i)=SNGL(xtry(i))
            if (tryy <= val) {
              RENUM(tryy, val, ix, j1, j2, j3, j4);
              for (int i = 1; i <= 5; i++)
                x[i] = double(xtry[i]);
            }
          }
        }
      }
    }
    //  203 continue
    //      if(val.eq.1.d+30) go to 230
    //      DO 204 I=1,4
    //      xhi(i)=xlo(i)+DBLE(ix(i)+1)*xstep(i)
    //  204 xlo(i)=xlo(i)+DBLE(ix(i)-3)*xstep(i)
    if (val == 1.0e+30) goto p230;
    for (int i = 1; i <= 4; i++) {
      xhi[i] = xlo[i] + double(ix[i] + 1) * xstep[i];
      xlo[i] = xlo[i] + double(ix[i] - 3) * xstep[i];
    }
  }
  //  208 CONTINUE

  //      DO 229 I=1,5
  //  229 XMEM(3,I)=X(I)
  //      VMEM(3)=SNGL(VAL)
  for (int i = 1; i <= 5; i++)
    xmem[3][i] = x[i];
  vmem[3] = val;

  //C     Search for X1,X3,X4,X5; X2 is calculated:
  //  230 val=1.d+30
  //      do 301 i=1,4
  //      xlo(i)=DBLE(-1.*10.**IEXP)
  //      xhi(i)=DBLE(10.**IEXP)
  //  301 ix(i)=0
  p230: val = 1.0e+30;
  for (int i = 1; i <= 4; i++) {
    xlo[i] = -1.0 * pow(10, IEXP);
    xhi[i] = pow(10, IEXP);
    ix[i] = 0;
  }

  //      DO 308 l=1,50
  for (int l = 1; l <= 50; l++) {
    PROGRESS(l + 250, 350);
    //      do 302 i=1,4
    //  302 xstep(i)=(xhi(i)-xlo(i))/SIX
    for (int i = 1; i <= 4; i++)
      xstep[i] = (xhi[i] - xlo[i]) / SIX;

    //      size=xhi(1)-xlo(1)
    //      iter=iter+1
    //      xtry(2)=ZERO
    //      CALL POSTEP(METH,JTER,J7)
    //SIZE = xhi[1] - xlo[1];
    iter++;
    xtry[2] = 0.0;
    //POSTEP(METH,jter,j7);

    //      do 303 j1=1,7
    for (int j1 = 1; j1 <= 7; j1++) {
      //      xtry(1)=xlo(1)+DBLE(j1-1)*xstep(1)
      xtry[1] = xlo[1] + double(j1 - 1) * xstep[1];
      //      do 103 j2=1,7
      for (int j2 = 1; j2 <= 7; j2++) {
        //      xtry(3)=xlo(2)+DBLE(j2-1)*xstep(2)
        xtry[3] = xlo[2] + double(j2 - 1) * xstep[2];
        //      do 303 j3=1,7
        for (int j3 = 1; j3 <= 7; j3++) {
          //      xtry(4)=xlo(3)+DBLE(j3-1)*xstep(3)
          xtry[4] = xlo[3] + double(j3 - 1) * xstep[3];
          //      do 303 j4=1,7
          for (int j4 = 1; j4 <= 7; j4++) {
            //      xtry(5)=xlo(4)+DBLE(j4-1)*xstep(4)
            xtry[5] = xlo[4] + double(j4 - 1) * xstep[4];
            //      if(dabs(xtry(4)+XTRY(1)).lt.1.d-6) go to 323
            if (fabs(xtry[4] + xtry[1]) < 1.0e-06) goto p323;

            //      do 321 i=1,5
            //  321 y(i)=xtry(i)*1.d-10
            for (int i = 1; i <= 5; i++)
              y[i] = xtry[i] * 1.0e-10;

            //      DEL=(TWO*Y(3)*Y(5))**2-FOUR*(Y(1)+Y(4))*(-y(1)**2*y(4)-Y(1)*
            //     $Y(4)**2-Y(5)**2*Y(1)-y(3)**2*Y(4))
            DEL = pow(TWO * y[3] * y[5], 2.0)
                - FOUR * (y[1] + y[4])
                    * (-y[1] * y[1] * y[4] - y[1] * y[4] * y[4]
                        - y[5] * y[5] * y[1] - y[3] * y[3] * y[4]);

            //      IF(DEL.LT.ZERO) GO TO 303
            if (DEL < ZERO) continue;

            //      DEL=DSQRT(DEL)
            //      Y(2)=(-TWO*Y(3)*Y(5)-DEL)/TWO/(Y(1)+Y(4))
            //      XTRY(2)=Y(2)*1.D+10
            //      call f2(xtry,try)
            DEL = sqrt(DEL);
            y[2] = (-TWO * y[3] * y[5] - DEL) / TWO / (y[1] + y[4]);
            xtry[2] = y[2] * 1.0e+10;
            f2(xtry, tryy);

            //      IF(TRY.GT.VAL) GO TO 322
            //      CALL RENUM(TRY,VAL,IX,J1,J2,J3,J4)
            //      DO 312 i=1,5
            //  312 x(i)=SNGL(XTRY(I))
            if (tryy <= val) /* ONE : poprawiony b��d 2006.10.05 */
            {
              RENUM(tryy, val, ix, j1, j2, j3, j4);
              for (int i = 1; i <= 5; i++)
                x[i] = xtry[i];
            }

            //  322 Y(2)=(-TWO*Y(3)*Y(5)+DEL)/TWO/(Y(1)+Y(4))
            //      XTRY(2)=Y(2)*1.D+10
            //      call f2(xtry,try)
            y[2] = (-TWO * y[3] * y[5] + DEL) / TWO / (y[1] + y[4]);
            xtry[2] = y[2] * 1.0e+10;
            f2(xtry, tryy);

            //      IF(TRY.GT.VAL) GO TO 303
            //      CALL RENUM(TRY,VAL,IX,J1,J2,J3,J4)
            if (tryy > val) continue;
            RENUM(tryy, val, ix, j1, j2, j3, j4);

            //      DO 325 i=1,5
            //  325 x(i)=SNGL(xtry(i))
            //      go to 303
            for (int i = 1; i <= 5; i++)
              x[i] = xtry[i];
            continue;

            //  323 IF((DABS(XTRY(3)).LT.1.d-6).OR.(DABS(XTRY(5)).LT.1.d-6))
            //     $ GO TO 303
            p323: if (fabs(xtry[3]) < 1.0e-06 || fabs(xtry[5]) < 1.0e-06)
              continue;

            //      DO 327 i=1,5
            //  327 y(i)=xtry(i)*1.d-10
            for (int i = 1; i <= 5; i++)
              y[i] = xtry[i] * 1.0e-10;

            //      DEL=-y(1)**2*y(4)-Y(1)*y(4)**2-Y(5)**2*Y(1)-Y(3)**2*Y(4)
            DEL = -y[1] * y[1] * y[4] - y[1] * y[4] * y[4] - y[5] * y[5] * y[1]
                - y[3] * y[3] * y[4];

            //      Y(2)=-DEL/TWO/Y(3)/Y(5)
            //      XTRY(2)=Y(2)*1.d+10
            //      call f2(xtry,try)
            y[2] = -DEL / TWO / y[3] / y[5];
            xtry[2] = y[2] * 1.0e+10;
            f2(xtry, tryy);

            //      if(try.gt.val) go to 303
            //      CALL RENUM(TRY,VAL,IX,J1,J2,J3,J4)
            //      do 328 i=1,5
            //  328 x(i)=SNGL(xtry(i))
            if (tryy <= val) {
              RENUM(tryy, val, ix, j1, j2, j3, j4);
              for (int i = 1; i <= 5; i++)
                x[i] = xtry[i];
            }
          }
        }
      }
    }
    //  303 continue
    //      if(val.eq.1.d+30) go to 330
    //      DO 304 I=1,4
    //      xhi(i)=xlo(i)+DBLE(ix(i)+1)*xstep(i)
    //  304 xlo(i)=xlo(i)+DBLE(ix(i)-3)*xstep(i)
    if (val == 1.0e+30) goto p330;
    for (int i = 1; i <= 4; i++) {
      xhi[i] = xlo[i] + double(ix[i] + 1) * xstep[i];
      xlo[i] = xlo[i] + double(ix[i] - 3) * xstep[i];
    }
  }
  //  308 CONTINUE
  //      DO 329 I=1,5
  //  329 XMEM(4,I)=X(I)
  //      VMEM(4)=SNGL(VAL)
  for (int i = 1; i <= 5; i++)
    xmem[4][i] = x[i];
  vmem[4] = val;

  //C     Search for X2,X2,X3,X5; X1 is calculated:
  //  330 val=1.d+30
  p330: val = 1.0e+30;
  //      do 401 i=1,4
  //      xlo(i)=DBLE(-1.*10.**IEXP)
  //      xhi(i)=DBLE(10.**IEXP)
  //  401 ix(i)=0
  for (int i = 1; i <= 4; i++) {
    xlo[i] = -1.0 * pow(10, IEXP);
    xhi[i] = pow(10, IEXP);
    ix[i] = 0;
  }

  //      DO 408 l=1,50
  for (int l = 1; l <= 50; l++) {
    PROGRESS(l + 300, 350);
    //      do 402 i=1,4
    //  402 xstep(i)=(xhi(i)-xlo(i))/SIX
    for (int i = 1; i <= 4; i++)
      xstep[i] = (xhi[i] - xlo[i]) / SIX;

    //      size=xhi(1)-xlo(1)
    //      iter=iter+1
    //      xtry(1)=ZERO
    //      CALL POSTEP(METH,JTER,J7)
    //SIZE = xhi[1] - xlo[1];
    iter++;
    xtry[1] = 0.0; /* ONE : poprawiony b��d 2006.10.05 */
    //POSTEP(METH,jter,j7);

    //      do 403 j1=1,7
    for (int j1 = 1; j1 <= 7; j1++) {
      //      xtry(2)=xlo(1)+DBLE(j1-1)*xstep(1)
      xtry[2] = xlo[1] + double(j1 - 1) * xstep[1];
      //      do 403 j2=1,7
      for (int j2 = 1; j2 <= 7; j2++) {
        //      xtry(2)=xlo(2)+DBLE(j2-1)*xstep(2)
        xtry[2] = xlo[2] + double(j2 - 1) * xstep[2];
        //      do 403 j3=1,7
        for (int j3 = 1; j3 <= 7; j3++) {
          //      xtry(3)=xlo(3)+DBLE(j3-1)*xstep(3)
          xtry[3] = xlo[3] + double(j3 - 1) * xstep[3];
          //      do 403 j4=1,7
          for (int j4 = 1; j4 <= 7; j4++) {
            //      xtry(5)=xlo(4)+DBLE(j4-1)*xstep(4)
            xtry[5] = xlo[4] + double(j4 - 1) * xstep[4];

            //      if(dabs(xtry(4)).lt.1.d-6) go to 423
            if (fabs(xtry[4]) < 1.0e-06) goto p423;

            //      do 421 i=1,5
            //  421 y(i)=xtry(i)*1.d-10
            for (int i = 1; i <= 5; i++)
              y[i] = xtry[i] * 1.0e-10;

            //      DEL=(-y(4)**2-y(5)**2+y(2)**2)**2+FOUR*Y(4)*(TWO*Y(2)*Y(3)*Y(5)
            //     $-Y(4)*Y(3)**2+Y(4)*Y(2)**2)
            DEL = pow(-y[4] * y[4] - y[5] * y[5] + y[2] * y[2], 2.0)
                + FOUR * y[4]
                    * (TWO * y[2] * y[3] * y[5] - y[4] * y[3] * y[3]
                        + y[4] * y[2] * y[2]);

            //      IF(DEL.LT.ZERO) GO TO 403
            if (DEL < ZERO) continue;

            //      DEL=DSQRT(DEL)
            //      Y(1)=(-Y(4)**2-Y(5)**2+Y(2)**2+DEL)/TWO/Y(4)
            //      XTRY(1)=Y(1)*1.D+10
            //      call f2(xtry,try)
            DEL = sqrt(DEL);
            y[1] = (-y[4] * y[4] - y[5] * y[5] + y[2] * y[2] + DEL) / TWO
                / y[4];
            xtry[1] = y[1] * 1.0e+10;
            f2(xtry, tryy);

            //      IF(TRY.GT.VAL) GO TO 422
            //      CALL RENUM(TRY,VAL,IX,J1,J2,J3,J4)
            //      DO 412 i=1,5
            //  412 x(i)=SNGL(XTRY(I))
            if (tryy <= val) {
              RENUM(tryy, val, ix, j1, j2, j3, j4);
              for (int i = 1; i <= 5; i++)
                x[i] = xtry[i];
            }

            //  422 Y(1)=(-Y(4)**2-Y(5)**2+Y(2)**2-DEL)/TWO/Y(4)
            y[1] = (-y[4] * y[4] - y[5] * y[5] + y[2] * y[2] - DEL) / TWO
                / y[4];

            //      XTRY(1)=Y(1)*1.D+10
            //      call f2(xtry,try)
            //      IF(TRY.GT.VAL) GO TO 403
            //      CALL RENUM(TRY,VAL,IX,J1,J2,J3,J4)
            xtry[1] = y[1] * 1.0e+10;
            f2(xtry, tryy);
            if (tryy > val) continue;
            RENUM(tryy, val, ix, j1, j2, j3, j4);

            //      DO 425 i=1,5
            //  425 x(i)=SNGL(xtry(i))
            //      go to 403
            for (int i = 1; i <= 5; i++)
              x[i] = xtry[i];
            continue;
            //  423 DO 427 i=1,5
            //  427 y(i)=xtry(i)*1.d-10
            p423: for (int i = 1; i <= 5; i++)
              y[i] = xtry[i] * 1.0e-10;

            //      HELP=-y(4)**2-y(5)**2+y(2)**2
            //      IF(DABS(HELP).LT.1.d-20) GO TO 403
            //      DEL=TWO*y(2)*y(3)*y(5)-y(3)**2*y(4)+y(2)**2*y(4)
            help = -pow(y[4], 2.0) - pow(y[5], 2.0) + pow(y[2], 2.0);
            if (fabs(help) < 1.0e-20) continue;
            DEL = TWO * y[2] * y[3] * y[5] - y[3] * y[3] * y[4]
                + y[2] * y[2] * y[4];

            //      Y(1)=-DEL/HELP
            //      XTRY(1)=Y(1)*1.d+10
            //      call f2(xtry,try)
            y[1] = -DEL / help;
            xtry[1] = y[1] * 1.0e+10;
            f2(xtry, tryy);
            //      if(try.gt.val) go to 403
            //      CALL RENUM(TRY,VAL,IX,J1,J2,J3,J4)
            //      do 428 i=1,5
            //  428 x(i)=SNGL(xtry(i))
            if (tryy <= val) {
              RENUM(tryy, val, ix, j1, j2, j3, j4);
              for (int i = 1; i <= 5; i++)
                x[i] = xtry[i];
            }
          }
        }
      }
    }

    //  403 continue
    //      if(val.eq.1.d+30) go to 430
    //      DO 404 I=1,4
    //      xhi(i)=xlo(i)+DBLE(ix(i)+1)*xstep(i)
    //  404 xlo(i)=xlo(i)+DBLE(ix(i)-3)*xstep(i)
    if (val == 1.0e+30) goto p430;
    for (int i = 1; i <= 4; i++) {
      xhi[i] = xlo[i] + double(ix[i] + 1) * xstep[i];
      xlo[i] = xlo[i] + double(ix[i] - 3) * xstep[i];
    }
  }
  //  408 CONTINUE
  //      DO 429 I=1,5
  //  429 XMEM(5,I)=X(I)
  //      VMEM(5)=SNGL(VAL)
  for (int i = 1; i <= 5; i++)
    xmem[5][i] = x[i];
  vmem[5] = val;

  //  430 DO 31 I=1,4
  p430: for (int i = 1; i <= 4; i++) {
    //      IF(VMEM(I).GT.VMEM(5)) GO TO 31
    //      VMEM(5)=VMEM(I)
    //      DO 32 J=1,5
    //   32 X(J)=XMEM(I,J)
    if (vmem[i] > vmem[5]) break;
    vmem[5] = vmem[i];
    for (int j = 1; j <= 5; j++)
      x[j] = xmem[i][j];
  }
  //   31 CONTINUE
  //      return
  //      END
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void Taquart::UsmtCore::f2(double x[], double &ffg) {
  //      SUBROUTINE F2(X,FFG)
  //      CHARACTER PS(80),TITLE*40
  //      REAL U(80),ARR(80),AZM(80),TKF(80)
  //      INTEGER RO(80),VEL(80),R(80)
  //      COMMON/MDATA/ PS,U,ARR,AZM,TKF,RO,VEL,R,TITLE,N,TROZ
  //      COMMON/PDATA/ A(80,6)
  //      DOUBLE PRECISION SUM,X(5),FFG

  //      FFG=DBLE(0.)
  ffg = 0.0;
  double SUM = 0.0;

  //      DO 1 I=1,N
  for (int i = 1; i <= N; i++) {
    //      SUM=DBLE(0.)
    SUM = 0.0;

    //      DO 2 J=1,5
    //    2 SUM=SUM+DBLE(A(I,J))*X(J)
    for (int j = 1; j <= 5; j++)
      SUM += (A[i][j] * x[j]);

    //      SUM=SUM-DBLE(A(I,6))*(X(1)+X(4))
    SUM -= (A[i][6] * (x[1] + x[4]));

    //      FFG=FFG+DABS(SUM-DBLE(U(I)))
    ffg += fabs(SUM - U[i]);

  }
  //    1 CONTINUE
  //      IF(DABS(FFG).GT.1.D+30) FFG=1.D+30
  if (fabs(ffg) > 1.0e+30) ffg = 1.0e+30;
  //      RETURN
  //      END
}

//-----------------------------------------------------------------------------
void Taquart::UsmtCore::POSTEP(int &METH, int &ITER, int &IND1) {
  //      SUBROUTINE POSTEP(METH,ITER,IND1)
  //      INTEGER*4 METH,ITER,IND1
  //      CHARACTER*42 TEX1,TEX2
  int NCALL = 0;

  //      IF(ITER.GT.1) GO TO 5
  if (ITER <= 1) {
    //      GO TO (10,20,30), METH
    //   10 R=FLOAT(IND1)/7.
    //      PERC=2.*(FLOAT(ITER-1)+R)
    //      GO TO 40
    //   20 R=FLOAT(IND1)/7.
    //      PERC=2.*(FLOAT(ITER-1)+R)
    //      GO TO 40
    //   30 R=FLOAT(IND1)/7.
    //      PERC=(FLOAT(ITER-1)+R)*.4
#ifdef USMTCORE_DEBUG
    double R = 0.0;
    double PERC = 0.0;
#endif
    switch (METH) {
      case 1:
#ifdef USMTCORE_DEBUG
        R = double(IND1 / 7.0);
        PERC = 2.0 * (double(ITER - 1) + R);
#endif
        break;
      case 2:
#ifdef USMTCORE_DEBUG
        R = double(IND1 / 7.0);
        PERC = 2.0 * (double(ITER - 1) + R);
#endif
        break;
      case 3:
#ifdef USMTCORE_DEBUG
        R = double(IND1 / 7.0);
        PERC = 0.4 * (double(ITER - 1) + R);
#endif
        break;
    }

    //   40 WRITE(75,41) METH,PERC
    //	WRITE(*,41) METH,PERC
    //   41 FORMAT(1X,/,' Method',I2,' (',F6.1,'% blks)',/,' .',$)

#ifdef USMTCORE_DEBUG
    std::cout << "Method " << METH << " ( " << PERC << "% blocks )" << std::endl;
#endif
    //      ITER=100
    //      NCALL=1
    ITER = 100;
    NCALL = 1;
    return;
    //      RETURN
  }

  //    5 IF(NCALL.EQ.73) GO TO 6
  if (NCALL != 73) {
    //      WRITE(75,42)
    //	WRITE(*,42)
    //   42 FORMAT(1H.,$)
    //      NCALL=NCALL+1
    //      RETURN
#ifdef USMTCORE_DEBUG
    std::cout << "|";
#endif
    NCALL++;
  }
  else {
    //    6 WRITE(75,43)
    //	WRITE(*,43)
    //   43 FORMAT(1H.,/,' .',$)
    //      NCALL=1
    //      RETURN
    //      END
#ifdef USMTCORE_DEBUG
    std::cout << "|!";
#endif
    NCALL = 1;
  }

}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
double Taquart::UsmtCore::DETR(double T[], double X) {
  //      FUNCTION DETR(T,X)
  //      DIMENSION S(6),T(6)

  //      SKAL=ABS(T(1))+ABS(T(2))+ABS(T(3))+ABS(T(4))+ABS(T(5))+ABS(T(6))
  //      IF(SKAL.LT.1.) SKAL=1.
  double SKAL = fabs(T[1]) + fabs(T[2]) + fabs(T[3]) + fabs(T[4]) + fabs(T[5])
      + fabs(T[6]);
  if (SKAL < 1.0) SKAL = 1.0;

  //      DO 1 I=1,6
  //    1 S(I)=T(I)/SKAL
  double S[6 + 1];
  for (int i = 1; i <= 6; i++)
    S[i] = T[i] / SKAL;

  //      Y=X/SKAL
  double Y = X / SKAL;

  //      DETR=(S(1)-Y)*(S(4)-Y)*(S(6)-Y)+2.*S(2)*S(3)*S(5)-S(2)*S(2)
  //     $*(S(6)-Y)-S(3)*S(3)*(S(4)-Y)-S(5)*S(5)*(S(1)-Y)
  return (S[1] - Y) * (S[4] - Y) * (S[6] - Y) + 2.0 * S[2] * S[3] * S[5]
      - S[2] * S[2] * (S[6] - Y) - S[3] * S[3] * (S[4] - Y)
      - S[5] * S[5] * (S[1] - Y);

  //      RETURN
  //      END
}

//-----------------------------------------------------------------------------
void Taquart::UsmtCore::RENUM(double &TRY, double &VAL, int ix[], int &j1,
    int &j2, int &j3, int &j4) {
  VAL = TRY;
  ix[1] = j1;
  ix[2] = j2;
  ix[3] = j3;
  ix[4] = j4;
}

//-----------------------------------------------------------------------------
void Taquart::UsmtCore::RDINP(Taquart::SMTInputData &InputData) {
  N = InputData.Count();
  TROZ = InputData.GetRuptureTime();
  Taquart::SMTInputLine InputLine;
  for (int i = 1; i <= N; i++) {
    //RPSTID[i] = i-1;
    //KNID[i] = i-1;
    //RPSTCP[i] = 'Z';
    //PS[i] = ' ';
    InputData.Get(i - 1, InputLine);
    U[i] = InputLine.Displacement;
    //ARR[i] = InputLine.Incidence; /* TODO -o3.1.19 : Code for SV and SH is switched off by default. */
    AZM[i] = InputLine.Azimuth;
    TKF[i] = InputLine.TakeOff;
    RO[i] = InputLine.Density;
    VEL[i] = InputLine.Velocity;
    R[i] = InputLine.Distance;
    //ACTIV[i] = 1;
  }
}

//-----------------------------------------------------------------------------
void Taquart::UsmtCore::SIZEMM(int &IEXP) {
  double X = 0.0;
  for (int i = 1; i <= 6; i++) {
    X = amax1(X, fabs(RM[i][2]));
  }
  X = 10.0 * X;
  IEXP = int(alog10(X)) + 1;
}

//-----------------------------------------------------------------------------

