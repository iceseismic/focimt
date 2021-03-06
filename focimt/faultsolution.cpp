//---------------------------------------------------------------------------
#include <trilib/fortranmath.h>
#include <trilib/georoutines.h>
#include <math.h>
#include "faultsolution.h"

//---------------------------------------------------------------------------
using namespace Taquart;

//---------------------------------------------------------------------------
FaultSolution::FaultSolution(void) {
  // Empty constructor
  DLA = 0.0;
  U_n = 0;
  UERR = 0.0;
  T0 = 0.0; /*!< Rupture time in seconds. */
  M0 = 0.0; /*!< Scalar moment tensor value in Nm. */
  MT = 0.0; /*!< Total seismic moment tensor value in Nm. */
  ERR = 0.0; /*!< Maximum error of the scalar seismic moment tensor
   value in Nm. This value is a square root of the maximum
   element for the moment tensor solution covariance
   matrix. */
  EXPL = 0.0; /*!< The size of the explosive component in the
   seismic moment tensor, in percents. */
  CLVD = 0.0; /*!< The size of the CLVD, compensated linear vector
   dipole component in the seismic moment tensor,
   in percents. */
  DBCP = 0.0; /*!< The size of the shear, double-couple component
   in the seismic moment tensor, in percents. */
  FIA = 0.0; /*!< Strike of the first fault plane in degrees. */
  DLA = 0.0; /*!< Dip of the first fault plane in degrees. */
  RAKEA = 0.0; /*!< Rake of the first fault plane in degrees. */
  FIB = 0.0; /*!< Strike of the second fault plane in degrees. */
  DLB = 0.0; /*!< Dip of the second fault plane in degrees. */
  RAKEB = 0.0; /*!< Rake of the second fault plane in degrees. */
  PXTR = 0.0; /*!< P-axis trend in degrees. */
  PXPL = 0.0; /*!< P-axis plunge in degrees.*/
  PXAM = 0.0; /*!< P-axis plunge in degrees.*/
  TXTR = 0.0; /*!< T-axis trend in degrees. */
  TXPL = 0.0; /*!< T-axis plunge in degrees.*/
  TXAM = 0.0; /*!< T-axis plunge in degrees.*/
  BXTR = 0.0; /*!< B-axis trend in degrees. */
  BXPL = 0.0; /*!< B-axis plunge in degrees.*/
  BXAM = 0.0; /*!< B-axis plunge in degrees.*/
  QI = 0.0;
  MAGN = 0.0;
}

//---------------------------------------------------------------------------
FaultSolution::~FaultSolution(void) {
  // Empty desctructor
}

//---------------------------------------------------------------------------
FaultSolution::FaultSolution(const FaultSolution &Source) {
  Assign(Source);
}

//---------------------------------------------------------------------------
/*
 XMLNode FaultSolution::xmlExport(XMLExporter &Exporter, String SectionName,
 String Variant, String Description)
 {
 XMLNode CurrentNode = Exporter.CurrentNode();
 XMLNode NewNode = Exporter.AddNode("fault_plane_solution");
 if(SectionName.Length() == 0)
 SectionName = "undefined";
 if(Variant.Length() == 0 && Description.Length() == 0)
 {
 Exporter.PutAttribute("version","1.0");
 Exporter.PutAttribute("type",SectionName);
 }
 else
 {
 Exporter.PutAttribute("version","1.1");
 Exporter.PutAttribute("type",SectionName);
 Exporter.PutAttribute("variant",Variant);
 Exporter.PutAttribute("description",Description);
 }

 Exporter.PutFloat("rupture_time",T0);
 Exporter.PutFloat("seismic_moment",M0);
 Exporter.PutFloat("total_moment",MT);
 Exporter.PutFloat("error",ERR);
 Exporter.PutFloat("expl",EXPL);
 Exporter.PutFloat("clvd",CLVD);
 Exporter.PutFloat("dbcp",DBCP);
 Exporter.PutFloat("a_strike",FIA);
 Exporter.PutFloat("b_strike",FIB);
 Exporter.PutFloat("a_dip",DLA);
 Exporter.PutFloat("b_dip",DLB);
 Exporter.PutFloat("a_rake",RAKEA);
 Exporter.PutFloat("b_rake",RAKEB);
 Exporter.PutFloat("moment_magnitude",MAGN);
 Exporter.PutFloat("quality_factor",QI);
 Exporter.PutFloat("p_trend",PXTR);
 Exporter.PutFloat("t_trend",TXTR);
 Exporter.PutFloat("b_trend",BXTR);
 Exporter.PutFloat("p_plunge",PXPL);
 Exporter.PutFloat("t_plunge",TXPL);
 Exporter.PutFloat("b_plunge",BXPL);
 Exporter.PutString("classification",Type);

 Exporter.PutFloat("m11",M[1][1]);
 Exporter.PutFloat("m21",M[2][1]);
 Exporter.PutFloat("m31",M[3][1]);
 Exporter.PutFloat("m22",M[2][2]);
 Exporter.PutFloat("m32",M[3][2]);
 Exporter.PutFloat("m33",M[3][3]);

 for(int i=1; i<=6; i++)
 for(int j=1; j<=6; j++)
 {
 Exporter.PutFloat("c"+FormatFloat("0",i)+FormatFloat("0",j),Covariance[i][j]);
 }

 Exporter.CurrentNode(CurrentNode);
 return NewNode;
 }
 */

//---------------------------------------------------------------------------
/*
 void FaultSolution::SaveIni(TMemIniFile *const File, String SectionName)
 {
 File -> WriteFloat(SectionName,"RuptureTime",T0);
 File -> WriteFloat(SectionName,"SeismicMoment",M0);
 File -> WriteFloat(SectionName,"TotalMoment",MT);
 File -> WriteFloat(SectionName,"Error",ERR);
 File -> WriteFloat(SectionName,"EXPL",EXPL);
 File -> WriteFloat(SectionName,"CLVD",CLVD);
 File -> WriteFloat(SectionName,"DBCP",DBCP);
 File -> WriteFloat(SectionName,"AAxisTrend",FIA);
 File -> WriteFloat(SectionName,"BAxisTrend",FIB);
 File -> WriteFloat(SectionName,"AAxisPlunge",DLA);
 File -> WriteFloat(SectionName,"BAxisPlunge",DLB);
 File -> WriteFloat(SectionName,"RakeA",RAKEA);
 File -> WriteFloat(SectionName,"RakeB",RAKEB);
 File -> WriteFloat(SectionName,"MomentMagnitude",MAGN);
 File -> WriteFloat(SectionName,"QualityFactor",QI);
 File -> WriteFloat(SectionName,"P-AxisTrend",PXTR);
 File -> WriteFloat(SectionName,"T-AxisTrend",TXTR);
 File -> WriteFloat(SectionName,"B-AxisTrend",BXTR);
 File -> WriteFloat(SectionName,"P-AxisPlunge",PXPL);
 File -> WriteFloat(SectionName,"T-AxisPlunge",TXPL);
 File -> WriteFloat(SectionName,"B-AxisPlunge",BXPL);
 File -> WriteString(SectionName,"Classification",Type);

 File -> WriteFloat(SectionName,"Tensor1",M[1][1]);
 File -> WriteFloat(SectionName,"Tensor2",M[2][1]);
 File -> WriteFloat(SectionName,"Tensor3",M[3][1]);
 File -> WriteFloat(SectionName,"Tensor4",M[2][2]);
 File -> WriteFloat(SectionName,"Tensor5",M[3][2]);
 File -> WriteFloat(SectionName,"Tensor6",M[3][3]);

 for(int i=1; i<=6; i++)
 for(int j=1; j<=6; j++)
 File -> WriteFloat(SectionName,
 "Covariance"+FormatFloat("0",i)+FormatFloat("0",j),Covariance[i][j]);
 }
 */

//---------------------------------------------------------------------------
void FaultSolution::sincos(double a, double *s, double *c) {
  *s = sin(a);
  *c = cos(a);
}

//---------------------------------------------------------------------------
/*
 double FaultSolution::computed_rake2(double str1, double dip1, double str2,
 double dip2, double fault) {


 double sinrake2;
 double sd, cd, cd2, ss, cs;

 sincos((str1 - str2) * DEG2RAD, &ss, &cs);
 sd = sind(dip1);
 cd = cosd(dip1);
 cd2 = cosd(dip2);

 if (fabs(dip2 - 90.) < EPSIL)
 sinrake2 = fault * cd;
 else
 sinrake2 = -fault * sd * cs / cd2;

 return Taquart::datan2(sinrake2, -fault * sd * ss);
 }
 */

//---------------------------------------------------------------------------
/*
 void FaultSolution::LoadIni(TMemIniFile *const File, String SectionName)
 {
 T0 = ::String2Float(File -> ReadString(SectionName,"RuptureTime","0"));
 M0 = ::String2Float(File -> ReadString(SectionName,"SeismicMoment","0"));
 MT = ::String2Float(File -> ReadString(SectionName,"TotalMoment","0"));
 ERR = ::String2Float(File -> ReadString(SectionName,"Error","0"));
 EXPL = ::String2Float(File -> ReadString(SectionName,"EXPL","0"));
 CLVD = ::String2Float(File -> ReadString(SectionName,"CLVD","0"));
 DBCP = ::String2Float(File -> ReadString(SectionName,"DBCP","0"));
 FIA = ::String2Float(File -> ReadString(SectionName,"AAxisTrend","0"));
 FIB = ::String2Float(File -> ReadString(SectionName,"BAxisTrend","0"));
 DLA = ::String2Float(File -> ReadString(SectionName,"AAxisPlunge","0"));
 DLB = ::String2Float(File -> ReadString(SectionName,"BAxisPlunge","0"));
 RAKEA = ::String2Float(File -> ReadString(SectionName,"RakeA","0"));
 RAKEB = ::String2Float(File -> ReadString(SectionName,"RakeB","0"));
 MAGN = ::String2Float(File -> ReadString(SectionName,"MomentMagnitude","0"));
 QI = ::String2Float(File -> ReadString(SectionName,"QualityFactor","0"));
 PXTR = ::String2Float(File -> ReadString(SectionName,"P-AxisTrend","0"));
 TXTR = ::String2Float(File -> ReadString(SectionName,"T-AxisTrend","0"));
 BXTR = ::String2Float(File -> ReadString(SectionName,"B-AxisTrend","0"));
 PXPL = ::String2Float(File -> ReadString(SectionName,"P-AxisPlunge","0"));
 TXPL = ::String2Float(File -> ReadString(SectionName,"T-AxisPlunge","0"));
 BXPL = ::String2Float(File -> ReadString(SectionName,"B-AxisPlunge","0"));
 Type = File -> ReadString(SectionName,"Classification","-");


 double Strike1 = FIA;
 double Strike2 = FIB;
 double Dip1 = DLA;
 double Dip2 = DLB;
 double FaultType = PXPL > TXPL ? -1.0 : 1.0;
 double Rake1 = computed_rake2(Strike2,Dip2,Strike1,Dip1,FaultType);
 double Rake2 = computed_rake2(Strike1,Dip1,Strike2,Dip2,FaultType);
 RAKEA = Rake1;
 RAKEB = Rake2;
 //------

 M[1][1] = ::String2Float(File -> ReadString(SectionName,"Tensor1","0"));
 M[2][1] = ::String2Float(File -> ReadString(SectionName,"Tensor2","0"));
 M[3][1] = ::String2Float(File -> ReadString(SectionName,"Tensor3","0"));
 M[1][2] = M[2][1];
 M[2][2] = ::String2Float(File -> ReadString(SectionName,"Tensor4","0"));
 M[3][2] = ::String2Float(File -> ReadString(SectionName,"Tensor5","0"));
 M[1][3] = M[3][1];
 M[2][3] = M[3][2];
 M[3][3] = ::String2Float(File -> ReadString(SectionName,"Tensor6","0"));

 for(int i=1; i<=6; i++)
 for(int j=1; j<=6; j++)
 Covariance[i][j] = ::String2Float(File -> ReadString(SectionName,
 "Covariance"+FormatFloat("0",i)+FormatFloat("0",j),-1.0));
 }
 */

//---------------------------------------------------------------------------
void FaultSolution::Assign(const FaultSolution &Source) {
  for (int i = 1; i < 4; i++)
    for (int j = 1; j < 4; j++)
      M[i][j] = Source.M[i][j];

  for (int i = 1; i <= 6; i++)
    for (int j = 1; j <= 6; j++)
      Covariance[i][j] = Source.Covariance[i][j];

  T0 = Source.T0;
  M0 = Source.M0;
  MT = Source.MT;
  ERR = Source.ERR;
  EXPL = Source.EXPL;
  CLVD = Source.CLVD;
  DBCP = Source.DBCP;
  FIA = Source.FIA;
  DLA = Source.DLA;
  RAKEA = Source.RAKEA;
  FIB = Source.FIB;
  DLB = Source.DLB;
  RAKEB = Source.RAKEB;
  PXTR = Source.PXTR;
  PXPL = Source.PXPL;
  PXAM = Source.PXAM;
  TXTR = Source.TXTR;
  TXPL = Source.TXPL;
  TXAM = Source.TXAM;
  BXTR = Source.BXTR;
  BXPL = Source.BXPL;
  BXAM = Source.BXAM;
  QI = Source.QI;
  MAGN = Source.MAGN;
  Type = Source.Type;
  for (unsigned int i = 0; i < MAXCHANNEL; i++) {
    U_th[i] = Source.U_th[i];
    U_measured[i] = Source.U_measured[i];
  }
  U_n = Source.U_n;
  //U_th = Source.U_th;
  //U_measured = Source.U_measured;
  UERR = Source.UERR;
}

//---------------------------------------------------------------------------
FaultSolution & FaultSolution::operator=(const FaultSolution &Source) {
  Assign(Source);
  return *this;
}

//---------------------------------------------------------------------------
Taquart::String FaultSolution::SubString(Taquart::String Line, int Start,
    int End) {
  return Line.SubString(Start, End - Start).Trim();
}

//---------------------------------------------------------------------------
/*
 void FaultSolution::FillCovarianceStringList(TStringList *Output)
 {
 String Line = "          ";;
 for(int i=1; i<=6; i++)
 {
 Line+="   C@si"+String(i)+"@o     "+"@b";
 }
 Output -> Add(Line);
 for(int j=1; j<=6; j++)
 {
 Line ="   C@s"+String(j)+"j@o   "+"@b";
 for(int i=1; i<=6; i++)
 {
 Line+=(" " + FormatFloat("' '0.00e+0;'-'0.00e+0",Covariance[i][j]) + " ");
 }
 Output -> Add(Line);
 }
 }
 */

//---------------------------------------------------------------------------
/*
 void FaultSolution::FillOutputStringList(TStringList *Output, bool PolishVersion)
 {
 if(PolishVersion)
 {
 for(int i=1; i<4; i++)
 Output -> Add(
 "M@s"+String(i)+"1@o = "+"@b"+(M[1][i]>=0 ? FormatFloat("+0.00E+00",M[1][i]) : FormatFloat("0.00E+00",M[1][i]))+"@n [Nm] "+
 "M@s"+String(i)+"2@o = "+"@b"+(M[2][i]>=0 ? FormatFloat("+0.00E+00",M[2][i]) : FormatFloat("0.00E+00",M[2][i]))+"@n [Nm] "+
 "M@s"+String(i)+"3@o = "+"@b"+(M[3][i]>=0 ? FormatFloat("+0.00E+00",M[3][i]) : FormatFloat("0.00E+00",M[3][i]))+"@n [Nm]");
 Output -> Add("");
 Output -> Add("Czas rozrywu: @b"+FormatFloat("0.0",T0*1000)+"@n [ms]");
 Output -> Add("Moment sejsmiczny: @b"+FormatFloat("0.00E+00",M0)+"@n [Nm]");
 Output -> Add("Ca�kowity moment sejsmiczny: @b"+FormatFloat("0.00E+00",MT)+"@n [Nm]");
 Output -> Add("B��d wyznaczenia: @b"+FormatFloat("0.00E+00",ERR)+"@n [Nm]");
 Output -> Add("Magnituda momentu sejsmicznego: @b"+FormatFloat("0.0",MAGN)+"@n");
 Output -> Add("");
 Output -> Add("Sk�adowe tensora:");
 Output -> Add(
 (  (EXPL!=0) ? (" EXPL: @b"+FormatFloat("0.0",EXPL)+"@n%") : String("")) +
 (  (CLVD!=0) ? ("  CLVD: @b"+FormatFloat("0.0",CLVD)+"@n%") : String("")) +
 (  (DBCP!=0) ? ("  DBCP: @b"+FormatFloat("0.0",DBCP)+"@n%") : String("")));

 Output -> Add("");
 Output -> Add("P�aszczyzny nodalne ( trend / zag��bienie / kierunek) :");
 Output -> Add(" P�aszczyzna A: @b"+
 FormatFloat("0.0",FIA)+"@n� / @b"+FormatFloat("0.0",DLA)+"@n� / @b"+FormatFloat("0.0",RAKEA)+"@n P�aszczyzna B: @b"+
 FormatFloat("0.0",FIB)+"@n� / @b"+FormatFloat("0.0",DLB)+"@n� / @b"+FormatFloat("0.0",RAKEB)+"@n�");
 Output -> Add("");
 Output -> Add("Trendy osi ( trend / zag��bienie ) :");
 Output -> Add(" O� P: @b"+FormatFloat("0",PXTR)+"@n� / @b"+FormatFloat("0",PXPL)
 +"@n�   O� T: @b"+FormatFloat("0",TXTR)+"@n� / @b"+FormatFloat("0",TXPL)
 +"@n�   O� B: @b"+FormatFloat("0",BXTR)+"@n� / @b"+FormatFloat("0",BXPL)+"@n�");
 Output -> Add("");
 Output -> Add("Wsp�czynnik jako�ci rozwi�zania: @b"+FormatFloat("0",QI)+"@n%");

 String PlType;
 if(Type.Pos("Normal fault"))
 PlType = "Uskok normalny";
 else if(Type.Pos("Reverse fault"))
 PlType = "Uskok odwr�cony";
 else if(Type.Pos("Strike fault"))
 PlType = "Uskok przesuwczy";

 Output -> Add("Klasyfikacja: @b"+PlType+"@n");
 }
 else
 {
 for(int i=1; i<4; i++)
 Output -> Add(
 "M@s"+String(i)+"1@o = "+"@b"+(M[1][i]>=0 ? FormatFloat("+0.00E+00",M[1][i]) : FormatFloat("0.00E+00",M[1][i]))+"@n [Nm] "+
 "M@s"+String(i)+"2@o = "+"@b"+(M[2][i]>=0 ? FormatFloat("+0.00E+00",M[2][i]) : FormatFloat("0.00E+00",M[2][i]))+"@n [Nm] "+
 "M@s"+String(i)+"3@o = "+"@b"+(M[3][i]>=0 ? FormatFloat("+0.00E+00",M[3][i]) : FormatFloat("0.00E+00",M[3][i]))+"@n [Nm]");
 Output -> Add("");
 Output -> Add("Rupture time: "+FormatFloat("0.0",T0*1000)+"[ms]");
 Output -> Add("Seismic moment: "+FormatFloat("0.00E+00",M0)+"[Nm]");
 Output -> Add("Total seismic moment: "+FormatFloat("0.00E+00",MT)+"[Nm]");
 Output -> Add("Error: "+FormatFloat("0.00E+00",ERR)+"[Nm]");
 Output -> Add("Seismic moment magnitude: "+FormatFloat("0.0",MAGN));
 Output -> Add("");
 Output -> Add("Decomposition:");
 Output -> Add(
 ((EXPL!=0) ? " EXPL: "+FormatFloat("0.0",EXPL) : String("")) +
 ((CLVD!=0) ? "  CLVD: "+FormatFloat("0.0",CLVD) : String("")) +
 ((DBCP!=0) ? "  DBCP: "+FormatFloat("0.0",DBCP) : String("")));

 Output -> Add("");
 Output -> Add("Nodal planes:  ");
 Output -> Add(" 1: "+
 FormatFloat("0.0",FIA)+"/"+FormatFloat("0.0",DLA)+"  2: "+
 FormatFloat("0.0",FIB)+"/"+FormatFloat("0.0",DLB));
 Output -> Add("");
 Output -> Add("Axes (trend/plunge):");
 Output -> Add(" P: "+FormatFloat("0",PXTR)+"/"+FormatFloat("0",PXPL)
 +"  T: "+FormatFloat("0",TXTR)+"/"+FormatFloat("0",TXPL)
 +"  B: "+FormatFloat("0",BXTR)+"/"+FormatFloat("0",BXPL));
 Output -> Add("");
 Output -> Add("Quality factor: "+FormatFloat("0",QI));
 Output -> Add("Classification: "+Type);
 }
 }
 */
//---------------------------------------------------------------------------
