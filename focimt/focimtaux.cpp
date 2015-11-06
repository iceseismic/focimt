//-----------------------------------------------------------------------------
// Source: focimtaux.cpp
// Module: focimt
// Auxiliary function.
//
// Copyright (c) 2013-2015, Grzegorz Kwiatek.
//
// Permission is hereby granted, free of charge, to any person or organization
// obtaining a copy of the software and accompanying documentation covered by
// this license (the "Software") to use, reproduce, display, distribute,
// execute, and transmit the Software, and to prepare derivative works of the
// Software, and to permit third-parties to whom the Software is furnished to
// do so, all subject to the following:
//
// The copyright notices in the Software and this entire statement, including
// the above license grant, this restriction and the following disclaimer,
// must be included in all copies of the Software, in whole or in part, and
// all derivative works of the Software, unless such copies or derivative
// works are solely in the form of machine-executable object code generated by
// a source language processor.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
// SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
// FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//-----------------------------------------------------------------------------
#include "focimtaux.h"
#include <iostream>
//-----------------------------------------------------------------------------

// Default values.
bool DrawStations = true;
bool DrawAxes = true;
bool DrawCross = true;
bool DrawDC = true;
bool WulffProjection = false;
bool LowerHemisphere = true;

//-----------------------------------------------------------------------------
void SetFaultSolution(Taquart::FaultSolution &fu, double M11, double M12,
    double M13, double M22, double M23, double M33, double strike, double dip,
    double rake) {
  Taquart::nodal_plane NP;
  fu.M[1][1] = M11;
  fu.M[1][2] = M12;
  fu.M[1][3] = M13;
  fu.M[2][1] = M12;
  fu.M[2][2] = M22;
  fu.M[2][3] = M23;
  fu.M[3][1] = M13;
  fu.M[3][2] = M23;
  fu.M[3][3] = M33;
  fu.FIA = strike;
  fu.DLA = dip;
  fu.RAKEA = rake;
  NP.str = strike;
  NP.dip = dip;
  NP.rake = rake;
  fu.FIB = Taquart::computed_strike1(NP);
  fu.DLB = Taquart::computed_dip1(NP);
  fu.RAKEB = Taquart::computed_rake1(NP);
}

//-----------------------------------------------------------------------------
unsigned int CountSlash(Taquart::String Input) {
  Taquart::String null;
  Dispatch(Input, null, ":"); // Remove following blocks separated with ":" if they exists
  Input = null;
  unsigned int n = 0;
  for (int i = 1; i <= Input.Length(); i++) {
    if (Input[i] == '/')
      n++;
  }
  return n;
}

//-----------------------------------------------------------------------------
bool Dispatch(Taquart::String &Input, Taquart::String &Chunk,
    Taquart::String delimiter) {
  if (Input.Pos(delimiter) == 0)
    return false;
  else {
    Chunk = Input.SubString(1, Input.Pos(delimiter) - 1);
    Input = Input.SubString(Input.Pos(delimiter) + 1, 10000);
    return true;
  }
}

//-----------------------------------------------------------------------------
void GenerateBallCairo(Taquart::TriCairo_Meca &Meca,
    std::vector<Taquart::FaultSolutions> &FSList, Taquart::SMTInputData &id,
    Taquart::String Type) {

  //if (FSList.size() == 0) return;

  Taquart::FaultSolution * s;

  if (FSList.size() > 0) {
    if (Type == Taquart::String("dc")) {
      s = &FSList[0].DoubleCoupleSolution;
    }
    else if (Type == "deviatoric") {
      s = &FSList[0].TraceNullSolution;
    }
    else if (Type == "full") {
      s = &FSList[0].FullSolution;
    }
  }

  // Setup solution properties.
  Meca.DrawAxis = DrawAxes;
  Meca.DrawStations = DrawStations;
  Meca.DrawCross = DrawCross;
  Meca.DrawDC = DrawDC;

  // Draw seismic moment tensor solution.
  Meca.Projection = WulffProjection ? Taquart::prWulff : Taquart::prSchmidt;
  Meca.Hemisphere = LowerHemisphere ? Taquart::heLower : Taquart::heUpper;

  double cmt[6];
  Taquart::TriCairo_Axis P, T, N;
  if (FSList.size() > 0) {
    cmt[0] = s->M[3][3];
    cmt[1] = s->M[1][1];
    cmt[2] = s->M[2][2];
    cmt[3] = s->M[1][3];
    cmt[4] = s->M[2][3] * -1.0;
    cmt[5] = s->M[1][2] * -1.0;

    Taquart::TriCairo_MomentTensor mt;
    for (int i = 0; i < 6; i++)
      mt.f[i] = cmt[i];

    const double scal =
        sqrt(
            FOCIMT_SQ(mt.f[0]) + FOCIMT_SQ(mt.f[1]) + FOCIMT_SQ(mt.f[2])
                + 2.
                    * (FOCIMT_SQ(mt.f[3]) + FOCIMT_SQ(mt.f[4])
                        + FOCIMT_SQ(mt.f[5]))) / M_SQRT2;
    for (int i = 0; i < 6; i++)
      mt.f[i] = mt.f[i] / scal;

    Meca.GMT_momten2axe(mt, &T, &N, &P);

    // Overwrite DC with those calculated from MT
    Taquart::nodal_plane A, B;
    Taquart::axe2dc(T, P, &A, &B);
    s->FIA = A.str;
    s->DLA = A.dip;
    s->FIB = B.str;
    s->DLB = B.dip;
    //std::cout << A.str << " " << A.dip << " " << A.rake;
    //std::cout << B.str << " " << B.dip << " " << B.rake;

    Meca.Tensor(T, N, P);
  }

  if (DrawCross)
    Meca.CenterCross();

  // If MORE than one solution on the list, plot additional DC lines.
  if (FSList.size() > 1) {
    for (unsigned int i = 1; i < FSList.size(); i++) {
      Meca.BDCColor = Taquart::TCColor(0.5, 0.5, 0.5, 0.5);

      Taquart::FaultSolution * s;

      if (Type == "dc") {
        s = &FSList[i].DoubleCoupleSolution;
      }
      if (Type == "deviatoric") {
        s = &FSList[i].TraceNullSolution;
      }
      if (Type == "full") {
        s = &FSList[i].FullSolution;
      }

      // Set color in response to the type of the fault.
      if (s->Type == "NF") {
        Meca.BDCColor = Taquart::TCColor(0.0, 0.0, 1.0, 0.5);
      }
      else if (s->Type == "TF") {
        Meca.BDCColor = Taquart::TCColor(1.0, 0.0, 0.0, 0.5);
      }
      else {
        Meca.BDCColor = Taquart::TCColor(0.0, 1.0, 0.0, 0.5);
      }

      //std::cout << s -> FIA << " " << s -> DLA << std::endl;

      Meca.DoubleCouple(s->FIA, s->DLA);
      Meca.DoubleCouple(s->FIB, s->DLB);
    }
  }

  // Draw double-couple lines.
  if (DrawDC && FSList.size() > 0) {
    Meca.BDCColor = Taquart::TCColor(0.0, 0.0, 0.0, 1.0);
    //std::cout << s.FIA << " " << s.DLA << std::endl;
    Meca.DoubleCouple(s->FIA, s->DLA);
    Meca.DoubleCouple(s->FIB, s->DLB);
  }

  // Draw P and T axes' directions.
  if (DrawAxes && FSList.size() > 0) {
    Meca.Axis(P, "P");
    Meca.Axis(T, "T");
  }

  // Draw stations.
  if (DrawStations) {
    for (unsigned int i = 0; i < id.Count(); i++) {
      Taquart::SMTInputLine il;
      id.Get(i, il);

      // Calculate Gamma and Displacement.
      double GA[5];

      double Tko = il.TakeOff;
      if (Tko == 90.0f)
        Tko = 89.75f;

      GA[3] = cos(Tko * DEG2RAD);
      double help = sqrt(1.0f - GA[3] * GA[3]);
      GA[1] = cos(il.Azimuth * DEG2RAD) * help;
      GA[2] = sin(il.Azimuth * DEG2RAD) * help;
      double U = il.Displacement;
      Taquart::String Name = il.Name;

      double mx, my;
      Meca.Station(GA, U, Name, mx, my);
    }
  }

}

//-----------------------------------------------------------------------------
void DispatchStations(Taquart::String &StationString,
    Taquart::SMTInputData &InputData) {

  Taquart::String temp;
  double azimuth, takeoff, polarity;
  Taquart::String name;

  while (StationString.Length()) {
    Dispatch(StationString, temp, "/");
    azimuth = temp.ToDouble();
    Dispatch(StationString, temp, "/");
    takeoff = temp.ToDouble();
    Dispatch(StationString, temp, "/");
    polarity = temp.ToDouble();
    if (StationString.Pos(":")) {
      name = StationString.SubString(1, StationString.Pos(":") - 1);
      Dispatch(StationString, temp, ":"); // Cut rake part first, as it was already interpreted.
    }
    else {
      name = StationString.Trim();
      StationString = "";
    }

    Taquart::SMTInputLine il;
    il.Name = name; /*!< Station name.*/
    il.Id = 0; /*!< Station id number.*/
    il.Component = "Z";
    il.MarkerType = "P";
    il.Start = 0.0;
    il.End = 0.0;
    il.Duration = 0.0;
    il.Displacement = polarity; // area below the first P wave pulse is divided by angle of incicence. (vertical sensor)
    il.Incidence = takeoff;
    il.Azimuth = azimuth;
    il.TakeOff = takeoff;
    il.Distance = 1000;
    il.Density = 1000;
    il.Velocity = 1000;
    il.PickActive = true;
    il.ChannelActive = true;
    InputData.Add(il);
  }
}

//-----------------------------------------------------------------------------
void Dispatch2(Taquart::String &Input, double &v1, double &v2) {
  Taquart::String temp;
  Dispatch(Input, temp, "/");
  v1 = temp.Trim().ToDouble();
  v2 = Input.Trim().ToDouble();
}

//-----------------------------------------------------------------------------
void String2MT(Taquart::String &FaultString, double &M11, double &M12,
    double &M13, double &M22, double &M23, double &M33) {

// Extract moment tensor components from a string.
  Taquart::String temp;
  Dispatch(FaultString, temp, "/");
  M11 = temp.ToDouble();
  Dispatch(FaultString, temp, "/");
  M12 = temp.ToDouble();
  Dispatch(FaultString, temp, "/");
  M13 = temp.ToDouble();
  Dispatch(FaultString, temp, "/");
  M22 = temp.ToDouble();
  Dispatch(FaultString, temp, "/");
  M23 = temp.ToDouble();
  if (FaultString.Pos(":")) {
    M33 = FaultString.SubString(1, FaultString.Pos(":") - 1).ToDouble();
    Dispatch(FaultString, temp, ":"); // Cut rake part first, as it was already interpreted.
  }
  else {
    M33 = FaultString.Trim().ToDouble();
    FaultString = "";
  }
}

//-----------------------------------------------------------------------------
void String2SDR(Taquart::String &FaultString, double &strike, double &dip,
    double &rake) {
// Extract strike/dip/rake triplet from a string.
  Taquart::String temp;
  Dispatch(FaultString, temp, "/");
  strike = temp.ToDouble();
  Dispatch(FaultString, temp, "/");
  dip = temp.ToDouble();
  if (FaultString.Pos(":")) {
    rake = FaultString.SubString(1, FaultString.Pos(":") - 1).ToDouble();
    Dispatch(FaultString, temp, ":"); // Cut rake part first, as it was already interpreted.
  }
  else {
    rake = FaultString.Trim().ToDouble();
    FaultString = "";
  }
}

//-----------------------------------------------------------------------------
void DispatchFaults(Taquart::String &FaultString,
    std::vector<Taquart::FaultSolutions> &FSList, bool onefault) {
  Taquart::FaultSolutions fs;
  Taquart::FaultSolution fu;

// Read strike, dip and rake for the first fault plane
  Taquart::String temp;
  double strike = 0.0, dip = 0.0, rake = 0.0;
  double M11 = 0.0, M22 = 0.0, M33 = 0.0, M12 = 0.0, M13 = 0.0, M23 = 0.0;

  unsigned int n = CountSlash(FaultString);
  if (n == 2) {
    // String is strike/dip/rake
    String2SDR(FaultString, strike, dip, rake);
    Taquart::StrikeDipRake2MT(strike * DEG2RAD, dip * DEG2RAD, rake * DEG2RAD,
        M11, M22, M33, M12, M13, M23);
  }
  else {
    // String is M11/M12/M13/M22/M23/M33
    String2MT(FaultString, M11, M12, M13, M22, M23, M33);
  }

// Transfer strike/dip/rake to moment tensor.

  fs.Type = 'N';
  fs.Channel = 0;
  SetFaultSolution(fu, M11, M12, M13, M22, M23, M33, strike, dip, rake);
  fs.FullSolution = fu;
  fs.TraceNullSolution = fu;
  fs.DoubleCoupleSolution = fu;
  FSList.push_back(fs);
//std::cout << "FPS: " << strike << " " << dip << " " << rake << std::endl;

  if (onefault)
    return;

  while (FaultString.Length()) {
    if (n == 2) {
      String2SDR(FaultString, strike, dip, rake);
      Taquart::StrikeDipRake2MT(strike * DEG2RAD, dip * DEG2RAD, rake * DEG2RAD,
          M11, M22, M33, M12, M13, M23);
    }
    else {
      String2MT(FaultString, M11, M12, M13, M22, M23, M33);
    }
    fs.Type = 'J';
    fs.Channel = 0;
    SetFaultSolution(fu, M11, M12, M13, M22, M23, M33, strike, dip, rake);
    fs.FullSolution = fu;
    fs.TraceNullSolution = fu;
    fs.DoubleCoupleSolution = fu;
    FSList.push_back(fs);
    //std::cout << "FPS: " << fu.FIA << " " << fu.DLA << " " << fu.RAKEA<< std::endl;
  }

}

//-----------------------------------------------------------------------------
void SplitFilename(Taquart::String& str2, Taquart::String &file,
    Taquart::String &path) {
  size_t found;
  std::string str(str2.c_str());
//cout << "Splitting: " << str << endl;
  found = str.find_last_of("/\\");
//cout << " folder: " << str.substr(0,found) << endl;
  path = Taquart::String(str.substr(0, found));
//cout << " file: " << str.substr(found+1) << endl;
  file = Taquart::String(str.substr(found + 1));
}

//-----------------------------------------------------------------------------
void PlotStations(Taquart::String FaultString, Taquart::String FilenameOut,
    unsigned int Size) {
  std::vector<Taquart::FaultSolutions> FSList;
  Taquart::SMTInputData InputData;

// Dispatch "faults" string.
  DispatchStations(FaultString, InputData);

  Taquart::String OutName = FilenameOut + ".png";
  Taquart::TriCairo_Meca Meca(Size, Size, Taquart::ctSurface);
  GenerateBallCairo(Meca, FSList, InputData, "dc");
  Meca.Save(OutName);
}

//-----------------------------------------------------------------------------
void DrawFaultsStations(Taquart::String FaultString,
    Taquart::String StationString, Taquart::String FilenameOut,
    unsigned int Size) {
  std::vector<Taquart::FaultSolutions> FSList;
  Taquart::SMTInputData InputData;

  DispatchStations(StationString, InputData);
  DispatchFaults(FaultString, FSList, false);

  Taquart::String OutName = FilenameOut + ".png";
  Taquart::TriCairo_Meca Meca(Size, Size, Taquart::ctSurface);
  GenerateBallCairo(Meca, FSList, InputData, "dc");
  Meca.Save(OutName);

}

//-----------------------------------------------------------------------------
void DrawFaults(Taquart::String FaultString, Taquart::String FilenameOut,
    unsigned int Size) {
  std::vector<Taquart::FaultSolutions> FSList;
  Taquart::SMTInputData InputData;

// Dispatch "faults" string.
  DispatchFaults(FaultString, FSList, false);

  Taquart::String OutName = FilenameOut + ".png";
  Taquart::TriCairo_Meca Meca(Size, Size, Taquart::ctSurface);
  GenerateBallCairo(Meca, FSList, InputData, "dc");
  Meca.Save(OutName);

}

//-----------------------------------------------------------------------------
void DrawFault(Taquart::String FaultString, Taquart::String FilenameOut,
    unsigned int Size) {

  std::vector<Taquart::FaultSolutions> FSList;
  Taquart::SMTInputData InputData;

// Dispatch "faults" string.
  DispatchFaults(FaultString, FSList, true);

  Taquart::String OutName = FilenameOut + ".png";
  Taquart::TriCairo_Meca Meca(Size, Size, Taquart::ctSurface);
  GenerateBallCairo(Meca, FSList, InputData, "dc");
  Meca.Save(OutName);
}

//-----------------------------------------------------------------------------
double rand_normal(double mean, double stddev) { //Box muller method
  static double n2 = 0.0;
  static int n2_cached = 0;
  if (!n2_cached) {
    double x, y, r;
    do {
      x = 2.0 * rand() / RAND_MAX - 1;
      y = 2.0 * rand() / RAND_MAX - 1;

      r = x * x + y * y;
    } while (r == 0.0 || r > 1.0);
    {
      double d = sqrt(-2.0 * log(r) / r);
      double n1 = x * d;
      n2 = y * d;
      double result = n1 * stddev + mean;
      n2_cached = 1;
      return result;
    }
  }
  else {
    n2_cached = 0;
    return n2 * stddev + mean;
  }
}

//-----------------------------------------------------------------------------
void PrepareHelp(Options &listOpts) {
// 0
  listOpts.addOption("i", "input", "Full path to the input file", true);
// 1
  listOpts.addOption("o", "output",
      "Output file name (without extension).                \n\n"
          "    If specified, the output solution data in ASCII format will be exported to \n"
          "    a single file. Otherwise, 'fileid' field from input file will be used      \n"
          "    instead and moment tensor solution data will be exported to multiple files.\n",
      true);
// 2
  listOpts.addOption("s", "solution",
      "Output solution type.                                \n\n"
          "    Arguments: [F][T][D] for the (F)ull, (T)race-null and (D)ouble-couple      \n"
          "    solutions. Defines which moment tensor inversion will be performed. The    \n"
          "    default option is '-s D'. Combine three options to export desired moment   \n"
          "    tensor solutions, e.g. '-s DFT' will produce all three solutions at once.  \n",
      true);
// 3
  listOpts.addOption("t", "type",
      "Output file type.                                    \n\n"
          "    Arguments: [NONE][PNG][SVG][PS][PDF] for different output file types.      \n"
          "    Produce graphical representation of the moment tensor solution in a form of\n"
          "    the beach ball. More than one output file format can be specified. The     \n"
          "    default value is '-t PNG'.                                                 \n",
      true);
// 4
  listOpts.addOption("n", "norm",
      "Norm type.                                           \n\n"
          "    Arguments: [L1|L2] for L1 and L2 norm, respectively. Defines norm used in  \n"
          "    seismic moment tensor inversion. The default option is '-n L2' (faster).   \n"
          "    When Jacknife method is used the option is ignored and L2 norm is always   \n"
          "    used.                                                                      \n",
      true);
// 5
  listOpts.addOption("p", "projection",
      "Projection type.                                     \n\n"
          "    Arguments: [W|S][U|L]: Defines projection for the graphical representation \n"
          "    of the seismic moment tensor. Choose either (W)ulff projection or (S)chmidt\n"
          "    projection. Then select (U)pper hemisphere or (L)ower hemispere projection \n"
          "    The default option is '-p SL' (Schmidt projection, Lower hemisphere).      \n",
      true);
// 6
  listOpts.addOption("b", "ball",
      "The details of the beach ball picture                \n\n"
          "    Arguments: [S][A][C][D]: Defines features of the graphical representation  \n"
          "    of seismic moment tensor. Plot (S)tations, (A)xes, (C)enter cross, best    \n"
          "    (D)ouble-couple lines. The default option is '-b SACD' (all features are   \n"
          "    displayed on the beach ball.                                               \n",
      true);
// 7
  listOpts.addOption("d", "dump",
      "Output data format and order.                        \n\n"
          "    Arguments: [M][C][F][D][A][W][Q][T][U][*].                                 \n"
          "    [M]: Moment tensor components in Aki's convention: M11,M12,M13,M22,M23,M33.\n"
          "         The moment tensor components are in [Nm]                              \n"
          "    [C]: Moment tensor components in CMT conventions: M33,M11,M22,M13,-M23,-M12\n"
          "         The moment tensor components are in [Nm]                              \n"
          "    [F]: Fault plane solutions in format: STRIKEA/DIPA/RAKEA/STRIKEB/DIPB/RAKEB\n"
          "         (all values are in degrees)                                           \n"
          "    [D]: Decomposition of the moment tensor into Isotropic, Compensated linear \n"
          "         vector dipole and double-couple in format: ISO/CLVD/DBCP. The numbers \n"
          "         are provided in percents and calculated according to Jost and Herrmann\n"
          "         (1989) approach.                                                      \n"
          "    [Y]: Decomposition of the moment tensor into Isotropic, Compensated linear \n"
          "         vector dipole and double-couple in format: ISO/CLVD/DBCP. The numbers \n"
          "         are provided in percents and calculated according to Vavrycuk (2001)  \n"
          "         approach.                                                             \n"
          "    [A]: P/T/B Axes orientations in format:                                    \n"
          "         PTREND/PPLUNGE/TTREND/TPLUNGE/BTREND/BPLUNGE                          \n"
          "         All values are in degrees.                                            \n"
          "    [W]: Seismic moment, total seismic moment, maximum error of the seismic    \n"
          "         moment tensor estimate and the moment magnitude calculated using      \n"
          "         Hanks & Kanamori formula. The first three values are in [Nm]          \n"
          "    [Q]: Quality index                                                         \n"
          "    [T]: Fault type. 'SS','NF' or 'TF' will be exported depending whether the  \n"
          "         faulting style is strike-slip, normal or thrust, respectively.        \n"
          "    [U]: Vector of synthetic displacements calculated (the number of exported  \n"
          "         numbers correspond to the number of amplitudes in the input file.     \n"
          "    [E]: Scaled RMS Error calculated between theoretical and measured seismic  \n"
          "         moments.                                                              \n"
          "    [V]: Diagonal elements of the MT covariance matrix in the following order: \n"
          "         C11, C22, C33, C44, C55, C66                                          \n"
          "    [*]: Export new line character                                             \n"
          "                                                                               \n"
          "    NOTE #1:                                                                   \n"
          "    The order of arguments determine to order of output. For example -d FAD    \n"
          "    exports firstly fault plane solutions, then P, T and B axes directions and \n"
          "    finally the moment tensor decomposition into ISO/CLVD/DBCP. The output file\n"
          "    will have the following structure:                                         \n"
          "    STRIKEA/DIPA/RAKEA/STRIKEB/DIPB/RAKEB/PTREND/PPLUNGE/TTREND/TPLUNGE/BTREND \n"
          "    /BPLUNGE/ISO/CLVD/DBCP                                                     \n"
          "                                                                               \n"
          "    NOTE #2:                                                                   \n"
          "    Use lowercase arguments in order to export data in eye-friendly format.    \n",
      true);
// 8
  listOpts.addOption("m", "model",
      "Velocity model file (with extension)                 \n\n"
          "    Velocity model in HYPO71 format. Imposes different ASCII input file format.\n",
      true);
// 9
  listOpts.addOption("j", "jacknife", "Performs station Jacknife test.\n");
// 10
  listOpts.addOption("a", "amplitude",
      "Perform amplitude resampling.                        \n\n"
          "    Arguments: x[/y] where x is a floating-point positive number that describes\n"
          "    the level of noise applied to each amplitude: A+x*A*N(0,1)/3 where N is a  \n"
          "    normal distribution with mean 0 and std 1. The default value of x is 1     \n"
          "    (i.e.amplitude vary by a max. factor of ~2). Optional parameter /y is      \n"
          "    a number of samples (default value is 100).                                \n",
      true);
// 11
  listOpts.addOption("f", "drawfault",
      "Generate picture with fault plane solution           \n\n"
          "    Arguments: strike/dip/rake                                                 \n"
          "                                                                               \n",
      true);
// 12
  listOpts.addOption("fj", "drawfaults",
      "Generate picture with jacknife solutions             \n\n"
          "    Arguments: strike/dip/rake[:s1/d1/r1][:s2/d2/r2]...      \n"
          "                                                                               \n",
      true);
// 13
  listOpts.addOption("fs", "drawstations",
      "Generate picture with station and polarity data      \n\n"
          "    Arguments: azimuth/takeoff/polarity/name[:a2/t2/p2/n2][:a3/t3/p3/n3]...    \n"
          "                                                                               \n",
      true);
// 14
  listOpts.addOption("z", "size",
      "Beach ball file size                                 \n\n"
          "    Size of the beach ball figure in pixels.                                   \n",
      true);
// 15
  listOpts.addOption("rp", "resampling_polarity",
      "Perform phase polarities resampling                  \n\n"
          "    Performs additional MT inversions on resampled input data with randomly    \n"
          "    toggled polarities.                                                        \n"
          "    Arguments: x/y where x is the number of resamplings of the original dataset\n"
          "    and y is the fraction of reversed amplitudes.                              \n",
      true);
  // 16
  listOpts.addOption("rr", "resampling_rejection",
      "Perform station rejection resampling               \n\n"
          "    Performs additional MT inversions on resampled input data with randomly    \n"
          "    rejected stations.                                                         \n"
          "    Arguments: x/y where x is the number of resamplings of the original        \n"
          "    dataset and y is the fraction of rejected stations                         \n",
      true);
  // 17
  listOpts.addOption("ra", "resampling_amplitude",
      "Perform amplitude resampling                       \n\n"
          "    Performs additional MT inversions on resampled input data with randomly    \n"
          "    modified input amplitude data.                                             \n"
          "    Arguments: x/y where x is the number of resamplings of the original        \n"
          "    dataset and y is the amplitude variation factor (see option -a for details)\n",
      true);
  // 18
  listOpts.addOption("mt", "modeltakeoff",
      "Export raytracing data                               \n\n"
          "    Procedure export raytracing data for specific set of epicentral distances  \n"
          "    and epicentral depths for 1D velocity model file specified with option -m  \n"
          "    Arguments: dstart/dstep/dend/estart/estep/eend in [km]                     \n",
      true);

// 19
  listOpts.addOption("v", "version", "Display focimt version info");
}
