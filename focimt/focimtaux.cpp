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

  Taquart::FaultSolution s;

  if (Type == Taquart::String("dbcp")) {
    s = FSList[0].DoubleCoupleSolution;
  }
  else if (Type == "clvd") {
    s = FSList[0].TraceNullSolution;
  }
  else if (Type == "full") {
    s = FSList[0].FullSolution;
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
  cmt[0] = s.M[3][3];
  cmt[1] = s.M[1][1];
  cmt[2] = s.M[2][2];
  cmt[3] = s.M[1][3];
  cmt[4] = s.M[2][3] * -1.0;
  cmt[5] = s.M[1][2] * -1.0;

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

  Taquart::TriCairo_Axis P, T, N;
  Meca.GMT_momten2axe(mt, &T, &N, &P);
  Meca.Tensor(T, N, P);

  // Draw stations.
  if (DrawStations) {
    for (unsigned int i = 0; i < id.Count(); i++) {
      Taquart::SMTInputLine il;
      id.Get(i, il);

      // Calculate Gamma and Displacement.
      double GA[5];

      double Tko = il.TakeOff;
      if (Tko == 90.0f) Tko = 89.75f;

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

  // Draw P and T axes' directions.
  if (DrawAxes) {
    Meca.Axis(P, "P");
    Meca.Axis(T, "T");
  }

  if (DrawCross) Meca.CenterCross();

  // Draw double-couple lines.
  if (DrawDC) {
    Meca.BDCColor = Taquart::TCColor(0.0, 0.0, 0.0, 1.0);
    Meca.DoubleCouple(s.FIA, s.DLA);
    Meca.DoubleCouple(s.FIB, s.DLB);
  }

  // If MORE than one solution on the list, plot additional DC lines.
  if (FSList.size() > 1) {
    for (unsigned int i = 1; i < FSList.size(); i++) {
      Meca.BDCColor = Taquart::TCColor(0.5, 0.5, 0.5, 0.7);

      Taquart::FaultSolution * s;

      if (Type == "dbcp") {
        s = &FSList[i].DoubleCoupleSolution;
      }
      if (Type == "clvd") {
        s = &FSList[i].TraceNullSolution;
      }
      if (Type == "full") {
        s = &FSList[i].FullSolution;
      }

      // Set color in response to the type of the fault.
      if (s->Type == "NF") {
        Meca.BDCColor = Taquart::TCColor(0.0, 0.0, 1.0, 0.7);
      }
      else if (s->Type == "TF") {
        Meca.BDCColor = Taquart::TCColor(1.0, 0.0, 0.0, 0.7);
      }
      else {
        Meca.BDCColor = Taquart::TCColor(0.0, 1.0, 0.0, 0.7);
      }

      Meca.DoubleCouple(s->FIA, s->DLA);
      Meca.DoubleCouple(s->FIB, s->DLB);
    }
  }

}

//-----------------------------------------------------------------------------
void DrawFaults(Taquart::String FaultString, Taquart::String FilenameOut,
    unsigned int Size) {
  // Read strike, dip and rake.
  Taquart::String temp;
  double strike = 0, dip = 0, rake = 0;
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

  // Transfer strike/dip/rake to tensor.
  double M11, M22, M33, M12, M13, M23;
  Taquart::StrikeDipRake2MT(strike * DEG2RAD, dip * DEG2RAD, rake * DEG2RAD,
      M11, M22, M33, M12, M13, M23);

  std::vector<Taquart::FaultSolutions> FSList;
  Taquart::SMTInputData InputData;
  Taquart::FaultSolutions fs;
  Taquart::FaultSolution fu;

  fs.Type = 'N';
  fs.Channel = 0;
  SetFaultSolution(fu, M11, M12, M13, M22, M23, M33, strike, dip, rake);
  fs.FullSolution = fu;
  fs.TraceNullSolution = fu;
  fs.DoubleCoupleSolution = fu;
  FSList.push_back(fs);

  while (FaultString.Length()) {
    Dispatch(FaultString, temp, "/");
    strike = temp.ToDouble();
    Dispatch(FaultString, temp, "/");
    dip = temp.ToDouble();
    if (FaultString.Pos(":")) {
      rake = FaultString.SubString(1, FaultString.Pos(":") - 1).ToDouble();
      // Dispatch station data...
      Dispatch(FaultString, temp, ":"); // Cut rake part first, as it was already interpreted.
    }
    else {
      rake = FaultString.Trim().ToDouble();
      FaultString = "";
    }

    double M11, M22, M33, M12, M13, M23;
    Taquart::StrikeDipRake2MT(strike * DEG2RAD, dip * DEG2RAD, rake * DEG2RAD,
        M11, M22, M33, M12, M13, M23);

    fs.Type = 'J';
    fs.Channel = 0;
    SetFaultSolution(fu, M11, M12, M13, M22, M23, M33, strike, dip, rake);
    fs.FullSolution = fu;
    fs.TraceNullSolution = fu;
    fs.DoubleCoupleSolution = fu;
    FSList.push_back(fs);
  }

  Taquart::String OutName = FilenameOut + ".png";
  Taquart::TriCairo_Meca Meca(Size, Size, Taquart::ctSurface);
  GenerateBallCairo(Meca, FSList, InputData, "dbcp");
  Meca.Save(OutName);

}

//-----------------------------------------------------------------------------
void DrawFault(Taquart::String FaultString, Taquart::String FilenameOut,
    unsigned int Size) {
  // Read strike, dip and rake.
  Taquart::String temp;
  double strike = 0, dip = 0, rake = 0;
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

  // Transfer strike/dip/rake to tensor.
  double M11, M22, M33, M12, M13, M23;
  Taquart::StrikeDipRake2MT(strike * DEG2RAD, dip * DEG2RAD, rake * DEG2RAD,
      M11, M22, M33, M12, M13, M23);

  std::vector<Taquart::FaultSolutions> FSList;
  Taquart::SMTInputData InputData;
  Taquart::FaultSolutions fs;
  Taquart::FaultSolution fu;

  fs.Type = 'N';
  fs.Channel = 0;
  SetFaultSolution(fu, M11, M12, M13, M22, M23, M33, strike, dip, rake);
  fs.FullSolution = fu;
  fs.TraceNullSolution = fu;
  fs.DoubleCoupleSolution = fu;
  FSList.push_back(fs);

  Taquart::String OutName = FilenameOut + ".png";
  Taquart::TriCairo_Meca Meca(Size, Size, Taquart::ctSurface);
  GenerateBallCairo(Meca, FSList, InputData, "dbcp");
  Meca.Save(OutName);
}

void PrepareHelp(Options &listOpts) {
  // 0
  listOpts.addOption("i", "input", "Full path to the input file", true);
  // 1
  listOpts.addOption("o", "output", "Output file name (without extension)",
      true);
  // 2
  listOpts.addOption("s", "solution",
      "Output solution type.                                \n\n"
          "    Arguments: [F][T][D] for the full, trace-null and double-couple solution.  \n"
          "    Default option is '-s D'. Combine options to get multiple solutions, e.g.  \n"
          "    '-s DFT' produces all three solutions at once.                             \n",
      true);
  // 3
  listOpts.addOption("t", "type",
      "Output file type.                                    \n\n"
          "    Arguments: [NONE|PNG|SVG|PS|PDF] for different file types (only one can) be\n"
          "    specified as an output. The default value is '-t PNG'                      \n",
      true);
  // 4
  listOpts.addOption("n", "norm",
      "Norm type.               \n\n"
          "    Arguments: [L1|L2] for L1 and L2 norm, respectively. The default option is \n"
          "    is '-n L2' (faster). When Jacknife method is used the option is ignored and\n"
          "    L2 norm is used.                                                           \n",
      true);
  // 5
  listOpts.addOption("p", "projection",
      "Projection type.                                     \n\n"
          "    Arguments: [W|S][U|L]: Choose either (W)ulff projection or (S)chmidt       \n"
          "    projection. Then select (U)pper hemisphere or (L)ower hemispere projection \n"
          "    The default option is '-p SL'.                                             \n",
      true);
  // 6
  listOpts.addOption("b", "ball",
      "The details of the beach ball picture                \n\n"
          "    Arguments: [S][A][C][D]: Plot (S)tations, (A)xes, (C)enter cross, best     \n"
          "    (D)ouble-couple lines. The default option is '-b SACD' (all features are   \n"
          "    displayed.                                                                 \n",
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
          "         (all values in degrees)                                               \n"
          "    [D]: Decomposition of the moment tensor into Isotropic, Compensated linear \n"
          "         vector dipole and double couple in format: ISO/CLVD/DBCP. The numbers \n"
          "         are provided in %.                                                    \n"
          "    [A]: P/T/B Axes orientations in format:                                    \n"
          "         PTREND/PPLUNGE/TTREND/TPLUNGE/BTREND/BPLUNGE                          \n"
          "    [W]: Seismic moment, total seismic moment, maximum error of the seismic    \n"
          "         moment tensor estiamte and the moment magnitude calculated using      \n"
          "         using Hanks&Kanamori formula. The first three values are in [Nm]      \n"
          "    [Q]: Quality index                                                         \n"
          "    [T]: Fault Type (strike slip/normal/inverse)                               \n"
          "    [U]: Vector of synthetic moments calculated (the number of exported numbers\n"
          "         correspond to the number of amplitudes in the input file (specified   \n"
          "         using '-l' option.                                                    \n"
          "    [E]: RMS Error calculated from theoretical and measured ground             \n"
          "         displacements.                                                        \n"
          "    [V]: Diagonal elements of covariance matrix in the following order:        \n"
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
          "    Use lowercase arguments in order to data in eye-friendly format.           \n",
      true);
  // 8
  listOpts.addOption("m", "model",
      "Velocity model file (with extension)                 \n\n"
          "    Velocity model in HYPO71 format. Forces different input file format.       \n",
      true);
  // 9
  listOpts.addOption("j", "jacknife", "Switches on/off Jacknife test.\n");
  // 10
  listOpts.addOption("a", "amplitude",
      "Perform amplitude test.                              \n\n"
          "    Arguments: x[/y] where x is a floating-point positive number that describes\n"
          "    the level of noise applied to each amplitude: A+x*A*N(0,1)/3 where N is a  \n"
          "    normal distribution with mean 0 and std 1. The default value of x is 1     \n"
          "    (i.e.amplitude vary by a max. factor of ~2). Optional parameter /y is      \n"
          "    a number of samples (default value is 100).                                \n",
      true);
  // 11
  listOpts.addOption("f", "fault",
      "Draw fault plane solution directly (and stations).   \n\n"
          "    Arguments: strike/dip/rake[:azimuth1/takeoff1][:azimuth2/takeoff2]...      \n"
          "                                                                               \n",
      true);
  // 12
  listOpts.addOption("g", "faults",
      "Draw fault plane solution and bootstrap solutions.   \n\n"
          "    Arguments: strike/dip/rake[:s1/d1/r1][:s2/d2/r2]...      \n"
          "                                                                               \n",
      true);
  // 13
  listOpts.addOption("z", "size",
      "Beach ball file size                                 \n\n"
          "    Size of the beach ball figure in pixels.                                   \n",
      true);
  // 14
  listOpts.addOption("v", "version", "Display version number");
}
