/*
 * focimtaux.cpp
 *
 *  Created on: Feb 12, 2015
 *      Author: kwiatek
 */

#include "focimtaux.h"

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

