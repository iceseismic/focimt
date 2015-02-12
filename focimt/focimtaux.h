//-----------------------------------------------------------------------------
#ifndef FOCIMTAUX_H_
#define FOCIMTAUX_H_
//-----------------------------------------------------------------------------
#include <vector>
#include <tricairo/tricairo_meca.h>
#include "faultsolution.h"
#include "inputdata.h"

extern bool DrawStations;
extern bool DrawAxes;
extern bool DrawCross;
extern bool DrawDC;
extern bool WulffProjection;
extern bool LowerHemisphere;

//-----------------------------------------------------------------------------
void SetFaultSolution(Taquart::FaultSolution &fu, double M11, double M12,
    double M13, double M22, double M23, double M33, double strike, double dip,
    double rake);

//-----------------------------------------------------------------------------
void GenerateBallCairo(Taquart::TriCairo_Meca &Meca,
    std::vector<Taquart::FaultSolutions> &FSList,
    Taquart::SMTInputData &InputData, Taquart::String Type);

//-----------------------------------------------------------------------------
bool Dispatch(Taquart::String &Input, Taquart::String &Chunk,
    Taquart::String delimiter);

//-----------------------------------------------------------------------------
#endif /* FOCIMTAUX_H_ */
