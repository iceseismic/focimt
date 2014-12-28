//---------------------------------------------------------------------------
#include <trilib/tristat.h>
#include <algorithm>
#include "inputdata.h"
#include "timedist.h"
#include "moment_tensor.h"

//============================================================================
/*
 XMLNode Foci::SMTInputLine::xmlExport(XMLExporter &Exporter)
 {
 XMLNode CurrentNode = Exporter.CurrentNode();
 XMLNode NewNode = Exporter.AddNode("smtinput");
 Exporter.PutAttribute("version","1.0");
 Exporter.PutAttribute("name",Name);
 Exporter.PutAttribute("id",Id);
 Exporter.PutAttribute("marker",MarkerType);
 Exporter.PutString("name",Name);
 Exporter.PutInteger("id",Id);
 Exporter.PutString("component",Component);
 Exporter.PutString("marker",MarkerType);
 Exporter.PutFloat("start",Start);
 Exporter.PutFloat("end",End);
 Exporter.PutFloat("duration",Duration);
 Exporter.PutFloat("displacement",Displacement);
 Exporter.PutFloat("incidence",Incidence);
 Exporter.PutFloat("azimuth",Azimuth);
 Exporter.PutFloat("takeoff",TakeOff);
 Exporter.PutFloat("distance",Distance);
 Exporter.PutFloat("density",Density);
 Exporter.PutFloat("velocity",Velocity);
 Exporter.PutBoolean("pickactive",PickActive);
 Exporter.PutBoolean("channelactive",ChannelActive);

 Exporter.CurrentNode(CurrentNode);
 return NewNode;
 }
 */
//============================================================================
Taquart::SMTInputData::SMTInputData(void) {
  Clear();
  Key = 0;
}

//---------------------------------------------------------------------------
Taquart::SMTInputData::~SMTInputData(void) {
  Clear();
}

//---------------------------------------------------------------------------
unsigned int Taquart::SMTInputData::Count(void) {
  return InputData.size();
}

//---------------------------------------------------------------------------
void Taquart::SMTInputData::AddRuptureTime(double ARuptureTime) {
  RuptureTime = ARuptureTime;
}

//---------------------------------------------------------------------------
void Taquart::SMTInputData::GetStatistics(double &AMeanDuration, double &AStdDuration,
    double &AMeanDisplacement, double &AStdDisplacement) {
  AMeanDuration = MeanDuration;
  AStdDuration = StdDuration;
  AMeanDisplacement = MeanDisplacement;
  AStdDisplacement = StdDisplacement;
}

//---------------------------------------------------------------------------
double Taquart::SMTInputData::CountRuptureTime(bool &Result) {
  Result = false;

  // Calculate average rupture time on the basis of 3 stations.
  if (InputData.size() >= MIN_ALLOWED_CHANNELS) {
    std::vector<Taquart::TimeDist> Lista;
    for (unsigned int i = 0; i < InputData.size(); i++) {
      /* TODO 5 -c3.1.12 : This constraint will be  possibly too weak in the future. */
      if (InputData[i].PickActive == true && InputData[i].ChannelActive == true) {
        Taquart::TimeDist Item;
        Item.Time = InputData[i].Duration;
        Item.Distance = InputData[i].Distance;

        Lista.push_back(Item);
      }
    }

    if (Lista.size() > 5) {
      Result = true;
      //std::sort(Lista.begin(), Lista.end());
      std::sort(Lista.begin(), Lista.end(), Taquart::TimeDistComparator());

      const double Average = (Lista[1].Time + Lista[2].Time + Lista[3].Time) / 3.0f;
      RuptureTime = Average;
      return Average;
    }
    else {
      Result = false;
      RuptureTime = 0.0f;
      return 0.0f;
    }
  }
  RuptureTime = 0.0f;
  return 0.0f;
}

//---------------------------------------------------------------------------
bool Taquart::SMTInputData::Recalculate(void) {
  if (InputData.size() > 0) {
    std::vector<double> Displacements;
    std::vector<double> Durations;
    Displacements.reserve(InputData.size());
    Durations.reserve(InputData.size());
    for (unsigned int i = 0; i < InputData.size(); i++) {
      Displacements.push_back(InputData[i].Displacement);
      Durations.push_back(InputData[i].Duration);
    }

    //double a, b;
    MeanDisplacement = Taquart::mean(&Displacements[0], Displacements.size());
    StdDisplacement = Taquart::std(&Displacements[0], Displacements.size());
    MeanDuration = Taquart::mean(&Durations[0], Durations.size());
    StdDuration = Taquart::std(&Durations[0], Durations.size());
    return true;
  }
  else {
    MeanDuration = 0.0;
    StdDuration = 0.0;
    MeanDisplacement = 0.0;
    StdDisplacement = 0.0;
    return false;
  }
}

//---------------------------------------------------------------------------
unsigned int Taquart::SMTInputData::Add(Taquart::SMTInputLine &InputLine) {
  InputLine.Key = Key++;
  InputData.push_back(InputLine);
  return InputData.size();
}

//---------------------------------------------------------------------------
void Taquart::SMTInputData::Get(unsigned int Index, Taquart::SMTInputLine &InputLine)
    throw (Taquart::TriEOutOfRange) {
  if (Index >= InputData.size())
    throw Taquart::TriEOutOfRange("Index out of range for "
        "Foci::SMTInputData class member: std::vector<SMTInputLine> InputData");
  else
    InputLine = InputData[Index];
}

void Taquart::SMTInputData::Set(unsigned int Index, Taquart::SMTInputLine &InputLine)
    throw (Taquart::TriEOutOfRange) {
  if (Index >= InputData.size())
    throw Taquart::TriEOutOfRange("Index out of range for "
        "Foci::SMTInputData class member: std::vector<SMTInputLine> InputData");
  else
    InputData[Index] = InputLine;
}

//---------------------------------------------------------------------------
double Taquart::SMTInputData::GetDisplacement(const unsigned int &Index)
    throw (Taquart::TriEOutOfRange) {
  if (Index >= InputData.size())
    throw Taquart::TriEOutOfRange("Index out of range for "
        "Foci::SMTInputData class member: std::vector<SMTInputLine> InputData");
  else
    return InputData[Index].Displacement;
}

//---------------------------------------------------------------------------
bool Taquart::SMTInputData::Find(Taquart::String &ChannelName, unsigned int &Index) {
  for (unsigned int i = 0; i < InputData.size(); i++) {
    if (ChannelName == InputData[i].Name) {
      Index = i;
      return true;
    }
  }
  return false;
}

//---------------------------------------------------------------------------
bool Taquart::SMTInputData::Find(int Key, unsigned int &Index) {
  for (unsigned int i = 0; i < InputData.size(); i++) {
    if (Key == InputData[i].Key) {
      Index = i;
      return true;
    }
  }
  return false;
}

//---------------------------------------------------------------------------

void Taquart::SMTInputData::Remove(unsigned int Index) throw (Taquart::TriEOutOfRange) {
  if (Index >= InputData.size())
    throw Taquart::TriEOutOfRange("Index out of range for "
        "Foci::SMTInputData class member: std::vector<SMTInputLine> InputData");
  else
    InputData.erase(InputData.begin() + Index);
}

//---------------------------------------------------------------------------
/*
 XMLNode Taquart::SMTInputData::xmlExport(XMLExporter &Exporter) {
 XMLNode CurrentNode = Exporter.CurrentNode();
 XMLNode NewNode = Exporter.AddNode("smt_input_data");
 Exporter.PutAttribute("version", "1.0");
 Exporter.PutAttribute("count", InputData.size());
 Exporter.PutFloat("rupture_time", RuptureTime);
 for (unsigned int i = 0; i < InputData.size(); i++) {
 InputData[i].xmlExport(Exporter);
 }
 Exporter.CurrentNode(CurrentNode);
 return NewNode;
 }
 */

//---------------------------------------------------------------------------
Taquart::SMTInputData::SMTInputData(const Taquart::SMTInputData &Source) {
  Assign(Source);
}

//---------------------------------------------------------------------------
void Taquart::SMTInputData::Clear(void) {
  InputData.clear();
  RuptureTime = 0.0;
  MeanDuration = 0.0;
  StdDuration = 0.0;
  MeanDisplacement = 0.0;
  StdDisplacement = 0.0;
}

//---------------------------------------------------------------------------
Taquart::SMTInputData & Taquart::SMTInputData::operator=(const Taquart::SMTInputData &Source) {
  if (this != &Source) Assign(Source);
  return *this;
}

//---------------------------------------------------------------------------
void Taquart::SMTInputData::Assign(const SMTInputData &Source) {
  RuptureTime = Source.RuptureTime;
  InputData = Source.InputData;
  MeanDuration = Source.MeanDuration;
  StdDuration = Source.StdDuration;
  MeanDisplacement = Source.MeanDisplacement;
  StdDisplacement = Source.StdDisplacement;
}

//---------------------------------------------------------------------------
double Taquart::SMTInputData::GetRuptureTime(void) {
  return RuptureTime;
}

