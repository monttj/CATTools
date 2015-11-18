#include "TFile.h"
#include "TCut.h"

class MySample{

  public:

  MySample(string label, string fileName, const double lumi){
    m_file = TFile::Open(fileName.c_str());
    m_label = label;
    m_lumi = lumi;
    m_scale = 1.0;
  }

  MySample(string label, string fileName, double xsec, double nEvents, Color_t color, const double lumi, bool doStack = true){
    m_file = TFile::Open(fileName.c_str());
    m_label = label;
    m_xsec = xsec;
    m_nEvents = nEvents;
    m_color = color;
    m_lumi = lumi;
    m_doStack = doStack;
    m_scale = m_lumi * m_xsec / m_nEvents; 
  }

  string label(){
    return m_label;
  }

  TFile * GetFile(){
    return m_file;
  }

  double GetXSection(){
    return m_xsec;
  }

  double GetNEvents(){
    return m_nEvents;
  }

  Color_t GetColor(){
    return m_color;
  }

  bool doStack(){
    return m_doStack;
  }

  double GetLumi(){
    return m_lumi;
  }

  double GetScale(){
    return m_scale;
  }

  private:
    TFile * m_file;
    double m_xsec;
    double m_nEvents;
    string m_label;
    Color_t m_color;
    bool m_doStack;
    double m_lumi;
    double m_scale;
};
