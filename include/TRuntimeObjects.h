#ifndef _RUNTIMEOBJECTS_H_
#define _RUNTIMEOBJECTS_H_

#include <string>

#include "TCutG.h"
#include "TList.h"

/// Object passed to the online histograms.
/**
   For each event, an instance of this type will be passed to the custom histogrammer.
   This class contains all detectors present, and all existing cuts and histograms.
 */
class TRuntimeObjects : public TObject {
public:
  /// Constructor
  TRuntimeObjects(TList* detectors,
                  TList* objects,
                  TList* variables);

  /// Returns a pointer to the detector of type T
  template<typename T>
  T* GetDetector(){
    TIter next(detectors);
    TObject* obj;
    while((obj = next())){
      if(obj->InheritsFrom(T::Class())){
        return (T*)obj;
      }
    }
    return NULL;
  }

  TCutG* GetCut(const std::string& name);

  TList& GetDetectors();
  TList& GetObjects();
  TList& GetVariables();

  void FillHistogram(std::string name,
                     int bins, double low, double high, double value);
  void FillHistogram(std::string name,
                     int Xbins, double Xlow, double Xhigh, double Xvalue,
                     int Ybins, double Ylow, double Yhigh, double Yvalue);

  double GetVariable(const char* name);

private:
  TList* detectors;
  TList* objects;
  TList* variables;


  ClassDef(TRuntimeObjects, 0);
};

#endif /* _RUNTIMEOBJECTS_H_ */
