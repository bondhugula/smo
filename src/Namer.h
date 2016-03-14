#ifndef _Namer_H
#define _Namer_H

#include "SMO.h"

namespace Smo{
  // !!! ideally,should be a singleton
  class Namer{
    static vector<int> _nameIndex;
  public:
    enum NameType {kParam,kInput,kOutput,kFarkasMult,kDecision,
		   kHyperplaneCoeff,kMisc,kEtaIntra,kEtaInter};
    static StringVec GetNames(int n,Namer::NameType type);
    static StringVec GetNames(int n,Namer::NameType type,string prefix);
  private:
    static StringVec _GenNames(string name,int start,int n);

    void operator=(const Namer&); // suppress
    Namer(const Namer&); // suppress
  };
}

#endif // _Namer_H

