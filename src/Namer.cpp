#include <sstream>
#include <cassert>

#include "Namer.h"

namespace Smo{
  // 7 kind of name types,counters initilized to 0
  vector<int> Namer::_nameIndex(9,0);

  StringVec Namer::GetNames(int n,Namer::NameType t){
    return GetNames(n,t,string(""));
  }

  StringVec Namer::GetNames(int n,Namer::NameType t,string prefix){
    StringVec nameVec;

    switch(t){
    case kParam: nameVec=_GenNames(string("N").append(prefix),_nameIndex[kParam],n);
      break;
    case kInput: nameVec=_GenNames(string("x").append(prefix),_nameIndex[kInput],n);
      break;
    case kOutput: nameVec=_GenNames(string("y").append(prefix),_nameIndex[kOutput],n);
      break;
    case kFarkasMult: nameVec=_GenNames(string("L").append(prefix),_nameIndex[kFarkasMult],n);
      break;
    case kDecision: nameVec=_GenNames(string("d").append(prefix),_nameIndex[kDecision],n);
      break;
    case kHyperplaneCoeff: nameVec=_GenNames(string("c").append(prefix),_nameIndex[kHyperplaneCoeff],n);
      break;
    case kMisc: nameVec=_GenNames(string("u").append(prefix),_nameIndex[kMisc],n);
      break;
    case kEtaIntra: nameVec=_GenNames(string("n").append(prefix),_nameIndex[kEtaIntra],n);
      break;
    case kEtaInter: nameVec=_GenNames(string("n'").append(prefix),_nameIndex[kEtaInter],n);
      break;
    default: assert(false);
    }
    _nameIndex[t] += n;

    return nameVec;
  }

  StringVec Namer::_GenNames(string name,int start,int n){
    StringVec nameVec;
    for(int i=0; i<n; ++i){
      stringstream ss;
      ss << start+i;
      nameVec.push_back(name + ss.str());
      ss.flush();
    }
    return nameVec;
  }
}

