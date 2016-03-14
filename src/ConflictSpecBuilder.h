#ifndef _ConflictSpecBuilder_H
#define _ConflictSpecBuilder_H

#include "SMO.h"
#include "ConflictSpec.h"

namespace Smo{
  class ConflictSpecBuilder{
    struct isl_ctx *_ctx;

    BMapPtrVec _rawdepPtrVec,_writePtrVec,_readPtrVec,_schedulePtrVec;
    BSetPtrVec _domPtrVec;
  
  public:
    ConflictSpecBuilder(string filename);

    bool FindMaxUtilitySpan(int stmtIndex,int tileStartDim,int tileEndDim,vector<int> &coefficients);
    void FindMaxUtilitySpanDeprecated(int stmtId,int tileStartDim,int tileEndDim,vector<int> &coefficients);
    ConflictSpec &Build(int stmtId,int tileStartDim,
			int tileEndDim,vector<int> &coefficients,
			bool inclusive,isl_union_set *liveOut);
    isl_union_set *InferLiveOut(int stmtId,int tileStart,int tileEnd);

    void PrintProgramInfo();

    struct isl_ctx *GetContextPtr() const { return _ctx;}
  private:
    BSetPtr mapBSetToParameterizedSpace(BSetPtr inbsetPtr,isl_space *newspace);
    BSetPtr parameterize(BSetPtr dom,StringVec &paramVec,StringVec &varVec);
    void parameterizeStmtDomains();

    int findInnerMostCommonTimeDim(BMapPtr rawDepPtr);
    int findInnerMostCommonTimeDim(BMapPtr rawDepPtr,int threshold);

    BSetPtr getDomain(int stmtId);

    ConflictSpecBuilder();
    void operator=(const Smo::ConflictSpec&); // suppress
    ConflictSpecBuilder(const Smo::ConflictSpec&); // suppress
  };

}

#endif // _ConflictSpecBuilder_H
