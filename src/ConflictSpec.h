#ifndef _ConflictSpec_H
#define _ConflictSpec_H

#include "SMO.h"

namespace Smo{
  // get the next line in this input file stream
  string NextLine(ifstream &ifs);

  class ConflictPoly;
  class ConflictSpec;

  typedef vector<ConflictPoly> ConflictPolyVec;
  typedef vector<ConflictSpec *> ConflictSpecPtrVec;

  class ConflictPoly{
    int srcStmtId;
    int destStmtId;
    BMapPtr bMapPtr;

  public:
    int GetSrcStmtId() const { return srcStmtId; }
    int GetDestStmtId() const { return destStmtId; }

    BMapPtr GetBMapPtr() const { return bMapPtr; }
    void SetBMapPtr(BMapPtr value) { bMapPtr=value; }

    ConflictPoly(int srcId, int destId, BMapPtr cPolyPtr)
      :srcStmtId(srcId),
       destStmtId(destId),
       bMapPtr(cPolyPtr)
    {}

    ConflictPoly(const ConflictPoly& cPolyRef)
      :srcStmtId(cPolyRef.GetSrcStmtId()),
       destStmtId(cPolyRef.GetDestStmtId()),
       bMapPtr(cPolyRef.GetBMapPtr())
    {}
  };

  class ConflictSpec{
    struct isl_ctx *_ctx;
    int _nStmts;
    int _nParams;
    int _nInputs;
    int _nOutputs;

    ConflictPolyVec &_cPolyVecRef;

    isl_printer *_printer;
  public:

    ConflictPolyVec &GetConflictPolyVecRef() const;
    bool IsEmpty();

    static ConflictSpec &Copy(const ConflictSpec &specToCopyRef);
    static ConflictSpec &Create(ifstream &ifs,int nStmts,int nConflictPoly,
				isl_ctx *ctx,BMapPtr orderPairPtr);
    static ConflictSpec &Create(__isl_keep isl_ctx *ctx,
				__isl_take BMapPtr map,
				vector<int> &uspan,
				bool inclusive,
				BSetPtrVec &liveOutBSetPtrVec,
				int tileStartDim,
		       	        int tileEndDim);
    void Destroy();

    static void Coalesce(ConflictPolyVec &cPolyVecRef,int srcStmtId,int destStmtId);

    struct isl_ctx *GetContextPtr() const { return _ctx;}
    int  GetNStmts() const { return _nStmts;}
    int  GetNParams() const { return _nParams;}
    int  GetNInputs() const { return _nInputs;}
    int  GetNOutputs() const { return _nOutputs;}
  private:
    static void formulateConflictSets(__isl_keep isl_ctx *ctx,
				      Smo::BMapPtr &rootBMapPtr,ConflictPolyVec &cPolyVecRef,
				      vector<int> &uspan,bool inclusive,
				      int currDim,int tileStartDim);

    ConflictSpec(struct isl_ctx *ctx,int nStmts,int nParams,int nInputs,int nOutputs,
		 ConflictPolyVec &cPolyVecRef);
    ~ConflictSpec();

    void operator=(const ConflictSpec&); // suppress
    ConflictSpec(const ConflictSpec&); // suppress
  };
}

#endif // _ConflictSpec_H
