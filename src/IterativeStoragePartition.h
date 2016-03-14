#ifndef _IterativeStoragePartition_H
#define _IterativeStoragePartition_H

#include "ConflictSpec.h"
#include "SMO.h"

namespace Smo{
  class IterativeStoragePartition{
    ConflictSpec &_cSpecRef;
  public:
    IterativeStoragePartition(ConflictSpec& cSpecRef);
    void FindStorageHyperplanes(bool enumerate);

    // getters and setters
    ConflictSpec &GetConflictSpecRef() const { return _cSpecRef; }
    isl_ctx *GetContextPtr() const { return _cSpecRef.GetContextPtr(); }

  private:
    void operator=(const IterativeStoragePartition&); // suppress
    IterativeStoragePartition(const IterativeStoragePartition&); // suppress
  };
}

#endif // _IterativeStoragePartition_H

