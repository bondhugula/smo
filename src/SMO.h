#ifndef _SMO_H
#define _SMO_H

#include <isl/space.h>
#include <isl/set.h>
#include <isl/map.h>
#include <isl/constraint.h>
#include <isl/ctx.h>

#include <vector>
#include <string>
#include <string.h>
#include <map>
#include <iostream>

#include <assert.h>

#define DEBUG 0

#if DEBUG
#define DPRINT(x) x
#else
#define DPRINT(x) 
#endif

#define BIGC 10
using namespace std;

namespace Smo{
  typedef isl_map * MapPtr;
  typedef vector<MapPtr> MapPtrVec;
  typedef MapPtrVec::iterator MapPtrVecIter;

  typedef isl_basic_map * BMapPtr;
  typedef vector<BMapPtr> BMapPtrVec;
  typedef BMapPtrVec::iterator BMapPtrVecIter;

  typedef isl_set * SetPtr;
  typedef vector<SetPtr> SetPtrVec;
  typedef SetPtrVec::iterator SetPtrVecIter;

  typedef isl_basic_set * BSetPtr;
  typedef vector<BSetPtr> BSetPtrVec;
  typedef BSetPtrVec::iterator BSetPtrVecIter;

  typedef isl_constraint * ConstraintPtr;
  typedef vector<ConstraintPtr> ConstraintPtrVec;
  typedef ConstraintPtrVec::iterator ConstraintPtrVecIter;

  typedef vector<string> StringVec;
  typedef vector<string>::iterator StringVecIter;
}

#endif // _SMO_H

