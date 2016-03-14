#ifndef _ISLUtils_H
#define _ISLUtils_H

#include <map>

#include "SMO.h"

namespace Smo{
  namespace ISL{
    void Print(isl_ctx *ctx, Smo::BMapPtrVec &bmapPtrVec);
    void Print(isl_ctx *ctx, Smo::BSetPtrVec &bsetPtrVec);

    typedef map<string,int> NameCoeffMap;
    typedef map<string,NameCoeffMap> NameCoeffMapByLiteral;

    __isl_give isl_space * isl_space_set_dim_names(__isl_take isl_space *space,
		enum isl_dim_type type,StringVec nameVec);

    void isl_space_get_dim_names(__isl_keep isl_space *space,
		enum isl_dim_type type,StringVec &nameVec);
    
    ConstraintPtr SetCoefficientVals(__isl_keep isl_space *spacePtr,
                 __isl_take ConstraintPtr cPtr,NameCoeffMap &nameCoeffMapRef);
    void GetCoefficientVals(__isl_keep ConstraintPtr c, __isl_keep isl_space *spacePtr,
			    StringVec names, isl_dim_type type,
			    StringVec mirrorNames, NameCoeffMap &nameCoeffMapRef);

    void GetILPSolution(__isl_keep BSetPtr bsetPtr, __isl_keep isl_space *spacePtr,
			int startIndex, int endIndex, vector<int> &coefficients);

    BSetPtr CreateSetFromConstraints(isl_ctx *ctx, isl_space *spacePtr,
				     StringVec &iParamVec, StringVec &oParamVec,
				     StringVec &iNameVec, StringVec &oNameVec,
				     ConstraintPtrVec &cstPtrVec);

    ConstraintPtr IneqFromNames(__isl_keep isl_space *spacePtr,
                                                 NameCoeffMap &nameCoeffMapRef);
    ConstraintPtr EqFromNames(__isl_keep isl_space *spacePtr,
                                                 NameCoeffMap &nameCoeffMapRef);

    int AddIfEqConstraint(ConstraintPtr c,void *cstVecPtr);
    int AddAnyConstraint(ConstraintPtr c,void *cstVecPtr);

    int AddAnyBasicSet(BSetPtr b,void *bsetVecPtr);
    int AddAnyBasicMap(BMapPtr b,void *bmapVecPtr);
    int AddAnySet(SetPtr b,void *setVecPtr);
    int AddAnyMap(MapPtr b,void *mapVecPtr);

    void ExtractAllBasicSet(isl_union_set *usetPtr,void *bsetVecPtr);
    void ExtractAllBasicMap(isl_union_map *umapPtr,void *bmapVecPtr);

    __isl_give isl_mat *AddRow(isl_ctx *ctx, __isl_take isl_mat *in, isl_mat *row);
    __isl_give isl_mat *AddRow(isl_ctx *ctx, __isl_take isl_mat *in, vector<int> row);
    __isl_give isl_mat *GetRowMatrix(isl_ctx *ctx, __isl_keep isl_mat *in, int row);
    __isl_give isl_mat *RemoveRow(isl_ctx *ctx, __isl_take isl_mat *in, int row);
    void PrintMatrix(__isl_give isl_mat *m);
  }
}

#endif // _ISLUtils_H

