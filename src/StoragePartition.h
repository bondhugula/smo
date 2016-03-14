#ifndef _StoragePartition_H
#define _StoragePartition_H

#include "ConflictSpec.h"
#include "ISLUtils.h"
#include "SMO.h"
#include <set>

#include <glpk.h>

namespace Smo{
  class StoragePartition{
    class FarkasLemma{
      isl_ctx *_ctx;

      isl_space *_searchSpacePtr;
      BMapPtr _polyPtr;

      ConstraintPtrVec _cstVec;
      ConstraintPtrVec _eqCstVec;
      StringVec _multNames;
    public:
      FarkasLemma(isl_ctx *contextPtr, isl_space *searchSpacePtr, BMapPtr polyPtr)
	:_ctx(contextPtr),
	 _searchSpacePtr(searchSpacePtr),
	 _polyPtr(isl_basic_map_copy(polyPtr))
      {
	isl_basic_map_foreach_constraint(polyPtr,&ISL::AddAnyConstraint,&_cstVec);
	isl_basic_map_foreach_constraint(polyPtr,&ISL::AddIfEqConstraint,&_eqCstVec);
      }

      BSetPtr Apply(ISL::NameCoeffMapByLiteral &nameCoeffByLit);
      isl_ctx *GetContextPtr() const { return _ctx; }

      ~FarkasLemma(){
	isl_basic_map_free(_polyPtr);
      }
    private:
      void equateCoeffsCore(ConstraintPtr cPtr,isl_dim_type dtype,
			     unsigned multIndex,bool isEq,
			     ISL::NameCoeffMapByLiteral &nameCoeffByLit);
      BSetPtr equateCoeffs(BSetPtr ilp,isl_dim_type dtype, 
			      ISL::NameCoeffMapByLiteral &nameCoeffByLit);
      void equateCoeffsConstCore(ConstraintPtr cPtr,
			      int multIndex, bool isEq,
			      ISL::NameCoeffMapByLiteral &nameCoeffByLit);
      BSetPtr equateCoeffsConst(BSetPtr ilp,
				 ISL::NameCoeffMapByLiteral &nameCoeffByLit);

      BSetPtr getEliminationSpace(int nVars);

      void operator=(const FarkasLemma&); // suppress
      FarkasLemma(const FarkasLemma&); // suppress
    };

    ConflictSpec &_cSpecRef;

    StringVec _decisionVarNameVec;
    StringVec _satDecisionVarNameVec;

    StringVec _expModNameVec;
    StringVec _expModInterNameVec;
    StringVec _etaIntraNameVec;
    StringVec _etaInterNameVec;
    StringVec _coeffNameVec;

    StringVec _paramNameVec;

    string _etaIntraMax;
    string _etaInterMax;
    StringVec _coeffMaxVec;
    /* StringVec _expModMaxNameVec; */

    isl_space *_searchSpacePtr;

    vector<isl_mat *> _transformMatPtrVec; // transformation matrix
#if DEBUG
    vector<isl_mat *> _intraStmtModuloMatPtrVec; // moduli matrix
    vector<isl_mat *> _interStmtModuloMatPtrVec; // moduli matrix
#endif
    vector<isl_mat *> _moduloMatPtrVec; // moduli matrix

    StringVec GetExpModNameVec(int stmt){
      int size=_cSpecRef.GetNParams()+1;
      int offset=stmt*size;
      StringVec sub(_expModNameVec.begin()+offset,_expModNameVec.begin()+offset+size);
      return sub;
    }

    StringVec GetExpModInterNameVec(int stmt){
      int size=_cSpecRef.GetNParams()+1;
      int offset=stmt*size;
      StringVec sub(_expModInterNameVec.begin()+offset,_expModInterNameVec.begin()+offset+size);
      return sub;
    }

    string GetEtaIntra(int stmtId){ return _etaIntraNameVec[stmtId]; }
    string GetEtaInter(int stmtId){ return _etaInterNameVec[stmtId]; }

    string GetCoeffName(int stmt,int pos){
      return _coeffNameVec[stmt*(_cSpecRef.GetNInputs()+1)+pos];
    }

  public:
    StoragePartition(__isl_keep ConflictSpec &cSpec);
    StoragePartition(__isl_keep ConflictSpec &cSpec,
		     __isl_keep vector<isl_mat *> transformMatPtrVec,
	      	     __isl_keep vector<isl_mat *> moduloMatPtrVec);

    // getters and setters
    isl_ctx *GetContextPtr() const { return _cSpecRef.GetContextPtr(); }

    __isl_give isl_mat *GetTransformMatPtr(int stmtId) const {
      return isl_mat_copy(_transformMatPtrVec[stmtId]);
    }

#if DEBUG
    __isl_give isl_mat *GetIntraStmtModuloMatPtr(int stmtId) const {
      return isl_mat_copy(_intraStmtModuloMatPtrVec[stmtId]);
    }

    __isl_give isl_mat *GetInterStmtModuloMatPtr(int stmtId) const {
      return isl_mat_copy(_interStmtModuloMatPtrVec[stmtId]);
    }
#endif

    __isl_give isl_mat *GetModuloMatPtr(int stmtId) const {
      return isl_mat_copy(_moduloMatPtrVec[stmtId]);
    }

    isl_mat *GetContractionMod(int stmtId);

    __isl_give ConflictSpec &GetConflictSpecRef() const {
      return ConflictSpec::Copy(_cSpecRef);
    }

    void FindStorageHyperplane(bool enumerate,vector<StoragePartition *> &altStoragePartitionVec);
    bool AreAllConflictsSatisfied();
    void PrintStorageMapping(int stmtId);

  private:
    void init();

    void findAltStorageHyperplanes(__isl_keep BSetPtr ilp,
				   set<int> satisfiedConflictPolySet,
				   vector<StoragePartition *> &altStoragePartitionVec,
				   vector<isl_mat *> hyperplanesFound);

    BSetPtr createBasicConstraints();
    BSetPtr createBasicConstraintsCore(ConflictPoly &cPolyRef,
     					       int nStmts,unsigned index);

    BSetPtr approxExpModConstraint(ConflictPoly &cPolyRef,int nStmts,
				   StringVec &expModNameVec,bool flag);
    BSetPtr expModConstraints(ConflictPoly &cPolyRef,bool flag,
	                  		      int nStmts,StringVec &expModNameVec);
    BSetPtr decisionConstraints(ConflictPoly &cPolyRef,unsigned index);

    BSetPtr interStmtConflictSatConstraints(ConflictPoly &cPolyRef,
					    int index,bool flag,BSetPtr ilp,int stmtId);
    BSetPtr miscModConstraints(BSetPtr ilp,int nStmts,
					  map<int,int> stmtPolyCountMap);

    __isl_give BSetPtr constrainDecisionVars(__isl_take BSetPtr ilp,StringVec &varNameVec);
    __isl_give BSetPtr constrainSatDecisionVars(__isl_take BSetPtr ilp);

    BSetPtr boundEtaIntra(BSetPtr ilp,int nStmts);
    BSetPtr boundEtaInter(BSetPtr ilp,int nStmts);
    BSetPtr boundEtaSums(BSetPtr ilp,int nStmts);
    /* BSetPtr boundExpMods(BSetPtr ilp); */
    BSetPtr boundHyperplaneCoeffs(BSetPtr ilp);

    BSetPtr AddAltHyperplaneConstraint(BSetPtr ilp,set<int> satisfiedConflictPolySet);

    void solve(BSetPtr ilp,set<int> &satisfiedConflictPolySet);
    void glpkSolve(BSetPtr ilp, int ndims, int nParams,
				   vector<vector<int> > &hyperplanes,
				   vector<vector<int> > &moduli,
		                   vector<vector<int> > &interStmtModuli,
		                   vector<bool> &interStmtConflictSatisfactionVec,
		                   set<int> &satisfiedConflictPolySet);
    double glpkMinExpModuli(glp_prob *prob,int startIdx,int nStmts,
			  int nParams,double startWeight,int reduceFactor);
    void getHyperplaneCoeffs(isl_set *lexmin,vector<vector<int> > &hyperplanes,
			     vector<vector<int> > &moduli);

    isl_constraint *reviseInterStmtConflictPoly(ISL::NameCoeffMap dict,isl_space *spacePtr,
      					  isl_mat *modMatPtr,bool flag);
    isl_constraint *reviseIntraStmtConflictPoly(ISL::NameCoeffMap &dict,ConflictPoly &cPolyRef,
						       isl_space *spacePtr);
    void reviseConflictSet();
    isl_mat *smallerOf(isl_mat *mod1,isl_mat *mod2);
    void getHyperplaneCoeffs(isl_set *lexmin,vector<int> &hyp,
			     vector<int> &modulo);

    BSetPtr createLinearIndepConstraints(__isl_take BSetPtr ilp,
				      __isl_take isl_mat *rowMat);
    BSetPtr createLinearIndepConstraintsCore(__isl_take BSetPtr ilp,
			                   __isl_keep isl_mat *rker,
					     int index, bool flag);
    BSetPtr constrainSkew(__isl_take BSetPtr ilp,
			  __isl_keep isl_mat *rowMat,
	                  __isl_keep isl_mat *rker,
			  int col, int index, bool flag);

    BSetPtr boundCoefficientsSum(__isl_take BSetPtr ilp,int nStmts);

    void operator=(const StoragePartition&); // suppress
    StoragePartition(const StoragePartition&); // suppress
  };
}

#endif // _StoragePartition_H

