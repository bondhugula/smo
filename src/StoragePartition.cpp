#include "StoragePartition.h"
#include "ISLUtils.h"
#include "Namer.h"

#include <algorithm> // for reverse
#include <math.h> //for pow

#include <isl/mat.h>


namespace Smo{

  StoragePartition::StoragePartition(__isl_keep ConflictSpec &cSpec)
    :_cSpecRef(cSpec)
  {
    _transformMatPtrVec.resize(_cSpecRef.GetNStmts());
#if DEBUG
    _intraStmtModuloMatPtrVec.resize(_cSpecRef.GetNStmts());
    _interStmtModuloMatPtrVec.resize(_cSpecRef.GetNStmts());
#endif
    _moduloMatPtrVec.resize(_cSpecRef.GetNStmts());
  }

  StoragePartition::StoragePartition(__isl_keep ConflictSpec &cSpec,
				     __isl_keep vector<isl_mat *> transformMatPtrVec,
				     __isl_keep vector<isl_mat *> moduloMatPtrVec)
    :_cSpecRef(cSpec),
     _transformMatPtrVec(transformMatPtrVec),
#if DEBUG
     _intraStmtModuloMatPtrVec(moduloMatPtrVec),
#endif
     _moduloMatPtrVec(moduloMatPtrVec)
  {
  }

  StoragePartition::StoragePartition(const StoragePartition &storagePartition)
    :_cSpecRef(ConflictSpec::Copy(storagePartition._cSpecRef)),
     _decisionVarNameVec(storagePartition._decisionVarNameVec),
     _satDecisionVarNameVec(storagePartition._satDecisionVarNameVec),
     _expModNameVec(storagePartition._expModNameVec),
     _expModInterNameVec(storagePartition._expModInterNameVec),
     _etaIntraNameVec(storagePartition._etaIntraNameVec),
     _etaInterNameVec(storagePartition._etaInterNameVec),
     _coeffNameVec(storagePartition._coeffNameVec),
     _paramNameVec(storagePartition._paramNameVec),
     _etaIntraMax(storagePartition._etaIntraMax),
     _etaInterMax(storagePartition._etaInterMax),
     _coeffMaxVec(storagePartition._coeffMaxVec),
     _searchSpacePtr(storagePartition._searchSpacePtr)
  {
    _transformMatPtrVec.resize(_cSpecRef.GetNStmts());
    _moduloMatPtrVec.resize(_cSpecRef.GetNStmts());
    for(int i=0; i<_cSpecRef.GetNStmts(); ++i){
      _transformMatPtrVec[i]=isl_mat_copy(storagePartition._transformMatPtrVec[i]);
      _moduloMatPtrVec[i]=isl_mat_copy(storagePartition._moduloMatPtrVec[i]);
    }
  }

  void StoragePartition::init(){
    ConflictPolyVec &conflictSet=_cSpecRef.GetConflictPolyVecRef();

// #if DEBUG
    // verify if any pair of conflict polytopes can be coalesced into one
    for(int i=0; i<conflictSet.size(); ++i){
      for(int j=i+1; j<conflictSet.size(); ++j){
        if(conflictSet[i].GetSrcStmtId()!=conflictSet[j].GetSrcStmtId()
	   || conflictSet[i].GetDestStmtId()!=conflictSet[j].GetDestStmtId())
	  continue;
        ConflictPoly first=ConflictPoly(conflictSet[i]);
	ConflictPoly second=ConflictPoly(conflictSet[j]);
	ConflictPolyVec polyVec;
	polyVec.push_back(first);
	polyVec.push_back(second);
	ConflictSpec::Coalesce(polyVec,first.GetSrcStmtId(),first.GetDestStmtId());
	if(polyVec.size()==1){
	  cout << "A pair of input conflict polytopes -- " << i+1 << " and " << j+1 << " can be coalesced!" << endl;
          isl_printer *p=isl_printer_to_file(GetContextPtr(),stdout);
	  isl_printer_print_basic_map(p,polyVec[0].GetBMapPtr());
	  isl_printer_free(p);
	  assert(false);
	}
      }
    }
// #endif

    int nStmts=_cSpecRef.GetNStmts();
    int nInputs=_cSpecRef.GetNInputs();
    int nParams=_cSpecRef.GetNParams();

    int nDecisions=2*conflictSet.size();
    _decisionVarNameVec=Namer::GetNames(nDecisions,Namer::kDecision);
    // _satDecisionVarNameVec=Namer::GetNames(nStmts,Namer::kDecision);

    _coeffMaxVec=Namer::GetNames(nStmts,Namer::kHyperplaneCoeff);

    // for each stmt, allocate names for etaIntra,etaInter,expansion modulo,hyp coeffs
    _expModNameVec.clear();
    _expModInterNameVec.clear();
#if DEBUG
    _etaIntraNameVec.clear();
    _etaInterNameVec.clear();
#endif
    _coeffNameVec.clear();
    // _expModMaxNameVec.clear();

    for(int i=0; i<nStmts; ++i){
      char str[4];
      sprintf(str,"%d_",i);
      StringVec names=Namer::GetNames(nParams+1,Namer::kMisc,str);
      _expModNameVec.insert(_expModNameVec.end(),names.begin(),names.end());

      names=Namer::GetNames(nParams+1,Namer::kMisc,str);
      _expModInterNameVec.insert(_expModInterNameVec.end(),names.begin(),names.end());

      names=Namer::GetNames(nInputs+1,Namer::kHyperplaneCoeff,str);
      _coeffNameVec.insert(_coeffNameVec.end(),names.begin(),names.end());
    }

#if DEBUG
    StringVec names=Namer::GetNames(nStmts,Namer::kEtaIntra);
    _etaIntraNameVec.insert(_etaIntraNameVec.end(),names.begin(),names.end());

    names=Namer::GetNames(nStmts,Namer::kEtaInter);
    _etaInterNameVec.insert(_etaInterNameVec.end(),names.begin(),names.end());
#endif

    _etaIntraMax=*Namer::GetNames(1,Namer::kEtaIntra).begin();
    _etaInterMax=*Namer::GetNames(1,Namer::kEtaInter).begin();

    // _expModMaxNameVec=Namer::GetNames(nParams+1,Namer::kMisc);

    const char *namePtr=NULL;
    isl_space *spacePtr=isl_basic_map_get_space(conflictSet.front().GetBMapPtr());

    // clear the param vector first and then populate it for this init
    _paramNameVec.clear();
    for(int i=0; i<nParams; ++i){
      namePtr=isl_space_get_dim_name(spacePtr,isl_dim_param,i);
      _paramNameVec.push_back(string(namePtr));
    }

    StringVec dimNameVec;
    dimNameVec.reserve(1 /* bound on the etaIntras */
		       + 1 /* bound on the etaInters */
		       // +_expModMaxNameVec.size()
		       + _coeffMaxVec.size() /* bounds on the coefficients */
		       +_expModNameVec.size()
		       +_expModInterNameVec.size()
#if DEBUG
		       +_etaIntraNameVec.size()
		       +_etaInterNameVec.size()
#endif 
		       +_decisionVarNameVec.size()
		       // +_satDecisionVarNameVec.size() 
		       +_coeffNameVec.size()
		       +_paramNameVec.size());
    dimNameVec.push_back(_etaIntraMax);
    dimNameVec.push_back(_etaInterMax);
    // dimNameVec.insert(dimNameVec.end(),_expModMaxNameVec.begin(),_expModMaxNameVec.end());
    dimNameVec.insert(dimNameVec.end(),_coeffMaxVec.begin(),_coeffMaxVec.end());
    dimNameVec.insert(dimNameVec.end(),_expModNameVec.begin(),_expModNameVec.end());
    dimNameVec.insert(dimNameVec.end(),_expModInterNameVec.begin(),_expModInterNameVec.end());
#if DEBUG
    dimNameVec.insert(dimNameVec.end(),_etaIntraNameVec.begin(),_etaIntraNameVec.end());
    dimNameVec.insert(dimNameVec.end(),_etaInterNameVec.begin(),_etaInterNameVec.end());
#endif 
    dimNameVec.insert(dimNameVec.end(),_decisionVarNameVec.begin(),_decisionVarNameVec.end());
    // dimNameVec.insert(dimNameVec.end(),_satDecisionVarNameVec.begin(),_satDecisionVarNameVec.end());
    dimNameVec.insert(dimNameVec.end(),_coeffNameVec.begin(),_coeffNameVec.end());

    _searchSpacePtr=isl_space_set_alloc(_cSpecRef.GetContextPtr(),_paramNameVec.size(),dimNameVec.size());
    _searchSpacePtr=ISL::isl_space_set_dim_names(_searchSpacePtr,isl_dim_param,_paramNameVec);
    _searchSpacePtr=ISL::isl_space_set_dim_names(_searchSpacePtr,isl_dim_set,dimNameVec);
  }

  static bool isRepeatedHyperplane(isl_mat * hyperplanesFound,isl_mat *newHyperplane){
    int nr=isl_mat_rows(hyperplanesFound);
    int nc=isl_mat_cols(hyperplanesFound);
    assert(isl_mat_rows(newHyperplane)==1);
    assert(isl_mat_cols(newHyperplane)==nc);

    for(int i=0; i<nr; ++i){
      bool same=true,negation=true;
      for(int j=0; j<nc; ++j){
        isl_val *val=isl_mat_get_element_val(newHyperplane,0,j);
	long newCoeff=isl_val_get_num_si(val);

        val=isl_mat_get_element_val(hyperplanesFound,i,j);
	long oldCoeff=isl_val_get_num_si(val);

	same=same && (newCoeff==oldCoeff);
	negation=negation && (newCoeff+oldCoeff==0);
      }
      if(same||negation){
	return true;
      }
    }
    return false;
  }

  void StoragePartition::findAltStorageHyperplanes(__isl_keep BSetPtr ilp,
						  set<int> satisfiedConflictPolySet,
						   vector<StoragePartition *> &altStoragePartitionVec,
						   vector<isl_mat *> hyperplanesFound){
    // keep a copy in case we need to find alternate storage hyperplanes
    StoragePartition *sp = new StoragePartition(*this); 
    BSetPtr ilpAlt=AddAltHyperplaneConstraint(isl_basic_set_copy(ilp),satisfiedConflictPolySet);

    satisfiedConflictPolySet.clear();
    solve(ilpAlt,satisfiedConflictPolySet);

    if(satisfiedConflictPolySet.size()>0){
      cout << "Satisfied " << satisfiedConflictPolySet.size() << " conflict polytopes" << endl;

      isl_ctx *ctx=GetContextPtr();
      vector<isl_mat *> hyperplanes;
      bool areHyperplanesRepeated=false;
      for(int i=0; i<_cSpecRef.GetNStmts(); ++i){
	  isl_mat *transformMatPtr=_transformMatPtrVec[i];
	  hyperplanes.push_back(ISL::GetRowMatrix(ctx,transformMatPtr,isl_mat_rows(transformMatPtr)-1));
	  areHyperplanesRepeated=areHyperplanesRepeated||isRepeatedHyperplane(hyperplanesFound[i],hyperplanes[i]);
      }

      if(!areHyperplanesRepeated){
	for(int i=0; i<_cSpecRef.GetNStmts(); ++i){
	  hyperplanesFound[i]=ISL::AddRow(ctx,hyperplanesFound[i],hyperplanes[i]);
	}
      }

      sp->findAltStorageHyperplanes(ilpAlt,satisfiedConflictPolySet,altStoragePartitionVec,hyperplanesFound);
      if(!areHyperplanesRepeated){
        altStoragePartitionVec.insert(altStoragePartitionVec.begin(),this);
      }
    }

    reviseConflictSet();
  }

  void StoragePartition::FindStorageHyperplane(bool enumerate,vector<StoragePartition *> &altStoragePartitionVec){
    init();
    BSetPtr ilp=createBasicConstraints();

    // keep a copy in case we need to find alternate storage hyperplanes
    StoragePartition *sp = new StoragePartition(*this); 

    set<int> satisfiedConflictPolySet;
    solve(ilp,satisfiedConflictPolySet);

    if(enumerate && satisfiedConflictPolySet.size()>0){
      cout << "Satisfied " << satisfiedConflictPolySet.size() << " conflict polytopes" << endl;
      
      vector<isl_mat *> hyperplanesFound; 
      hyperplanesFound.resize(_cSpecRef.GetNStmts());
      for(int i=0; i<_cSpecRef.GetNStmts(); ++i){
	  isl_ctx *ctx=GetContextPtr();
	  isl_mat *transformMatPtr=_transformMatPtrVec[i];
	  isl_mat *hyperplane=ISL::GetRowMatrix(ctx,transformMatPtr,isl_mat_rows(transformMatPtr)-1);
	  hyperplanesFound[i]=ISL::AddRow(ctx,hyperplanesFound[i],hyperplane);
      }

      sp->findAltStorageHyperplanes(ilp,satisfiedConflictPolySet,altStoragePartitionVec,hyperplanesFound);
    }
    reviseConflictSet();
  }

  BSetPtr StoragePartition::AddAltHyperplaneConstraint(BSetPtr ilp,set<int> satisfiedConflictPolySet){
    isl_printer *p=isl_printer_to_file(GetContextPtr(),stdout);
    assert(p!=NULL);

    ISL::NameCoeffMap dict;
    dict["1"]=satisfiedConflictPolySet.size()-1;
    for(set<int>::iterator i=satisfiedConflictPolySet.begin();
	i!=satisfiedConflictPolySet.end(); ++i){ 
      int conflictPolyIndex=*i;
      dict[_decisionVarNameVec[2*conflictPolyIndex]]=-1;
      dict[_decisionVarNameVec[2*conflictPolyIndex+1]]=-1;
    }

    isl_constraint *cstPtr=ISL::IneqFromNames(isl_basic_set_get_space(ilp),dict);

    DPRINT(cout << endl << endl);
    DPRINT(p=isl_printer_print_constraint(p,cstPtr));
    DPRINT(cout << endl);

    ilp=isl_basic_set_add_constraint(ilp,cstPtr);

    isl_printer_free(p);
    return ilp;
  }

  BSetPtr StoragePartition::createBasicConstraints(){
    isl_printer *p=isl_printer_to_file(GetContextPtr(),stdout);
    assert(p!=NULL);

    ConflictPolyVec &conflictSet=_cSpecRef.GetConflictPolyVecRef();
    int nStmts=_cSpecRef.GetNStmts();

    BSetPtr ilp=isl_basic_set_universe(isl_space_copy(_searchSpacePtr));
    unsigned index=0;

    DPRINT(cout << "The conflict sets are as follows..." << endl);
    BMapPtrVec bMapPtrVec;
    for(ConflictPolyVec::iterator i=conflictSet.begin();
	i!=conflictSet.end(); ++i,++index){
      bMapPtrVec.push_back(i->GetBMapPtr());
    }
    DPRINT(ISL::Print(GetContextPtr(),bMapPtrVec));

    index=0;
    map<int,int> stmtPolyCountMap;
    for(int i=0; i<nStmts; ++i){
      stmtPolyCountMap[i]=0;
    }
    for(ConflictPolyVec::iterator i=conflictSet.begin();
	i!=conflictSet.end(); ++i,++index){
      if(i->GetSrcStmtId()==i->GetDestStmtId()){
        stmtPolyCountMap[i->GetSrcStmtId()]+=1;
      }

      cout << "=========== TREATING A NEW CONFLICT POLYHEDRON ("
	   << i->GetSrcStmtId() << " " << i->GetDestStmtId() << ") =============\n";
      ilp=isl_basic_set_intersect(ilp,createBasicConstraintsCore(*i,nStmts,index));
      // ilp=isl_basic_set_remove_divs(ilp);
      //ilp=isl_basic_set_intersect(ilp,boundExpMods(cPolyRef));
    }

    ilp=constrainDecisionVars(ilp,_decisionVarNameVec); // add constraints: 0 <= d <= 1
    // ilp=constrainDecisionVars(ilp,_satDecisionVarNameVec); // add constraints: 0 <= d <= 1
    // ilp=constrainSatDecisionVars(ilp);
#if DEBUG
    ilp=boundEtaIntra(ilp,nStmts); 
    ilp=boundEtaInter(ilp,nStmts); 
#endif 
    ilp=boundEtaSums(ilp,nStmts); 
    ilp=boundHyperplaneCoeffs(ilp); // add constraints: -4 <= c <= 4
    ilp=boundCoefficientsSum(ilp,nStmts); // constraints:|c0|+|c1|+|c2|<=c3

    ilp=miscModConstraints(ilp,nStmts,stmtPolyCountMap); 

    isl_printer_free(p);
    return ilp;
  }

  void StoragePartition::solve(BSetPtr ilp,set<int> &satisfiedConflictPolySet){
    assert(satisfiedConflictPolySet.size()==0);
    isl_printer *p=isl_printer_to_file(GetContextPtr(),stdout);
    assert(p!=NULL);

    int ndims=isl_space_dim(_searchSpacePtr,isl_dim_set);
    int nParams=isl_space_dim(_searchSpacePtr,isl_dim_param);

    vector<vector<int> > hyperplanes,intraStmtModuli,interStmtModuli;
    vector<bool> interStmtConflictSatisfactionVec; // true if any inter stmt conflict satisfied

    int nStmts=_cSpecRef.GetNStmts();
    hyperplanes.resize(nStmts);
    intraStmtModuli.resize(nStmts);
    interStmtModuli.resize(nStmts);
    interStmtConflictSatisfactionVec.resize(nStmts);

    for(int i; i<nStmts; ++i){
      hyperplanes[i].resize(_coeffNameVec.size()/nStmts);
      intraStmtModuli[i].resize(nParams+1);
      interStmtModuli[i].resize(nParams+1);
      interStmtConflictSatisfactionVec[i]=false;
    }

    DPRINT(cout << endl);
    DPRINT(isl_printer_print_basic_set(p,ilp));
    DPRINT(cout << "\n\n======== LEXMIN ========"<< endl);

    if(true){ // use glpk for solving the ilp always
      ilp=isl_basic_set_remove_divs(ilp);
      glpkSolve(ilp,_coeffNameVec.size()/nStmts,nParams,
		hyperplanes,intraStmtModuli,interStmtModuli,
		interStmtConflictSatisfactionVec,satisfiedConflictPolySet);
    }
    else{
      isl_set *lexmin=isl_basic_set_lexmin(ilp);
      DPRINT(isl_printer_print_set(p,lexmin));
      DPRINT(cout << endl);

      getHyperplaneCoeffs(lexmin,hyperplanes,intraStmtModuli);
    }

    for(int i=0; i<_cSpecRef.GetNStmts(); ++i){
      _transformMatPtrVec[i]=ISL::AddRow(GetContextPtr(),_transformMatPtrVec[i],hyperplanes[i]);
#if DEBUG
      _intraStmtModuloMatPtrVec[i]=ISL::AddRow(GetContextPtr(),_intraStmtModuloMatPtrVec[i],intraStmtModuli[i]);
      _interStmtModuloMatPtrVec[i]=ISL::AddRow(GetContextPtr(),_interStmtModuloMatPtrVec[i],interStmtModuli[i]);
#endif
      _moduloMatPtrVec[i]=ISL::AddRow(GetContextPtr(),_moduloMatPtrVec[i],
				      interStmtConflictSatisfactionVec[i]?interStmtModuli[i]:intraStmtModuli[i]);
      cout << "Stmt " << i << ": any inter-statement conflict polytope satisfied? ";
      cout << (interStmtConflictSatisfactionVec[i] ? "Yes" : "No") << endl;
      // ISL::PrintMatrix(_transformMatPtrVec[i]);
      // ISL::PrintMatrix(_moduliMatPtrVec[i]);
    }
  }

  bool StoragePartition::AreAllConflictsSatisfied(){
    return _cSpecRef.IsEmpty();
  }

  void StoragePartition::PrintStorageMapping(int stmtId){
    int ndims=_cSpecRef.GetNInputs();
    isl_mat *transformMatPtr=GetTransformMatPtr(stmtId);
    isl_mat *moduloMatPtr=GetModuloMatPtr(stmtId);
    int nrows=isl_mat_rows(transformMatPtr);

    cout << "S" << stmtId << ": A" << stmtId << "[";
    for(int i=0; i<ndims; ++i){
      cout << "i" << i << ((i==ndims-1)?"":",");
    }
    cout << "] --> A" << stmtId << "[";
    for(int i=0; i<nrows; ++i){
      bool firstCoeffPrinted=false;
      int ncols=isl_mat_cols(transformMatPtr);
      for(int col=0; col<ncols; ++col){
	bool isConstantDim=(col==ncols-1);

        isl_val *val=isl_mat_get_element_val(transformMatPtr,i,col);
	long coeff=isl_val_get_num_si(val);

	if(coeff==0)
	  continue;
	if(firstCoeffPrinted || coeff<0)
	  cout << ((coeff>0)?"+":"-");
	coeff=abs(coeff);
	if(coeff!=1 || isConstantDim)
	  cout << coeff;
	if(!isConstantDim)// non-constant dimension
	  cout << "i" << col;
	firstCoeffPrinted=true;
      }
      if(!firstCoeffPrinted) // can we get into this scenario?
	cout << "0";
      cout << " mod ";

      firstCoeffPrinted=false;
      ncols=isl_mat_cols(moduloMatPtr);
      for(int col=0; col<ncols; ++col){
	bool isConstantDim=(col==ncols-1);

        isl_val *val=isl_mat_get_element_val(moduloMatPtr,i,col);
	long coeff=isl_val_get_num_si(val);

	if(isConstantDim)
	  coeff++;
	if(coeff==0)
	  continue;
	if(firstCoeffPrinted || coeff<0)
	  cout << ((coeff>0)?"+":"-");
	coeff=abs(coeff);
	if(coeff!=1 || isConstantDim)
	  cout << coeff;
	if(!isConstantDim)// non-constant dimension
	  cout << "N" << col;
	firstCoeffPrinted=true;
      }
      cout << ((i==nrows-1)?"":",");
    }
    cout << "]" << endl;
  }

  double StoragePartition::glpkMinExpModuli(glp_prob *prob,int startIdx,int nStmts,
					  int nParams,double startWeight,int reduceFactor){
    for(int offset=0; offset<(nParams+1); ++offset){
      for(int stmtId=0; stmtId<nStmts; ++stmtId){
        glp_set_obj_coef(prob,startIdx+offset+stmtId*(nParams+1)+1,startWeight); // |u0_off| + |u1_off| + .. + |u2_off|
      }
      startWeight/=reduceFactor;
    }
    return startWeight;
  }

  void StoragePartition::glpkSolve(BSetPtr ilp, int ndims, int nParams,
				   vector<vector<int> > &hyperplanes,
				   vector<vector<int> > &moduli,
				   vector<vector<int> > &interStmtModuli,
				   vector<bool> &interStmtConflictSatisfactionVec,
				   set<int> &satisfiedConflictPolySet){
    assert(satisfiedConflictPolySet.size()==0);

    ConstraintPtrVec cstPtrVec;
    isl_basic_set_foreach_constraint(ilp,&ISL::AddAnyConstraint,&cstPtrVec);

    int size=100000+1;
    int ia[100001];
    int ja[100001];
    double ar[100001];

    glp_prob *prob=glp_create_prob();
    glp_set_prob_name(prob,"problem");
    glp_set_obj_dir(prob,GLP_MIN);

    glp_add_rows(prob,cstPtrVec.size());
    int i=0;
    for(ConstraintPtrVecIter it=cstPtrVec.begin();
	it!=cstPtrVec.end();++it,++i){
      char buf[5];
      sprintf(buf,"%d",i+1);
      glp_set_row_name(prob,i+1,buf);

      long v=isl_val_get_num_si(isl_constraint_get_constant_val(*it));
      if(isl_constraint_is_equality(*it)){
	glp_set_row_bnds(prob,i+1,GLP_FX,-v,0.0);
      }
      else{
	glp_set_row_bnds(prob,i+1,GLP_LO,-v,0.0);
      }
    }

    isl_space *spacePtr=isl_basic_set_get_space(ilp);
    glp_add_cols(prob,isl_space_dim(spacePtr,isl_dim_set));
    for(int i=0; i<isl_space_dim(spacePtr,isl_dim_set); ++i){
      const char *dimName=isl_space_get_dim_name(spacePtr,isl_dim_set,i);
      glp_set_col_name(prob,i+1,dimName);
      if(find(_decisionVarNameVec.begin(),_decisionVarNameVec.end(),string(dimName))==_decisionVarNameVec.end()){
	 // || find(_satDecisionVarNameVec.begin(),_satDecisionVarNameVec.end(),string(dimName))==_satDecisionVarNameVec.end()){
	glp_set_col_kind(prob,i+1,GLP_IV);
	glp_set_col_bnds(prob,i+1,GLP_FR,0.0,0.0);
      }
      else{
	glp_set_col_kind(prob,i+1,GLP_BV);
      }
    }

    // Set the linear objective function coefficients
    // glp_set_obj_coef(prob,1,1000.0);
    // glp_set_obj_coef(prob,2,10.0);
    // glp_set_obj_coef(prob,3,0.1);

    int nStmts=_cSpecRef.GetNStmts();
    int factor= min(BIGC,50);
    double coeff=1000000.0;

    // 1. Maximize intra-stmt conflict poly satisfaction
    glp_set_obj_coef(prob,1,coeff); // etaIntra
    coeff/=factor;

    // 2. Minimize the intra-stmt modulo for each statement
    coeff=glpkMinExpModuli(prob,2+/*nParams+1+*/nStmts,nStmts,nParams,coeff,factor);

    // 3. Maximize inter-stmt conflict poly satisfaction
    glp_set_obj_coef(prob,2,coeff); // etaInter
    coeff/=factor;

    // 4. Minimize the inter-stmt modulo for each statement
    coeff=glpkMinExpModuli(prob,2+/*nParams+1+*/nStmts+nStmts*(nParams+1),nStmts,nParams,coeff,factor);

    for(i=2; i<2+/*nParams+1+*/nStmts; ++i){
      glp_set_obj_coef(prob,i+1,coeff); // sum of (bound on sum of coefficients) 
    }

    i=0;
    int count=1;
    for(ConstraintPtrVecIter it=cstPtrVec.begin();
	it!=cstPtrVec.end(); ++it,++i){
      isl_constraint *cstPtr=*it;
      for(int j=0; j<isl_space_dim(spacePtr,isl_dim_set); ++j){
	isl_val *val=isl_constraint_get_coefficient_val(cstPtr,isl_dim_set,j);
        long v=isl_val_get_num_si(val);
	if(v!=0){
	  ia[count]=i+1;
	  ja[count]=j+1;
	  ar[count]=v;
	  count=count+1;
	}
      }
    }

    DPRINT(cout << count-1 << endl);
    glp_load_matrix(prob,count-1,ia,ja,ar);

    glp_write_lp(prob,0,"./lp.lp");
    glp_adv_basis(prob,0);
    glp_simplex(prob,NULL);
    glp_intopt(prob,NULL);

    /* index of the zero'th dim of the zero'th modulo */
    int offset=1+1+/*_expModMaxNameVec.size()+*/nStmts; 
    /* index of the zero'th dim of the zero'th modulo */
    int interStmtModuloOffset=1+1+/*_expModMaxNameVec.size()+*/nStmts+(nParams+1)*nStmts;
    int decisionVarOffset=isl_space_dim(spacePtr,isl_dim_set)-_decisionVarNameVec.size()-_coeffNameVec.size()+1;
    ConflictPolyVec &conflictPolyVecRef=_cSpecRef.GetConflictPolyVecRef();

    for(i=1; i<=isl_space_dim(spacePtr,isl_dim_set); ++i){
      int val=glp_mip_col_val(prob,i);
      DPRINT(cout << isl_space_get_dim_name(spacePtr,isl_dim_set,i-1));
      DPRINT(cout << "=" << val << endl);

      int moduloLen=nParams+1;
      int moduloId=((i-offset-1)/moduloLen);
      int hyperplaneId=(nStmts-1)-((isl_space_dim(spacePtr,isl_dim_set)-i)/ndims);
      int interStmtModuloId=((i-interStmtModuloOffset-1)/moduloLen);

      if(i>=(offset+1) && moduloId<nStmts){
      	moduli[moduloId][moduloLen-1-(i-offset)%moduloLen]=val;
      }
      else if(hyperplaneId>=0){
      	hyperplanes[hyperplaneId][ndims-1-(isl_space_dim(spacePtr,isl_dim_set)-i)%ndims]=val;
      }
      else if(i>=(interStmtModuloOffset+1) && interStmtModuloId<nStmts){
      	interStmtModuli[interStmtModuloId][moduloLen-1-(i-interStmtModuloOffset)%moduloLen]=val;
      }
      else if(i>=decisionVarOffset && i<decisionVarOffset+_decisionVarNameVec.size()){
	int srcStmtId=conflictPolyVecRef[(i-decisionVarOffset)/2].GetSrcStmtId();
	int destStmtId=conflictPolyVecRef[(i-decisionVarOffset)/2].GetDestStmtId();
	if(srcStmtId!=destStmtId && val==1){
	  interStmtConflictSatisfactionVec[srcStmtId]=true;
	  interStmtConflictSatisfactionVec[destStmtId]=true;
	}
	if(val==1){
	  satisfiedConflictPolySet.insert((i-decisionVarOffset)/2);
	}
      }
    }
    glp_delete_prob(prob);
  }


  // flag ? (conflict difference) >= mod + 1
  //      : (conflict difference) <= -(mod + 1)
  isl_constraint *StoragePartition::reviseInterStmtConflictPoly(ISL::NameCoeffMap dict,isl_space *spacePtr,
						    isl_mat *modMatPtr,bool flag){

    int nModuli=isl_mat_rows(modMatPtr);
    if(!flag){
      for(ISL::NameCoeffMap::iterator it=dict.begin();
	  it!=dict.end(); ++it){
	it->second=-it->second;
      }
    }

    int nParams=isl_space_dim(spacePtr,isl_dim_param);
    for(int i=0; i<nParams; ++i){
      isl_val *val=isl_mat_get_element_val(modMatPtr,nModuli-1,i);
      int num_si=isl_val_get_num_si(val);

      dict[string(isl_space_get_dim_name(spacePtr,isl_dim_param,i))]=-num_si;
    }
    isl_val *val=isl_mat_get_element_val(modMatPtr,nModuli-1,nParams);
    int num_si=isl_val_get_num_si(val);// the constant dimension
    dict["1"]-=(num_si+1);

    isl_constraint *cstPtr=ISL::IneqFromNames(spacePtr,dict);
    return cstPtr;
  }

  isl_constraint *StoragePartition::reviseIntraStmtConflictPoly(ISL::NameCoeffMap &dict,ConflictPoly &cPolyRef,
						     isl_space *spacePtr){
    isl_mat *srcTransformMatPtr=GetTransformMatPtr(cPolyRef.GetSrcStmtId());
    isl_mat *destTransformMatPtr=GetTransformMatPtr(cPolyRef.GetDestStmtId());
    int nHyperplanes=isl_mat_rows(srcTransformMatPtr);
    unsigned i=0;
    for(; i<isl_space_dim(spacePtr,isl_dim_in); ++i){
      isl_val *val=isl_mat_get_element_val(srcTransformMatPtr,nHyperplanes-1,i);
      int src_num_si=isl_val_get_num_si(val);
      val=isl_mat_get_element_val(destTransformMatPtr,nHyperplanes-1,i);
      int dest_num_si=isl_val_get_num_si(val);

      dict[string(isl_space_get_dim_name(spacePtr,isl_dim_in,i))]=src_num_si;
      dict[string(isl_space_get_dim_name(spacePtr,isl_dim_out,i))]=-dest_num_si;
    }

    // the constant dimension
    isl_val *val=isl_mat_get_element_val(srcTransformMatPtr,nHyperplanes-1,i);
    int src_num_si=isl_val_get_num_si(val);
    val=isl_mat_get_element_val(destTransformMatPtr,nHyperplanes-1,i);
    int dest_num_si=isl_val_get_num_si(val);
    dict["1"]=src_num_si;
    dict["1"]-=dest_num_si;

    if(cPolyRef.GetSrcStmtId()==cPolyRef.GetDestStmtId())
      assert(dict["1"]==0);

    return ISL::EqFromNames(spacePtr,dict);
  }

  isl_mat *StoragePartition::GetContractionMod(int stmtId){
    return GetModuloMatPtr(stmtId);
  }

  isl_mat *StoragePartition::smallerOf(isl_mat *mod1,isl_mat *mod2){
    int c=isl_mat_cols(mod1);
    assert(c==isl_mat_cols(mod2));

    int a,b;
    for(int i=0; i<c; ++i){
      isl_val *val=isl_mat_get_element_val(mod1,0,i);
      a=isl_val_get_num_si(val);
      val=isl_mat_get_element_val(mod2,0,i);
      b=isl_val_get_num_si(val);
      if(a<b)
	return mod1;
      else if(a>b)
	return mod2;
    }
    return mod1;
  }

  void StoragePartition::reviseConflictSet(){
    // DPRINT(cout << "REVISING CONFLICT SET" << endl);
    isl_printer *p=isl_printer_to_file(GetContextPtr(),stdout);

    ConflictPolyVec &cPolyVec=_cSpecRef.GetConflictPolyVecRef();
    isl_space *spacePtr=isl_basic_map_get_space(cPolyVec.front().GetBMapPtr());

    ConflictPolyVec mapsToRemove;
    ConflictPolyVec mapsToAdd;
    for(ConflictPolyVec::iterator iter=cPolyVec.begin(); iter!=cPolyVec.end(); ++iter){
      // DPRINT(cout << endl << srcStmtId << " " << destStmtId << endl);
      ISL::NameCoeffMap m;
      BMapPtr origBMapPtr=iter->GetBMapPtr();
      int srcStmtId=iter->GetSrcStmtId();
      int destStmtId=iter->GetDestStmtId();

      isl_constraint *cstPtr=reviseIntraStmtConflictPoly(m,*iter,spacePtr);
      // DPRINT(p=isl_printer_prit_constraint(p,cstPtr); cout << endl);

      iter->SetBMapPtr(isl_basic_map_add_constraint(isl_basic_map_copy(origBMapPtr),isl_constraint_copy(cstPtr)));
      bool satisfied=isl_basic_map_is_empty(iter->GetBMapPtr());
      if(satisfied)
	mapsToRemove.push_back(*iter);

      // revising inter-statement conflict polyhedron
      if(srcStmtId!=destStmtId){
	DPRINT(cout << endl);
        ConflictPolyVec toAdd;
	isl_mat *srcModMatPtr=GetContractionMod(srcStmtId);
        isl_mat *destModMatPtr=GetContractionMod(destStmtId);
	isl_mat *smallerModMatPtr=smallerOf(srcModMatPtr,destModMatPtr);

	// (conflict difference) >= mod(src) + 1
        cstPtr=reviseInterStmtConflictPoly(m,spacePtr,smallerModMatPtr,true);
        // DPRINT(p=isl_printer_print_constraint(p,cstPtr); cout << endl);
	BMapPtr bMapPtr=isl_basic_map_add_constraint(isl_basic_map_copy(origBMapPtr),cstPtr);
	if(!isl_basic_map_is_empty(bMapPtr)){
          // DPRINT(p=isl_printer_print_basic_map(p,bMapPtr); cout << endl);
          mapsToAdd.push_back(ConflictPoly(srcStmtId,destStmtId,bMapPtr));
	}

	// // (conflict difference) >= mod(dest) + 1
        // cstPtr=reviseInterStmtConflictPoly(m,spacePtr,destModMatPtr,true);
	// bMapPtr=isl_basic_map_add_constraint(isl_basic_map_copy(origBMapPtr),cstPtr);
	// if(!isl_basic_map_is_empty(bMapPtr)){
        //   // DPRINT(p=isl_printer_print_basic_map(p,bMapPtr); cout << endl);
        //   mapsToAdd.push_back(ConflictPoly(srcStmtId,destStmtId,bMapPtr));
	// }

	// (conflict difference) <= -(mod(src) + 1)
        cstPtr=reviseInterStmtConflictPoly(m,spacePtr,smallerModMatPtr,false);
	bMapPtr=isl_basic_map_add_constraint(isl_basic_map_copy(origBMapPtr),cstPtr);
	if(!isl_basic_map_is_empty(bMapPtr)){
          // DPRINT(p=isl_printer_print_basic_map(p,bMapPtr); cout << endl);
          mapsToAdd.push_back(ConflictPoly(srcStmtId,destStmtId,bMapPtr));
	}

	// // (conflict difference) <= -(mod(dest) + 1)
        // cstPtr=reviseInterStmtConflictPoly(m,spacePtr,destModMatPtr,false);
	// bMapPtr=isl_basic_map_add_constraint(isl_basic_map_copy(origBMapPtr),cstPtr);
	// if(!isl_basic_map_is_empty(bMapPtr)){
        //   // DPRINT(p=isl_printer_print_basic_map(p,bMapPtr); cout << endl);
        //   mapsToAdd.push_back(ConflictPoly(srcStmtId,destStmtId,bMapPtr));
	// }
      }
    }

    // remove empty conflict polyhedra
    for(ConflictPolyVec::iterator r=mapsToRemove.begin(); r!=mapsToRemove.end(); ++r){
      for(ConflictPolyVec::iterator i=cPolyVec.begin(); i!=cPolyVec.end(); ++i){
	if(r->GetBMapPtr()==i->GetBMapPtr()){
	  cPolyVec.erase(i);
	  break;
	}
      }
    }
    // add non-empty conflict polyhedra
    for(ConflictPolyVec::iterator r=mapsToAdd.begin(); r!=mapsToAdd.end(); ++r){
      if(!isl_basic_map_is_empty(r->GetBMapPtr())){
        // DPRINT(p=isl_printer_print_basic_map(p,r->GetBMapPtr()); cout << endl);
	cPolyVec.push_back(*r);
      }
    }
    for(int i=0; i<_cSpecRef.GetNStmts(); ++i){
      for(int j=0; j<_cSpecRef.GetNStmts(); ++j){
	ConflictSpec::Coalesce(cPolyVec,i,j);
      }
    }
    isl_printer_free(p);
  }

  void StoragePartition::getHyperplaneCoeffs(isl_set *lexmin,
			        vector<vector<int> > &hyperplanes,
			        vector<vector<int> > &moduli){
    int nStmts=_cSpecRef.GetNStmts();
    vector<BSetPtr> bsetVec;
    isl_set_foreach_basic_set(lexmin,&ISL::AddAnyBasicSet,&bsetVec);
    assert(bsetVec.size()==1);
    ConstraintPtrVec cstPtrVec;
    isl_basic_set_foreach_constraint(bsetVec.front(),
				     &ISL::AddAnyConstraint,&cstPtrVec);

    int ndims=isl_space_dim(_searchSpacePtr,isl_dim_set);

    // determine the hyperplane coefficients
    int index=0;
    for(int stmtId=0; stmtId<nStmts; ++stmtId){
      for(unsigned i=0; i<_coeffNameVec.size()/nStmts; ++i,++index){
	for(ConstraintPtrVecIter iter=cstPtrVec.begin();
	    iter!=cstPtrVec.end(); ++iter){
	  isl_val *v=isl_constraint_get_coefficient_val(*iter,isl_dim_set,ndims-index-1);
	  long num=isl_val_get_num_si(v);
	  if(num==1){
	    hyperplanes[stmtId].push_back(-isl_val_get_num_si(isl_constraint_get_constant_val(*iter)));
	    break;
	  }
	}
      }
      std:reverse(hyperplanes[stmtId].begin(),hyperplanes[stmtId].end());
    }

    // determine the modulo coefficients
    index=1;
    for(int stmtId=0; stmtId<nStmts; ++stmtId){
      for(unsigned i=0; i<_paramNameVec.size()+1; ++i,++index){
	for(ConstraintPtrVecIter iter=cstPtrVec.begin();
	    iter!=cstPtrVec.end(); ++iter){
	  isl_val *v=isl_constraint_get_coefficient_val(*iter,isl_dim_set,index);
	  long num=isl_val_get_num_si(v);
	  if(num==1){
	    moduli[stmtId].push_back(-isl_val_get_num_si(isl_constraint_get_constant_val(*iter)));
	    break;
	  }
	}
      }
    }
  }

  BSetPtr StoragePartition::boundEtaIntra(BSetPtr ilp,int nStmts){
    isl_printer *p=isl_printer_to_file(GetContextPtr(),stdout);
    isl_constraint *cstPtr=NULL;
    isl_space *spacePtr=isl_basic_set_get_space(ilp);

    for(int stmtId=0; stmtId<nStmts; ++stmtId){ 
      // create the etaIntra bound for each stmt
      ISL::NameCoeffMap m;
      int index=0;
      m["1"]=0;
      for(ConflictPolyVec::iterator it=_cSpecRef.GetConflictPolyVecRef().begin();
	  it!=_cSpecRef.GetConflictPolyVecRef().end(); ++it,++index){

	if(it->GetSrcStmtId()!=stmtId || it->GetDestStmtId()!=stmtId)
	  continue; // ignore if this polytope is not related

	m["1"]+=1;
	m[GetEtaIntra(stmtId)]=-1;

	m[_decisionVarNameVec[2*index]]=-1;
	m[_decisionVarNameVec[2*index+1]]=-1;
      }
      cstPtr=ISL::EqFromNames(spacePtr,m);

      DPRINT(cout << endl << endl);
      DPRINT(p=isl_printer_print_constraint(p,cstPtr));
      DPRINT(cout << endl);

      // n_i = |K|-(sum of decision variables)
      ilp=isl_basic_set_add_constraint(ilp,cstPtr);

      m.clear();
      m[GetEtaIntra(stmtId)]=1;
      cstPtr=ISL::IneqFromNames(spacePtr,m); // n_i >= 0
      ilp=isl_basic_set_add_constraint(ilp,cstPtr);

      m.clear();
      m[GetEtaIntra(stmtId)]=-1;
      m[_etaIntraMax]=1;
      cstPtr=ISL::IneqFromNames(spacePtr,m); // n_max >= n_i

      DPRINT(p=isl_printer_print_constraint(p,cstPtr));
      DPRINT(cout << endl << endl);

      ilp=isl_basic_set_add_constraint(ilp,cstPtr);

    }

    isl_printer_free(p);

    return ilp;
  }

  // !!! try to merget his with boundEtaIntra
  BSetPtr StoragePartition::boundEtaInter(BSetPtr ilp,int nStmts){
    isl_printer *p=isl_printer_to_file(GetContextPtr(),stdout);
    isl_constraint *cstPtr=NULL;
    isl_space *spacePtr=isl_basic_set_get_space(ilp);

    for(int stmtId=0; stmtId<nStmts; ++stmtId){ 
      // create the etaIntra bound for each stmt
      ISL::NameCoeffMap m;
      int index=0;
      m["1"]=0;
      for(ConflictPolyVec::iterator it=_cSpecRef.GetConflictPolyVecRef().begin();
	  it!=_cSpecRef.GetConflictPolyVecRef().end(); ++it,++index){

	if((it->GetSrcStmtId()!=stmtId && it->GetDestStmtId()!=stmtId)
	   || it->GetSrcStmtId()==it->GetDestStmtId())
	  continue; // ignore if this polytope is not related

	m["1"]+=1;
	m[GetEtaInter(stmtId)]=-1;

	m[_decisionVarNameVec[2*index]]=-1;
	m[_decisionVarNameVec[2*index+1]]=-1;
      }
      isl_constraint *cstPtr=ISL::EqFromNames(spacePtr,m);

      DPRINT(cout << endl << endl);
      DPRINT(p=isl_printer_print_constraint(p,cstPtr));
      DPRINT(cout << endl);

      // n_i = |K|-(sum of decision variables)
      ilp=isl_basic_set_add_constraint(ilp,cstPtr);

      m.clear();
      m[GetEtaInter(stmtId)]=1;
      cstPtr=ISL::IneqFromNames(spacePtr,m); // n_i >= 0
      ilp=isl_basic_set_add_constraint(ilp,cstPtr);

      m.clear();
      m[GetEtaInter(stmtId)]=-1;
      m[_etaInterMax]=1;
      cstPtr=ISL::IneqFromNames(spacePtr,m); // n_max >= n_i

      DPRINT(p=isl_printer_print_constraint(p,cstPtr));
      DPRINT(cout << endl << endl);

      ilp=isl_basic_set_add_constraint(ilp,cstPtr);

    }

    isl_printer_free(p);
    return ilp;
  }

  BSetPtr StoragePartition::boundEtaSums(BSetPtr ilp,int nStmts){
    isl_printer *p=isl_printer_to_file(GetContextPtr(),stdout);
    isl_constraint *cstPtr=NULL;
    isl_space *spacePtr=isl_basic_set_get_space(ilp);

    ISL::NameCoeffMap dict;

    // First, set up the bound for the eta_intra
    dict["1"]=0;
    dict[_etaIntraMax]=1;

    int index=0;
    // n_max = Sigma(d_1i+d_2i) for all K_i which are intra-statement polyhedra
    for(ConflictPolyVec::iterator it=_cSpecRef.GetConflictPolyVecRef().begin();
	  it!=_cSpecRef.GetConflictPolyVecRef().end(); ++it,++index){

	if(it->GetSrcStmtId() != it->GetDestStmtId())
	  continue; // ignore inter-statement polyhedra

	dict["1"]+=-1;
	dict[_decisionVarNameVec[2*index]]=1;
	dict[_decisionVarNameVec[2*index+1]]=1;
    }


    DPRINT(cout << "Bound on sum of eta_intra values:" << endl);
    cstPtr=ISL::EqFromNames(spacePtr,dict); 

    DPRINT(p=isl_printer_print_constraint(p,cstPtr));
    DPRINT(cout << endl << endl);

    ilp=isl_basic_set_add_constraint(ilp,cstPtr);

    // Next, set up the bound for the eta_inter
    dict.clear();
    dict["1"]=0;
    dict[_etaInterMax]=1;

    // n_max = Sigma(d_1i+d_2i) for all K_i which are intra-statement polyhedra
    index=0;
    for(ConflictPolyVec::iterator it=_cSpecRef.GetConflictPolyVecRef().begin();
	  it!=_cSpecRef.GetConflictPolyVecRef().end(); ++it,++index){

	if(it->GetSrcStmtId()==it->GetDestStmtId())
	  continue; // ignore inter-statement polyhedra

	dict["1"]+=-1;
	dict[_decisionVarNameVec[2*index]]=1;
	dict[_decisionVarNameVec[2*index+1]]=1;
    }

    DPRINT(cout << "Bound on sum of eta_inter values:" << endl);
    cstPtr=ISL::EqFromNames(spacePtr,dict); // n_max = Sigma(n_i)

    DPRINT(p=isl_printer_print_constraint(p,cstPtr));
    DPRINT(cout << endl << endl);

    ilp=isl_basic_set_add_constraint(ilp,cstPtr);

    isl_printer_free(p);
    return ilp;
  }

  // BSetPtr StoragePartition::boundExpMods(BSetPtr ilp){
  //   isl_printer *p=isl_printer_to_file(GetContextPtr(),stdout);

  //   for(int stmtId=0; stmtId<nStmts; ++stmtId){ 
  //     for(unsigned i=0; i<_expModMaxNameVec.size(); ++i){
  //       ISL::NameCoeffMap dict;
  // 	dict[GetExpModName(stmtId,i)]=-1;
  // 	dict[_expModMaxNameVec[i]]=1;
  //       isl_constraint *cstPtr=ISL::EqFromNames(spacePtr,m);

  // 	DPRINT(cout << endl << endl);
  // 	DPRINT(p=isl_printer_print_constraint(p,cstPtr));
  // 	DPRINT(cout << endl);

  //       ilp=isl_basic_set_add_constraint(ilp,cstPtr);
  //     }
  //   }

  //   isl_printer_free(p);
  //   return ilp;
  // }

  BSetPtr StoragePartition::boundHyperplaneCoeffs(BSetPtr ilp){
    ISL::NameCoeffMap m;
    for(unsigned i=0; i<_coeffNameVec.size(); ++i){
      isl_space *spacePtr=isl_basic_set_get_space(ilp);
      m[_coeffNameVec[i]]=1; m["1"]=8;
      ilp=isl_basic_set_add_constraint(ilp,ISL::IneqFromNames(spacePtr,m));
      
      m.clear();
      m[_coeffNameVec[i]]=-1; m["1"]=8;
      ilp=isl_basic_set_add_constraint(ilp,ISL::IneqFromNames(spacePtr,m));
    }
    return ilp;
  }

  // flag ? (u.P + w) <= cP + c
  //      : (u.P + w) <= u_max.P + w_max
  BSetPtr StoragePartition::approxExpModConstraint(ConflictPoly &cPolyRef,
						   int nStmts,StringVec &expModNameVec,bool flag){
    BMapPtr polyPtr=cPolyRef.GetBMapPtr();
    DPRINT(cout << "\nAPPROX UBOUND CONSTRAINT\n");
    ISL::NameCoeffMap dict;
    ISL::NameCoeffMapByLiteral nameCoeffByLit;

    int nExpModVars=expModNameVec.size();
    DPRINT(cout << nExpModVars << endl);
    string lastName=expModNameVec[nExpModVars-1];
    // !!! Somashekar...what about destStmtId?

    if(flag)
      dict["1"]=BIGC;
    // else
    //   dict[_expModMaxNameVec[nExpModVars-1]]=1;
    dict[lastName]=-1;
    nameCoeffByLit["1"]=dict; 

    isl_space *polySpacePtr=isl_basic_map_get_space(polyPtr);
    ISL::NameCoeffMap map2;
    for(unsigned i=0; i<isl_space_dim(polySpacePtr,isl_dim_param); ++i){
      map2.clear();
      map2[expModNameVec[i]]=-1;
      if(flag)
        map2["1"]=BIGC;
      // else
      //   map2[_expModMaxNameVec[i]]=1;
      string name(isl_space_get_dim_name(polySpacePtr,isl_dim_param,i));
      nameCoeffByLit[name]=map2;
    }

    FarkasLemma lemma(GetContextPtr(),_searchSpacePtr,polyPtr);
    BSetPtr ilp=lemma.Apply(nameCoeffByLit);
    return ilp;
  }

  // flag ? (h.s - h.t) <= u.P + w
  //      :-(h.s - h.t) <= u.P + w
BSetPtr StoragePartition::expModConstraints(ConflictPoly &cPolyRef,bool flag,
					    int nStmts,StringVec &expModNameVec){
    BMapPtr polyPtr=cPolyRef.GetBMapPtr();
    DPRINT(cout << "\nUBOUND CONSTRAINT\n");
    ISL::NameCoeffMap dict;
    ISL::NameCoeffMapByLiteral nameCoeffByLit;

    int nExpModVars=expModNameVec.size();

    isl_space *polySpacePtr=isl_basic_map_get_space(polyPtr);
    for(unsigned i=0; i<isl_space_dim(polySpacePtr,isl_dim_param); ++i){
      dict.clear();
      dict[expModNameVec[i]]=1;
      string name(isl_space_get_dim_name(polySpacePtr,isl_dim_param,i));
      nameCoeffByLit[name]=dict;
    }

    int ndims=isl_space_dim(polySpacePtr,isl_dim_in);
    // the hyperplane coefficients for the dimensions
    for(unsigned i=0; i<ndims; ++i){
      dict.clear();
      dict[GetCoeffName(cPolyRef.GetSrcStmtId(),i)]=(flag ? -1 : 1);
      string name1(isl_space_get_dim_name(polySpacePtr,isl_dim_in,i));
      nameCoeffByLit[name1]=dict;

      dict.clear();
      dict[GetCoeffName(cPolyRef.GetDestStmtId(),i)]=(flag ? 1 : -1);
      string name2(isl_space_get_dim_name(polySpacePtr,isl_dim_out,i));
      nameCoeffByLit[name2]=dict;
    }

    // the hyperplane coefficient for the constant dimensions
    {
      dict.clear();
      if(cPolyRef.GetSrcStmtId()!=cPolyRef.GetDestStmtId()){
        dict[GetCoeffName(cPolyRef.GetSrcStmtId(),ndims)]=(flag ? -1 : 1); // from the lhs
        dict[GetCoeffName(cPolyRef.GetDestStmtId(),ndims)]=(flag ? 1 : -1); // fromthe lhs
      }

      dict[expModNameVec[nExpModVars-1]]=1;

      nameCoeffByLit["1"]=dict;
    }

    FarkasLemma lemma(GetContextPtr(),_searchSpacePtr,polyPtr);
    BSetPtr ilp=lemma.Apply(nameCoeffByLit);

    return ilp;
  }

  // (index % 2 == 0) ? (h.s - h.t) >= 1 - (1 - d1)(cP + c + 1)
  //                  : (h.s - h.t) <= -1 + (1 -d2)(cP + c + 1) 
  // that is...
  // (index % 2 == 0) ? (h.s - h.t) +cP +c -cPd1 -cd1 -d1 >= 0
  //                  : -(h.s - h.t) +cP +c -cPd2 -cd2 -d2 >= 0
  BSetPtr StoragePartition::decisionConstraints(ConflictPoly &cPolyRef,
                                                    unsigned index){
    BMapPtr polyPtr=cPolyRef.GetBMapPtr();
    DPRINT(cout << "\nDECISION CONSTRAINT\n");
    ISL::NameCoeffMap dict;
    ISL::NameCoeffMapByLiteral nameCoeffByLit;

    isl_space *polySpacePtr=isl_basic_map_get_space(polyPtr);
    for(unsigned i=0; i<isl_space_dim(polySpacePtr,isl_dim_param); ++i){
      dict.clear();
      dict["1"]=BIGC; dict[_decisionVarNameVec[index]]=-BIGC;
      string name(isl_space_get_dim_name(polySpacePtr,isl_dim_param,i));
      nameCoeffByLit[name]=dict;
    }

    int ndims=isl_space_dim(polySpacePtr,isl_dim_in);

    // the hyperplane coefficients for the dimensions
    for(unsigned i=0; i<ndims; ++i){
      dict.clear();
      dict[GetCoeffName(cPolyRef.GetSrcStmtId(),i)]=(((index % 2) == 0) ? 1 : -1);
      string name1(isl_space_get_dim_name(polySpacePtr,isl_dim_in,i));
      nameCoeffByLit[name1]=dict;

      dict.clear();
      dict[GetCoeffName(cPolyRef.GetDestStmtId(),i)]=(((index % 2) == 0) ? -1 : 1);
      string name2(isl_space_get_dim_name(polySpacePtr,isl_dim_out,i));
      nameCoeffByLit[name2]=dict;
    }

    // the hyperplane coefficient for the constant dimensions
    {
      dict.clear();
      if(cPolyRef.GetSrcStmtId()!=cPolyRef.GetDestStmtId()){
	dict[GetCoeffName(cPolyRef.GetSrcStmtId(),ndims)]=(((index % 2) == 0) ? 1 : -1);   // from the lhs
	dict[GetCoeffName(cPolyRef.GetDestStmtId(),ndims)]=(((index % 2) == 0) ? -1 : 1); // from the lhs
      }

      dict["1"]=-1 + BIGC + 1;
      dict[_decisionVarNameVec[index]]=-BIGC - 1;

      nameCoeffByLit["1"]=dict; 
    }

    FarkasLemma lemma(GetContextPtr(),_searchSpacePtr,polyPtr);
    BSetPtr ilp=lemma.Apply(nameCoeffByLit);
    return ilp;
  }

  // // flag ? u.P - u'P >= (1 - d1 - d2)
  // //      : u.P - u'P <= (1 - d1 - d2)(cP + c)
  // // that is...
  // // P(u - u') - (1 - d1 - d2) >= 0
  // // P(-u + u' + c -c*d1 - c*d2) + c(1 - d1 - d2) >= 0
  // BSetPtr StoragePartition::interStmtConflictSatConstraints(ConflictPoly &cPolyRef,int index,bool flag){
  //   BMapPtr polyPtr=cPolyRef.GetBMapPtr();
  //   DPRINT(cout << "\nINTER STATEMENT CONFLICT SATISFACTION CONSTRAINT\n");
  //   ISL::NameCoeffMap dict;
  //   ISL::NameCoeffMapByLiteral nameCoeffByLit;

  //   StringVec srcExpModNameVec=GetExpModInterNameVec(cPolyRef.GetSrcStmtId());
  //   StringVec destExpModNameVec=GetExpModInterNameVec(cPolyRef.GetDestStmtId());

  //   isl_space *polySpacePtr=isl_basic_map_get_space(polyPtr);

  //   // equate the multipliers of the parameters
  //   for(unsigned i=0; i<isl_space_dim(polySpacePtr,isl_dim_param); ++i){
  //     dict.clear();
  //     dict[srcExpModNameVec[i]]=(flag ? 1 : -1);
  //     dict[destExpModNameVec[i]]=(flag ? -1 : 1);

  //     if(!flag){
  //       dict["1"]=BIGC;
  //       dict[_decisionVarNameVec[2*index]]=-BIGC;
  //       dict[_decisionVarNameVec[2*index+1]]=-BIGC;
  //     }

  //     string name(isl_space_get_dim_name(polySpacePtr,isl_dim_param,i));
  //     nameCoeffByLit[name]=dict;
  //   }

  //   // equate the multipliers in the constant dimension
  //   {
  //     dict.clear();
  //     dict["1"]=1;
  //     dict[_decisionVarNameVec[2*index]]=-1;
  //     dict[_decisionVarNameVec[2*index+1]]=-1;

  //     dict["1"]*=(flag ? -1 : BIGC);
  //     dict[_decisionVarNameVec[2*index]]*=(flag ? -1 : BIGC);
  //     dict[_decisionVarNameVec[2*index+1]]*=(flag ? -1 : BIGC);

  //     nameCoeffByLit["1"]=dict;
  //   }

  //   FarkasLemma lemma(GetContextPtr(),_searchSpacePtr,polyPtr);
  //   BSetPtr ilp=lemma.Apply(nameCoeffByLit);

  //   return ilp;
  // }


  // flag ? Sigma(u'_i - u_i) >= (1 - d1 - d2)
  //      : Sigma(u'_i - u_i) <= (1 - d1 - d2)(c*nParams)
  BSetPtr StoragePartition::interStmtConflictSatConstraints(ConflictPoly &cPolyRef,
						  int index,bool flag,BSetPtr ilp,int stmtId){
    assert(stmtId==cPolyRef.GetSrcStmtId()
	   ||stmtId==cPolyRef.GetDestStmtId());
    isl_printer *p=isl_printer_to_file(GetContextPtr(),stdout);
    DPRINT(cout << "\nINTER STATEMENT CONFLICT SATISFACTION CONSTRAINT\n");

    BMapPtr polyPtr=cPolyRef.GetBMapPtr();
    isl_space *polySpacePtr=isl_basic_map_get_space(polyPtr);

    StringVec srcExpModNameVec=GetExpModInterNameVec(stmtId);
    StringVec destExpModNameVec=GetExpModNameVec(stmtId);

    ISL::NameCoeffMap dict;
    // the multipliers of the parameters
    int nParams=isl_space_dim(polySpacePtr,isl_dim_param);
    for(unsigned i=0; i<nParams; ++i){
      dict[srcExpModNameVec[i]]=(flag ? 1 : -1);
      dict[destExpModNameVec[i]]=(flag ? -1 : 1);

      if(!flag){
	// !!! this if block seems redundant
        dict["1"]=BIGC;
        dict[_decisionVarNameVec[2*index]]=-BIGC;
        dict[_decisionVarNameVec[2*index+1]]=-BIGC;
      }
    }

    dict["1"]=1;
    dict[_decisionVarNameVec[2*index]]=-1;
    dict[_decisionVarNameVec[2*index+1]]=-1;

    dict["1"]*=(flag ? -1 : BIGC);
    dict[_decisionVarNameVec[2*index]]*=(flag ? -1 : BIGC*nParams);
    dict[_decisionVarNameVec[2*index+1]]*=(flag ? -1 : BIGC*nParams);

    isl_constraint *cstPtr=ISL::IneqFromNames(isl_basic_set_get_space(ilp),dict);
    DPRINT(cout << endl; p=isl_printer_print_constraint(p,cstPtr); cout << endl);
    ilp=isl_basic_set_add_constraint(ilp,cstPtr);

    isl_printer_free(p);
    return ilp;
  }

  BSetPtr StoragePartition::miscModConstraints(BSetPtr ilp,int nStmts,
					       map<int,int> stmtPolyCountMap){

    isl_printer *p=isl_printer_to_file(GetContextPtr(),stdout);
    // let contraction modulo of stmts without polyhedra be equal to 1
    for(int i=0; i<nStmts; ++i){
      if(stmtPolyCountMap[i]==0){
	ISL::NameCoeffMap dict;
	StringVec expModNameVec=GetExpModNameVec(i);
	for(int index=0; index<expModNameVec.size(); ++index){
	  dict.clear();
	  // if(index!=expModNameVec.size()-1)
	    dict["1"]=0;
	  // else
	  //   dict["1"]=-1;
	  dict[expModNameVec[index]]=1;

          isl_constraint *cstPtr=ISL::IneqFromNames(isl_basic_set_get_space(ilp),dict);
          DPRINT(cout << endl; p=isl_printer_print_constraint(p,cstPtr); cout << endl);
          ilp=isl_basic_set_add_constraint(ilp,cstPtr);
	}
      }
    }

    // u.P - u'P >= 0 for each statement 
    for(int stmtId=0; stmtId<nStmts; ++stmtId){
      ISL::NameCoeffMap dict;
      ISL::NameCoeffMapByLiteral nameCoeffByLit;

      StringVec expModNameVec=GetExpModNameVec(stmtId);
      StringVec expModInterNameVec=GetExpModInterNameVec(stmtId);

      // >= constraint, dimension wise
      for(unsigned i=0; i<expModNameVec.size(); ++i){
	dict.clear();
	dict[expModInterNameVec[i]]=1;
	dict[expModNameVec[i]]=-1;

	isl_constraint *cstPtr=ISL::IneqFromNames(isl_basic_set_get_space(ilp),dict);
	DPRINT(cout << endl; p=isl_printer_print_constraint(p,cstPtr); cout << endl);
	ilp=isl_basic_set_add_constraint(ilp,cstPtr);
      }
    }
    
    isl_printer_free(p);
    return ilp;
  }

  BSetPtr StoragePartition::createBasicConstraintsCore(ConflictPoly &cPolyRef,
						       int nStmts,unsigned index){
    BSetPtr ilp=isl_basic_set_universe(isl_space_copy(_searchSpacePtr));


    if(cPolyRef.GetSrcStmtId()==cPolyRef.GetDestStmtId()){ // use the ubound of GetDestStmtId()
      StringVec destExpModNameVec=GetExpModNameVec(cPolyRef.GetDestStmtId());
      ilp=isl_basic_set_intersect(ilp,approxExpModConstraint(cPolyRef,nStmts,destExpModNameVec,true));
      ilp=isl_basic_set_intersect(ilp,expModConstraints(cPolyRef,true,nStmts,destExpModNameVec));
      ilp=isl_basic_set_intersect(ilp,expModConstraints(cPolyRef,false,nStmts,destExpModNameVec));
    }
    else{
      StringVec srcExpModNameVec=GetExpModInterNameVec(cPolyRef.GetSrcStmtId());
      StringVec destExpModNameVec=GetExpModInterNameVec(cPolyRef.GetDestStmtId());

      ilp=isl_basic_set_intersect(ilp,approxExpModConstraint(cPolyRef,nStmts,srcExpModNameVec,true));
      ilp=isl_basic_set_intersect(ilp,approxExpModConstraint(cPolyRef,nStmts,destExpModNameVec,true));

      ilp=isl_basic_set_intersect(ilp,expModConstraints(cPolyRef,true,nStmts,srcExpModNameVec));
      ilp=isl_basic_set_intersect(ilp,expModConstraints(cPolyRef,false,nStmts,srcExpModNameVec));
      ilp=isl_basic_set_intersect(ilp,expModConstraints(cPolyRef,true,nStmts,destExpModNameVec));
      ilp=isl_basic_set_intersect(ilp,expModConstraints(cPolyRef,false,nStmts,destExpModNameVec));

      ilp=interStmtConflictSatConstraints(cPolyRef,index,false,ilp,cPolyRef.GetSrcStmtId());
      ilp=interStmtConflictSatConstraints(cPolyRef,index,false,ilp,cPolyRef.GetDestStmtId());
    }

    ilp=isl_basic_set_intersect(ilp,decisionConstraints(cPolyRef,2*index));
    ilp=isl_basic_set_intersect(ilp,decisionConstraints(cPolyRef,2*index+1));

    isl_printer *p=isl_printer_to_file(GetContextPtr(),stdout);
    assert(ilp!=NULL);
    DPRINT(isl_printer_print_basic_set(p,ilp));
    DPRINT(cout << endl);

    return ilp;
  }

  __isl_give BSetPtr StoragePartition::constrainDecisionVars(__isl_take BSetPtr ilp,StringVec &varNameVec){
    ISL::NameCoeffMap dict;

    for(int index=0; index<varNameVec.size(); ++index){
      // d >= 0
      dict.clear();
      dict[varNameVec[index]]=1;
      ilp=isl_basic_set_add_constraint(ilp,
	       		     ISL::IneqFromNames(_searchSpacePtr,dict));

      // d <= 1
      dict.clear();
      dict[varNameVec[index]]=-1;
      dict["1"]=1;
      ilp=isl_basic_set_add_constraint(ilp,
				 ISL::IneqFromNames(_searchSpacePtr,dict));
    }
    return ilp;
  }

  __isl_give BSetPtr StoragePartition::constrainSatDecisionVars(__isl_take BSetPtr ilp){
    map<int,int> stmtToIntraPolyCountMap;
    ConflictPolyVec &conflictSet=_cSpecRef.GetConflictPolyVecRef();
    int nStmts=_cSpecRef.GetNStmts();
    isl_printer *p=isl_printer_to_file(GetContextPtr(),stdout);
    isl_constraint *cstPtr=NULL;

    for(int stmtId=0; stmtId<nStmts; ++stmtId){
      stmtToIntraPolyCountMap[stmtId]=0;
    }

    for(ConflictPolyVec::iterator i=conflictSet.begin(); i!=conflictSet.end(); ++i){
      if(i->GetSrcStmtId()==i->GetDestStmtId())
        stmtToIntraPolyCountMap[i->GetSrcStmtId()]+=1;
    }

    DPRINT(cout << "Inter-statement conflict satisfaction constraints" << endl);
    ISL::NameCoeffMap dict;
    for(int stmtId=0; stmtId<nStmts; ++stmtId){
      // n_intra >= d
      dict.clear();
      dict["1"]=0;
      dict[_satDecisionVarNameVec[stmtId]]=-1;
      dict[_etaIntraNameVec[stmtId]]=1;

      cstPtr=ISL::IneqFromNames(_searchSpacePtr,dict);
      DPRINT(p=isl_printer_print_constraint(p,cstPtr);cout << endl);
      ilp=isl_basic_set_add_constraint(ilp,cstPtr);

      // n_intra <= |CS_stmtId| * d
      dict.clear();
      dict["1"]=0;
      dict[_satDecisionVarNameVec[stmtId]]=stmtToIntraPolyCountMap[stmtId];
      dict[_etaIntraNameVec[stmtId]]=-1;

      cstPtr=ISL::IneqFromNames(_searchSpacePtr,dict);
      DPRINT(p=isl_printer_print_constraint(p,cstPtr);cout << endl);
      ilp=isl_basic_set_add_constraint(ilp,cstPtr);
    }

    // consciously satisfy inter-statement conflicts only when
    // intra-statement conflicts are already satisfied
    int index=0;
    for(ConflictPolyVec::iterator i=conflictSet.begin();
	i!=conflictSet.end(); ++i,++index){
      if(i->GetSrcStmtId()==i->GetDestStmtId())
	continue;

      dict.clear();
      dict["1"]=1;
      dict[_satDecisionVarNameVec[i->GetSrcStmtId()]]=-1;
      dict[_satDecisionVarNameVec[i->GetDestStmtId()]]=-1;
      dict[_decisionVarNameVec[2*index]]=-1;
      dict[_decisionVarNameVec[2*index+1]]=-1;

      cstPtr=ISL::IneqFromNames(_searchSpacePtr,dict);
      DPRINT(p=isl_printer_print_constraint(p,cstPtr);cout << endl);
      ilp=isl_basic_set_add_constraint(ilp,cstPtr);
    }
    return ilp;
  }

  void StoragePartition::FarkasLemma::equateCoeffsCore(ConstraintPtr cPtr,
			      isl_dim_type dtype, unsigned multIndex, bool isEq,
			      ISL::NameCoeffMapByLiteral &nameCoeffByLit){
    isl_space *spacePtr=isl_basic_map_get_space(_polyPtr);
    for(unsigned pos=0; pos<isl_space_dim(spacePtr,dtype); ++pos){
      string name=string(isl_space_get_dim_name(spacePtr,dtype,pos));

      ISL::NameCoeffMap m1;
      if(nameCoeffByLit.find(name)!=nameCoeffByLit.end()){
        m1=nameCoeffByLit[name];
      }
      isl_val *v=isl_constraint_get_coefficient_val(cPtr,dtype,pos);
      long num=isl_val_get_num_si(v);
      m1[_multNames[multIndex]]=(isEq?num:-num);
      // cout << name << " : " << "[" << multIndex << "]" << _multNames[multIndex] << ", " << (isEq?num:-num) << endl;
      nameCoeffByLit[name]=m1;
    }
  }

  BSetPtr StoragePartition::FarkasLemma::equateCoeffs(BSetPtr ilp,
			           isl_dim_type dtype, 
			           ISL::NameCoeffMapByLiteral &nameCoeffByLit){
    int index=0;
    for(ConstraintPtrVecIter i=_cstVec.begin(); i!=_cstVec.end(); ++i,++index){
      // Get coefficient of specific isl_dim_in dimension in the face
      // Some of these constraints might be equalities...still
      equateCoeffsCore(*i,dtype,index,false,nameCoeffByLit);
    }

    for(ConstraintPtrVecIter i=_eqCstVec.begin();
	i!=_eqCstVec.end(); ++i,++index){
      // get coefficient of specific isl_dim_in dimension in the face
      equateCoeffsCore(*i,dtype,index,true,nameCoeffByLit);
    }

    // create constraints for each name
    isl_space *spacePtr=isl_basic_map_get_space(_polyPtr);
    for(unsigned pos=0; pos<isl_space_dim(spacePtr,dtype); ++pos){
      string name=string(isl_space_get_dim_name(spacePtr,dtype,pos));
      if(nameCoeffByLit.find(name)!=nameCoeffByLit.end()){
	ISL::NameCoeffMap m=nameCoeffByLit.find(name)->second;
	ConstraintPtr cPtr=ISL::EqFromNames(isl_basic_set_get_space(ilp),m);

        isl_printer *p=isl_printer_to_file(GetContextPtr(),stdout);
        DPRINT(cout << endl << "Adding constraint due to " << name << endl);
	DPRINT(isl_printer_print_constraint(p,cPtr));
	DPRINT(cout << endl << endl);

	ilp=isl_basic_set_add_constraint(ilp,cPtr);
      }
    }

    return ilp;
  }

  void StoragePartition::FarkasLemma::equateCoeffsConstCore(ConstraintPtr cPtr,
			      int multIndex, bool isEq,
			      ISL::NameCoeffMapByLiteral &nameCoeffByLit){
    // constant values due to the face
    ISL::NameCoeffMap m2;
    if(nameCoeffByLit.find("1")!=nameCoeffByLit.end()){
      m2=nameCoeffByLit["1"];
    }
    isl_val *v=isl_constraint_get_constant_val(cPtr);
    long num=isl_val_get_num_si(v);
    m2[_multNames[multIndex]]=(isEq?num:-num);
    nameCoeffByLit["1"]=m2;
  }


  BSetPtr StoragePartition::FarkasLemma::equateCoeffsConst(BSetPtr ilp,
			      ISL::NameCoeffMapByLiteral &nameCoeffByLit){
    int index=0;
    for(ConstraintPtrVecIter i=_cstVec.begin()
	  ; i!=_cstVec.end(); ++i,++index){
      // Get coefficient of specific isl_dim_in dimension in the face
      // Some of these constraints might be equalities...still
      equateCoeffsConstCore(*i,index,false,nameCoeffByLit);
    }

    for(ConstraintPtrVecIter i=_eqCstVec.begin();
	i!=_eqCstVec.end(); ++i,++index){
      // get coefficient of specific isl_dim_in dimension in the face
      equateCoeffsConstCore(*i,index,true,nameCoeffByLit);
    }

    // constant values due to the constant Farkas multiplier
    ISL::NameCoeffMap m=nameCoeffByLit["1"];
    m[_multNames.back()]=-1;
    nameCoeffByLit["1"]=m;

    // create constraints due to the constants
    m=nameCoeffByLit.find("1")->second;
    ConstraintPtr cPtr=ISL::EqFromNames(isl_basic_set_get_space(ilp),m);

        isl_printer *p=isl_printer_to_file(GetContextPtr(),stdout);
        DPRINT(cout << endl << "Adding constraint due to equating constants" << endl);
	DPRINT(isl_printer_print_constraint(p,cPtr));
	DPRINT(cout << endl << endl);

    ilp=isl_basic_set_add_constraint(ilp,cPtr);

    return ilp;
  }

  //!!! this can be done once, before formulating the various constraints
  BSetPtr StoragePartition::FarkasLemma::getEliminationSpace(int nVars){
    int nMults=_cstVec.size() + _eqCstVec.size() + 1;
    _multNames=Namer::GetNames(nMults,Namer::kFarkasMult);

    BSetPtr ilp=isl_basic_set_universe(isl_space_copy(_searchSpacePtr));
    ilp=isl_basic_set_add_dims(ilp,isl_dim_set,nMults);

    int index=0;
    // a Farkas multiplier for each face
    for(ConstraintPtrVecIter i=_cstVec.begin();
	i!=_cstVec.end(); ++i,++index){
      ilp=isl_basic_set_set_dim_name(ilp,isl_dim_set,
                                     nVars + index,_multNames[index].c_str());
    }

    index=0;
    // a Farkas multiplier for each 'equality' face
    for(ConstraintPtrVecIter i=_eqCstVec.begin();
	i!=_eqCstVec.end(); ++i,++index){
      ilp=isl_basic_set_set_dim_name(ilp,isl_dim_set,
		        nVars + index + _cstVec.size(),
			_multNames[_cstVec.size()+index].c_str());
    }

    // a Farkas multiplier for the constant part in the affine 
    // combination of the faces
    ilp=isl_basic_set_set_dim_name(ilp,isl_dim_set,nVars + nMults - 1,
                                   _multNames[nMults-1].c_str());

    // all the added Farkas multipliers are non-negative
    ISL::NameCoeffMap dict;
    index=0;
    for(ConstraintPtrVecIter i=_cstVec.begin();
	i!=_cstVec.end(); ++i,++index){
      dict.clear(); dict[_multNames[index]]=1;
      ConstraintPtr c=ISL::IneqFromNames(isl_basic_set_get_space(ilp),dict);
      ilp=isl_basic_set_add_constraint(ilp,c);
    }
    index=0;
    for(ConstraintPtrVecIter i=_eqCstVec.begin();
	i!=_eqCstVec.end(); ++i,++index){
      dict.clear(); dict[_multNames[_cstVec.size() + index]]=1;
      ConstraintPtr c=ISL::IneqFromNames(isl_basic_set_get_space(ilp),dict);
      ilp=isl_basic_set_add_constraint(ilp,c);
    }
    dict.clear(); dict[_multNames[nMults-1]]=1;
    ConstraintPtr c=ISL::IneqFromNames(isl_basic_set_get_space(ilp),dict);
    ilp=isl_basic_set_add_constraint(ilp,c);

    DPRINT(cout << "Elimination Space: ");
    isl_printer *p=isl_printer_to_file(GetContextPtr(),stdout);
    DPRINT(isl_printer_print_basic_set(p,ilp));
    DPRINT(cout << endl);

    return ilp;
  }

  BSetPtr StoragePartition::FarkasLemma::Apply(ISL::NameCoeffMapByLiteral &nameCoeffByLit){
    int nVars=isl_space_dim(_searchSpacePtr,isl_dim_set);
    BSetPtr ilp=getEliminationSpace(nVars);
    ilp=equateCoeffs(ilp,isl_dim_param,nameCoeffByLit);
    ilp=equateCoeffs(ilp,isl_dim_out,nameCoeffByLit);
    ilp=equateCoeffs(ilp,isl_dim_in,nameCoeffByLit);
    ilp=equateCoeffsConst(ilp,nameCoeffByLit);

    DPRINT(cout << "Before Elimination: " << endl);
    isl_printer *p=isl_printer_to_file(GetContextPtr(),stdout);
    DPRINT(isl_printer_print_basic_set(p,ilp));
    DPRINT(cout << endl);

    unsigned nMults=_cstVec.size() + _eqCstVec.size() + 1;
    for(unsigned i=0; i<nMults; ++i){
      ilp=isl_basic_set_project_out(ilp,isl_dim_set,nVars,1);
    }

    DPRINT(cout << "After Elimination: " << endl);
    DPRINT(isl_printer_print_basic_set(p,ilp));
    DPRINT(cout << endl);

    ilp=isl_basic_set_remove_divs(ilp);
    return ilp;
  }


  // B=8nb^2
  // (flag) ? H.l + (1-d1)B - H.h >= 0
  //          H.l + (1-d1)B + H.h >=0
  //        : -H.l + (1-d2)B -H.h >= 0
  //          -H.l + (1-d2)B + H.h >=0
  // where l is a basis of the rkernel
  // and h is the hyperplane which we want to ignore
  BSetPtr StoragePartition::constrainSkew(__isl_take BSetPtr ilp,
				     __isl_keep isl_mat *rowMat,
			             __isl_keep isl_mat *rker,
				     int col, int index, bool flag){
    ISL::NameCoeffMap dict;

    int n=isl_mat_rows(rker);
    int ubound=8*n*(int)pow(BIGC,2);

    dict["1"]=ubound;
    dict[_decisionVarNameVec[index]]=-ubound;

    for(int i=0; i<2; ++i){
      for(int r=0; r<n; ++r){ // iterate over the rows
        // coefficient due to the column matrix in the kernel
        isl_val *val=isl_mat_get_element_val(rker,r,col);
        int num1=isl_val_get_num_si(val);
	num1=(flag ? num1 : -num1);

        // coefficient due to the given hyperplane (row matrix)
        val=isl_mat_get_element_val(rowMat,0,r);
        int num2=isl_val_get_num_si(val);

	if(0==i)
          dict[_coeffNameVec[r]]=num1+num2;
	else
          dict[_coeffNameVec[r]]=num1-num2;
      }

      ConstraintPtr c=ISL::IneqFromNames(isl_basic_set_get_space(ilp),dict);

      isl_printer *p=isl_printer_to_file(GetContextPtr(),stdout);
      DPRINT(cout << endl << "Constraining the skew" << endl);
      DPRINT(isl_printer_print_constraint(p,c));
      DPRINT(cout << endl << endl);

      ilp=isl_basic_set_add_constraint(ilp,c);
    }

    return ilp;
  }

  // 
  // (flag) ? -SUM((B^i)h.c) >= 1 - (1 - d1)(B^(n))
  //        :  SUM((B^i)h.c) >= 1 - d1(B^(n))
  BSetPtr StoragePartition::createLinearIndepConstraintsCore(
				     __isl_take BSetPtr ilp,
			             __isl_keep isl_mat *rker,
				     int index,
                                     bool flag){
    DPRINT(cout << "\nLINEAR INDEPENDENCE DECISION CONSTRAINT\n");
    ISL::NameCoeffMap dict;

    int n=isl_mat_rows(rker);
    int base=n*BIGC+1;

    dict["1"]=(flag ? (-1+(int)pow(base,n)) : -1);
    dict[_decisionVarNameVec[index]]=(flag ? (-(int)pow(base,n))
				           : (int)pow(base,n));

    for(int c=0; c<n-1; ++c){ // iterate over the columns
      int multiplier=(int)pow(base,c);
      for(int r=0; r<=n-1; ++r){ // over the rows, ignore the row for constant dimension
	isl_val *val=isl_mat_get_element_val(rker,r,c);
	int num=isl_val_get_num_si(val);
	num=num*multiplier;
	num=(flag ? -num : num);
	dict[_coeffNameVec[r]]=num;
      }
    }

    ConstraintPtr c=ISL::IneqFromNames(isl_basic_set_get_space(ilp),dict);

    isl_printer *p=isl_printer_to_file(GetContextPtr(),stdout);
    DPRINT(cout << endl << "Adding linear independence constraint" << endl);
    DPRINT(isl_printer_print_constraint(p,c));
    DPRINT(cout << endl << endl);

    ilp=isl_basic_set_add_constraint(ilp,c);

    return ilp;
  }

  BSetPtr StoragePartition::createLinearIndepConstraints(__isl_take BSetPtr ilp,
							 __isl_take isl_mat *rowMat){
    isl_mat *rker=isl_mat_right_kernel(isl_mat_copy(rowMat));
    DPRINT(cout << "Right kernel:" << endl);
    ISL::PrintMatrix(rker);

    ConflictPolyVec &conflictSet=_cSpecRef.GetConflictPolyVecRef();
    int starti=2*conflictSet.size();

    // !!! Somashekar...just use the number of columns here
    int ncols=isl_mat_cols(rker);
    for(int i=0; i<ncols; ++i){
      ilp=createLinearIndepConstraintsCore(ilp,rker,starti+i,true);
      ilp=createLinearIndepConstraintsCore(ilp,rker,starti+i,false);
    }

    // 1.constrain the skew of the linear indep hyperplane
    starti+=ncols;

    // 1a. first add the constraints on the components
    for(int i=0; i<ncols; ++i){
      ilp=constrainSkew(ilp,rowMat,rker,i,starti+2*i,true);
      ilp=constrainSkew(ilp,rowMat,rker,i,starti+2*i+1,false);
    }

    // 1b. add constraints on these decision variables
    ISL::NameCoeffMap dict;
    for(int i=0; i<ncols; ++i){
      dict[_decisionVarNameVec[starti+2*i]]=1;
      dict[_decisionVarNameVec[starti+2*i+1]]=1;
    }
    dict["1"]=-1;
    ConstraintPtr c=ISL::IneqFromNames(isl_basic_set_get_space(ilp),dict);

    isl_printer *p=isl_printer_to_file(GetContextPtr(),stdout);
    DPRINT(cout << endl << "Adding linear independence constraint" << endl);
    DPRINT(isl_printer_print_constraint(p,c));
    DPRINT(cout << endl << endl);

    ilp=isl_basic_set_add_constraint(ilp,c);

    return ilp;
  }

  BSetPtr StoragePartition::boundCoefficientsSum(__isl_take BSetPtr ilp, int nStmts){
    ISL::NameCoeffMap dict;
    isl_printer *p=isl_printer_to_file(GetContextPtr(),stdout);

    int ncoeffs=_coeffNameVec.size()/nStmts;
    for(int stmtId=0; stmtId<nStmts; ++stmtId){// foreach stmt, 
      DPRINT(cout << endl << "Bound on sum of coefficients, stmtId=" << stmtId << endl);

      // 1. first place a lower bound of 0 on the sum
      dict.clear();
      dict[_coeffMaxVec[stmtId]]=1;
      ConstraintPtr c=ISL::IneqFromNames(isl_basic_set_get_space(ilp),dict);

      DPRINT(isl_printer_print_constraint(p,c));
      DPRINT(cout << endl << endl);

      ilp=isl_basic_set_add_constraint(ilp,c);
    
      // 2. add constraints :|c0|+|c1|+|c2|.. |c_n|<=c3
      dict.clear();
      dict[_coeffMaxVec[stmtId]]=1;
      unsigned max=(unsigned)pow(2,(unsigned)ncoeffs)-1;
      for(unsigned bpattern=0; bpattern<=max; ++bpattern){
	unsigned bits=bpattern;
	for(int i=0; i<ncoeffs; ++i){
	  unsigned lastBit=1&bits;
	  if(lastBit)
	    dict[GetCoeffName(stmtId,i)]=1;
	  else
	    dict[GetCoeffName(stmtId,i)]=-1;
	  bits=bits>>1;
	}
	c=ISL::IneqFromNames(isl_basic_set_get_space(ilp),dict);

	DPRINT(cout << bpattern << " : ");
	DPRINT(isl_printer_print_constraint(p,c));
	DPRINT(cout << endl << endl);

	ilp=isl_basic_set_add_constraint(ilp,c);
      }
      // 3. finally, place an upper bound on the sum
      dict.clear();
      dict[_coeffMaxVec[stmtId]]=-1;
      dict["1"]=ncoeffs*BIGC;

      c=ISL::IneqFromNames(isl_basic_set_get_space(ilp),dict);
      DPRINT(isl_printer_print_constraint(p,c));
      DPRINT(cout << endl << endl);

      ilp=isl_basic_set_add_constraint(ilp,c);
    }
    return ilp;
  }
  
}
