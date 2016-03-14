#include "ConflictSpec.h"
#include "ISLUtils.h"
#include "Namer.h"

#include <fstream>
#include <sstream>

namespace Smo{
  // get the next line in this input file stream
  string NextLine(ifstream &ifs){
    char str[1024];
    ifs.getline(str,1024);
    assert((ifs.rdstate() & ifstream::failbit)==0);

    while(str[0] == '#' || str[0] == '\n' || string(str).empty()){
      DPRINT(cout << string(str) << endl);
      ifs.getline(str,1024);
      assert((ifs.rdstate() & ifstream::failbit)==0);
    }

    DPRINT(cout << string(str) << endl);
    return string(str);
  }

  ConflictSpec::ConflictSpec(struct isl_ctx *ctx,int nStmts,int nParams,int nInputs,
                             int nOutputs,ConflictPolyVec &cPolyVecRef)
    :_ctx(ctx),
     _nStmts(nStmts),
     _nParams(nParams),
     _nInputs(nInputs),
     _nOutputs(nOutputs),
     _cPolyVecRef(cPolyVecRef)
  {
    _printer=isl_printer_to_file(ctx,stdout);
  }

  ConflictSpec::~ConflictSpec(){
    isl_ctx_free(_ctx);
    ConflictPolyVec::iterator i=_cPolyVecRef.begin();
    for(; i!=_cPolyVecRef.end(); ++i){
      isl_basic_map_free(i->GetBMapPtr());
    }
    _cPolyVecRef.clear();
    delete &_cPolyVecRef;
  }

  ConflictPolyVec &ConflictSpec::GetConflictPolyVecRef() const{
    return _cPolyVecRef;
  }

  bool ConflictSpec::IsEmpty(){
    ConflictPolyVec::iterator i=GetConflictPolyVecRef().begin();
    for(; i!=GetConflictPolyVecRef().end(); ++i){
      if(!isl_basic_map_is_empty(i->GetBMapPtr()))
	return false;
    }
    return true;
  }

  ConflictSpec &ConflictSpec::Copy(const ConflictSpec &specToCopyRef){
    struct isl_ctx *ctx=specToCopyRef.GetContextPtr();
    int nStmts=specToCopyRef.GetNStmts();
    int nParams=specToCopyRef.GetNParams();
    int nInputs=specToCopyRef.GetNInputs();
    int nOutputs=specToCopyRef.GetNOutputs();
    ConflictPolyVec *cPolyVecPtr=new ConflictPolyVec();
    for(ConflictPolyVec::iterator i=specToCopyRef.GetConflictPolyVecRef().begin();
	i!=specToCopyRef.GetConflictPolyVecRef().end(); ++i){
      BMapPtr bMapPtr=i->GetBMapPtr();
      assert(bMapPtr!=NULL);
      ConflictPoly cPoly(i->GetSrcStmtId(),i->GetDestStmtId(),
			     isl_basic_map_copy(bMapPtr));
      cPolyVecPtr->push_back(cPoly);
    }
    return *(new ConflictSpec(ctx,nStmts,nParams,nInputs,
                                 nOutputs,*cPolyVecPtr));
  }

  ConflictSpec &ConflictSpec::Create(ifstream &ifs,int nStmts,int nConflictPoly,
				     isl_ctx *ctx,BMapPtr orderPairPtr){
    // struct isl_ctx *ctx=isl_ctx_alloc();
    // assert(ctx!=NULL);

    // // read the domain of the conflict set relations
    // ifstream ifs;
    // ifs.open(filename.c_str(),ios::in);
    // assert(ifs.is_open() && (ifs.rdstate() & ifstream::failbit) == 0);

    // BSetPtr domain=isl_basic_set_read_from_str(ctx,NextLine(ifs).c_str());
    
    // int nParams=isl_space_dim(isl_basic_set_get_space(domain),isl_dim_param);
    // int nVars=isl_space_dim(isl_basic_set_get_space(domain),isl_dim_set);
    // StringVec paramNameVec=Namer::GetNames(nParams,Namer::kParam);
    // StringVec varNameVec=Namer::GetNames(nVars,Namer::kInput);
    // isl_space *space=isl_space_set_alloc(ctx,nParams,nVars);
    // space=ISL::isl_space_set_dim_names(space,isl_dim_param,paramNameVec);
    // space=ISL::isl_space_set_dim_names(space,isl_dim_set,varNameVec);

    // isl_space *copySpacePtr=isl_space_copy(space);
    // domain=(BSetPtr)isl_set_reset_space((isl_set *)domain,copySpacePtr);

    // // read the image of the conflict set relations
    // BSetPtr image=isl_basic_set_read_from_str(ctx,NextLine(ifs).c_str());
    // image=(BSetPtr )isl_set_reset_space((isl_set *)image,space);

    // BMapPtr orderPairPtr=isl_basic_map_from_domain_and_range(domain,image);

    // int nConflictPoly=atoi(NextLine(ifs).c_str());

    // // !!! ideally,the order pair should be constructed using a range
    // // with these output names
    // string outputNames=NextLine(ifs);
    
    ConflictPolyVec *cPolyVecPtr=new ConflictPolyVec;

    int nParams=isl_space_dim(isl_basic_map_get_space(orderPairPtr),isl_dim_param);
    int nVars=isl_space_dim(isl_basic_map_get_space(orderPairPtr),isl_dim_in);

    for(int i=0; i<nConflictPoly; ++i){
      stringstream stream(Smo::NextLine(ifs));
      int srcStmtId,destStmtId;
      stream >> srcStmtId;
      stream >> destStmtId;

      BMapPtr poly=isl_basic_map_read_from_str(ctx,NextLine(ifs).c_str());
      assert(poly!=NULL);

      poly=isl_basic_map_intersect(poly,isl_basic_map_copy(orderPairPtr));
      cPolyVecPtr->push_back(ConflictPoly(srcStmtId,destStmtId,poly));
    }

    // print all conflict polyhedra
    int index=1;
    for(ConflictPolyVec::iterator iter=cPolyVecPtr->begin();
	iter!=cPolyVecPtr->end(); ++iter){
      DPRINT(cout << "Conflict Polyhedron #" << index++ << endl);
      isl_printer *p=isl_printer_to_file(ctx,stdout);
      assert(p!=NULL);
      DPRINT(isl_printer_print_basic_map(p,iter->GetBMapPtr()));
      DPRINT(cout << endl);
    }

    // isl_basic_map_free(orderPairPtr);
    // ifs.close();

    return *(new ConflictSpec(ctx,nStmts,nParams,nVars,nVars,*cPolyVecPtr));
  }

  ConflictSpec &ConflictSpec::Create(__isl_keep isl_ctx *ctx,
				     __isl_take BMapPtr bMapPtr,
				     vector<int> &uspan,
				     bool inclusive,
				     BSetPtrVec &liveOutBSetPtrVec,
				     int tileStartDim,
				     int tileEndDim){
    isl_printer *p=isl_printer_to_file(ctx,stdout);

    int nParams=isl_space_dim(isl_basic_map_get_space(bMapPtr),isl_dim_param);
    int nVars=isl_space_dim(isl_basic_map_get_space(bMapPtr),isl_dim_set);

    ConflictPolyVec *cPolyVecPtr=new ConflictPolyVec();

    // the conflicts due to the max utility span
    DPRINT(cout << endl << "Creating conflicts due to the dependences" << endl << endl);
    DPRINT(p=isl_printer_print_basic_map(p,bMapPtr));
    DPRINT(cout << endl);
    BMapPtr bMapCopyPtr=isl_basic_map_copy(bMapPtr);
    formulateConflictSets(ctx,bMapCopyPtr,*cPolyVecPtr,uspan,inclusive,tileStartDim,tileStartDim);

    // the conflicts due to the live out
    for(BSetPtrVecIter it=liveOutBSetPtrVec.begin();
    	it!=liveOutBSetPtrVec.end(); ++it){
      DPRINT(cout << endl << "Creating conflicts due to the live out" << endl << endl);
      DPRINT(p=isl_printer_print_basic_set(p,*it));
      DPRINT(cout << endl << endl);

      // create the basic relation for the live out conflicts
      BSetPtr range=isl_basic_map_range(isl_basic_map_copy(bMapPtr));
      BMapPtr baseRelation=isl_basic_map_from_domain_and_range(*it,range);
      isl_space *space=isl_basic_map_get_space(baseRelation);

      // add constraints to consider only the intra-tile conflicts
      for(int i=0; i<tileStartDim; ++i){
    	ISL::NameCoeffMap dict;
    	dict[string(isl_space_get_dim_name(space,isl_dim_out,i))]=1;
    	dict[string(isl_space_get_dim_name(space,isl_dim_in,i))]=-1;
    	dict["1"]=0;
    	baseRelation=isl_basic_map_add_constraint(baseRelation,ISL::EqFromNames(space,dict));
      }

      // create the lexicogrpahically less than relation for this space
      space=isl_map_get_space(isl_map_from_basic_map(baseRelation));
      isl_map *ltmapPtr=isl_map_lex_lt_first(space,nVars);
      BMapPtrVec bMapPtrVec;
      isl_map_foreach_basic_map(ltmapPtr,&ISL::AddAnyBasicMap,&bMapPtrVec);

      // intersect with the basic relation for liveout conflicts
      for(BMapPtrVecIter mapIter=bMapPtrVec.begin();
    	  mapIter!=bMapPtrVec.end(); ++mapIter){
    	*mapIter=isl_basic_map_intersect(isl_basic_map_copy(baseRelation),*mapIter);

	// isl_map *mapPtr=isl_map_from_basic_map(isl_basic_map_copy(*mapIter));

	// // subtract any conflicts which are already represented
	// // in the conflict polytopes inferred due to the dependences
        // for(BMapPtrVecIter cPolyIt=cPolyVecPtr->begin();
	//     cPolyIt!=cPolyVecPtr->end(); ++cPolyIt){
	//   mapPtr=isl_map_subtract(mapPtr,isl_map_from_basic_map(isl_basic_map_copy(*cPolyIt)));

	//   DPRINT(cout << "Before subtraction..." << endl);
    	//   DPRINT(p=isl_printer_print_basic_map(p,*mapIter));
    	//   DPRINT(cout << endl);
	//   DPRINT(cout << "After subtraction..." << endl);
    	//   p=isl_printer_print_basic_map(p,*mapIter);
    	//   DPRINT(cout << endl << endl);
	// }

    	if(!isl_basic_map_is_empty(*mapIter)){
	  // !!! Somashekar... add arguments for srcId and destId
    	  cPolyVecPtr->push_back(ConflictPoly(0,0,*mapIter));
    	  DPRINT(p=isl_printer_print_basic_map(p,*mapIter));
	  DPRINT(cout << endl);
    	}
      }

      isl_basic_map_free(baseRelation);
    }

    // !!! Somashekar...must coalesce the conflict polytopes here
    Coalesce(*cPolyVecPtr,1,1);

    // cout << endl << "The conflict sets are as follows:" << endl;
    // ISL::Print(ctx,*cPolyVecPtr);

    isl_printer_free(p);
    isl_basic_map_free(bMapPtr);

    return *(new ConflictSpec(ctx,1,nParams,nVars,nVars,*cPolyVecPtr));
  }

  void ConflictSpec::Destroy(){
    this->~ConflictSpec();
  }

// the given map is split into at most three different conflict sets (some may be empty)
// for example, given the space [y,x]->[y',x'] and utility span (2,1),
// we first obtain three different conflict sets by adding the constraints
// y-y'=2; y-y'>0 and y-y' <2; y-y'=0 ; The first of these is then split
// again based on the remaining components of the utility span
  void ConflictSpec::formulateConflictSets(__isl_keep isl_ctx *ctx,
					   Smo::BMapPtr &rootBMapPtr, ConflictPolyVec &cPolyVecRef,
					   vector<int> &uspan,bool inclusive,
					   int currDim,int tileStartDim){

    isl_printer *p=isl_printer_to_file(ctx,stdout);
    
    isl_space *spacePtr=isl_basic_map_get_space(rootBMapPtr);
    int ndims=isl_space_dim(spacePtr,isl_dim_set);

    if(currDim>=ndims)
      return;

    // the first conflict set
    Smo::BMapPtr cset1=isl_basic_map_copy(rootBMapPtr);
    Smo::ISL::NameCoeffMap dict;

    
    if(currDim>=uspan.size()+tileStartDim){
      // when inner most common time dimension is not at full depth
      dict[isl_space_get_dim_name(spacePtr,isl_dim_set,currDim)]=1;
      dict[isl_space_get_dim_name(spacePtr,isl_dim_in,currDim)]=-1;
      dict["1"]=-1;
      cset1=isl_basic_map_add_constraint(cset1,
			Smo::ISL::IneqFromNames(isl_basic_map_get_space(cset1),dict));
      cPolyVecRef.push_back(ConflictPoly(0,0,cset1)); // !!! intra-array conflicts

      DPRINT(p=isl_printer_print_basic_map(p,cset1));
      DPRINT(cout << endl);

      if(currDim<ndims-1){
	// the second conflict set
	Smo::BMapPtr cset2=rootBMapPtr;
	dict.clear();
	dict[isl_space_get_dim_name(spacePtr,isl_dim_set,currDim)]=1;
	dict[isl_space_get_dim_name(spacePtr,isl_dim_in,currDim)]=-1;
	dict["1"]=0;
	cset2=isl_basic_map_add_constraint(cset2,
					   Smo::ISL::EqFromNames(isl_basic_map_get_space(cset2),dict));
        cPolyVecRef.push_back(ConflictPoly(0,0,cset2));

        DPRINT(p=isl_printer_print_basic_map(p,cset2));
        DPRINT(cout << endl);
      }
    }
    //else if(uspan[currDim-tileStartDim]!=0){
    else {
      // !!! must assert somewhere that uspan is not a zero vector
      // otherwise, we will get a basic map universe as the conflic set
      dict[isl_space_get_dim_name(spacePtr,isl_dim_out,currDim)]=1;
      dict[isl_space_get_dim_name(spacePtr,isl_dim_in,currDim)]=-1;
      dict["1"]=-uspan[currDim-tileStartDim];
	
      // if(currDim+1==uspan.size()+tileStartDim){ 
      // 	isl_basic_map_free(cset1);
      // 	rootBMapPtr=isl_basic_map_add_constraint(rootBMapPtr,
      // 			Smo::ISL::IneqFromNames(isl_basic_map_get_space(rootBMapPtr),dict));
      // 	if(currDim==0){
      // 	  bMapPtrVec.push_back(rootBMapPtr);
      //     DPRINT(p=isl_printer_print_basic_map(p,cset1));
      //     DPRINT(cout << endl);
      // 	}
      // }
      // else
	{ // if there are more components remaining in the uspan
	cset1=isl_basic_map_add_constraint(cset1,
			Smo::ISL::EqFromNames(isl_basic_map_get_space(cset1),dict));
	// the second conflict set
	Smo::BMapPtr cset2=isl_basic_map_copy(rootBMapPtr);

	// upper bound on the difference
	dict.clear();
	dict[isl_space_get_dim_name(spacePtr,isl_dim_out,currDim)]=-1;
	dict[isl_space_get_dim_name(spacePtr,isl_dim_in,currDim)]=1;
	dict["1"]=uspan[currDim-tileStartDim]-1;
	cset2=isl_basic_map_add_constraint(cset2,
			Smo::ISL::IneqFromNames(isl_basic_map_get_space(cset2),dict));

	if(currDim==tileStartDim){
	  // lower bound on the difference is needed only for the tile start dim
	  dict.clear();
	  dict[isl_space_get_dim_name(spacePtr,isl_dim_out,currDim)]=1;
	  dict[isl_space_get_dim_name(spacePtr,isl_dim_in,currDim)]=-1;
	  dict["1"]=-1;
	  cset2=isl_basic_map_add_constraint(cset2,
					     Smo::ISL::EqFromNames(isl_basic_map_get_space(cset2),dict));
	  // the third conflict set - also needed only for the tile start dim
	  Smo::BMapPtr cset3=rootBMapPtr;
	  dict.clear();
	  dict[isl_space_get_dim_name(spacePtr,isl_dim_out,currDim)]=1;
	  dict[isl_space_get_dim_name(spacePtr,isl_dim_in,currDim)]=-1;
	  dict["1"]=0;
	  cset3=isl_basic_map_add_constraint(cset3,
					     Smo::ISL::EqFromNames(isl_basic_map_get_space(cset3),dict));
	  dict.clear();
	  dict[isl_space_get_dim_name(spacePtr,isl_dim_out,currDim+1)]=1;
	  dict[isl_space_get_dim_name(spacePtr,isl_dim_in,currDim+1)]=-1;
	  dict["1"]=-1;
	  cset3=isl_basic_map_add_constraint(cset3,
					     Smo::ISL::IneqFromNames(isl_basic_map_get_space(cset3),dict));

	  cPolyVecRef.push_back(ConflictPoly(0,0,cset3));
	  DPRINT(p=isl_printer_print_basic_map(p,cset3));
	  DPRINT(cout << " here" << endl);
	}

	cPolyVecRef.push_back(ConflictPoly(0,0,cset2));
	DPRINT(p=isl_printer_print_basic_map(p,cset2));
	DPRINT(cout << " here" << endl);

	formulateConflictSets(ctx,cset1,cPolyVecRef,uspan,inclusive,currDim+1,tileStartDim);
	if(inclusive && currDim==ndims-1){ // only push this for the last dimension
	  cPolyVecRef.push_back(ConflictPoly(0,0,cset1));
	  DPRINT(p=isl_printer_print_basic_map(p,cset1));
	  DPRINT(cout << " here" << endl);
	}
      }
    }
    // else{
    //   cout << "here" << endl;
    //   // add an equality constraint involving the curr dimensions
    //   dict[isl_space_get_dim_name(spacePtr,isl_dim_set,currDim)]=1;
    //   dict[isl_space_get_dim_name(spacePtr,isl_dim_in,currDim)]=-1;
    //   dict["1"]=0;
    //   rootBMapPtr=isl_basic_map_add_constraint(rootBMapPtr,
    // 					   Smo::ISL::EqFromNames(isl_basic_map_get_space(rootBMapPtr),dict));
    //   formulateConflictSets(ctx,rootBMapPtr,bMapPtrVec,uspan,currDim+1,tileStartDim);
    // }
  }

  void ConflictSpec::Coalesce(ConflictPolyVec &cPolyVecRef,int srcStmtId,int destStmtId){
    if(cPolyVecRef.empty())
      return;

    ConflictPolyVec excludedPolyVec;
    isl_union_map *unionMap=NULL;

    // only consider the polytopes with the given srcStmtId and destStmtId for coalescing.
    for(ConflictPolyVec::iterator it=cPolyVecRef.begin();
	it!=cPolyVecRef.end(); ++it){
      if(srcStmtId==it->GetSrcStmtId() && destStmtId==it->GetDestStmtId()){
	if(unionMap==NULL)
	  unionMap=isl_union_map_from_basic_map(it->GetBMapPtr());
	else
	  unionMap=isl_union_map_union(unionMap,isl_union_map_from_basic_map(it->GetBMapPtr()));
      }
      else
	excludedPolyVec.push_back(*it);
    }
    unionMap=isl_union_map_coalesce(unionMap);

    BMapPtrVec bMapPtrVec;
    ISL::ExtractAllBasicMap(unionMap,&bMapPtrVec);

    cPolyVecRef.clear();
    for(BMapPtrVecIter i=bMapPtrVec.begin();i!=bMapPtrVec.end(); ++i){
      ConflictPoly cPoly=ConflictPoly(srcStmtId,destStmtId,*i);
      cPolyVecRef.push_back(cPoly);
    }
    // add back the excluded polytopes
    for(ConflictPolyVec::iterator i=excludedPolyVec.begin(); i!=excludedPolyVec.end(); ++i){
      cPolyVecRef.push_back(*i);
    }
  }

}

