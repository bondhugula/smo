#include "ConflictSpecBuilder.h"
#include "ISLUtils.h"
#include "Namer.h"

#include <pet.h>
#include <isl/flow.h>
#include <isl/constraint.h>

#include <algorithm>

// !!! find a better place for this global
// additional parameters introduced in place of huge constants
map<long,string> gParams; 

namespace Smo{
  ConflictSpecBuilder::ConflictSpecBuilder(string filename){
    struct pet_scop *scop;

    // extract the scop using pet
    _ctx=isl_ctx_alloc();
    scop=pet_scop_extract_from_C_source(_ctx,filename.c_str(),NULL);

    if(!scop){
      cout << "No scop extracted or some error while extracting the scop" << endl;
      isl_ctx_free(_ctx);
    }

    // compute the dependences using isl
    isl_space *space=isl_set_get_space(scop->context);

    isl_union_set *domains=pet_scop_collect_domains(scop);
    isl_union_map *writes=pet_scop_collect_may_writes(scop);
    isl_union_map *reads=pet_scop_collect_may_reads(scop);
    isl_union_map *schedule=pet_scop_collect_schedule(scop);
    isl_union_map *empty=isl_union_map_empty(isl_space_copy(space));

    isl_union_map *rawdeps=NULL;
    isl_union_map_compute_flow(isl_union_map_copy(reads),
			       isl_union_map_copy(writes),
			       isl_union_map_copy(empty),
			       isl_union_map_copy(schedule),
			       &rawdeps, NULL, NULL, NULL);
    isl_union_map_coalesce(rawdeps);

    ISL::ExtractAllBasicSet(domains,&_domPtrVec);
    ISL::ExtractAllBasicMap(rawdeps,&_rawdepPtrVec);
    ISL::ExtractAllBasicMap(writes,&_writePtrVec);
    ISL::ExtractAllBasicMap(reads,&_readPtrVec);
    ISL::ExtractAllBasicMap(schedule,&_schedulePtrVec);

    // // ensure that the raw deps are for the given iteration spaces
    // for(BMapPtrVecIter it=_rawdepPtrVec.begin();
    // 	it!=_rawdepPtrVec.end(); ++it){
    //   BMapPtr rawdepPtr=*it;
    //   int src=atoi(isl_basic_map_get_tuple_name(rawdepPtr,isl_dim_in)+2);
    //   int dest=atoi(isl_basic_map_get_tuple_name(rawdepPtr,isl_dim_out)+2);
    //   BSetPtr dom=isl_basic_set_copy(getDomain(src));
    //   BSetPtr range=isl_basic_set_copy(getDomain(dest));
    //   cout << dest << src << endl;
    //   BMapPtr orderPairPtr=isl_basic_map_from_domain_and_range(dom,range);
    //   isl_printer *p=isl_printer_to_file(_ctx,stdout);
    //   p=isl_printer_print_basic_map(p,rawdepPtr);
    //   cout << endl;
    //   p=isl_printer_print_basic_map(p,orderPairPtr);
    //   cout << endl << endl;
    //   *it=isl_basic_map_intersect(rawdepPtr,orderPairPtr);
    //   isl_printer_free(p);
    // }

    parameterizeStmtDomains();

    PrintProgramInfo();
  }

  void ConflictSpecBuilder::PrintProgramInfo(){
    cout << endl << "The domains:" << endl;
    ISL::Print(_ctx, _domPtrVec);

    cout << endl << "The writes:" << endl;
    ISL::Print(_ctx, _writePtrVec);

    cout << endl << "The reads:" << endl;
    ISL::Print(_ctx, _readPtrVec);

    cout << endl << "The schedule:" << endl;
    ISL::Print(_ctx, _schedulePtrVec);

    cout << endl << "The RAW dependences:" << endl;
    ISL::Print(_ctx, _rawdepPtrVec);
  }


  int ConflictSpecBuilder::findInnerMostCommonTimeDim(BMapPtr rawDepPtr){
    // !!! assume here that we won't be dealing with a 100-d loop nest!
    return findInnerMostCommonTimeDim(rawDepPtr,100);
  }

  //!!! using the schedule, we need to find the appropriate depth up
  // to which we need to find the utility span but for now we just
  // assume that it is for the full depth....the following is only a
  // temporary solution
  int ConflictSpecBuilder::findInnerMostCommonTimeDim(BMapPtr rawDepPtr,int threshold){
    int src=atoi(isl_basic_map_get_tuple_name(rawDepPtr,isl_dim_in)+2);
    int dest=atoi(isl_basic_map_get_tuple_name(rawDepPtr,isl_dim_out)+2);

    BMapPtr bmap1Ptr=_schedulePtrVec[src];
    BMapPtr bmap2Ptr=_schedulePtrVec[dest];

    ConstraintPtrVec cstPtrVec1,cstPtrVec2;
    isl_basic_map_foreach_constraint(_schedulePtrVec[src],&ISL::AddAnyConstraint,&cstPtrVec1);
    isl_basic_map_foreach_constraint(_schedulePtrVec[dest],&ISL::AddAnyConstraint,&cstPtrVec2);

    isl_space *space1=isl_basic_map_get_space(_schedulePtrVec[src]);
    isl_space *space2=isl_basic_map_get_space(_schedulePtrVec[dest]);

    StringVec inVec1,outVec1,paramVec1;
    StringVec inVec2,outVec2,paramVec2;
    ISL::isl_space_get_dim_names(space1,isl_dim_in,inVec1);
    ISL::isl_space_get_dim_names(space1,isl_dim_out,outVec1);
    ISL::isl_space_get_dim_names(space1,isl_dim_param,paramVec1);
    ISL::isl_space_get_dim_names(space2,isl_dim_in,inVec2);
    ISL::isl_space_get_dim_names(space2,isl_dim_out,outVec2);
    ISL::isl_space_get_dim_names(space2,isl_dim_param,paramVec2);

    int depth=0;
    for(ConstraintPtrVec::reverse_iterator i1=cstPtrVec1.rbegin(),i2=cstPtrVec2.rbegin();
	i1!=cstPtrVec1.rend() && i2!=cstPtrVec2.rend() && depth<threshold; ++i1,++i2){
      ISL::NameCoeffMap dict1,dict2;

      ISL::GetCoefficientVals(*i1,space1,inVec1,isl_dim_in,inVec1,dict1);
      ISL::GetCoefficientVals(*i1,space1,outVec1,isl_dim_out,outVec1,dict1);
      ISL::GetCoefficientVals(*i1,space1,paramVec1,isl_dim_param,paramVec1,dict1);

      ISL::GetCoefficientVals(*i2,space2,inVec2,isl_dim_in,inVec2,dict2);
      ISL::GetCoefficientVals(*i2,space2,outVec2,isl_dim_out,outVec2,dict2);
      ISL::GetCoefficientVals(*i2,space2,paramVec2,isl_dim_param,paramVec2,dict2);

      if(equal(dict1.begin(),dict1.end(),dict2.begin())){
	bool increment=false;
	for(StringVecIter s1=inVec1.begin();
	    !increment && s1!=inVec1.end() ; ++s1){
	  increment=(dict1[*s1]!=0);
	}
	for(StringVecIter s1=paramVec1.begin();
	    !increment && s1!=paramVec1.end() ; ++s1){
	  increment=(dict1[*s1]!=0);
	}
	if(increment)
	  ++depth;
      }
      else 
	break;
    }

    cout << "The innermost common time dimension is at depth = " << depth << endl;
    return depth;
  }

  bool ConflictSpecBuilder::FindMaxUtilitySpan(int stmtId,int tileStartDim,
					       int tileEndDim,vector<int> &coefficients){
    isl_set *lexmax=NULL;
    isl_printer *p=isl_printer_to_file(_ctx,stdout);
    bool inclusive=true;

    cout << endl << "FINDING MAX UTILITY SPAN" << endl << endl;
    for(BMapPtrVecIter depIter=_rawdepPtrVec.begin();
	depIter!=_rawdepPtrVec.end(); ++depIter){
      BMapPtr rawDepPtr=isl_basic_map_copy(*depIter);

      int src=atoi(isl_basic_map_get_tuple_name(rawDepPtr,isl_dim_in)+2);
      int dest=atoi(isl_basic_map_get_tuple_name(rawDepPtr,isl_dim_out)+2);
      if(src!=stmtId) continue;

      cout << endl << "Examining the following raw dependence" << endl;
      isl_printer_print_basic_map(p,rawDepPtr);
      cout << endl;

      int depth=findInnerMostCommonTimeDim(rawDepPtr);

      // project out the inner variables which are at a depth that
      // does not matter for finding the utility span
      isl_space *space=isl_basic_map_get_space(rawDepPtr);
      rawDepPtr=isl_basic_map_project_out(rawDepPtr,isl_dim_in,
					  depth,isl_space_dim(space,isl_dim_in)-depth);
      rawDepPtr=isl_basic_map_project_out(rawDepPtr,isl_dim_out,
					  depth,isl_space_dim(space,isl_dim_out)-depth);

      BSetPtr deltaBSetPtr=isl_basic_map_deltas(rawDepPtr);
      isl_space *deltaSpacePtr=isl_basic_set_get_space(deltaBSetPtr);

      // we only need to consider the intra-tile deps
      // so zero out the deltas for the inter-tile iterators
      for(int i=0; i<min(depth,tileStartDim); ++i){
	ISL::NameCoeffMap nameCoeffMap;
	nameCoeffMap[isl_space_get_dim_name(deltaSpacePtr,isl_dim_set,i)]=1;
	nameCoeffMap["1"]=0;
	deltaBSetPtr=isl_basic_set_add_constraint(deltaBSetPtr,ISL::EqFromNames(deltaSpacePtr,nameCoeffMap));
      }

      deltaBSetPtr=isl_basic_set_remove_divs(deltaBSetPtr);
      p=isl_printer_print_basic_set(p,deltaBSetPtr);
      cout << endl;

      cout << "The utility span is..." << endl;
      isl_set *setPtr=isl_basic_set_lexmax(deltaBSetPtr);
      p=isl_printer_print_set(p,setPtr);
      cout << endl << endl;

      // !!! can't we find lexmax over all these bsets (actually, they
      // are sets!) once?
      if(NULL==lexmax){
	lexmax=setPtr;
        inclusive=(depth<tileEndDim || src!=dest);
      }
      else if(isl_set_plain_cmp(lexmax,setPtr)<0){
	lexmax=setPtr;
        inclusive=(depth<tileEndDim || src!=dest);
      }
    }

    BSetPtrVec bsetPtrVec;
    isl_set_foreach_basic_set(lexmax,&ISL::AddAnyBasicSet,&bsetPtrVec);
    assert(!bsetPtrVec.empty());

    cout << "THE MAX UTILITY SPAN:" << endl;
    // !!! we should consider each basic set in the lexmax only one of
    // them will be for the average case, the rest of them handle
    // corner cases
    p=isl_printer_print_basic_set(p,bsetPtrVec.front());
    cout << endl << "Inclusive? " << inclusive << endl << endl;

    // !!! must pay attention to the other constraints in the basic
    // set during codegen
    ISL::GetILPSolution(bsetPtrVec.front(),isl_basic_set_get_space(bsetPtrVec.front()),
			tileStartDim,tileEndDim,coefficients);

    isl_printer_free(p);

    return inclusive; // inclusive or exclusive max utility span
  }

  // !!! check if isl_basic_map_deltas_map can be used for this
  void ConflictSpecBuilder::FindMaxUtilitySpanDeprecated(int stmtId,int tileStartDim,
					       int tileEndDim,vector<int> &coefficients){
    int tileDims=tileEndDim-tileStartDim;
    isl_set *lexmax=NULL;
    isl_printer *p=isl_printer_to_file(_ctx,stdout);

    cout << endl << "FINDING MAX UTILITY SPAN" << endl << endl;
    for(BMapPtrVecIter depIter=_rawdepPtrVec.begin();
	depIter!=_rawdepPtrVec.end(); ++depIter){
      BMapPtr rawdepPtr=*depIter;

      int depth=findInnerMostCommonTimeDim(rawdepPtr);

      isl_space *spacePtr=isl_basic_map_get_space(rawdepPtr);
      int ndims=isl_space_dim(spacePtr,isl_dim_set);

      StringVec uspanVec=Namer::GetNames(ndims,Namer::kMisc);
      StringVec paramVec,inVec,outVec,mirrorInVec,mirrorOutVec;

      // get the param names
      ISL::isl_space_get_dim_names(spacePtr,isl_dim_param,paramVec);

      for(int i=0; i<ndims; ++i){
	const char *namePtr=isl_space_get_dim_name(spacePtr,isl_dim_out,i);
	uspanVec.push_back(string(namePtr).append("o"));
	outVec.push_back(string(namePtr));
	mirrorOutVec.push_back(string(namePtr).append("o"));

	namePtr=isl_space_get_dim_name(spacePtr,isl_dim_in,i);
	uspanVec.push_back(string(namePtr).append("i"));
	inVec.push_back(string(namePtr));
	mirrorInVec.push_back(string(namePtr).append("i"));
      }

      isl_space *setSpacePtr=isl_space_set_alloc(_ctx,paramVec.size(),uspanVec.size());
      setSpacePtr=ISL::isl_space_set_dim_names(setSpacePtr,isl_dim_set,uspanVec);
      setSpacePtr=ISL::isl_space_set_dim_names(setSpacePtr,isl_dim_param,paramVec);

      BSetPtr bsetPtr=isl_basic_set_universe(setSpacePtr);

      // populate the set with the constraints in the RAW dependence map
      ConstraintPtrVec cstPtrVec;
      isl_basic_map_foreach_constraint(rawdepPtr,&ISL::AddAnyConstraint,&cstPtrVec);
      for(ConstraintPtrVecIter i=cstPtrVec.begin();
	  i!=cstPtrVec.end(); ++i){
	// !!! isn't there a simpler way to add the map's constraints to a set?
	ISL::NameCoeffMap nameCoeffMap;
	ISL::GetCoefficientVals(*i,spacePtr,inVec,isl_dim_in,mirrorInVec,nameCoeffMap);
	ISL::GetCoefficientVals(*i,spacePtr,inVec,isl_dim_out,mirrorOutVec,nameCoeffMap);
	ISL::GetCoefficientVals(*i,spacePtr,paramVec,isl_dim_param,paramVec,nameCoeffMap);
	if(isl_constraint_is_equality(*i))
	  bsetPtr=isl_basic_set_add_constraint(bsetPtr,ISL::EqFromNames(setSpacePtr,nameCoeffMap));
	else
	  bsetPtr=isl_basic_set_add_constraint(bsetPtr,ISL::IneqFromNames(setSpacePtr,nameCoeffMap));
      }

      for(int i=0; i<ndims; ++i){
	ISL::NameCoeffMap nameCoeffMap;
	nameCoeffMap[mirrorOutVec[i]]=1;
	nameCoeffMap[mirrorInVec[i]]=-1;
	nameCoeffMap[uspanVec[i]]=-1;
	bsetPtr=isl_basic_set_add_constraint(bsetPtr,ISL::EqFromNames(setSpacePtr,nameCoeffMap));
      }

      // we only need to consider the intra-tile dependences
      // !!! each raw dep could result in a different depth value!
      if(depth<tileEndDim){
      	tileEndDim=depth;
      	tileDims=tileEndDim-tileStartDim;
      }
      // if(tileStartDim>0){
      // 	bsetPtr=isl_basic_set_project_out(bsetPtr,isl_dim_set,ndims,2*(tileStartDim));
      // 	bsetPtr=isl_basic_set_project_out(bsetPtr,isl_dim_set,0,tileStartDim);
      // }
      // if(ndims>tileEndDim){
      //   bsetPtr=isl_basic_set_project_out(bsetPtr,isl_dim_set,tileDims,ndims-tileEndDim);
      //   bsetPtr=isl_basic_set_project_out(bsetPtr,isl_dim_set,3*tileDims,2*(ndims-tileEndDim));
      // }

      // add equality constraint for the inter-tile iterators
      for(int i=0; i<tileStartDim; ++i){
	ISL::NameCoeffMap nameCoeffMap;
	nameCoeffMap[mirrorOutVec[i]]=1;
	nameCoeffMap[mirrorInVec[i]]=-1;
	nameCoeffMap["1"]=0;
	bsetPtr=isl_basic_set_add_constraint(bsetPtr,ISL::EqFromNames(setSpacePtr,nameCoeffMap));
      }

      cout << "Finding max utility span..." << endl;
      cout << "The ILP is as follows..." << endl;
      p=isl_printer_print_basic_set(p,bsetPtr);
      cout << endl << "The solution is..." << endl;

      isl_set *setPtr=isl_basic_set_lexmax(bsetPtr);

      // !!! can't we find lexmax over all these bsets once?
      if(NULL==lexmax)
	lexmax=setPtr;
      else if(isl_set_plain_cmp(lexmax,setPtr)<0)
	lexmax=setPtr;

      p=isl_printer_print_set(p,setPtr);
      cout << endl << endl;
    }

    BSetPtrVec bsetPtrVec;
    isl_set_foreach_basic_set(lexmax,&ISL::AddAnyBasicSet,&bsetPtrVec);
    assert(!bsetPtrVec.empty());

    p=isl_printer_print_basic_set(p,bsetPtrVec.front());
    cout << endl;

    ISL::GetILPSolution(bsetPtrVec.front(),isl_basic_set_get_space(bsetPtrVec.front()),
			0,tileDims,coefficients);

    isl_printer_free(p);
  }

  // inbsetPtr constraints will be transferred to another basic set
  // constructed with the given new space which differs the old space
  // of the constraints only in the sense that it has more parameters
  BSetPtr ConflictSpecBuilder::mapBSetToParameterizedSpace(BSetPtr inbsetPtr,
							     isl_space *newspace){
    isl_printer *p=isl_printer_to_file(_ctx,stdout);

    ConstraintPtrVec cstPtrVec;
    isl_basic_set_foreach_constraint(inbsetPtr,&ISL::AddAnyConstraint,&cstPtrVec);

    StringVec paramVec,varVec,originalParamVec;
    ISL::isl_space_get_dim_names(isl_basic_set_get_space(inbsetPtr),
				        isl_dim_param,originalParamVec);
    ISL::isl_space_get_dim_names(newspace,isl_dim_param,paramVec);
    ISL::isl_space_get_dim_names(newspace,isl_dim_set,varVec);

    BSetPtr outbsetPtr=isl_basic_set_universe(newspace);
    outbsetPtr=isl_basic_set_set_tuple_name(outbsetPtr,
					    isl_basic_set_get_tuple_name(inbsetPtr));

    for(ConstraintPtrVecIter i=cstPtrVec.begin(); i!=cstPtrVec.end(); ++i){
      isl_val *v=isl_constraint_get_constant_val(*i);
      long num=isl_val_get_num_si(v);
      long absnum=abs(num);
      ISL::NameCoeffMap dict;

      ISL::GetCoefficientVals(*i,newspace,originalParamVec,
			      isl_dim_param,paramVec,dict);
      ISL::GetCoefficientVals(*i,newspace,varVec,
			      isl_dim_set,varVec,dict);
      if(absnum>BIGC){
	if(gParams.find(absnum)!=gParams.end()){
	  dict[gParams[absnum]]=((absnum!=num) ? -1 : 1);
	  dict["1"]=0;
	}
	else{
	  bool found=false;
	  for(map<long,string>::iterator it=gParams.begin();
	      it!=gParams.end() && !found; ++it){
	    long key=it->first;
	    long absdiff=abs(absnum-key);
	    if(absdiff<BIGC){
	      int sign=((absnum!=num) ? -1 : 1);
	      dict[gParams[key]]=sign;
	      dict["1"]= sign*((key>absnum) ? (absnum-key) : (key-absnum));
	      found=true;
	    }
	  }
	  // cout << absnum << endl;
	  //assert(found);
	}
      }

      if(isl_constraint_is_equality(*i))
	outbsetPtr=isl_basic_set_add_constraint(outbsetPtr,ISL::EqFromNames(newspace,dict));
      else
	outbsetPtr=isl_basic_set_add_constraint(outbsetPtr,ISL::IneqFromNames(newspace,dict));
    }

    p=isl_printer_print_basic_set(p,outbsetPtr);
    cout << endl;
    isl_printer_free(p);
    return outbsetPtr;
  }

  void ConflictSpecBuilder::parameterizeStmtDomains(){
    isl_printer *p=isl_printer_to_file(_ctx,stdout);

    cout << endl << "PARAMETERIZATION" << endl << endl;
    for(BSetPtrVecIter it=_domPtrVec.begin();
	it!=_domPtrVec.end(); ++it){
      BSetPtr domain=isl_basic_set_remove_redundancies(*it);
      isl_space *spacePtr=isl_basic_set_get_space(domain);

      StringVec paramVec,varVec;
      ISL::isl_space_get_dim_names(spacePtr,isl_dim_param,paramVec);
      ISL::isl_space_get_dim_names(spacePtr,isl_dim_set,varVec);

      cout << "The given iteration domain is as follows" << endl;
      p=isl_printer_print_basic_set(p,domain);
      cout << endl;
      *it=parameterize(domain,paramVec,varVec);
    }
    isl_printer_free(p);
  }

  BSetPtr ConflictSpecBuilder::parameterize(BSetPtr dom, StringVec &paramVec,
					    StringVec &varVec){
    StringVec originalParamVec=paramVec;
    ConstraintPtrVec cstPtrVec;
    isl_printer *p=isl_printer_to_file(_ctx,stdout);

    isl_basic_set_foreach_constraint(dom,&ISL::AddAnyConstraint,&cstPtrVec);

    // first introduce as many new parameters as required
    bool needsParameterization=false;
    for(ConstraintPtrVecIter i=cstPtrVec.begin(); i!=cstPtrVec.end(); ++i){
      isl_val *v=isl_constraint_get_constant_val(*i);
      long num=abs(isl_val_get_num_si(v));
      if(num>BIGC){
	needsParameterization=true;
	if(gParams.find(num)==gParams.end()){
	  gParams[num]=Namer::GetNames(1,Namer::kParam).front();
	  paramVec.push_back(gParams[num]);
	  cout << "Introducing parameter " << gParams[num]
	       << " for " << num << " due to the constraint" << endl;
	  p=isl_printer_print_constraint(p,*i);
	  cout << endl;
	}
	else if(find(paramVec.begin(),paramVec.end(),gParams[num])==paramVec.end()){
	  paramVec.push_back(gParams[num]);
	}
      }
    }
    cout << endl;

    if(!needsParameterization)
      return dom; 

    // create another space with to accomodate any new names
    isl_space *space=isl_space_set_alloc(_ctx,paramVec.size(),varVec.size());
    space=ISL::isl_space_set_dim_names(space,isl_dim_param,paramVec);
    space=ISL::isl_space_set_dim_names(space,isl_dim_set,varVec);

    // replace big constants with the corresponding new parameters
    // BSetPtr bsetPtr=isl_basic_set_universe(space);
    // for(ConstraintPtrVecIter i=cstPtrVec.begin(); i!=cstPtrVec.end(); ++i){
    //   isl_val *v=isl_constraint_get_constant_val(*i);
    //   long num=isl_val_get_num_si(v);
    //   long absnum=abs(num);
    //   ISL::NameCoeffMap map1;

    //   ISL::GetCoefficientVals(*i,space,originalParamVec,
    // 			      isl_dim_param,paramVec,map1);
    //   ISL::GetCoefficientVals(*i,space,varVec,
    // 			      isl_dim_set,varVec,map1);
    //   if(num>BIGC){
    // 	map1[gParams[num]]=((absnum!=num) ? -1 : 1);
    // 	map1["1"]=0;
    //   }

    //   if(isl_constraint_is_equality(*i))
    // 	bsetPtr=isl_basic_set_add_constraint(bsetPtr,ISL::EqFromNames(space,map1));
    //   else
    // 	bsetPtr=isl_basic_set_add_constraint(bsetPtr,ISL::IneqFromNames(space,map1));
    // }

    BSetPtr bsetPtr=mapBSetToParameterizedSpace(dom,space);

    p=isl_printer_print_basic_set(p,bsetPtr);
    cout << endl;
    isl_printer_free(p);
    return bsetPtr;
  }

  ConflictSpec &ConflictSpecBuilder::Build(int stmtId,int tileStartDim,
			     int tileEndDim,vector<int> &coefficients,
					   bool inclusive,isl_union_set *liveOut){
    isl_printer *p=isl_printer_to_file(_ctx,stdout);

    BSetPtr domain=isl_basic_set_copy(getDomain(stmtId));
    isl_space *spacePtr=isl_basic_set_get_space(domain);

    StringVec paramVec,varVec;
    ISL::isl_space_get_dim_names(spacePtr,isl_dim_param,paramVec);
    ISL::isl_space_get_dim_names(spacePtr,isl_dim_set,varVec);

    // cout << endl << "PARAMETERIZATION" << endl << endl;

    // cout << "The given iteration domain is as follows" << endl;
    // p=isl_printer_print_basic_set(p,domain);
    // cout << endl;
    // domain=parameterize(domain,paramVec,varVec);
    // spacePtr=isl_basic_set_get_space(domain);

    // create a copy of the domain with just the names (no statement name)
    ConstraintPtrVec cstPtrVec;
    isl_basic_set_foreach_constraint(domain,&ISL::AddAnyConstraint,&cstPtrVec);

    StringVec mirParamVec,mirVarVec;
    for(StringVecIter i=varVec.begin(); i!=varVec.end(); ++i){
      string s=*i;
      mirVarVec.push_back(s.append("'"));
    }

    BSetPtr range=ISL::CreateSetFromConstraints(_ctx,spacePtr,
						paramVec,paramVec,
						varVec,mirVarVec,cstPtrVec);

    // project out the dims we are not interested in contracting
    cout << endl << "BUILDING THE CONFLICT SPEC " << endl << endl;
    cout << "The conflict relation is between tthe following spaces" << endl;
    int ndims=isl_space_dim(spacePtr,isl_dim_set);
    int tileDims=tileEndDim-tileStartDim;

    //domain=isl_basic_set_project_out(domain,isl_dim_set,0,ndims-tileDims);
    domain=isl_basic_set_remove_divs(domain);
    //domain=isl_basic_set_remove_redundancies(domain);
    //range=isl_basic_set_project_out(range,isl_dim_set,0,ndims-tileDims);
    range=isl_basic_set_remove_divs(range);
    //range=isl_basic_set_remove_redundancies(range);

    p=isl_printer_print_basic_set(p,domain);
    cout << endl;
    p=isl_printer_print_basic_set(p,range);
    cout << endl << endl;

    isl_space *domainSpacePtr=isl_basic_set_get_space(domain);
    BMapPtr orderedPairPtr=isl_basic_map_from_domain_and_range(domain,range);

    for(int i=0; i<tileStartDim; ++i){
      ISL::NameCoeffMap dict;
      isl_space *orderedPairSpace=isl_basic_map_get_space(orderedPairPtr);
      dict[string(isl_space_get_dim_name(orderedPairSpace,isl_dim_out,i))]=1;
      dict[string(isl_space_get_dim_name(orderedPairSpace,isl_dim_in,i))]=-1;
      dict["1"]=0;
      orderedPairPtr=isl_basic_map_add_constraint(orderedPairPtr,ISL::EqFromNames(orderedPairSpace,dict));
    }

    cout << "Parameterizing the live out" << endl;
    BSetPtrVec liveOutBSetPtrVec;
    ISL::ExtractAllBasicSet(liveOut,&liveOutBSetPtrVec);
    for(BSetPtrVecIter it=liveOutBSetPtrVec.begin();
	it!=liveOutBSetPtrVec.end(); ++it){
      // !!! what if more parameters need to be added?
      *it=mapBSetToParameterizedSpace(*it,domainSpacePtr);
    }

    ConflictSpec &cSpec = ConflictSpec::Create(_ctx,orderedPairPtr,coefficients,
					       inclusive,liveOutBSetPtrVec,
					       tileStartDim,tileEndDim);
    isl_printer_free(p);
    return cSpec;
  }

  isl_union_set *ConflictSpecBuilder::InferLiveOut(int stmtId,int tileStart,int tileEnd){
    assert(tileStart>0);
    isl_printer *p=isl_printer_to_file(_ctx,stdout);

    BMapPtrVec tileDepPtrVec;
    cout << "INFERRING LIVE OUT" << endl ;
    for(BMapPtrVecIter depIter=_rawdepPtrVec.begin();
	             depIter!=_rawdepPtrVec.end(); ++depIter){
      // work on a copy as we modify the out dim names
      BMapPtr rawDepPtr=isl_basic_map_copy(*depIter);

      int src=atoi(isl_basic_map_get_tuple_name(rawDepPtr,isl_dim_in)+2);
      if(src!=stmtId) continue;

      cout << endl << "Examining the following raw dependence" << endl;
      isl_printer_print_basic_map(p,rawDepPtr);
      cout << endl;

      int depth=findInnerMostCommonTimeDim(rawDepPtr,tileStart);

      isl_space *spacePtr=isl_basic_map_get_space(rawDepPtr);

      // !!! sort this out in a better way
      // first set the dim names of the out variables to make sure
      // they are different from those of the in variables
      for(int i=0; i<isl_space_dim(spacePtr,isl_dim_out); ++i){
	const char *oldname=isl_space_get_dim_name(spacePtr,isl_dim_out,i);
	char newname[4];
	sprintf(newname,"%s'",oldname);
	spacePtr=isl_space_set_dim_name(spacePtr,isl_dim_out,i,newname);
      }

      BMapPtr rootBMapPtr=isl_basic_map_universe(spacePtr);
      rootBMapPtr=isl_basic_map_intersect(rootBMapPtr,rawDepPtr);

      for(int i=0; i<depth; ++i){ // segregate all the tile dependences
        Smo::ISL::NameCoeffMap dict;

	BMapPtr tileDepPtr=isl_basic_map_copy(rootBMapPtr);
	dict[isl_space_get_dim_name(spacePtr,isl_dim_out,i)]=1;
	dict[isl_space_get_dim_name(spacePtr,isl_dim_in,i)]=-1;
	dict["1"]=-1;

	tileDepPtr=isl_basic_map_add_constraint(tileDepPtr,
		     Smo::ISL::IneqFromNames(isl_basic_map_get_space(tileDepPtr),dict));
	if(isl_basic_map_is_empty(tileDepPtr)==0){
	  tileDepPtrVec.push_back(tileDepPtr);
	  cout << endl << "Found an inter-tile dependence " << endl;
	  p=isl_printer_print_basic_map(p,tileDepPtr);
	  cout << endl;
	}

	dict.clear();
	dict[isl_space_get_dim_name(spacePtr,isl_dim_out,i)]=1;
	dict[isl_space_get_dim_name(spacePtr,isl_dim_in,i)]=-1;
	dict["1"]=0;
	rootBMapPtr=isl_basic_map_add_constraint(rootBMapPtr,
					   Smo::ISL::EqFromNames(isl_basic_map_get_space(rootBMapPtr),dict));
      }
      isl_basic_map_free(rootBMapPtr);
    }

    isl_union_set *liveOut=NULL;
    cout << endl << "Live out sets inferred are as follows" << endl;
    for(BMapPtrVecIter iter=tileDepPtrVec.begin();
	iter!=tileDepPtrVec.end(); ++iter){
      // project out the tile iterators
      BMapPtr bmapPtr=*iter;
      // !!! this creates problem for the kind of tiled code we work with
      // bmapPtr=isl_basic_map_project_out(bmapPtr,isl_dim_in,0,tileStart);
      // bmapPtr=isl_basic_map_project_out(bmapPtr,isl_dim_out,0,tileStart);

      isl_basic_set *liveOutBSetPtr=isl_basic_set_remove_redundancies(isl_basic_map_domain(bmapPtr));
      isl_set *liveOutSetPtr=isl_set_from_basic_set(liveOutBSetPtr);

      // !!! is this too early to remove divs?
      liveOutSetPtr=isl_set_remove_divs(liveOutSetPtr);
      // p=isl_printer_print_set(p,liveOutSetPtr);
      // cout << endl;

      if(NULL==liveOut){
	liveOut=isl_union_set_from_set(liveOutSetPtr);
      }
      else{
	liveOut=isl_union_set_union(isl_union_set_from_set(liveOutSetPtr),liveOut);
      }
    }

    liveOut=isl_union_set_coalesce(liveOut);

    p=isl_printer_print_union_set(p,liveOut);
    cout << endl ;
    isl_printer_free(p);

    return liveOut;
  }

  BSetPtr ConflictSpecBuilder::getDomain(int stmtId){
    isl_basic_set *domain=NULL;
    for(BSetPtrVecIter it=_domPtrVec.begin();
	it!=_domPtrVec.end(); ++it){
      domain=*it;
      int src=atoi(isl_basic_set_get_tuple_name(domain)+2);
      if(src==stmtId) break;
    }
    return domain;
  }
} // namespece Smo

