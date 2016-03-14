#include "ISLUtils.h"

namespace Smo{
  namespace ISL{
    void Print(isl_ctx *ctx, Smo::BMapPtrVec &bmapPtrVec){
      isl_printer *p=isl_printer_to_file(ctx,stdout);
      for(Smo::BMapPtrVecIter i=bmapPtrVec.begin();
	  i!=bmapPtrVec.end(); ++i){
	p=isl_printer_print_basic_map(p,*i);
	cout << endl;
      }
      cout << endl;
      isl_printer_free(p);
    }

    void Print(isl_ctx *ctx, Smo::BSetPtrVec &bsetPtrVec){
      isl_printer *p=isl_printer_to_file(ctx,stdout);
      for(Smo::BSetPtrVecIter i=bsetPtrVec.begin();
	  i!=bsetPtrVec.end(); ++i){
	p=isl_printer_print_basic_set(p,*i);
	cout << endl;
      }
      cout << endl;
      isl_printer_free(p);
    }

    __isl_give isl_space * isl_space_set_dim_names(__isl_take isl_space *space,
		enum isl_dim_type type,StringVec nameVec){
      for(unsigned i=0; i<isl_space_dim(space,type) ; ++i){
	space=isl_space_set_dim_name(space,type,i,nameVec[i].c_str());
      }
      return space;
    }

    void isl_space_get_dim_names(__isl_keep isl_space *space,
	        enum isl_dim_type type,StringVec &nameVec){
      for(unsigned i=0; i<isl_space_dim(space,type) ; ++i){
	const char *name=isl_space_get_dim_name(space,type,i);
	if(name==NULL){
	  char str[4];
	  sprintf(str,"o%d",i);
	  nameVec.push_back(string(str));
	}
	else{
	  nameVec.push_back(string(name));
	}
      }
    }

    ConstraintPtr SetCoefficientVals(__isl_keep isl_space *spacePtr,
                 __isl_take ConstraintPtr cPtr,NameCoeffMap &nameCoeffMapRef){
      for(NameCoeffMap::iterator mapIter=nameCoeffMapRef.begin(); 
                                   mapIter != nameCoeffMapRef.end(); ++mapIter){
	string name=mapIter->first;

	if(strcmp("1",name.c_str()) == 0){
	  int val=mapIter->second;
	  cPtr=isl_constraint_set_constant_si(cPtr,val);
	}
	else{
	  int pos=isl_space_find_dim_by_name(spacePtr,isl_dim_set,name.c_str());
	  isl_dim_type type=isl_dim_set;
	  if(pos == -1){
	    pos=isl_space_find_dim_by_name(spacePtr,isl_dim_in,name.c_str());
	    type=isl_dim_in;
	    if(pos == -1){
	      pos=isl_space_find_dim_by_name(spacePtr,isl_dim_param,name.c_str());
	      type=isl_dim_param;
	    }
	  }
	  // assert(pos>=0);

	  int val=mapIter->second;
	  cPtr=isl_constraint_set_coefficient_si(cPtr,type,pos,val);
	}
      }
      return cPtr;
    }

    // gets the coefficient values for the names in <names>
    // maps the correspond names in <mirrorNames> to the value obtained
    void GetCoefficientVals(__isl_keep ConstraintPtr c, __isl_keep isl_space *spacePtr,
			    StringVec names, isl_dim_type type,
			    StringVec mirrorNames, NameCoeffMap &nameCoeffMapRef){
      int pos=0;
      for(StringVecIter s=names.begin(), ms=mirrorNames.begin();
	  s!=names.end(); ++s,++ms,++pos){
	//int pos=isl_space_find_dim_by_name(spacePtr,type,s->c_str()); // problem when dims are not named!
	isl_val *v=isl_constraint_get_coefficient_val(c,type,pos);
	nameCoeffMapRef[*ms]=isl_val_get_num_si(v);
      }
      isl_val *v=isl_constraint_get_constant_val(c);
      nameCoeffMapRef["1"]=isl_val_get_num_si(v);
    }

    // this method can be used to obtain the solution to an ilp
    // i.e., it populates the vector <coefficients>
    void GetILPSolution(__isl_keep BSetPtr bsetPtr, __isl_keep isl_space *spacePtr,
			int startIndex, int endIndex, vector<int> &coefficients){
      assert(bsetPtr!=NULL && spacePtr!=NULL);
      ConstraintPtrVec cstPtrVec;
      isl_basic_set_foreach_constraint(bsetPtr,&ISL::AddAnyConstraint,&cstPtrVec);
      int ndims=isl_space_dim(spacePtr,isl_dim_set);
    
      // determine the coefficients
      for(unsigned i=startIndex; i<endIndex; ++i){
        for(ConstraintPtrVecIter iter=cstPtrVec.begin();
	    iter!=cstPtrVec.end(); ++iter){
	  isl_val *v=isl_constraint_get_coefficient_val(*iter,isl_dim_set,i);
	  long num=isl_val_get_num_si(v);
	  if(num==1){
	    coefficients.push_back(-isl_val_get_num_si(isl_constraint_get_constant_val(*iter)));
	    break;
	  }
        }
      }
    }

    BSetPtr CreateSetFromConstraints(isl_ctx *ctx, isl_space *spacePtr,
				     StringVec &iParamVec, StringVec &oParamVec,
				     StringVec &iNameVec, StringVec &oNameVec,
				     ConstraintPtrVec &cstPtrVec){
      isl_space *setSpacePtr=isl_space_set_alloc(ctx,oParamVec.size(),oNameVec.size());
      setSpacePtr=Smo::ISL::isl_space_set_dim_names(setSpacePtr,isl_dim_param,oParamVec);
      setSpacePtr=Smo::ISL::isl_space_set_dim_names(setSpacePtr,isl_dim_set,oNameVec);
      Smo::BSetPtr bsetPtr=isl_basic_set_universe(setSpacePtr);

      // populate the set with the constraints in the RAW dependence map
      for(Smo::ConstraintPtrVecIter i=cstPtrVec.begin();
	    i!=cstPtrVec.end(); ++i){
	// !!! isn't there a simpler way to add the map's constraints to a set?
	Smo::ISL::NameCoeffMap nameCoeffMap;
	Smo::ISL::GetCoefficientVals(*i,spacePtr,iNameVec,isl_dim_set,oNameVec,nameCoeffMap);
	Smo::ISL::GetCoefficientVals(*i,spacePtr,iParamVec,isl_dim_param,oParamVec,nameCoeffMap);
	if(isl_constraint_is_equality(*i))
	  bsetPtr=isl_basic_set_add_constraint(bsetPtr,Smo::ISL::EqFromNames(setSpacePtr,nameCoeffMap));
	else
	  bsetPtr=isl_basic_set_add_constraint(bsetPtr,Smo::ISL::IneqFromNames(setSpacePtr,nameCoeffMap));

      // isl_printer *p=isl_printer_to_file(ctx,stdout);
      // p=isl_printer_print_constraint(p,*i); cout << "   ";
      // p=isl_printer_print_basic_set(p,bsetPtr);
      // cout << endl;
      // isl_printer_free(p);
      }
      return bsetPtr;
    }

    ConstraintPtr IneqFromNames(__isl_keep isl_space *spacePtr,
                                                 NameCoeffMap &nameCoeffMapRef){
      isl_space *spaceCopyPtr=isl_space_copy(spacePtr);
      isl_local_space *lspacePtr=isl_local_space_from_space(spaceCopyPtr);
      ConstraintPtr cPtr=isl_inequality_alloc(lspacePtr);
      return SetCoefficientVals(spacePtr,cPtr,nameCoeffMapRef);
    }

    ConstraintPtr EqFromNames(__isl_keep isl_space *spacePtr,
                                                 NameCoeffMap &nameCoeffMapRef){
      isl_space *spaceCopyPtr=isl_space_copy(spacePtr);
      isl_local_space *lspacePtr=isl_local_space_from_space(spaceCopyPtr);
      ConstraintPtr cPtr=isl_equality_alloc(lspacePtr);
      return SetCoefficientVals(spacePtr,cPtr,nameCoeffMapRef);
    }

    int AddIfEqConstraint(ConstraintPtr c,void *cstVecPtr){
      if(isl_constraint_is_equality(c)){
	static_cast<ConstraintPtrVec *>(cstVecPtr)->push_back(c);
      }
      return 0;
    }

    int AddAnyConstraint(ConstraintPtr c,void *cstVecPtr){
      static_cast<ConstraintPtrVec *>(cstVecPtr)->push_back(c);
      return 0;
    }

    int AddAnyBasicSet(BSetPtr b,void *bsetVecPtr){
      static_cast<BSetPtrVec *>(bsetVecPtr)->push_back(b);
      return 0;
    }

    int AddAnyBasicMap(BMapPtr b,void *bmapVecPtr){
      static_cast<BMapPtrVec *>(bmapVecPtr)->push_back(b);
      return 0;
    }

    int AddAnySet(SetPtr b,void *setVecPtr){
      static_cast<SetPtrVec *>(setVecPtr)->push_back(b);
      return 0;
    }

    int AddAnyMap(MapPtr b,void *mapVecPtr){
      static_cast<MapPtrVec *>(mapVecPtr)->push_back(b);
      return 0;
    }

    void ExtractAllBasicSet(isl_union_set *usetPtr,void *bsetVecPtr){
      SetPtrVec setPtrVec;
      isl_union_set_foreach_set(usetPtr,AddAnySet,&setPtrVec);
      for(SetPtrVecIter i=setPtrVec.begin(); i!=setPtrVec.end(); ++i){
	isl_set_foreach_basic_set(*i,AddAnyBasicSet,bsetVecPtr);
      }
    }

    void ExtractAllBasicMap(isl_union_map *umapPtr,void *bmapVecPtr){
      MapPtrVec mapPtrVec;
      isl_union_map_foreach_map(umapPtr,AddAnyMap,&mapPtrVec);
      for(MapPtrVecIter i=mapPtrVec.begin(); i!=mapPtrVec.end(); ++i){
	isl_map_foreach_basic_map(*i,AddAnyBasicMap,bmapVecPtr);
      }
    }

    __isl_give isl_mat *AddRow(isl_ctx *ctx, __isl_take isl_mat *in, isl_mat *row){
      int nr=isl_mat_rows(row);
      int nc=isl_mat_cols(row);
      assert(nr==1);

      vector<int> rowVec;
      for(int j=0; j<nc; ++j){
        isl_val *val=isl_mat_get_element_val(row,0,j);
	rowVec.push_back(isl_val_get_num_si(val));
      }
      return AddRow(ctx,in,rowVec);
    }

    __isl_give isl_mat *AddRow(isl_ctx *ctx, __isl_take isl_mat *in, vector<int> row){
      int nr=1;
      int nc=row.size();
      isl_mat *out=NULL;

      if(NULL!=in){
	nr=isl_mat_rows(in)+1;
	nc=isl_mat_cols(in);

	assert(nc==row.size());

	// copy the input matrix rows to the new matrix
	out=isl_mat_alloc(ctx,nr,nc);
	for(int i=0; i<nr-1; ++i)
	  for(int j=0; j<nc; ++j){
	    isl_val *val=isl_mat_get_element_val(in,i,j);
	    out=isl_mat_set_element_si(out,i,j,isl_val_get_num_si(val));
	  }
      }
      else{
	out=isl_mat_alloc(ctx,1,nc);
      }

      // make the given row the last row of the matrix
      for(int i=0; i<nc; ++i)
	out=isl_mat_set_element_si(out,nr-1,i,row[i]);

      isl_mat_free(in);
      return out;
    }

    __isl_give isl_mat *GetRowMatrix(isl_ctx *ctx, __isl_keep isl_mat *in, int row){
      assert(row<isl_mat_rows(in));

      int nc=isl_mat_cols(in);
      isl_mat *out=isl_mat_alloc(ctx,1,nc);

      // copy the specified row to the row matrix
      for(int i=0; i<nc; ++i){
	isl_val *val=isl_mat_get_element_val(in,row,i);
	out=isl_mat_set_element_si(out,0,i,isl_val_get_num_si(val));
      }
      return out;
    }

    __isl_give isl_mat *RemoveRow(isl_ctx *ctx, __isl_take isl_mat *in, int row){
      if(NULL==in || 1==isl_mat_rows(in))
	return NULL;
      assert(row<isl_mat_rows(in));

      int nr=isl_mat_rows(in);
      int nc=isl_mat_cols(in);
      isl_mat *out=isl_mat_alloc(ctx,nr-1,nc);

      int r=0,i=0;
      for(; r<nr;){
        // do not copy the specified row to the row matrix
	if(r==row){
	  ++r;
          continue;
	}
	// copy the others
        for(int c=0; c<nc; ++c){
	  isl_val *val=isl_mat_get_element_val(in,r,c);
	  out=isl_mat_set_element_si(out,i,c,isl_val_get_num_si(val));
        }
	++i;
      }
      return out;
    }

    void PrintMatrix(__isl_give isl_mat *m){
      int nr=isl_mat_rows(m);
      int nc=isl_mat_cols(m);

      for(int i=0; i<nr; ++i){
	for(int j=0; j<nc; ++j){
	  isl_val *val=isl_mat_get_element_val(m,i,j);
	  cout << isl_val_get_num_si(val) << " ";
	}
	cout << endl;
      }
    }
  } // namespace ISL
}

