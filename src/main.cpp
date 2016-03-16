/*
 *
 * SMO: A memory/storage optimizer for affine loop nests
 *
 * The MIT License (MIT)

 * Copyright (c) 2015-2016 Somashekaracharya Bhaskaracharya

 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */
#include <stdio.h>
#include <getopt.h>
#include <fstream>
#include <sstream>

#include "IterativeStoragePartition.h"
#include "ConflictSpecBuilder.h"

#include "Namer.h"
#include "ISLUtils.h"

struct SMOOptions{
  int enableDebug;
  int glpkSolve;
  int enumerate;
  int extractScop;
} SMOOpts;

void printHelpMessage()
{
    fprintf(stdout, "Usage: ./smo <xyz.smo.input>|[-i scopfile] [options]\n");
    fprintf(stdout, "\nOptions:\n");
    fprintf(stdout, "       --enable-debug|-d        Verbose debug prints\n");
    fprintf(stdout, "       --extract-scop|-s        Extract scop from given C code\n");
    // fprintf(stdout, "       --glpk-solve  |-g          Use GLPK as the ilp solver\n");
    fprintf(stdout, "       --enumerate|-e          Enumerate various feasible storage mappings\n");
    fprintf(stdout, "       --help        | -h         Print this help menu\n");
}

void storageOptimizeScopInput(string inputFileName){
  Smo::ConflictSpecBuilder cSpecBuilder(inputFileName);

  // do array optimization for each stmt in the tilespaces file
  ifstream tilespacesFile((inputFileName+string(".tilespaces")).c_str());
  string line;
  while(!tilespacesFile.fail() && getline(tilespacesFile,line)){
    // every line should be of the form <stmtId,tileStartIdx,tileEndIdx>
    int stmtId,tileStartIdx,tileEndIdx;
    istringstream iss(line);
    if(!(iss >> stmtId >> tileStartIdx >> tileEndIdx))
      assert("Specification in wrong format");

    // !!! need to handle the scenario when there are many RAW deps
    // !!! also need to handle this statement wise perhaps?
    vector<int> coefficients;
    bool inclusive=cSpecBuilder.FindMaxUtilitySpan(stmtId,tileStartIdx,
				     tileEndIdx,coefficients);
    isl_union_set *liveOut=NULL;
    if(tileStartIdx>0){
      liveOut=cSpecBuilder.InferLiveOut(stmtId,tileStartIdx,tileEndIdx);
    }

    Smo::ConflictSpec &cSpec=cSpecBuilder.Build(stmtId,tileStartIdx,
				     tileEndIdx,coefficients,inclusive,liveOut);

    Smo::IterativeStoragePartition iterativePartition(cSpec);
    iterativePartition.FindStorageHyperplanes(SMOOpts.enumerate);
  }
}

void storageOptimizeSmoInput(string filename){

  struct isl_ctx *ctx=isl_ctx_alloc();
  assert(ctx!=NULL);

  // read the domain of the conflict set relations
  ifstream ifs;
  ifs.open(filename.c_str(),ios::in);
  assert(ifs.is_open() && (ifs.rdstate() & ifstream::failbit) == 0);

  //  the domain
  Smo::BSetPtr domain=isl_basic_set_read_from_str(ctx,Smo::NextLine(ifs).c_str());
  
  // reset the dim names
  int nParams=isl_space_dim(isl_basic_set_get_space(domain),isl_dim_param);
  int nVars=isl_space_dim(isl_basic_set_get_space(domain),isl_dim_set);
  Smo::StringVec paramNameVec=Smo::Namer::GetNames(nParams,Smo::Namer::kParam);
  Smo::StringVec varNameVec=Smo::Namer::GetNames(nVars,Smo::Namer::kInput);
  isl_space *space=isl_space_set_alloc(ctx,nParams,nVars);
  space=Smo::ISL::isl_space_set_dim_names(space,isl_dim_param,paramNameVec);
  space=Smo::ISL::isl_space_set_dim_names(space,isl_dim_set,varNameVec);

  // reset the space of the domain
  isl_space *copySpacePtr=isl_space_copy(space);
  domain=(Smo::BSetPtr)isl_set_reset_space((isl_set *)domain,copySpacePtr);

  // read the image of the conflict set relations
  Smo::BSetPtr image=isl_basic_set_read_from_str(ctx,Smo::NextLine(ifs).c_str());
  image=(Smo::BSetPtr )isl_set_reset_space((isl_set *)image,space);

  // create an ordered pair from the given domain and range
  Smo::BMapPtr orderPairPtr=isl_basic_map_from_domain_and_range(domain,image);

  // !!! ideally,the order pair should be constructed using a range
  // with these output names
  string outputNames=Smo::NextLine(ifs);

  // read the number of statements
  int nStmnts=atoi(Smo::NextLine(ifs).c_str());

  int nConflictPoly=atoi(Smo::NextLine(ifs).c_str());

  // read the conflict polytopes
  // Smo::ConflictSpec *intraCSpecPtr=NULL,*interCSpecPtr=NULL;
  // if(nIntraArrPoly>0)
  // intraCSpecPtr=&Smo::ConflictSpec::reate(ifs,nIntraArrPoly,ctx,orderPairPtr);
  // if(nInterArrPoly>0)
  // interCSpecPtr=&Smo::ConflictSpec::Create(ifs,nInterArrPoly,ctx,orderPairPtr);

  Smo::ConflictSpec &cSpecRef=Smo::ConflictSpec::Create(ifs,nStmnts,nConflictPoly,ctx,orderPairPtr);

  // clean-up
  isl_basic_map_free(orderPairPtr);
  ifs.close();

  Smo::IterativeStoragePartition iterativePartition(cSpecRef);
  iterativePartition.FindStorageHyperplanes(SMOOpts.enumerate);
}

int main(int argc, char **argv){
  if(argc<=1){
    printHelpMessage();
    return 1;
  }

  int optIndex=0;
  const struct option longOptions[] =
    {
      {"enable-debug", no_argument, &SMOOpts.enableDebug, 1},
      {"glpk-solve", no_argument, &SMOOpts.glpkSolve, 1},
      {"extract-scop", no_argument, &SMOOpts.extractScop, 1},
      {"enumerate", no_argument, &SMOOpts.enumerate, 1},
      {"help", no_argument, 0, 'h'},
      {0, 0, 0, 0}
    };

  while(true){
    int opt=getopt_long(argc, argv, "dgehs", longOptions,
                &optIndex);
    if(opt==-1)
      break;

    switch(opt){
      case 'd': SMOOpts.enableDebug=true;
                break;
      // case 'g': SMOOpts.glpkSolve=true;
      //           break;
      case 'e': SMOOpts.enumerate=true;
                break;
      case 's': SMOOpts.extractScop=true;
	        break;
      case 'h': 
      case '?': 
      default : printHelpMessage();
	        return 1;
    }
  }

  if(optind<=argc-1){
    string inputFileName=string(argv[optind]);

    if(SMOOpts.extractScop){ // extract the scop
      storageOptimizeScopInput(inputFileName);
    }
    else{ // read the given conflict specification and solve
      storageOptimizeSmoInput(inputFileName);
    }
  }
  else{ // no input file was specified
    printHelpMessage();
    return 1;
  }

  return 0;
}
