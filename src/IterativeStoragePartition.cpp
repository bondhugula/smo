#include "IterativeStoragePartition.h"
#include "ISLUtils.h"
#include "Namer.h"
#include "StoragePartition.h"

#include <sys/time.h>

namespace Smo{

/* Subtract the `struct timeval' values X and Y,
 * storing the result in RESULT.
 *
 * Return 1 if the difference is negative, otherwise 0.
 */
  static int timeval_subtract(struct timeval *result, struct timeval *x, struct timeval *y) {
    /* Perform the carry for the later subtraction by updating y. */
    if (x->tv_usec < y->tv_usec) {
      int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;

      y->tv_usec -= 1000000 * nsec;
      y->tv_sec += nsec;
    }

    if (x->tv_usec - y->tv_usec > 1000000) {
      int nsec = (x->tv_usec - y->tv_usec) / 1000000;

      y->tv_usec += 1000000 * nsec;
      y->tv_sec -= nsec;
    }

    /* Compute the time remaining to wait.
     * tv_usec is certainly positive.
     */
    result->tv_sec = x->tv_sec - y->tv_sec;
    result->tv_usec = x->tv_usec - y->tv_usec;

    /* Return 1 if result is negative. */
    return x->tv_sec < y->tv_sec;
  }


  IterativeStoragePartition::IterativeStoragePartition(ConflictSpec &cSpecRef)
    :_cSpecRef(cSpecRef)
  {}

  void IterativeStoragePartition::FindStorageHyperplanes(bool enumerate){
    // queue of partitioning problems
    vector<StoragePartition *> inqueue, outqueue;

    // create the seed problem and enqueue it
    StoragePartition *sPtr = new StoragePartition(GetConflictSpecRef());
    inqueue.insert(inqueue.begin(),sPtr);

    StoragePartition *elem=inqueue.back();
    inqueue.pop_back();

    int storagePartitionCount=0;
    do{ // while the  queue is not empty
      cout << "--------------- STORAGE PARTITION #" << storagePartitionCount++ << "-------------------" << endl;

      StoragePartition &s=*elem;

      // then we iteratively satisfy all the conflicts
      int iterCount=0;
      struct timeval start,end;

      cout << "FINDING STORAGE HYPERPLANES" << endl;
      gettimeofday (&start,NULL);

      while(!s.AreAllConflictsSatisfied()){
	ConflictSpec &cSpecRef = s.GetConflictSpecRef();

	vector<StoragePartition *> altStoragePartitionVec;
	s.FindStorageHyperplane(enumerate,altStoragePartitionVec);

	for(int stmtId=0; stmtId<cSpecRef.GetNStmts(); ++stmtId){
	  cout << "Partial transformation matrix #" << stmtId << " :\n";
	  ISL::PrintMatrix(s.GetTransformMatPtr(stmtId));
	  cout << "Partial modulo matrix #" << stmtId << " :\n";
	  ISL::PrintMatrix(s.GetModuloMatPtr(stmtId));
	}

	inqueue.insert(inqueue.begin(),altStoragePartitionVec.begin(),altStoragePartitionVec.end());

	//cSpecRef.Destroy(); !!! Somashekar...ctx will get freed!
	++iterCount;
      }

      gettimeofday (&end,NULL);
      struct timeval elapsedTime;
      timeval_subtract(&elapsedTime,&end,&start);
      cout << "Start Time: \t" << start.tv_sec << "s " << start.tv_usec << "ms " << endl;
      cout << "End Time: \t" << end.tv_sec << "s " << end.tv_usec << "ms " << endl;
      cout << "Elapsed Time: \t" << elapsedTime.tv_sec << "s " << elapsedTime.tv_usec << "ms " << endl;
      cout << "Time Taken: \t" << (elapsedTime.tv_sec*1000000.0+elapsedTime.tv_usec)/1000000.0 << "s" <<endl;
	            // << 10e-6*elapsedTime.tv_usec << " seconds" << endl;
      
      outqueue.push_back(elem);

      if(inqueue.size()!=0){
	elem=inqueue.back();
	inqueue.pop_back();
      }
      else
	elem=NULL;

    } while(NULL!=elem);
    
    cout << endl << "Storage mapping(s) found:" << endl;
    cout << "---------------------------" << endl;
    for(vector<StoragePartition *>::iterator iter=outqueue.begin(); iter!=outqueue.end(); ++iter){
      StoragePartition &s=**iter;
      for(int stmtId=0; stmtId<s.GetConflictSpecRef().GetNStmts(); ++stmtId){
	s.PrintStorageMapping(stmtId);
        // cout << "Partial transformation matrix #" << stmtId << " :\n";
        // ISL::PrintMatrix(s.GetTransformMatPtr(stmtId));
        // cout << "Partial modulo matrix #" << stmtId << " :\n";
        // ISL::PrintMatrix(s.GetModuloMatPtr(stmtId));
      }
      if(s.GetConflictSpecRef().GetNStmts()>1)
        cout << endl;
    }
  }
}

