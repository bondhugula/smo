#include <stdio.h>

#include <unistd.h>
#include <sys/time.h>
#include <math.h>

#ifdef TIME
#define IF_TIME(foo) foo;
#else
#define IF_TIME(foo)
#endif


#ifndef M
#define M 256
#endif

#ifndef B
#define B 32
#endif

#define in(k,x,y) IN[k][tx*B+x][ty*B+y]
#define blury(k,x,y) BLURY[k][tx*B+x][ty*B+y]
#define sharpen(k,x,y) SHARPEN[k][tx*B+x][ty*B+y]
#define mask(k,x,y) MASK[k][tx*B+x][ty*B+y]

// As P is enclosed by 4 loops, after total expansion of the 
// array blurx, the write access to blurx becomes a 4-d access
#define blurx(k,i,j) ((j>=0)? A[k][tx][ty % 2][i][j] : A[k][tx][(ty-1) % 2][i][B+j])
#define blurx_pos(k,i,j) A[k][tx][ty % 2][i][j]

#ifdef VERIFY
#define blurx_verify(k,x,y) BLURXV[k][x][y]
#define blury_verify(k,x,y) BLURYV[k][x][y]
#endif

// As P is enclosed by 4 loops, after total expansion of the 
// array blurx, the write access to blurx becomes a 4-d access
// #define blurx(x,y) (y>0)? A[tx,ty,x,y] : A[tx,ty-1,x,B-y]


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

//#include "decls.h"

#ifdef PERFCTR
#include "papiStdEventDefs.h"
#include <papi.h>
#include "papi_defs.h"
#endif

// #include "util.h"

float SHARPEN[3][M][M];
float MASK[3][M][M];

double IN[3][M+4][M];
double A[3][M/B][2][B][B];
double BLURY[3][M][M];

#ifdef VERIFY
float sharpen_verify[3][M][M];
#define in_verify(k,x,y) IN[k][x][y]
double BLURXV[3][M][M];
double BLURYV[3][M][M];
#endif

double t_start, t_end;

void init_array()
{
  int i, j;

  for (i=0; i<M+2; i++) {
    for (j=0; j<M; j++) {
      IN[0][i][j] = (i + j);
      IN[1][i][j] = (i + j);
      IN[2][i][j] = (i + j);
    }
  }
}

double rtclock()
{
  struct timezone Tzp;
  struct timeval Tp;
  int stat;
  stat = gettimeofday (&Tp, &Tzp);
  if (stat != 0) printf("Error return from gettimeofday: %d",stat);
  return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}

int main()
{
  int tx, ty, x, y;
  int k, i, j;
  int trial;
  double thresh = 0.23432f;
  double weight = 0.23432f;
  double _ct3, _ct4, _ct5;

  init_array();

#ifdef PERFCTR
  PERF_INIT; 
#endif

  IF_TIME(t_start = rtclock());

  for (trial = 0; trial < 10 ; ++trial)
    {
#pragma scop
      for (k = 0; k <= 2; ++k)
	{
	  for(ty = 0; ty <= (M-1)/B; ++ty)
	    // #pragma omp parallel for private(tx,x,y)
	    for(tx = 0; tx <= (M-1)/B; ++tx){
	      for(x = 0; x <= B-1; ++x){	
		for(y = 0; y <= B-1; ++y)
		  blurx_pos(k,x,y) = (in(k,x,y) * 0.0625f) +
		    (in(k,(x+1),y) * 0.25f) +
		    (in(k,(x+2),y) * 0.375f) +
		    (in(k,(x+3),y) * 0.25f) +
		    (in(k,(x+4),y) * 0.0625f);
		for(y = 0; y <= B-1; ++y)
		  if(ty*B+y>=4)
		    {
		      blury(k,x,y) = (blurx(k,x,y) * 0.0625f) +
			(blurx(k,x,(y-1)) * 0.25f) +
			(blurx(k,x,(y-2)) * 0.375f) +
			(blurx(k,x,(y-3)) * 0.25f) +
			(blurx(k,x,(y-4)) * 0.0625f);
		      sharpen(k,x,y) = ((in(k,(x+2),(y-2)) * (1 + weight)) + (blury(k,x,y) * -(weight)));
		      _ct3 = in(k,(x+2),(y-2));
		      _ct4 = sharpen(k,x,y);
		      _ct5 = ((abs((in(k,(x+2),(y-2)) - blury(k,x,y))) < thresh)? _ct3: _ct4);
		      mask(k,x,y) = _ct5;
		    }

	      }
	    }
	}
#pragma endscop
    }

  IF_TIME(t_end = rtclock());
  IF_TIME(fprintf(stderr, "File:%s \t\t M=%d,T=%d \t Runtime=%0.6lfs\n", __FILE__, M, B, (t_end - t_start)/10));

#ifdef PERFCTR
  PERF_EXIT; 
#endif


#ifdef VERIFY
    
  for (k = 0; k <= 2; ++k)
    {
      for (i = 0; i <= M-1; ++i)
	{
	  for (j = 0; j <= M-1; j = ++j)
	    {
	      blurx_verify(k,i,j) = (in_verify(k,i,j) * 0.0625f) +
		(in_verify(k,(i+1),j) * 0.25f) +
		(in_verify(k,(i+2),j) * 0.375f) +
		(in_verify(k,(i+3),j) * 0.25f) +
		(in_verify(k,(i+4),j) * 0.0625f);
	    }
	}
    }
  for (k = 0; k <= 2; ++k)
    {
      for (i = 0; i <= M-1; ++i)
	{
	  for (j = 4; j <= M-1; j = ++j)
	    {
	      blury_verify(k,i,j) = (blurx_verify(k,i,j) * 0.0625f) +
		(blurx_verify(k,i,(j-1)) * 0.25f) +
		(blurx_verify(k,i,(j-2)) * 0.375f) +
		(blurx_verify(k,i,(j-3)) * 0.25f) +
		(blurx_verify(k,i,(j-4)) * 0.0625f);
	      if(blury_verify(k,i,j) != BLURY[k][i][j])
		{
		  printf("Blury Difference at %d,%d,%d %f != %f\n", k, i, j, blury_verify(k,i,j), BLURY[k][i][j]);
		  break;
		}
	    }
	}
    }
  for (k = 0; k <= 2; ++k)
    {
      for (i = 0; i <= M-1; ++i)
	{
	  for (j = 4; j <= M-1; j = ++j)
	    {
	      sharpen_verify[k][i][j] = ((in_verify(k,(i+2),(j-2)) * (1 + weight)) + (blury_verify(k,i,j) * -(weight)));
	      if(sharpen_verify[k][i][j] != SHARPEN[k][i][j])
		{
		  printf("Sharpen Difference at %d,%d,%d %f != %f\n", k, i, j, sharpen_verify[k][i][j], SHARPEN[k][i][j]);
		  break;
		}
	    }
	}
    }
  for (k = 0; k <= 2; ++k)
    {
      for (i = 0; i <= M-1; ++i)
	{
	  for (j = 4; j <= M-1; j = ++j)
	    {
	      _ct3 = in_verify(k,(i+2),(j-2));
	      _ct4 = sharpen_verify[k][i][j];
	      _ct5 = ((abs((in_verify(k,(i+2),(j-2)) - blury_verify(k,i,j))) < thresh)? _ct3: _ct4);
	      if(_ct5 != MASK[k][i][j])
		{
		  printf("Difference at %d,%d,%d \n", k, i, j);
		  break;
		}
	    }
	}
    }
#endif

  /*
    if (fopen(".test", "r")) {
    // print_array();
    }
  */

  return 0;
}
