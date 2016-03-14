#include <stdio.h>

#include <unistd.h>
#include <sys/time.h>


/* #ifndef N */
/* #define N 2048 */
/* #endif //N */
unsigned N;

#ifndef B
#define B 64
#endif //B

#ifdef TIME
#define IF_TIME(foo) foo;
#else
#define IF_TIME(foo)
#endif

#define in(x,y) IN[tx*B+x][ty*B+y]
#define blurx(x,y) BLURX[tx*B+x][ty*B+y]
#define out(x,y) OUT[tx*B+x][ty*B+y]

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

/* double IN[N+2][N]; */
/* double BLURX[N][N]; */
/* double OUT[N][N]; */
double **IN;
double **BLURX;
double **OUT;


double t_start, t_end;

void init_array()
{
  int i, j;

  for (i=0; i<N+2; i++) {
    for (j=0; j<N; j++) {
      IN[i][j] = (i + j);
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
  int trial;

  IN = (double **)malloc((N+2) * sizeof(double *));
  for (int i=0; i<N+2; i++)
    IN[i] = (double *)malloc((N) * sizeof(double));

  BLURX = (double **)malloc((N) * sizeof(double *));
  for (int i=0; i<N; i++)
    BLURX[i] = (double *)malloc((N) * sizeof(double));

  OUT = (double **)malloc((N) * sizeof(double *));
  for (int i=0; i<N; i++)
    OUT[i] = (double *)malloc((N) * sizeof(double));

  init_array();

#ifdef PERFCTR
  PERF_INIT; 
#endif

  IF_TIME(t_start = rtclock());

  for (trial = 0; trial < 10 ; ++trial)
    {

#pragma scop
      for(ty = 0; ty <= floord((N-1),B); ++ty)
#pragma omp parallel for private(tx,x,y)
	for(tx = 0; tx <= floord((N-1),B); ++tx){
	  for(x = 0; x <= B-1; ++x){	
	    for(y = 0; y <= B-1; ++y)
	      /* P */     blurx(x,y)=in(x,y)+in((x+1),y)+in((x+2),y);
	    for(y = 0; y <= B-1; ++y)
	      if(ty*B+y>=2)
		/* Q */        out(x,y)=blurx(x,y)+blurx(x,(y-1))+blurx(x,(y-2));
	  }
	}
#pragma endscop
    }

  IF_TIME(t_end = rtclock());
  IF_TIME(fprintf(stderr, "File:%s \t\t N=%d,T=%d \t Runtime=%0.6lfs\n", __FILE__, N, B, (t_end - t_start)/10));

#ifdef PERFCTR
  PERF_EXIT; 
#endif


#ifdef VERIFY
  for(x = 0; x <= N-1; ++x)
    for(y = 0; y <= N-1; ++y)
      BLURX[x][y]=IN[x][y]+IN[x+1][y]+IN[x+2][y];
  // Stage 2: vertical blur
  for(x = 0; x <= N-1; ++x)
    for(y = 2; y <= N-1; ++y)
      {
	if(OUT[x][y] != BLURX[x][y]+BLURX[x][y-1]+BLURX[x][y-2])
	  {
	    printf("Difference at (%d, %d) : %f versus %f\n", x, y, OUT[x][y], BLURX[x][y]+BLURX[x][y-1]+BLURX[x][y-2]);
	    break;
	  }
      }
#endif

  if (fopen(".test", "r")) {
    // print_array();
  }

  return 0;
}
