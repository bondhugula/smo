/*
 * blur filter
 *
 * somashekaracharya <gbs@csa.iisc.ernet.in>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>

//#define DEBUG

/*
 * M is the height of the image
 * N is the width of the image
 */
/*
*/


#ifndef M
#define M 2048
#endif

#ifndef N
#define N 3072
#endif

#define TIME

#include <unistd.h>
#include <sys/time.h>

#ifdef TIME
#define IF_TIME(foo) foo;
#else
#define IF_TIME(foo)
#endif

double t_start, t_end;

double rtclock()
{
    struct timezone Tzp;
    struct timeval Tp;
    int stat;
    stat = gettimeofday (&Tp, &Tzp);
    if (stat != 0) printf("Error return from gettimeofday: %d",stat);
    return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}

/* Subtract the `struct timeval' values X and Y,
 * storing the result in RESULT.
 *
 * Return 1 if the difference is negative, otherwise 0.
 */
int timeval_subtract(struct timeval *result, struct timeval *x, struct timeval *y) {
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

double in[N+2][M];
double blurx[N][M];
double  out[N][M-2];

#if VERIFY
double A[N][N];
#endif

void print_array()
{
    int i, j;

    for(i = 0; i < M-2; ++i){
	    for (j= 0; j < N; ++j) {
		fprintf(stderr, "%f ", out[j][i]);
		if (j%80 == 20) fprintf(stderr, "\n");
	    }
	    fprintf(stderr, "\n");
    }
}

void print_inarray()
{
    int i, j;

    for(i = 0; i < M; ++i){
	    for (j= 0; j < N+2; ++j) {
		fprintf(stderr, "%f ", in[j][i]);
		if (j%80 == 20) fprintf(stderr, "\n");
	    }
	    fprintf(stderr, "\n");
    }
}


void foo(){
	int y,x,trial;
    IF_TIME(t_start = rtclock());

for (trial=0;trial<10;++trial)
{
#pragma scop
	for (y = 0; y <= M-1; ++y)
   		for(x = 0; x <= N-1; ++x) {
	         	blurx[y][x]=in[x][y]+in[x+1][y]+in[x+2][y];
			if (y >= 2)
		 		out[x][y-2]=blurx[y-2][x]+blurx[y-1][x]+blurx[y][x];
   		}
#pragma endscop
}


    IF_TIME(t_end = rtclock());
    IF_TIME(fprintf(stdout, "%s\t\t(M=%d,N=%d) \t %0.6lfs\n", __FILE__, M, N, (t_end - t_start)/trial));

#ifdef VERIFY
for(x = 0; x <= N-1; ++x)
	for(y = 0; y <= M-1; ++y)
		A[x][y]=in[x][y]+in[x+1][y]+in[x+2][y];
// Stage 2: vertical blur
for(x = 0; x <= N-1; ++x)
	for(y = 2; y <= M-1; ++y)
	{
		if(out[x][y-2] != A[x][y]+A[x][y-1]+A[x][y-2])
		{
			printf("blur-smo.c: Difference at (%d, %d) : %f versus %f\n", x, y, out[x][y-2], A[x][y]+A[x][y-1]+A[x][y-2]);
		}
	}
#endif
}

int main(int argc, char * argv[]) {
    int y, x, t, i;
    const int BASE = 1024;
    //  double point[2][N+1];// declaration here overflows the activation record stack. Declare globally

    // for timekeeping
    int ts_return = -1;
    struct timeval start, end, result;
    double tdiff = 0.0;

    // printf("Number of points = %ld\nNumber of timesteps = %ld\n", N, T);

    // seed with a constant value to verify results
    srand(57);

    for(y = 0; y < M; ++y ){
    	for (x = 0; x < N+2; ++x) {
       		in[x][y] = y*x; 
    	}
    }

    if (fopen(".test", "r")) {
#ifdef MPI
        if (my_rank == 0) {
            print_array();
        }
#else
        print_inarray();
#endif
    }

    foo();

    if (fopen(".test", "r")) {
#ifdef MPI
        if (my_rank == 0) {
            print_array();
        }
#else
        print_array();
#endif
    }

    return 0;
}
