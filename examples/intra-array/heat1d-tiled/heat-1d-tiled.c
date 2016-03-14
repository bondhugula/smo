//#include <omp.h>
/* #define ceild(n,d)  ceil(((double)(n))/((double)(d))) */
/* #define floord(n,d) floor(((double)(n))/((double)(d))) */
/* #define max(x,y)    ((x) > (y)? (x) : (y)) */
/* #define min(x,y)    ((x) < (y)? (x) : (y)) */

/*
 * Discretized 1D heat equation stencil
 * Adapted from Pochoir test bench example
 *
 * Irshad Pananilath: irshad@csa.iisc.ernet.in
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>

//#define DEBUG

/*
 * N is the number of points
 * T is the number of timesteps
 */

#define B 32

/* #define N 8192 */
/* #define T 8192 */
unsigned N, T;

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


//double POINT[T+1][N+1];
//double point[T+1][N+1];
double **point;

/* double seed[N+2]; */
/* double output[N+2]; */
/* double POINT[3][((T+N)/B)+1][B][B]; */

/* #define point(a,b,c,d) (((c)==0)||((d)==0)||((d)==(N+1)) ? seed[d] : (((c)==T) ? output[d] : ((((c)<lb2)&&((d)+(c)<(B*t2))) ? map(a-2,b-1,B-lb2+c,B*t2-lbv+(d)) : (((c)<lb2) ? map(a-1,b,B-lb2+c,(d)) : (((d)+(c)<(B*t2)) ? map(a-1,b-1,c,B*t2-lbv+(d)) : map(a,b,(c),(d))))))) */
/* #define map(a,b,c,d) POINT[(a)%3][b][(c)%B][(d)%B] */

/* #define andpoint(a,b,c,d) (((c)==0)||((d)==0)||((d)==(N+1)) ? &seed[d] : (((c)==T) ? &output[d] : ((((c)<lb2)&&((d)+(c)<(B*t2))) ? MAP(a-2,b-1,B-lb2+c,B*t2-lbv+(d)) : (((c)<lb2) ? MAP(a-1,b,B-lb2+c,(d)) : (((d)+(c)<(B*t2)) ? MAP(a-1,b-1,c,B*t2-lbv+(d)) : MAP(a,b,(c),(d))))))) */
/* #define MAP(a,b,c,d) &POINT[(a)%3][b][(c)%B][(d)%B] */

void print_array()
{
  int j;

  for (j=1; j<=N; j++) {
    //fprintf(stderr, "%lf ", output[j]);
    fprintf(stderr, "%lf ", point[T][j]);
    if (j%80 == 20) fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");
}




int main(int argc, char * argv[]) {
  long int t, i;

  point= (double **)malloc((T+1) * sizeof(double *));
  for (i=0; i<T+1; i++)
    point[i] = (double *)malloc((N+1) * sizeof(double));


  const int BASE = 1024;
  //  double point[2][N+1];// declaration here overflows the activation record stack. Declare globally

  // for timekeeping
  int trial;
  int ts_return = -1;
  struct timeval start, end, result;
  double tdiff = 0.0;

  // printf("Number of points = %ld\nNumber of timesteps = %ld\n", N, T);

  // seed with a constant value to verify results
  srand(57);

  for (i = 0; i < N+2; i++) {
    point[0][i] = (double)(i+1) / 100; 
  }

  IF_TIME(t_start = rtclock());

  for (trial = 0; trial < 10 ; ++trial)
    {



      int t1, t2, t3, t4;
      int lb1,ub1,lb2,ub2,lbv,ubv;
	  /* lb1=max(ceild(t1,2),ceild(B*t1-T,B)); */
	  /* ub1=min(min(floord(T+N,B),floord(B*t1+N+(B-1),(2*B))),t1); */
	    /* lb2=max(max(1,B*t1-B*t2),B*t2-N); */
	    /* ub2=min(min(T,B*t2+(B-2)),B*t1-B*t2+(B-1)); */
		/* lbv=max(B*t2,t3+1); ubv=min(B*t2+(B-1),t3+N); */

      /* Generated from PLUTO-produced CLooG file by CLooG 0.16.2-20-ga0c4004 gmp bits in 0.01s. */
#pragma scop
     // if ((N >= 1) && (T >= 1)) {
     if ((N >= 1) && (N <= 8192) && (T <= 8192) && (T >= 1)) {
	for (t1=0;t1<=floord(2*T+N,B);t1++){
	  /* lb1=max(ceild(t1,2),ceild(B*t1-T,B)); */
	  ub1=min(min(floord(T+N,B),floord(B*t1+N+(B-1),(2*B))),t1);
	  for (t2=max(ceild(t1,2),ceild(B*t1-T,B));t2<=ub1;t2++) {
	    for (t3=max(max(1,B*t1-B*t2),B*t2-N);t3<=min(min(T,B*t2+(B-2)),B*t1-B*t2+(B-1));t3++) {
	      for (t4=max(B*t2,t3+1); t4<=min(B*t2+(B-1),t3+N); t4++) {
		point[t3][-t3+t4]=point[t3-1][-t3+t4+1]-2.0*point[t3-1][-t3+t4]+point[t3-1][-t3+t4-1];
	      }
	    }
	  }
	}
      }
#pragma endscop
      /* End of CLooG code */

    }

  IF_TIME(t_end = rtclock());
  IF_TIME(fprintf(stdout, "File : %s \t N=%d,B=%d \t Runtime = %0.6lfs\n", __FILE__, N, B, (t_end - t_start)/10));


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
