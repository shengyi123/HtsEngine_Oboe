#ifndef INCLUDE_OEC_NTPU
#define INCLUDE_OEC_NTPU

/* allocate matrix */
void mtxalc(double ***a, int m, int n);

/* free matrix */
void mtxfree(double **a, int m, int n);

/*------------------------------------*/
/*      Matrix Inverse Algorithm      */
/*------------------------------------*/
void mtxinv(double **a, double **d, int ordmtx);

/*---------------------------------------------*/
/*      Matrix Multiplication Algorithm        */
/*  c         =     a        *       b         */
/*  m1 by n2        m1 by n1         m2 by n2  */
/*          n1 must be the same as m2          */
/*---------------------------------------------*/
void mtxmtp(double **c, double **a, double **b, int m1, int n1, int m2, int n2);

/*---------------------------------------------*/
/*      Matrix Substraction Algorithm          */
/*  c         =     a        -       b         */
/*  m by n        m by n         m by n        */
/*---------------------------------------------*/
void mtxsub(double **c, double **a, double **b, int m, int n);

/* print matrix element to stdout */
void mtxprt(double **c, int m, int n);

/* matrix copy */
/*  a <-- b    */
void mtxcpy(double **a, double **b, int m, int n);

/* matrix transpose */
void mtxtrsp(double **a, double **at, int m, int n);

/* Smooth_OEC */
void Smooth_OEC(double **a, double ***Rn, int *T, int N, int d, int *n_k, int K, double **x);

/* Smooth_Dyn_OEC */
void Smooth_Dyn_OEC(double *w, double *mu, double *b, double **a, double ***Rn, int *T, int N, int d, int *n_k, int K, double **x);

#endif