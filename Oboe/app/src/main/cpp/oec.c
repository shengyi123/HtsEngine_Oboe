#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "oec.h"


/* allocate matrix */
void mtxalc(double ***a, int m, int n)
{
	double **p;
	int i;
	// int j;
	(*a) = (double **)calloc(sizeof(double *), m);
	p = (*a);
	if( p==NULL )
	{
		fprintf(stderr, "Cannot allocate a %d by %d matrix\n", m, n);
		exit(1);
	}
	for(i=0;i<m;i++)
	{
		p[i] = (double *)calloc(sizeof(double), n);
		if( p[i]==NULL )
		{
			fprintf(stderr, "Cannot allocate a %d by %d matrix\n", m, n);
			exit(1);
		}
	}
}

/* free matrix */
void mtxfree(double **a, int m, int n)
{
	int i;
	for(i=0;i<m;i++)
	{
		free(a[i]);
	}
	free(a);
}


/*------------------------------------*/
/*      Matrix Inverse Algorithm      */
/*------------------------------------*/
void mtxinv(double **a, double **d, int ordmtx)
{
        int i,j,k;
        double   b, c, tmp;
		double **cp;
		mtxalc(&cp, ordmtx, ordmtx);
		mtxcpy(cp, a, ordmtx, ordmtx);
		
        for (i=0;i<ordmtx;i++)
            for (j=0;j<ordmtx;j++)
                if (i==j) d[i][i]=1.0;
                else d[i][j]=0.0;
        for (i=0;i<ordmtx;i++)
        {
            label1: if (a[i][i]!=0.0) /* 1 */
                    {
                        b=a[i][i];
                        for (j=0;j<ordmtx;j++)
                        {
                            a[i][j]/=b;
                            d[i][j]/=b;
                        } /* for j */
                        for (k=0;k<ordmtx;k++)
                            if (k!=i) /* 2 */
                            {
                                c=a[k][i]/a[i][i];
                                for (j=0;j<ordmtx;j++)
                                {
                                    a[k][j]-=c*a[i][j];
                                    d[k][j]-=c*d[i][j];
									//printf("**%d %d %f\n",k,j,d[k][j]);

                                } /* for j */
                            } /* if 2 */
                    } /* if 1 */
                    else
                    {
                        for (k=i;k<ordmtx;k++)
                        {
                            if (a[k][i]!=0.0) /* 3 */
                            {
                                for (j=0;j<ordmtx;j++)
                                {
                                    tmp=a[i][j];
                                    a[i][j]=a[k][j];
                                    a[k][j]=tmp;
                                    tmp=d[i][j];
                                    d[i][j]=d[k][j];
                                    d[k][j]=tmp;
									//printf("****%d %d %f\n",k,j,d[k][j]);
                                } /* for j */
                                goto label1;
                            } /* if 3 */
                        } /* for k */
                    } /* end else */
        } /* for i */

		   /*     for (i=0;i<ordmtx;i++)
            for (j=0;j<ordmtx;j++)
				if(i>=412)
				printf("****%d %d %f\n",i,j,d[i][j]);
*/
		mtxcpy(a, cp, ordmtx, ordmtx);
		mtxfree(cp, ordmtx, ordmtx);
} /* end proc mtxinv */


/*---------------------------------------------*/
/*      Matrix Multiplication Algorithm        */
/*  c         =     a        *       b         */
/*  m1 by n2        m1 by n1         m2 by n2  */
/*          n1 must be the same as m2          */
/*---------------------------------------------*/
void mtxmtp(double **c, double **a, double **b, int m1, int n1, int m2, int n2)
{
	int i, j, k;	
	if( n1 != m2 )
	{
		return;
	}
	/* matrix initialization*/
	for(i=0;i<m1;i++)
	{
		for(j=0;j<n2;j++)
		{
			c[i][j] = 0;
		}
	}
	/* element obtaining */
	for(i=0;i<m1;i++)
	{
		for(j=0;j<n2;j++)
		{
			for(k=0;k<n1;k++)
			{
				c[i][j] = c[i][j] + a[i][k] * b[k][j];
			//	printf("%e\n", c[i][j]);
			}
		}		
	}
}


/*---------------------------------------------*/
/*      Matrix Substraction Algorithm          */
/*  c         =     a        -       b         */
/*  m by n        m by n         m by n        */
/*---------------------------------------------*/
void mtxsub(double **c, double **a, double **b, int m, int n)
{
	int i, j;
	for(i=0;i<m;i++) {
		for(j=0;j<n;j++) {
			c[i][j] = a[i][j] - b[i][j];
		}
	}
}

/* print matrix element to stdout */
void mtxprt(double **c, int m, int n)
{
	int i, j;
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			fprintf(stdout, "%e ", c[i][j]);
		}
		fprintf(stdout, "\n");
	}
}
/* matrix copy */
/*  a <-- b    */
void mtxcpy(double **a, double **b, int m, int n)
{
	int i, j;
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			a[i][j] = b[i][j];
		}
	}
}

/* matrix transpose */
void mtxtrsp(double **a, double **at, int m, int n)
{
	int i, j;
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			at[j][i] = a[i][j];
		}
	}
}


double OEC_phi(int j, int i, int N)
{
	double x, y, z, a, b;
	switch(j)
	{
	case 0:
		return 1;
		break;

	case 1:
		if( N<1 )
			return 0;
		x = pow(12,0.5);
		y = pow(((double)N), 0.5);
		z = pow(((double)(N+2)), 0.5);
		a = ( ((double)i) / ((double)N)) - 0.5;
		b = x * y / z * a;
		return b;

		return ( pow( (double)12 * (double)(N) / (double)(N + 2) , (double)0.5 ) 
		* ( (double)i / (double)N - (double)0.5 ) );

		break;

	case 2:
		if( N<2 )
			return 0;
		return pow( ((double)180 * pow( (double)(N), (double)3 ) 
		/ (double)(N - 1) / (double)(N + 2) / (double)(N + 3))     , (double)0.5 )
		* ( pow((double)i/(double)N,(double)2 ) - (double)i / (double)N + (double)(N-1) / (double)(6*N));

		break;

	case 3:
		{
			double N1,N2,N3,N5;
			double I1,I2,I3;

			N1 = (double) N;
			N2 = pow(N,2);
			N5 = pow(N,5);
			N3 = pow(N,3);
			I1 = (double) i;
			I2 = pow(I1,2);
			I3 = pow(I1,3);

			if( N<3 )
				return 0;
			return 

				pow(2800*N5/(N1-1)/(N1-2)/(N1+2)/(N1+3)/(N1+4),0.5) *
				( I3/N3 - (double)3/(double)2*I2/N2 + (6*N2-3*N1+2)/(double)10/N3*I1 - (N1-1)*(N1-2)/(double)20/N2);



			
			break;
		}
	default:
		break;
	}
	return 0;
}

void OEC(int start, int end, double *fea, double *result)
{
	int N;
	int i,j;
	N = end - start ;

	for(j=0;j<=3;j++)
	{
		result[j] = 0;

		if( j > N )
		{
			continue;
		}
		for(i=0;i<=N;i++)
		{
			result[j] += fea[start+i] * OEC_phi(j,i,N);
		}
		result[j] = result[j] / (double)(N + 1);
	}

}


double OEC_MSE(int start, int end, double *original, double *coe)
{
	int N;
	int i;
	double error;
	double orthexp;
	N = end - start;

	error = 0;
	for(i=0;i<=N;i++)
	{
		orthexp = coe[0] * OEC_phi(0,i,N) + coe[1] * OEC_phi(1,i,N) + coe[2] * OEC_phi(2,i,N) + coe[3] * OEC_phi(3,i,N);
		error += ( original[start+i] - orthexp ) 
			* ( original[start+i] - orthexp );
	}

	return error / (double)(N+1);

}

void OEC_Expansion_Pitch(int n_smp, double *pitch, double *a)
{
	int N;
	int i;
	N = n_smp-1;
	for(i=0;i<=N;i++)
	{
		pitch[i] = a[0] * OEC_phi(0,i,N) + a[1] * OEC_phi(1,i,N) + a[2] * OEC_phi(2,i,N) + a[3] * OEC_phi(3,i,N);		
//		pitch[i] = a[0] * phi(0,i,N) + a[1] * phi(1,i,N);	
	}
}

/* Smooth_OEC */
void Smooth_OEC(double **a, double ***Rn, int *T, int N, int d, int *n_k, int K, double **x)
{
	
	double **R;
	double **R_trsp;
	double **R_trsp_inv;
	double **Rn_inv;
	double **Z;
	double **Z_trsp;
	double **R_trsp_inv_Z_trsp;
	double **Z_R_trsp_inv_Z_trsp;
	double **Z_R_trsp_inv_Z_trsp_inv;
	double **Za;
	double **Z_R_trsp_inv_Z_trsp_inv_Za;
	double **R_trsp_inv_Z_trsp_Z_R_trsp_inv_Z_trsp_inv_Za;
	int n, i ,j, k;
	mtxalc(&R, N*d, N*d);
	mtxalc(&R_trsp, N*d, N*d);
	mtxalc(&R_trsp_inv, N*d, N*d);
	mtxalc(&Rn_inv, d, d);
	mtxalc(&Z, 2*K, d*N);
	mtxalc(&Z_trsp, d*N, 2*K);
	mtxalc(&R_trsp_inv_Z_trsp, d*N, 2*K);
	mtxalc(&Z_R_trsp_inv_Z_trsp, 2*K, 2*K);
	mtxalc(&Z_R_trsp_inv_Z_trsp_inv, 2*K, 2*K);
	mtxalc(&Za, 2*K, 1);
	mtxalc(&Z_R_trsp_inv_Z_trsp_inv_Za, 2*K, 1);
	mtxalc(&R_trsp_inv_Z_trsp_Z_R_trsp_inv_Z_trsp_inv_Za, d*N, 1);
	// make R
	for(n=0;n<N;n++) {
		mtxinv(Rn[n], Rn_inv, d);
		for(i=0;i<d;i++) {
			for(j=0;j<d;j++) {
				R[n*d+i][n*d+j] = Rn_inv[i][j];
			}
		}
	}
	mtxtrsp(R, R_trsp, N*d, N*d);
	mtxinv(R_trsp, R_trsp_inv, N*d);// R_trsp_inv: (R')^(-1)

	// make Z
	for(k=1;k<=K;k++) {
		if( (n_k[k])*d >= d*N ) {
			continue;
		}
		for(i=0;i<d;i++) {
			Z[ k-1 ][ (n_k[k]-1)*d + i ] += OEC_phi(i, T[n_k[k]]+1, T[n_k[k]]);
		}
		for(i=0;i<d;i++) {
			Z[ k-1 ][ (n_k[k])*d + i ] -= OEC_phi(i, 0, T[n_k[k]]);
		}
		for(i=0;i<d;i++) {
			Z[ k-1 + K ][ (n_k[k]-1)*d + i ] += OEC_phi(i, T[n_k[k]], T[n_k[k]]);
		}
		for(i=0;i<d;i++) {
			Z[ k-1 + K ][ (n_k[k])*d + i ] -= OEC_phi(i, -1, T[n_k[k]]);
		}
	}
	mtxtrsp(Z, Z_trsp, 2*K, d*N);// Z_trsp: Z'
	mtxmtp(R_trsp_inv_Z_trsp, R_trsp_inv, Z_trsp, d*N, d*N, d*N, 2*K);// (R')^(-1) * Z'
	mtxmtp(Z_R_trsp_inv_Z_trsp, Z, R_trsp_inv_Z_trsp, 2*K, d*N, d*N, 2*K);// Z * (R')^(-1) * Z'
	mtxinv(Z_R_trsp_inv_Z_trsp, Z_R_trsp_inv_Z_trsp_inv, 2*K);// (Z * (R')^(-1) * Z')^(-1)
	mtxmtp(Za, Z, a, 2*K, d*N, d*N , 1);// Z * a 
	mtxmtp(Z_R_trsp_inv_Z_trsp_inv_Za,Z_R_trsp_inv_Z_trsp_inv, Za, 2*K, 2*K, 2*K, 1);
	mtxmtp(R_trsp_inv_Z_trsp_Z_R_trsp_inv_Z_trsp_inv_Za, R_trsp_inv_Z_trsp, Z_R_trsp_inv_Z_trsp_inv_Za, d*N, 2*K, 2*K, 1);
	mtxsub(x, a, R_trsp_inv_Z_trsp_Z_R_trsp_inv_Z_trsp_inv_Za, d*N, 1);

	mtxfree(R, N*d, N*d);
	mtxfree(R_trsp, N*d, N*d);
	mtxfree(R_trsp_inv, N*d, N*d);
	mtxfree(Rn_inv, d, d);
	mtxfree(Z, 2*K, d*N);
	mtxfree(Z_trsp, d*N, 2*K);
	mtxfree(R_trsp_inv_Z_trsp, d*N, 2*K);
	mtxfree(Z_R_trsp_inv_Z_trsp, 2*K, 2*K);
	mtxfree(Z_R_trsp_inv_Z_trsp_inv, 2*K, 2*K);
	mtxfree(Za, 2*K, 1);
	mtxfree(Z_R_trsp_inv_Z_trsp_inv_Za, 2*K, 1);
	mtxfree(R_trsp_inv_Z_trsp_Z_R_trsp_inv_Z_trsp_inv_Za, d*N, 1);
}

/* Smooth_Dyn_OEC */
void Smooth_Dyn_OEC(double *w, double *mu, double *b, double **a, double ***Rn, int *T, int N, int d, int *n_k, int K, double **x)
{
	int MAX_ITER_NUM = 50000;
	double **R;
	double **R_trsp;
	double **R_trsp_inv;
	double **Rn_inv;
	double **Z;
	double **Z_trsp;
	double **R_trsp_inv_Z_trsp;
	double **Z_R_trsp_inv_Z_trsp;
	double **Z_R_trsp_inv_Z_trsp_inv;
	double **Za;
	double **Z_R_trsp_inv_Z_trsp_inv_Za;
	double **R_trsp_inv_Z_trsp_Z_R_trsp_inv_Z_trsp_inv_Za;


	double **Cj[4];
	double **Cj_trsp[4];
//	double **I;
//	double **e;
//	double **e_trsp;
//	double **ee_trsp;
	double **I_minus_ee_trsp;
	double **Dj[4];
	double **Dj_trsp[4];
	double **I_minus_ee_trsp_Cj[4];
	double sj[4];// sj[j] = 0 if sj=-1, sj[j] = 1 if sj=1, 
	double **vj[4];
	double **Dj_x;
	double **x_trsp;
	double **Q[16]; // Q[j]: j = sj[0] + sj[1] * 2 + sj[2] * 4 + sj[3] * 8
	double **Q_inv[16]; 
	double **neg_Lambda[16];
	double **neg_Lambda_trsp[16];
	double **Z_Q_inv_Z_trsp[16];
	double **Z_Q_inv_Z_trsp_inv[16];
	double **Q_inv_Z_trsp;
	double **R_trsp_a;
	double **Z_trsp_neg_Lambda[16];
	double **Qx;
	double **g_U;
	double **Q_now;
	double **Z_trsp_neg_Lambda_now;
	int n, i ,j, k, u, v, idx;
	double g_U_norm;
	double eps;// learing rate;
	double rmd;
	double **x_minus_a;
	double **x_minus_a_trsp;
	double **R_x_minus_a;
	double **loglike;
	double log_likelihood, old_log_likelihood;
	double log_var_likelihood, old_log_var_likelihood;
	double lambda_likelihood, old_lambda_likelihood;
	double **Zx;
	double **neg_Lambda_trsp_now;
	double **Lambda_trsp_Zx;
	double Zx_2norm;
	mtxalc(&R, N*d, N*d);
	mtxalc(&R_trsp, N*d, N*d);
	mtxalc(&R_trsp_inv, N*d, N*d);
	mtxalc(&Rn_inv, d, d);
	mtxalc(&Z, 2*K, d*N);
	mtxalc(&Z_trsp, d*N, 2*K);
	mtxalc(&R_trsp_inv_Z_trsp, d*N, 2*K);
	mtxalc(&Z_R_trsp_inv_Z_trsp, 2*K, 2*K);
	mtxalc(&Z_R_trsp_inv_Z_trsp_inv, 2*K, 2*K);
	mtxalc(&Za, 2*K, 1);
	mtxalc(&Z_R_trsp_inv_Z_trsp_inv_Za, 2*K, 1);
	mtxalc(&R_trsp_inv_Z_trsp_Z_R_trsp_inv_Z_trsp_inv_Za, d*N, 1);
	// make R
	for(n=0;n<N;n++) {
		mtxinv(Rn[n], Rn_inv, d);
		for(i=0;i<d;i++) {
			for(j=0;j<d;j++) {
				R[n*d+i][n*d+j] = Rn_inv[i][j];
			}
		}
	}
	mtxtrsp(R, R_trsp, N*d, N*d);
	mtxinv(R_trsp, R_trsp_inv, N*d);// R_trsp_inv: (R')^(-1)

	// make Z
	for(k=1;k<=K;k++) {
		if( (n_k[k])*d >= d*N ) {
			continue;
		}
		for(i=0;i<d;i++) {
			Z[ k-1 ][ (n_k[k]-1)*d + i ] += OEC_phi(i, T[n_k[k]]+1, T[n_k[k]]);
		}
		for(i=0;i<d;i++) {
			Z[ k-1 ][ (n_k[k])*d + i ] -= OEC_phi(i, 0, T[n_k[k]]);
		}
		for(i=0;i<d;i++) {
			Z[ k-1 + K ][ (n_k[k]-1)*d + i ] += OEC_phi(i, T[n_k[k]], T[n_k[k]]);
		}
		for(i=0;i<d;i++) {
			Z[ k-1 + K ][ (n_k[k])*d + i ] -= OEC_phi(i, -1, T[n_k[k]]);
		}
	}
	mtxtrsp(Z, Z_trsp, 2*K, d*N);// Z_trsp: Z'
	mtxmtp(R_trsp_inv_Z_trsp, R_trsp_inv, Z_trsp, d*N, d*N, d*N, 2*K);// (R')^(-1) * Z'
	mtxmtp(Z_R_trsp_inv_Z_trsp, Z, R_trsp_inv_Z_trsp, 2*K, d*N, d*N, 2*K);// Z * (R')^(-1) * Z'
	mtxinv(Z_R_trsp_inv_Z_trsp, Z_R_trsp_inv_Z_trsp_inv, 2*K);// (Z * (R')^(-1) * Z')^(-1)
	mtxmtp(Za, Z, a, 2*K, d*N, d*N , 1);// Z * a 
	mtxmtp(Z_R_trsp_inv_Z_trsp_inv_Za,Z_R_trsp_inv_Z_trsp_inv, Za, 2*K, 2*K, 2*K, 1);
	mtxmtp(R_trsp_inv_Z_trsp_Z_R_trsp_inv_Z_trsp_inv_Za, R_trsp_inv_Z_trsp, Z_R_trsp_inv_Z_trsp_inv_Za, d*N, 2*K, 2*K, 1);
	mtxsub(x, a, R_trsp_inv_Z_trsp_Z_R_trsp_inv_Z_trsp_inv_Za, d*N, 1);




	/*-------------------- start to make dynamic ---------------------*/
	mtxalc(&R_trsp_a, d*N, 1);
	mtxalc(&Dj_x, d*N, 1); 
	mtxalc(&x_trsp, 1, d*N); 
	mtxalc(&Q_inv_Z_trsp, d*N, 2*K);
	mtxalc(&Qx, d*N, 1);
	mtxalc(&g_U, d*N, 1);
	mtxalc(&x_minus_a, d*N, 1);
	mtxalc(&x_minus_a_trsp, 1, d*N);
	mtxalc(&R_x_minus_a, d*N, 1);
	mtxalc(&loglike, 1, 1);
	mtxalc(&Zx, 2*K, d*N);
	mtxalc(&Lambda_trsp_Zx, 1,1);
	// make Cj
	for(j=0;j<d;j++) {
		mtxalc(&(Cj[j]), d*N, d*N);
		mtxalc(&(Cj_trsp[j]), d*N, d*N);
		mtxalc(&(I_minus_ee_trsp_Cj[j]), d*N, d*N);
		mtxalc(&(Dj[j]), d*N, d*N);
		mtxalc(&(Dj_trsp[j]), d*N, d*N);
		mtxalc(&(vj[j]), 1, 1);
		
		for(u=0;u<N;u++) {
			Cj[j][u*4+j][u*4+j] = 1.0 / sqrt((double)N);
		}
		mtxtrsp(Cj[j], Cj_trsp[j], N*d, N*d);
	}
	// make I_minus_ee_trsp
	mtxalc(&I_minus_ee_trsp, d*N, d*N);
	for(u=0;u<d*N;u++) {
		for(v=0;v<d*N;v++) {
			if( u==v ) {
				I_minus_ee_trsp[u][v] = 1.0 - 1.0/((double)N);
			}
			else {
				I_minus_ee_trsp[u][v] = - 1.0/((double)N);
			}
		}
	}
	// make I_minus_ee_trsp_Cj & Dj
	for(j=0;j<d;j++) {
		mtxmtp(I_minus_ee_trsp_Cj[j], I_minus_ee_trsp, Cj[j], d*N, d*N, d*N, d*N);
		mtxmtp(Dj[j], Cj_trsp[j], I_minus_ee_trsp_Cj[j], d*N, d*N, d*N, d*N);
		mtxtrsp(Dj[j], Dj_trsp[j], d*N, d*N);
	}

	for(j=0;j<16;j++) {
		mtxalc(&(Q[j]), d*N, d*N);
		mtxalc(&(Q_inv[j]), d*N, d*N);
		mtxalc(&(neg_Lambda[j]), 2*K, 1);
		mtxalc(&(neg_Lambda_trsp[j]), 1, 2*K);
		mtxalc(&(Z_Q_inv_Z_trsp[j]), 2*K, 2*K);
		mtxalc(&(Z_Q_inv_Z_trsp_inv[j]), 2*K, 2*K);
		mtxalc(&(Z_trsp_neg_Lambda[j]), d*N, 1);
	}
	
	

	// make Q
	for(j=0;j<16;j++) {
		idx = j;
		sj[3] = floor(((double)j) / 8.0);
		rmd = j - sj[3] * 8;
		sj[2] = floor(((double)rmd) / 4.0);
		rmd = rmd - sj[2] * 4;
		sj[1] = floor(((double)rmd) / 2.0);
		rmd = rmd - sj[1] * 2;
		sj[0] = floor(((double)rmd) / 1.0);
		for(u=0;u<d;u++) {
			if( sj[u]==0 ) {
				sj[u] = -1.0;
			}
		}
		// make Q = R' - sigma_j=0~3(w_j*s_j/b_j*Dj')
		for(u=0;u<d*N;u++) {
			for(v=0;v<d*N;v++) {
				Q[idx][u][v] = R_trsp[u][v];
				for(k=0;k<d;k++) {
					Q[idx][u][v] = Q[idx][u][v] + w[k] * sj[k] / b[k] * Dj[k][u][v];
				}
			}
		}	
	}


	// make Q^(-1), ZQ^(-1)Z', and -Lambda
	for(j=0;j<16;j++) {
		mtxinv(Q[j], Q_inv[j], N*d);
		mtxmtp(Q_inv_Z_trsp, Q_inv[j], Z_trsp, N*d, N*d, N*d, 2*K);
		mtxmtp(Z_Q_inv_Z_trsp[j], Z, Q_inv_Z_trsp, 2*K, d*N, d*N, 2*K); 
		mtxmtp(neg_Lambda[j], Z_Q_inv_Z_trsp[j], Za, 2*K, 2*K, 2*K, 1);
		mtxtrsp(neg_Lambda[i], neg_Lambda_trsp[i], 2*K, 1);
	}

	// make R'a
	mtxmtp(R_trsp_a, R_trsp, a, d*N, d*N, d*N, 1);

	// make -Z'Lambda
	for(j=0;j<16;j++) {
		mtxmtp(Z_trsp_neg_Lambda[j], Z_trsp, neg_Lambda[j], d*N, 2*K, 2*K, 1);
	}

	

	for(k=0;k<MAX_ITER_NUM;k++) {
		// calculate -0.5 * (x-a)'R(x-a)
		mtxsub(x_minus_a, x, a, d*N, 1);
		mtxtrsp(x_minus_a, x_minus_a_trsp, d*N, 1);
		mtxmtp(R_x_minus_a, R, x_minus_a, d*N, d*N, d*N, 1);
		mtxmtp(loglike, x_minus_a_trsp, R_x_minus_a, 1, d*N, d*N, 1);
		log_likelihood = loglike[0][0] * (-0.5);
		// make x_trsp & calculate wj/2/b(j)*|x'Djx-mu(j)|
		log_var_likelihood = 0;
		mtxtrsp(x, x_trsp, d*N, 1);
		// make sj
		for(j=0;j<d;j++) {
			// make Dj_x
			mtxmtp(Dj_x, Dj[j], x, d*N, d*N, d*N, 1);
			mtxmtp(vj[j], x_trsp, Dj_x, 1, d*N, d*N, 1);
			sj[j] = ( ( vj[j][0][0] - mu[j] ) > 0.0 ) * 1.0 - ( ( vj[j][0][0] - mu[j] ) <= 0.0 ) * 1.0;
			log_var_likelihood = log_var_likelihood - w[j] * 0.5 / b[j] * sj[j] * ( vj[j][0][0] - mu[j] );
		}
		// calculate -Lambda'Zx
		mtxmtp(Zx, Z, x, 2*K, d*N, d*N, 1);
		// neg_Lambda_now
		idx = (sj[0]>0)*1 + (sj[1]>0)*2 + (sj[2]>0)*4 + (sj[3]>0)*8;
		neg_Lambda_trsp_now = neg_Lambda_trsp[idx];
		mtxmtp(Lambda_trsp_Zx, neg_Lambda_trsp_now, Zx, 1, 2*K, 2*K, 1);
		lambda_likelihood = Lambda_trsp_Zx[0][0] * (-1.0);
		Zx_2norm = 0;
		for(j=0;j<2*K;j++) {
			Zx_2norm = Zx_2norm + Zx[j][0] * Zx[j][0];
		}
		


		// make Qx
		idx = (sj[0]>0)*1 + (sj[1]>0)*2 + (sj[2]>0)*4 + (sj[3]>0)*8;
		Q_now = Q[idx];
		Z_trsp_neg_Lambda_now = Z_trsp_neg_Lambda[idx];
		mtxmtp(Qx, Q_now, x, d*N, d*N, d*N, 1);

		// make gradient U
		g_U_norm = 0;
		for(j=0;j<d*N;j++) {
			g_U[j][0] = -Qx[j][0] + R_trsp_a[j][0] - Z_trsp_neg_Lambda_now[j][0];
		//	g_U[j][0] = -Qx[j][0] + R_trsp_a[j][0];
			g_U_norm = g_U_norm + g_U[j][0] * g_U[j][0];
		}

		printf("-------------------------------------\n");
		printf("iteration %d: gradient norm^2 = %lf\n", k , g_U_norm);
		printf("log_likelihood = %e\tlambda_likelihood = %e\tlog_var_likelihood = %e\n",
			log_likelihood, lambda_likelihood, log_var_likelihood);
		printf("total_log_likelihood = %e\n", log_likelihood + lambda_likelihood + log_var_likelihood);
		printf("Zx_2norm = %e\n", Zx_2norm);
		printf("%lf//%lf %lf//%lf %lf//%lf %lf//%lf\n", vj[0][0][0], mu[0], vj[1][0][0], mu[1], vj[2][0][0], mu[2], vj[3][0][0], mu[3]);
		printf("-------------------------------------\n");
		//system("pause");
		eps = 1 / g_U_norm;
		//if( g_U_norm < 10 )
		//	break;

		// update x
		for(j=0;j<d*N;j++) {
			x[j][0] = x[j][0] + eps * g_U[j][0];
		}

		


		

		old_log_likelihood = log_likelihood;
		old_log_var_likelihood = log_var_likelihood;
		old_lambda_likelihood = lambda_likelihood;
	}



	mtxfree(R, N*d, N*d);
	mtxfree(R_trsp, N*d, N*d);
	mtxfree(R_trsp_inv, N*d, N*d);
	mtxfree(Rn_inv, d, d);
	mtxfree(Z, 2*K, d*N);
	mtxfree(Z_trsp, d*N, 2*K);
	mtxfree(R_trsp_inv_Z_trsp, d*N, 2*K);
	mtxfree(Z_R_trsp_inv_Z_trsp, 2*K, 2*K);
	mtxfree(Z_R_trsp_inv_Z_trsp_inv, 2*K, 2*K);
	mtxfree(Za, 2*K, 1);
	mtxfree(Z_R_trsp_inv_Z_trsp_inv_Za, 2*K, 1);
	mtxfree(R_trsp_inv_Z_trsp_Z_R_trsp_inv_Z_trsp_inv_Za, d*N, 1);

	for(j=0;j<d;j++) {
		mtxfree(Cj[j], d*N, d*N);
		mtxfree(Cj_trsp[j], d*N, d*N);
		mtxfree(I_minus_ee_trsp_Cj[j], d*N, d*N);
		mtxfree(Dj[j], d*N, d*N);
		mtxfree(Dj_trsp[j], d*N, d*N);
		mtxfree(vj[j], 1, 1);		
	}
	for(j=0;j<16;j++) {
		mtxfree((Q[j]), d*N, d*N);
		mtxfree((Q_inv[j]), d*N, d*N);
		mtxfree((neg_Lambda[j]), 2*K, 1);
		mtxfree((Z_Q_inv_Z_trsp[j]), 2*K, 2*K);
		mtxfree((Z_Q_inv_Z_trsp_inv[j]), 2*K, 2*K);
		mtxfree((Z_trsp_neg_Lambda[j]), d*N, 1);
	}
	mtxfree(I_minus_ee_trsp, d*N, d*N);
	mtxfree(Dj_x, d*N, 1);
	mtxfree(x_trsp, 1, d*N); 
	mtxfree(R_trsp_a, d*N, 1);
	mtxfree(Qx, d*N, 1);
	mtxfree(g_U, d*N, 1);
	mtxfree(x_minus_a, d*N, 1);
	mtxfree(x_minus_a_trsp, 1, d*N);
	mtxfree(R_x_minus_a, d*N, 1);
	mtxfree(loglike, 1, 1);
	mtxfree(Zx, 2*K, d*N);
	mtxfree(Lambda_trsp_Zx, 1,1);
}