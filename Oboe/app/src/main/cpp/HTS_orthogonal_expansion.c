/* ----------------------------------------------------------------- */
/*           The Othogonal Expansion                                 */
/*           developed by Speech Processing Labratorying             */
/*           www.speech.cm.nctu.edu.tw                               */
/* ----------------------------------------------------------------- */
/*                                                                   */
/*  Copyright (c) 2001-2009  National Chiao Tung University, Taiwan  */
/*                           Department of Communication Engineering */
/*                                                                   */
/*                                                                   */
/* All rights reserved.                                              */
/*                                                                   */
/* ----------------------------------------------------------------- */
#include <math.h>
double phi(int j, int i, int N)
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

void Orthogonal_Expansion(int start, int end, double *fea, double *result)
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
			result[j] += fea[start+i] * phi(j,i,N);
		}
		result[j] = result[j] / (double)(N + 1);
	}

}


double MSE(int start, int end, double *original, double *coe)
{
	int N;
	int i;
	double error;
	double orthexp;
	N = end - start;

	error = 0;
	for(i=0;i<=N;i++)
	{
		orthexp = coe[0] * phi(0,i,N) + coe[1] * phi(1,i,N) + coe[2] * phi(2,i,N) + coe[3] * phi(3,i,N);
		error += ( original[start+i] - orthexp ) 
			* ( original[start+i] - orthexp );
	}

	return error / (double)(N+1);

}

void Expansion_Pitch(int n_smp, double *pitch, double *a)
{
	int N;
	int i;
	N = n_smp-1;
	for(i=0;i<=N;i++)
	{
		pitch[i] = a[0] * phi(0,i,N) + a[1] * phi(1,i,N) + a[2] * phi(2,i,N) + a[3] * phi(3,i,N);		
//		pitch[i] = a[0] * phi(0,i,N) + a[1] * phi(1,i,N);	
	}
}