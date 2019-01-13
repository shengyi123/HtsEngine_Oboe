/* ----------------------------------------------------------------- */
/*           The HMM-Based Speech Synthesis Engine "hts_engine API"  */
/*           developed by HTS Working Group                          */
/*           http://hts-engine.sourceforge.net/                      */
/* ----------------------------------------------------------------- */
/*                                                                   */
/*  Copyright (c) 2001-2011  Nagoya Institute of Technology          */
/*                           Department of Computer Science          */
/*                                                                   */
/*                2001-2008  Tokyo Institute of Technology           */
/*                           Interdisciplinary Graduate School of    */
/*                           Science and Engineering                 */
/*                                                                   */
/* All rights reserved.                                              */
/*                                                                   */
/* Redistribution and use in source and binary forms, with or        */
/* without modification, are permitted provided that the following   */
/* conditions are met:                                               */
/*                                                                   */
/* - Redistributions of source code must retain the above copyright  */
/*   notice, this list of conditions and the following disclaimer.   */
/* - Redistributions in binary form must reproduce the above         */
/*   copyright notice, this list of conditions and the following     */
/*   disclaimer in the documentation and/or other materials provided */
/*   with the distribution.                                          */
/* - Neither the name of the HTS working group nor the names of its  */
/*   contributors may be used to endorse or promote products derived */
/*   from this software without specific prior written permission.   */
/*                                                                   */
/* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND            */
/* CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,       */
/* INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF          */
/* MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE          */
/* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS */
/* BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,          */
/* EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED   */
/* TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,     */
/* DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON */
/* ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,   */
/* OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY    */
/* OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE           */
/* POSSIBILITY OF SUCH DAMAGE.                                       */
/* ----------------------------------------------------------------- */

#ifndef HTS_GSTREAM_C
#define HTS_GSTREAM_C

#ifdef __cplusplus
#define HTS_GSTREAM_C_START extern "C" {
#define HTS_GSTREAM_C_END   }
#else
#define HTS_GSTREAM_C_START
#define HTS_GSTREAM_C_END
#endif                          /* __CPLUSPLUS */

HTS_GSTREAM_C_START;

/* hts_engine libraries */

#include <time.h>
#include "HTS_hidden.h"
#
#include "LF0NormFct.h"
#include "oec.h"
#include "HTS_orthogonal_expansion.h"
#include "LOG.h"

#define R0 0.000442432880592355
#define R1 0.003277922825216
#define R2 0.001061948970553
#define R3 0.000358158046903697
/* HTS_GStreamSet_initialize: initialize generated parameter stream set */
void HTS_GStreamSet_initialize(HTS_GStreamSet * gss)
{
   gss->nstream = 0;
   gss->total_frame = 0;
   gss->total_nsample = 0;
   gss->gstream = NULL;
   gss->gspeech = NULL;
}

/* HTS_GStreamSet_create: generate speech */
/* (stream[0] == spectrum && stream[1] == lf0) */
HTS_Boolean HTS_GStreamSet_create(HTS_GStreamSet * gss, HTS_PStreamSet * pss, int stage, HTS_Boolean use_log_gain, int sampling_rate, int fperiod, double alpha, double beta, HTS_Boolean * stop, double volume, HTS_Audio * audio)
{
   int i, j, k;
   int msd_frame;
   HTS_Vocoder v;
   int nlpf = 0;
   double *lpf = NULL;

   /* check */
   if (gss->gstream || gss->gspeech) {
      HTS_error(1, "HTS_GStreamSet_create: HTS_GStreamSet is not initialized.\n");
      return FALSE;
   }

   /* initialize */
   gss->nstream = HTS_PStreamSet_get_nstream(pss);
   gss->total_frame = HTS_PStreamSet_get_total_frame(pss);
   gss->total_nsample = fperiod * gss->total_frame;
   gss->gstream = (HTS_GStream *) HTS_calloc(gss->nstream, sizeof(HTS_GStream));
   for (i = 0; i < gss->nstream; i++) {
      gss->gstream[i].static_length = HTS_PStreamSet_get_static_length(pss, i);
      gss->gstream[i].par = (double **) HTS_calloc(gss->total_frame, sizeof(double *));
      for (j = 0; j < gss->total_frame; j++)
         gss->gstream[i].par[j] = (double *) HTS_calloc(gss->gstream[i].static_length, sizeof(double));
   }
   HTS_wavebuffer = gss->gspeech = (short *) HTS_calloc(gss->total_nsample, sizeof(short));
   buffersize = gss->total_nsample;
   /* copy generated parameter */
   for (i = 0; i < gss->nstream; i++) {
      if (HTS_PStreamSet_is_msd(pss, i)) {      /* for MSD */
         for (j = 0, msd_frame = 0; j < gss->total_frame; j++)
            if (HTS_PStreamSet_get_msd_flag(pss, i, j)) {
               for (k = 0; k < gss->gstream[i].static_length; k++)
                  gss->gstream[i].par[j][k] = HTS_PStreamSet_get_parameter(pss, i, msd_frame, k);
               msd_frame++;
            } else
               for (k = 0; k < gss->gstream[i].static_length; k++)
                  gss->gstream[i].par[j][k] = LZERO;
      } else {                  /* for non MSD */
         for (j = 0; j < gss->total_frame; j++)
            for (k = 0; k < gss->gstream[i].static_length; k++)
               gss->gstream[i].par[j][k] = HTS_PStreamSet_get_parameter(pss, i, j, k);
      }
   }

   /* check */
   if (gss->nstream != 2 && gss->nstream != 3) {
      HTS_error(1, "HTS_GStreamSet_create: The number of streams should be 2 or 3.\n");
      HTS_GStreamSet_clear(gss);
      return FALSE;
   }
   if (HTS_PStreamSet_get_static_length(pss, 1) != 1) {
      HTS_error(1, "HTS_GStreamSet_create: The size of lf0 static vector should be 1.\n");
      HTS_GStreamSet_clear(gss);
      return FALSE;
   }
   if (gss->nstream >= 3 && gss->gstream[2].static_length % 2 == 0) {
      HTS_error(1, "HTS_GStreamSet_create: The number of low-pass filter coefficient should be odd numbers.");
      HTS_GStreamSet_clear(gss);
      return FALSE;
   }

   /* synthesize speech waveform */
   HTS_Vocoder_initialize(&v, gss->gstream[0].static_length - 1, stage, use_log_gain, sampling_rate, fperiod);
   if (gss->nstream >= 3)
      nlpf = (gss->gstream[2].static_length - 1) / 2;
   for (i = 0; i < gss->total_frame && (*stop) == FALSE; i++) {
      if (gss->nstream >= 3)
         lpf = &gss->gstream[2].par[i][0];
      HTS_Vocoder_synthesize(&v, gss->gstream[0].static_length - 1, gss->gstream[1].par[i][0], &gss->gstream[0].par[i][0], nlpf, lpf, alpha, beta, volume, &gss->gspeech[i * fperiod], audio);
   }
   HTS_Vocoder_clear(&v);
   if (audio)
      HTS_Audio_flush(audio);

   return TRUE;
}



/* HTS_GStreamSet_create_by_xlabel: generate speech given prthogonal expansion, added by Chne-Yu Chiang, 20140701 */
/* (stream[0] == spectrum && stream[1] == lf0) */
HTS_Boolean HTS_GStreamSet_create_by_xlabel(HTS_XLabel * xlabel, HTS_GStreamSet * gss, HTS_PStreamSet * pss, int stage, HTS_Boolean use_log_gain, int sampling_rate, int fperiod, double alpha, double beta, HTS_Boolean * stop, double volume, HTS_Audio * audio)
{

   int i, j, k, l;
   int start_idx, end_idx;
   int syl_idx;
   int lf0_len;
//   int pre_lf0_len;
   double *frame_lf0, *xxx_frame_lf0;
   HTS_XLabelString *tmp_xlabel;
   int msd_frame;
   HTS_Vocoder v;
   int nlpf = 0;
   double *lpf = NULL;
   /*20190113*/
  // On Smoothing and Increasing Dynamics of Discrete Orthogonal Polynomial Generated Pitch Contour
  // int N; // number of syllable
   //int K; // number of smooth juncture
   //int T[1024];// voiced frame number
   //   int n_k[1024]; // smooth juncture index array  /*20190113*/
   //int last_end_idx;
   // double **a; // matrix/vector for OEC /*20190113*/
   //   double ***Rn; // covariance matrix for each syllable logF0 /*20190113*/
   //   int d; // dimension for OEC   /*20190113*/
   //double **x;// smooth OEC
   //double w[4];
   //double mu[4];
   //double b[4];
   //int syl_map[1024];
   int clock1,clock2;



   /* check */
   if (gss->gstream || gss->gspeech) {
      HTS_error(1, "HTS_GStreamSet_create: HTS_GStreamSet is not initialized.\n");
      return FALSE;
   }

   /* initialize */
   gss->nstream = HTS_PStreamSet_get_nstream(pss);
   gss->total_frame = HTS_PStreamSet_get_total_frame(pss);
   gss->total_nsample = fperiod * gss->total_frame;
   gss->gstream = (HTS_GStream *) HTS_calloc(gss->nstream, sizeof(HTS_GStream));
   for (i = 0; i < gss->nstream; i++) {
      gss->gstream[i].static_length = HTS_PStreamSet_get_static_length(pss, i);
      gss->gstream[i].par = (double **) HTS_calloc(gss->total_frame, sizeof(double *));
      for (j = 0; j < gss->total_frame; j++)
         gss->gstream[i].par[j] = (double *) HTS_calloc(gss->gstream[i].static_length, sizeof(double));
   }
   gss->gspeech = (short *) HTS_calloc(gss->total_nsample, sizeof(short));

   /* copy generated parameter */
   for (i = 0; i < gss->nstream; i++) {
      if (HTS_PStreamSet_is_msd(pss, i)) {      /* for MSD */
         for (j = 0, msd_frame = 0; j < gss->total_frame; j++)
            if (HTS_PStreamSet_get_msd_flag(pss, i, j)) {
               for (k = 0; k < gss->gstream[i].static_length; k++)
                  gss->gstream[i].par[j][k] = HTS_PStreamSet_get_parameter(pss, i, msd_frame, k);
               msd_frame++;
            } else
               for (k = 0; k < gss->gstream[i].static_length; k++)
                  gss->gstream[i].par[j][k] = LZERO;
      } else {                  /* for non MSD */
         for (j = 0; j < gss->total_frame; j++)
            for (k = 0; k < gss->gstream[i].static_length; k++)
               gss->gstream[i].par[j][k] = HTS_PStreamSet_get_parameter(pss, i, j, k);
      }
   }


   /* orthogonal expansion lf0 */
   for (i = 1; i < 2; i++) {
      if (HTS_PStreamSet_is_msd(pss, i)) {      /* for MSD */
#ifdef OEC_SMOOTH
		  // count number of syllable & count number of smooth junctutes
		  N = 0;
		  K = 0;
		  j = 0;
		  last_end_idx = 0;
		  while( j < gss->total_frame ) {
            if (HTS_PStreamSet_get_msd_flag(pss, i, j)) {
			   start_idx = j;
			   syl_idx = pss->pstream[1].syllable_idx[j];
			   while ( j < gss->total_frame ){
				   if ( pss->pstream[1].syllable_idx[j] != syl_idx ){
					  break;
				   }
				   j ++;
			   }
			   end_idx = j - 1;
			   lf0_len = end_idx - start_idx + 1;
			   tmp_xlabel = HTS_XLabel_get_label(xlabel, syl_idx);
			   if( strcmp(tmp_xlabel->syllable, "sil")==0 ||
				   strcmp(tmp_xlabel->syllable, "sp")==0  ) {
			   }
			   else {
				   if( N > 0 && syl_idx == syl_map[N-1] ) {
					   if( lf0_len > pre_lf0_len ) {
						   pre_lf0_len = lf0_len;
					   }
					   T[N] = lf0_len;
					   //printf("re:%d %d\n", start_idx, end_idx);
					   last_end_idx = end_idx;
				   }
				   else {
					   syl_map[N] = syl_idx;
					   pre_lf0_len = lf0_len;
					   N++;
					   T[N] = lf0_len;
					   if( (start_idx/2) == ((last_end_idx)/2+1) && start_idx>0 ) {
						   K ++;
						   n_k[K] = N - 1;
						   tmp_xlabel = HTS_XLabel_get_label(xlabel, syl_map[N-1]);
						   printf("%d, %s\n",  syl_map[N-1], tmp_xlabel->syllable);
					   }
					   //printf("%d %d\n", start_idx, end_idx);
					   last_end_idx = end_idx;
				   }
			   }
			}
			j ++;
		  }
		  
		 printf("K=%d\n", K);
		 //system("pause");
		  if( K>0 ) {
			  // make OEC array
			  d = 4;// dimension for OEC
			  mtxalc(&a, N*d, 1);
			  for(l=0;l<N;l++) {
				  tmp_xlabel = HTS_XLabel_get_label(xlabel, syl_map[l]);
				  a[l*4]  [0] = tmp_xlabel->orth_coef[0];
				  a[l*4+1][0] = tmp_xlabel->orth_coef[1];
				  a[l*4+2][0] = tmp_xlabel->orth_coef[2];
				  a[l*4+3][0] = tmp_xlabel->orth_coef[3];
			  }
			  // make covariance matrix	
			  Rn = (double ***) HTS_calloc(N, sizeof(double **));
			  for(l=0;l<N;l++) {
				  mtxalc(Rn+l, d, d);
				  Rn[l][0][0] = R0;
				  Rn[l][1][1] = R1;
				  Rn[l][2][2] = R2;
				  Rn[l][3][3] = R3;
			  }
			  mtxalc(&x, N*d, 1);
			  // n_k: a syllable index array saving indices of smooth junctures
			  // n_k[k]: a syllable index of k-th smooth juncture, k = 1 ~ K
			  // T: an array to save length of voiced frame
			  // T[n]: the length of voiced frames for n-th syllable, n = 1 ~ N
			  // N: number of syllable
			  // K: number of smooth juncture
			  // a: matrix/vector for OEC
			  // Rn: covariance matrix for each syllable logF0
			  // d: dimension for OEC
			  // x: resultant smooth OEC
			  
			  
			  Smooth_OEC(a, Rn, T, N, d, n_k, K, x);
			  /*
			  //w[0] = 1.0; w[1] = 1.0; w[2] = 1.0; w[3] = 1.0;
			  for(l=0;l<4;l++)
				  w[l] = ((double) N);
			  mu[0] = 0.059657191514347; mu[1] = 0.006583952643522; mu[2] = 0.001440595136382; mu[3] = 0.000414419309673529;
			  //mu[0] = 0.189657191514347; mu[1] = 0.018583952643522; mu[2] = 0.013440595136382; mu[3] = 0.001914419309673529;
			  b[0] = 0.01;               b[1] = 0.001;              b[2] = 0.0002;             b[3] = 0.00007;
			  //b[0] = 0.001;               b[1] = 0.0001;              b[2] = 0.00002;             b[3] = 0.000007;
			  Smooth_Dyn_OEC(w, mu, b, a, Rn, T, N, d, n_k, K, x);
			  */
			  for(l=0;l<N;l++) {
				  tmp_xlabel = HTS_XLabel_get_label(xlabel, syl_map[l]);
				  printf("%s: %lf %lf %lf %lf\t%lf %lf %lf %lf\n", tmp_xlabel->syllable, tmp_xlabel->orth_coef[0], tmp_xlabel->orth_coef[1], tmp_xlabel->orth_coef[2], tmp_xlabel->orth_coef[3],
					  x[ l*4 + 0 ][0], x[ l*4 + 1 ][0], x[ l*4 + 2 ][0], x[ l*4 + 3 ][0]);
				  tmp_xlabel->orth_coef[0] = x[ l*4 + 0 ][0];
				  tmp_xlabel->orth_coef[1] = x[ l*4 + 1 ][0];
				  tmp_xlabel->orth_coef[2] = x[ l*4 + 2 ][0];
				  tmp_xlabel->orth_coef[3] = x[ l*4 + 3 ][0];
			  }
			  
			//  system("pause");
			  
			  
		  }
#endif
		  /* count total voiced frame of a syllable, added by Chen-Yu Chiang */


		  for (j = 0, msd_frame = 0; j < gss->total_frame; j++){
            if (HTS_PStreamSet_get_msd_flag(pss, i, j)) {
			   start_idx = j;
			   syl_idx = pss->pstream[1].syllable_idx[j];
			   for (; j < gss->total_frame; j++){
				   if ( pss->pstream[1].syllable_idx[j] != syl_idx ){
				      //end_idx = j - 1;
					  break;
				   }
			   }
			   j = j - 1;
			   end_idx=j; // modified by Liang, avoiding the sentence end is voiced part.
			   lf0_len = end_idx - start_idx + 1;
			   frame_lf0 = (double *) HTS_calloc((const size_t) lf0_len, sizeof(double));
               /* get frame lfo from orthogonal expansion coefficient */
               tmp_xlabel = HTS_XLabel_get_label(xlabel, syl_idx);
			   //Expansion_Pitch(lf0_len, frame_lf0, tmp_xlabel->orth_coef);

			   //because original 10ms,now 5 ms
			   xxx_frame_lf0 = (double *) HTS_calloc((const size_t) (int)(lf0_len / 2), sizeof(double));
			   Expansion_Pitch((int)(lf0_len/2), xxx_frame_lf0, tmp_xlabel->orth_coef);
			   for (l = 0; l < ((int)(lf0_len/2)); l ++ ){
				   frame_lf0[l*2+1] = (xxx_frame_lf0[l] - GMEAN) / GSTD * TSTD + TMEAN;
				   frame_lf0[l*2] = frame_lf0[l*2+1];
			   }
		       frame_lf0[lf0_len-1] = (xxx_frame_lf0[((int)(lf0_len/2))-1] - GMEAN) / GSTD * TSTD + TMEAN; // avoiding the length is odd.
			   HTS_free(xxx_frame_lf0);
			   for (l = 0; l < lf0_len; l ++){
			      for (k = 0; k < gss->gstream[i].static_length; k++)
                     gss->gstream[i].par[start_idx+l][k] = frame_lf0[l];
			   }
			   HTS_free(frame_lf0);
               msd_frame++;
            } else
               for (k = 0; k < gss->gstream[i].static_length; k++)
                  gss->gstream[i].par[j][k] = LZERO;
		  }

      } else {                  /* for non MSD */
         for (j = 0; j < gss->total_frame; j++)
            for (k = 0; k < gss->gstream[i].static_length; k++)
               gss->gstream[i].par[j][k] =
                   HTS_PStreamSet_get_parameter(pss, i, j, k);
      }
   }




   /* check */
   if (gss->nstream != 2 && gss->nstream != 3) {
      HTS_error(1, "HTS_GStreamSet_create: The number of streams should be 2 or 3.\n");
      HTS_GStreamSet_clear(gss);
      return FALSE;
   }
   if (HTS_PStreamSet_get_static_length(pss, 1) != 1) {
      HTS_error(1, "HTS_GStreamSet_create: The size of lf0 static vector should be 1.\n");
      HTS_GStreamSet_clear(gss);
      return FALSE;
   }
   if (gss->nstream >= 3 && gss->gstream[2].static_length % 2 == 0) {
      HTS_error(1, "HTS_GStreamSet_create: The number of low-pass filter coefficient should be odd numbers.");
      HTS_GStreamSet_clear(gss);
      return FALSE;
   }

   /* synthesize speech waveform */
   HTS_Vocoder_initialize(&v, gss->gstream[0].static_length - 1, stage, use_log_gain, sampling_rate, fperiod);
   if (gss->nstream >= 3)
      nlpf = (gss->gstream[2].static_length - 1) / 2;
    clock1 = clock();
   for (i = 0; i < gss->total_frame && (*stop) == FALSE; i++) {
      if (gss->nstream >= 3)
         lpf = &gss->gstream[2].par[i][0];

         HTS_Vocoder_synthesize(&v, gss->gstream[0].static_length - 1, gss->gstream[1].par[i][0], &gss->gstream[0].par[i][0], nlpf, lpf, alpha, beta, volume, &gss->gspeech[i * fperiod], audio);
   }
    clock2 = clock();
    LOGD("HTS_Vocoder_synthesize running time : %d", clock2-clock1);
   HTS_Vocoder_clear(&v);
   if (audio)
        HTS_Audio_flush(audio);

   return TRUE;
}

/* HTS_GStreamSet_get_total_nsample: get total number of sample */
int HTS_GStreamSet_get_total_nsample(HTS_GStreamSet * gss)
{
   return gss->total_nsample;
}

/* HTS_GStreamSet_get_total_frame: get total number of frame */
int HTS_GStreamSet_get_total_frame(HTS_GStreamSet * gss)
{
   return gss->total_frame;
}

/* HTS_GStreamSet_get_static_length: get static features length */
int HTS_GStreamSet_get_static_length(HTS_GStreamSet * gss, int stream_index)
{
   return gss->gstream[stream_index].static_length;
}

/* HTS_GStreamSet_get_speech: get synthesized speech parameter */
short HTS_GStreamSet_get_speech(HTS_GStreamSet * gss, int sample_index)
{
   return gss->gspeech[sample_index];
}

/* HTS_GStreamSet_get_parameter: get generated parameter */
double HTS_GStreamSet_get_parameter(HTS_GStreamSet * gss, int stream_index, int frame_index, int vector_index)
{
   return gss->gstream[stream_index].par[frame_index][vector_index];
}

/* HTS_GStreamSet_clear: free generated parameter stream set */
void HTS_GStreamSet_clear(HTS_GStreamSet * gss)
{
   int i, j;

   if (gss->gstream) {
      for (i = 0; i < gss->nstream; i++) {
         for (j = 0; j < gss->total_frame; j++)
            HTS_free(gss->gstream[i].par[j]);
         HTS_free(gss->gstream[i].par);
      }
      HTS_free(gss->gstream);
   }
   //if (gss->gspeech)
     // HTS_free(gss->gspeech); 20190110
   HTS_GStreamSet_initialize(gss);
}

HTS_GSTREAM_C_END;

#endif                          /* !HTS_GSTREAM_C */
