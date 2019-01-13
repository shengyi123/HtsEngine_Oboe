/* ----------------------------------------------------------------- */
/*           The HMM-Based Speech Synthesis System (HTS)             */
/*           developed by HTS Working Group                          */
/*           http://hts.sp.nitech.ac.jp/                             */
/* ----------------------------------------------------------------- */
/*                                                                   */
/*  Copyright (c) 2001-2008  Nagoya Institute of Technology          */
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

/* $id: HTS_label.c,v 1.1.2.1 2008/02/12 10:47:43 uratec Exp $ */

#include <stdlib.h>             /* for atoi() */
#include <ctype.h>              /* for isgraph(),isdigit() */
#include <string.h>

/* hts_engine libraries */
#include "HTS_hidden.h"
#include "LOG.h"

/* HTS_Label_initialize: initialize label */
void HTS_XLabel_initialize(HTS_XLabel * label)
{
   label->head = NULL;
   label->size = 0;
   label->hmm_size = 0;
   //label->version_str[0] = NULL;  /*20181027  clearec by shengyi*/
}

/* HTS_Label_load_from_fn: load label from file name */
void HTS_XLabel_load_from_fn(HTS_XLabel * label, int sampling_rate, int fperiod, char *fn)
{
   FILE *fp = HTS_fopen(fn, "r");
   HTS_XLabel_load_from_fp(label, sampling_rate, fperiod, fp);
   fclose(fp);
}

/* HTS_Label_load_from_fp: load label from file pointer */
void HTS_XLabel_load_from_fp(HTS_XLabel * label, int sampling_rate, int fperiod, FILE * fp)
{
   char buff[HTS_MAXBUFLEN];
   char syllable[HTS_MAXBUFLEN];
   int hmm_model_num = 0;
   char model_name[16][1024];
   int start = 0, end = 0;
   double dur = 0;
   double a[4];
   double el;
   int i;
   int zxc;
   HTS_XLabelString *lstring = NULL;
   const double rate = sampling_rate / (fperiod * 1e+7);


   if (label->head || label->size != 0)
      HTS_error(1, "XHTS_Label_load_from_fp: Xlabel is not initialized.\n");
      HTS_get_token(fp, buff);
       zxc = strcmp(buff, "ver1.0");
       LOGD("strcmp(buff, \"ver1.0\") = %d",zxc);
     if( strcmp(buff, "ver1.0")!=0) {/*20181116 add '!'*/
       // check header
	   /* parse label file (duration in sec) */
	   while (HTS_get_token(fp, buff)) {   

		  label->size++;
		  strcpy(syllable, buff);
		  HTS_get_token(fp, buff);
		  hmm_model_num = atoi(buff);
		  for (i=0; i<hmm_model_num; i++){
			  HTS_get_token(fp, buff);
			  strcpy(model_name[i], buff);
			  label->hmm_size ++;
		  }
		  HTS_get_token(fp, buff);
		  dur = atof(buff);
		  end = start + ((int) (dur * 10000000.0));
		  //start = atoi(buff);
		  //HTS_get_token(fp, buff);
		  //end = atoi(buff);
		  for (i=0; i<4; i++)
		  {
			HTS_get_token(fp, buff);
			a[i] = atof(buff);
		  }
		  HTS_get_token(fp, buff);
		  el = atof(buff);
      
		  if (lstring) {
			 lstring->next =
				 (HTS_XLabelString *) HTS_calloc(1, sizeof(HTS_XLabelString));
			 lstring = lstring->next;
		  } else {                  /* first time */
			 lstring = (HTS_XLabelString *) HTS_calloc(1, sizeof(HTS_XLabelString));
			 label->head = lstring;
		  }  
		  lstring->syllable = (char *) HTS_calloc(strlen(syllable)+1, sizeof(char));
		  strcpy(lstring->syllable, syllable);
		  lstring->hmm_model_num = hmm_model_num;
		  lstring->hmm_model_name = (char **) HTS_calloc(hmm_model_num, sizeof(char *));
		  for (i=0; i<hmm_model_num; i++){
			 lstring->hmm_model_name[i] = (char *) HTS_calloc(strlen(model_name[i])+1, sizeof(char));//modified by Liang, original size is wrong  
			 strcpy(lstring->hmm_model_name[i], model_name[i]);
		  }
		  lstring->syllable_start = ((double)start) * rate;
		  lstring->syllable_end = ((double)end) * rate;
		  for (i=0; i<4; i++){	  
			lstring->orth_coef[i] = a[i];
		  }
		  lstring->energy_level = el;
		  start = end;
	   }
   }
   else {
	   rewind(fp);
	   /* parse label file */
	   while (HTS_get_token(fp, buff)) {   
	   
		  label->size++;
		  strcpy(syllable, buff);
		  HTS_get_token(fp, buff);
		  hmm_model_num = atoi(buff);
		  for (i=0; i<hmm_model_num; i++){
			  HTS_get_token(fp, buff);
			  strcpy(model_name[i], buff);
			  label->hmm_size ++;
		  }
		  HTS_get_token(fp, buff);
		  start = atoi(buff);
		  HTS_get_token(fp, buff);
		  end = atoi(buff);
		  for (i=0; i<4; i++)
		  {
			HTS_get_token(fp, buff);
			a[i] = atof(buff);
		  }
		  HTS_get_token(fp, buff);
		  el = atof(buff);
      
		  if (lstring) {
			 lstring->next =
				 (HTS_XLabelString *) HTS_calloc(1, sizeof(HTS_XLabelString));
			 lstring = lstring->next;
		  } else {                  /* first time */
			 lstring = (HTS_XLabelString *) HTS_calloc(1, sizeof(HTS_XLabelString));
			 label->head = lstring;
		  }  
		  lstring->syllable = (char *) HTS_calloc(strlen(syllable)+1, sizeof(char));
		  strcpy(lstring->syllable, syllable);
		  lstring->hmm_model_num = hmm_model_num;
		  lstring->hmm_model_name = (char **) HTS_calloc(hmm_model_num, sizeof(char *));
		  for (i=0; i<hmm_model_num; i++){
			 lstring->hmm_model_name[i] = (char *) HTS_calloc(strlen(model_name[i])+1, sizeof(char));//modified by Liang, original size is wrong  
			 strcpy(lstring->hmm_model_name[i], model_name[i]);
		  }
		  lstring->syllable_start = ((double)start) * rate;
		  lstring->syllable_end = ((double)end) * rate;
		  for (i=0; i<4; i++){	  
			lstring->orth_coef[i] = a[i];
		  }
		  lstring->energy_level = el;
	   }
   }
}



/* HTS_Label_get_size: get number of label string */
int HTS_XLabel_get_size(HTS_XLabel * label)
{
   return label->size;
}

/* HTS_XLabel_get_hmm_size: get number of hmm label string */
int HTS_XLabel_get_hmm_size(HTS_XLabel * label)
{
   return label->hmm_size;
}

/* HTS_Label_get_label: get label */
HTS_XLabelString *HTS_XLabel_get_label(HTS_XLabel * label, int string_index)
{
   HTS_XLabelString *lstring = label->head;

   while (string_index-- && lstring)
      lstring = lstring->next;
   if (!lstring)
      return NULL;
   return lstring;
}

/* XHTS_Label_clear: free label */
void HTS_XLabel_clear(HTS_XLabel * label)
{
   HTS_XLabelString *lstring, *next_lstring;
   int i;

   for (lstring = label->head; lstring; lstring = next_lstring) {
      next_lstring = lstring->next;
	  HTS_free(lstring->syllable);
	  for(i=0;i<lstring->hmm_model_num;i++)
	  {
		  HTS_free(lstring->hmm_model_name[i]);
	  }
	  HTS_free(lstring->hmm_model_name);
      HTS_free(lstring);
   }
   HTS_XLabel_initialize(label);
}
