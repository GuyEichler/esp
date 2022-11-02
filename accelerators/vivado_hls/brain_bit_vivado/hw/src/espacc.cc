// Copyright (c) 2011-2022 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#include "../inc/espacc_config.h"
#include "../inc/espacc.h"
#include "hls_stream.h"
#include "hls_math.h"
#include <cstring>

void load(word_t _inbuff[SIZE_IN_CHUNK_DATA], dma_word_t *in1,
          /* <<--compute-params-->> */
	 const word_t avg,
	 const unsigned key_length,
	 const word_t std,
	 const word_t R,
	 const unsigned L,
	 const unsigned key_batch,
	  dma_info_t &load_ctrl, int chunk, int batch)
{
load_data:

    const unsigned length = round_up(key_length, VALUES_PER_WORD) / 1;
    const unsigned index = length * (batch * 1 + chunk);

    unsigned dma_length = length / VALUES_PER_WORD;
    unsigned dma_index = index / VALUES_PER_WORD;

    load_ctrl.index = dma_index;
    load_ctrl.length = dma_length;
    load_ctrl.size = SIZE_WORD_T;

    for (unsigned i = 0; i < dma_length; i++) {
    load_label0:for(unsigned j = 0; j < VALUES_PER_WORD; j++) {
	    _inbuff[i * VALUES_PER_WORD + j] = in1[dma_index + i].word[j];
    	}
    }
}

void store(word_t _outbuff[SIZE_OUT_CHUNK_DATA], dma_word_t *out,
          /* <<--compute-params-->> */
	 const word_t avg,
	 const unsigned key_length,
	 const word_t std,
	 const word_t R,
	 const unsigned L,
	 const unsigned key_batch,
	   dma_info_t &store_ctrl, int chunk, int batch)
{
store_data:

    const unsigned length = round_up(key_length, VALUES_PER_WORD) / 1;
    const unsigned store_offset = round_up(key_length, VALUES_PER_WORD) * key_batch;
    const unsigned out_offset = store_offset;
    const unsigned index = out_offset + length * (batch * 1 + chunk);

    unsigned dma_length = length / VALUES_PER_WORD;
    unsigned dma_index = index / VALUES_PER_WORD;

    store_ctrl.index = dma_index;
    store_ctrl.length = dma_length;
    store_ctrl.size = SIZE_WORD_T;

    for (unsigned i = 0; i < dma_length; i++) {
    store_label1:for(unsigned j = 0; j < VALUES_PER_WORD; j++) {
	    out[dma_index + i].word[j] = _outbuff[i * VALUES_PER_WORD + j];
	}
    }
}


void compute(word_t _inbuff[SIZE_IN_CHUNK_DATA],
             /* <<--compute-params-->> */
	 const word_t avg,
	 const unsigned key_length,
	 const word_t std,
	 const word_t R,
	 const unsigned L,
	 const unsigned key_batch,
             word_t _outbuff[SIZE_OUT_CHUNK_DATA])
{

    // TODO implement compute functionality
    const unsigned length = round_up(key_length, VALUES_PER_WORD) / 1;

    word_t Rs = (word_t)R * (word_t)std;
    unsigned result;

    for (int i = 0; i < length; i++){
        // _outbuff[i] = _inbuff[i];
        word_t val = _inbuff[i];
        bool filter;// = (fabs(val - avg) >= Rs);

#ifndef __SYNTHESIS__
        filter = (fabs((float)val - (float)avg) >= (float)Rs);
#endif

#ifdef __SYNTHESIS__
        filter = (fabs(val - avg) >= Rs);
#endif

        // if(val >= avg) filter = ((val - avg) >= Rs);
        // else filter = ((avg - val) >= Rs);

        if(!filter){
#ifndef __SYNTHESIS__
            result = floor((float)(((val - (avg - Rs)) / (2*Rs)) * L));
#endif
#ifdef __SYNTHESIS__
            result = floor(((val - (avg - Rs)) / (2*Rs)) * L);
#endif
            result = result % 2;
            _outbuff[i] = result;
        }
        else
            _outbuff[i] = 3;
    }
}


void top(dma_word_t *out, dma_word_t *in1,
         /* <<--params-->> */
	 const word_t conf_info_avg,
	 const unsigned conf_info_key_length,
	 const word_t conf_info_std,
	 const word_t conf_info_R,
	 const unsigned conf_info_L,
	 const unsigned conf_info_key_batch,
	 dma_info_t &load_ctrl, dma_info_t &store_ctrl)
{

    /* <<--local-params-->> */
	 const word_t avg = conf_info_avg;
	 const unsigned key_length = conf_info_key_length;
	 const word_t std = conf_info_std;
	 const word_t R = conf_info_R;
	 const unsigned L = conf_info_L;
	 const unsigned key_batch = conf_info_key_batch;

    // Batching
batching:
    for (unsigned b = 0; b < key_batch; b++)
    {
        // Chunking
    go:
        for (int c = 0; c < 1; c++)
        {
            word_t _inbuff[SIZE_IN_CHUNK_DATA];
            word_t _outbuff[SIZE_OUT_CHUNK_DATA];

            load(_inbuff, in1,
                 /* <<--args-->> */
	 	 avg,
	 	 key_length,
	 	 std,
	 	 R,
	 	 L,
	 	 key_batch,
                 load_ctrl, c, b);
            compute(_inbuff,
                    /* <<--args-->> */
	 	 avg,
	 	 key_length,
	 	 std,
	 	 R,
	 	 L,
	 	 key_batch,
                    _outbuff);
            store(_outbuff, out,
                  /* <<--args-->> */
	 	 avg,
	 	 key_length,
	 	 std,
	 	 R,
	 	 L,
	 	 key_batch,
                  store_ctrl, c, b);
        }
    }
}
