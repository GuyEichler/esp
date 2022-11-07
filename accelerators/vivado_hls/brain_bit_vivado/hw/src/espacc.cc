// Copyright (c) 2011-2022 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#include "../inc/espacc_config.h"
#include "../inc/espacc.h"
#include "hls_stream.h"
#include "hls_math.h"
#include <cstring>

void load(word_t _inbuff[SIZE_IN_CHUNK_DATA], dma_word_t *in1,
          /* <<--load-params-->> */
	 const word_t avg,
	 const unsigned key_length,
	 const word_t std,
	 const word_t R,
	 const unsigned L,
	 const unsigned key_batch,
          bool &is_load_finished,
          bool &load_values,
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

    if(batch == key_batch - 1){
        is_load_finished = true;

#ifndef __SYNTHESIS__
        std::cout << "LOAD : About to load last batch " << std::endl;
#endif

    }

    if(load_values){

#ifndef __SYNTHESIS__
        std::cout << "LOAD : Loading new values " << std::endl;
#endif

        for (unsigned i = 0; i < dma_length; i++) {
        load_label0:for(unsigned j = 0; j < VALUES_PER_WORD; j++) {
                _inbuff[i * VALUES_PER_WORD + j] = in1[dma_index + i].word[j];
            }
        }
    }
}

void store(word_t _outbuff[SIZE_OUT_CHUNK_DATA], dma_word_t *out,
          /* <<--store-params-->> */
	 const word_t avg,
	 const unsigned key_length,
	 const word_t std,
	 const word_t R,
	 const unsigned L,
	 const unsigned key_batch,
           bool &is_output_ready,
           bool &is_output_sent,
	   dma_info_t &store_ctrl, int chunk, int batch)
{
store_data:

    const unsigned length = round_up(key_length, VALUES_PER_WORD) / 1;
    const unsigned store_offset = round_up(key_length, VALUES_PER_WORD) * key_batch;
    const unsigned out_offset = store_offset;
    // const unsigned index = out_offset + length * (batch * 1 + chunk);
    static unsigned index = out_offset;

    unsigned dma_length = length / VALUES_PER_WORD;
    unsigned dma_index = index / VALUES_PER_WORD;

    store_ctrl.index = dma_index;
    store_ctrl.length = dma_length;
    store_ctrl.size = SIZE_WORD_T;

    //is_output_sent = false;

    if(is_output_ready){
        for (unsigned i = 0; i < dma_length; i++) {
        store_label1:for(unsigned j = 0; j < VALUES_PER_WORD; j++) {
                out[dma_index + i].word[j] = _outbuff[i * VALUES_PER_WORD + j];
            }
        }
        is_output_sent = true;
        index += length;

#ifndef __SYNTHESIS__
        std::cout << "STORE : Output was sent " << std::endl;
#endif

    }
    else
        is_output_sent = false;
}

void compute(word_t _inbuff[SIZE_IN_CHUNK_DATA],
             /* <<--compute-params-->> */
	 const word_t avg,
	 const unsigned key_length,
	 const word_t std,
	 const word_t R,
	 const unsigned L,
	 const unsigned key_batch,
         bool &is_output_ready,
         // bool &is_output_sent,
         // unsigned& output_idx,
         // bool &is_load_finished,
         bool &load_values,
         unsigned &add,
         word_t _outbuff[SIZE_OUT_CHUNK_DATA])
{

    // TODO implement compute functionality
    const unsigned in_length = round_up(key_length, VALUES_PER_WORD) / 1;

    word_t Rs = R * std;
    unsigned result;
    static unsigned output_idx = 0;
    static unsigned input_offset = 0;
    unsigned i;

    is_output_ready = false;
    //load_values = true;

#ifndef __SYNTHESIS__
    std::cout << "COMPUTE : Input offset is " << input_offset << std::endl;
#endif

    for (i = 0 + input_offset; i < in_length; i++){
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

        //while(is_output_ready && !is_output_sent){/*WAIT*/};

        if(!filter){
#ifndef __SYNTHESIS__
            result = floor((float)(((val - (avg - Rs)) / (2*Rs)) * L));
            // std::cout << "Output idx " << output_idx << " input idx " << i << " result " << result % 2 << std::endl;
#endif
#ifdef __SYNTHESIS__
            result = floor(((val - (avg - Rs)) / (2*Rs)) * L);
#endif
            result = result % 2;
            _outbuff[output_idx] = result;
            output_idx++;
        }
        // else{
        //     _outbuff[output_idx] = 3;
        //     output_idx++;
        // }

        if(output_idx == in_length){

#ifndef __SYNTHESIS__
            std::cout << "COMPUTE : Output idx equals in_length " << std::endl;
#endif

            is_output_ready = true;
            output_idx = 0;
        }
//         else if(is_load_finished && i == in_length - 1){
//             // is_output_ready = true;

// #ifndef __SYNTHESIS__
//             std::cout << "COMPUTE : Load was done and signaled to compute " << std::endl;
// #endif

//         }
        else{
            is_output_ready = false;
        }

        if(is_output_ready){
            if(i != in_length - 1){
                input_offset = i+1;

#ifndef __SYNTHESIS__
                std::cout << "COMPUTE : Input offset is set to " << input_offset << std::endl;
#endif

                break;
            }
            else{
                input_offset = 0;

#ifndef __SYNTHESIS__
                std::cout << "COMPUTE : Input offset is initialized " << std::endl;
#endif

            }

        }
    }

    if(i == in_length){
        load_values = true;
        input_offset = 0;

#ifndef __SYNTHESIS__
        std::cout << "COMPUTE : Load can load new values. Input offset initialized. Output index is " << output_idx << std::endl;
#endif

    }
    else{
        load_values = false;
        add = add + 1;

#ifndef __SYNTHESIS__
        std::cout << "COMPUTE : Load should be blocked " << std::endl;
#endif
#ifndef __SYNTHESIS__
        std::cout << "COMPUTE : Adding another batch. Total is " << key_batch + add << std::endl;
#endif


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
	 const unsigned conf_info_key_num,
	 dma_info_t &load_ctrl, dma_info_t &store_ctrl)
{

    /* <<--local-params-->> */
	 const word_t avg = conf_info_avg;
	 const unsigned key_length = conf_info_key_length;
	 const word_t std = conf_info_std;
	 const word_t R = conf_info_R;
	 const unsigned L = conf_info_L;
	 const unsigned key_batch = conf_info_key_batch;
	 const unsigned key_num = conf_info_key_num;

         static bool is_output_ready = false;
         static bool is_output_sent = false;
         static bool is_load_finished = false;
         static bool load_values = true;
         static unsigned add = 0;
         // static bool first = true;
         static unsigned keys_done = 0;
         // unsigned output_idx = 0;

    // Batching
batching:
    for (unsigned b = 0; b < key_batch + add; b++)
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
                is_load_finished,
                load_values,
                 load_ctrl, c, b-add);

            compute(_inbuff,
                    /* <<--args-->> */
	 	 avg,
	 	 key_length,
	 	 std,
	 	 R,
	 	 L,
	 	 key_batch,
                 is_output_ready,
                 // is_output_sent,
                 // is_load_finished,
                 // output_idx,
                 load_values,
                 add,
                _outbuff);

            store(_outbuff, out,
                  /* <<--args-->> */
	 	 avg,
	 	 key_length,
	 	 std,
	 	 R,
	 	 L,
	 	 key_batch,
                is_output_ready,
                is_output_sent,
                  store_ctrl, c, b);

//             if(!load_values){
//                 add = add + 1;

// #ifndef __SYNTHESIS__
//                std::cout << "TOP : Adding another batch. Total is " << key_batch + add << std::endl;
//                    // " b-add is "<< b-add << std::endl;
// #endif

//             }

            if(is_output_sent){
                keys_done = keys_done + 1;

#ifndef __SYNTHESIS__
                std::cout << "TOP : Keys generated " << keys_done << std::endl;
#endif

                if(keys_done == key_num){
                    b = key_batch + add - 1;

#ifndef __SYNTHESIS__
               std::cout << "TOP : Enough keys were generated " << std::endl;
#endif

                }
            }
        }
    }
}
