// Copyright (c) 2011-2022 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#include "../inc/espacc_config.h"
#include "../inc/espacc.h"
#include "hls_stream.h"
#include "hls_math.h"
#include <cstring>
#include <bitset>

void load(word_t _inbuff[SIZE_IN_CHUNK_DATA], dma_word_t *in1,
          /* <<--load-params-->> */
	 const word_t avg,
	 const unsigned key_length,
	 const word_t std,
	 const word_t R,
	 const unsigned L,
	 const unsigned key_batch,
          const unsigned key_num,
          bool &load_values,
	  dma_info_t &load_ctrl, int batch)
{
load_data:

    const unsigned length = round_up(key_length, VALUES_PER_WORD);
    const unsigned index = length * batch;

    unsigned dma_length = load_values ? length / VALUES_PER_WORD : 0;
    unsigned dma_index = index / VALUES_PER_WORD;

    load_ctrl.index = dma_index;
    load_ctrl.length = dma_length;
    load_ctrl.size = SIZE_WORD_T;

#ifndef __SYNTHESIS__
    if(!load_values){
        std::cout << "LOAD : Loading new values " << std::endl;
    }
#endif

    for (unsigned i = 0; i < dma_length; i++) {

#pragma HLS loop_tripcount max=128

    load_label0:for(unsigned j = 0; j < VALUES_PER_WORD; j++) {
            _inbuff[i * VALUES_PER_WORD + j] = in1[dma_index + i].word[j];
        }
    }
}

void store(out_dma_word_t *out,
          /* <<--store-params-->> */
	 const word_t avg,
	 const unsigned key_length,
	 const word_t std,
	 const word_t R,
	 const unsigned L,
	 const unsigned key_batch,
           const unsigned key_num,
         bool &is_output_ready,
         unsigned &keys_done,
	 dma_info_t &store_ctrl, int batch,
         ap_uint<32> _outbuff_bit[SIZE_OUT_BIT_DATA])
{
store_data:

    const unsigned length = round_up(key_length >> DATA_BITWIDTH_LOG, VALUES_PER_WORD);
    const unsigned store_offset = round_up(key_length, VALUES_PER_WORD) * key_batch;
    const unsigned out_offset = store_offset;

    const unsigned offset = keys_done * length;
    const unsigned index = out_offset + offset;

    unsigned dma_length = is_output_ready ? length / VALUES_PER_WORD : 0;
    unsigned dma_index = index / VALUES_PER_WORD;

    store_ctrl.index = dma_index;
    store_ctrl.length = dma_length;
    store_ctrl.size = SIZE_WORD_T;

    for (unsigned i = 0; i < dma_length; i++) {

#pragma HLS loop_tripcount max=4

    store_label1:for(unsigned j = 0; j < VALUES_PER_WORD; j++) {
            out[dma_index + i].word[j] = _outbuff_bit[i * VALUES_PER_WORD + j];

#ifndef __SYNTHESIS__
            std::cout << "STORE : Sent "
                      << std::bitset<32>(_outbuff_bit[i * VALUES_PER_WORD + j])
                      << " in memory "
                      << std::bitset<32>(out[dma_index + i].word[j]) << std::endl;
#endif

        }
    }

    if(is_output_ready){

#ifndef __SYNTHESIS__
        std::cout << "STORE : Output was sent " << std::endl;
#endif

        keys_done = keys_done + 1;

#ifndef __SYNTHESIS__
            std::cout << "STORE : Keys generated " << keys_done << std::endl;
#endif

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
            const unsigned key_num,
         bool &is_output_ready,
         bool &load_values,
         unsigned &add,
         ap_uint<32> _outbuff_bit[SIZE_OUT_BIT_DATA])
{

    const unsigned in_length = round_up(key_length, VALUES_PER_WORD);

    word_t Rs = R * std;
    unsigned result;
    ap_uint<1> result_b;
    static unsigned output_idx = 0;
    static unsigned input_offset = 0;
    unsigned i;

#ifndef __SYNTHESIS__
    std::cout << "COMPUTE : Input offset is " << input_offset << std::endl;
#endif

    for (i = 0 + input_offset; i < in_length; i++){

#pragma HLS loop_tripcount max=128

        word_t val = _inbuff[i];
        bool filter;// = (fabs(val - avg) >= Rs);

#ifndef __SYNTHESIS__
        filter = (fabs((float)val - (float)avg) >= (float)Rs);
#endif

#ifdef __SYNTHESIS__
        filter = (fabs(val - avg) >= Rs);
#endif

        if(!filter){
#ifndef __SYNTHESIS__
            result = floor((float)(((val - (avg - Rs)) / (2*Rs)) * L));
            // std::cout << "Output idx " << output_idx << " input idx " << i << " result " << result % 2 << std::endl;
#endif
#ifdef __SYNTHESIS__
            result = floor(((val - (avg - Rs)) / (2*Rs)) * L);
#endif

            result_b = result % 2;
            unsigned bit = output_idx % DATA_BITWIDTH;
            unsigned word = output_idx >> DATA_BITWIDTH_LOG;

            _outbuff_bit[word][bit] = result_b;

#ifndef __SYNTHESIS__
            // std::cout << "COMPUTE : bit equals " << result_b
            //           << " word is " << word
            //           << " bit idx is " << bit
            //           << " outbuff_bit is " << std::bitset<32>(_outbuff_bit[word]) << std::endl;
#endif

            output_idx++;
        }

        if(output_idx == in_length){

#ifndef __SYNTHESIS__
            std::cout << "COMPUTE : Output idx equals in_length " << std::endl;
#endif

            is_output_ready = true;
            output_idx = 0;
        }
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


void top(out_dma_word_t *out, dma_word_t *in1,
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

        bool is_output_ready = false;
        bool load_values = true;
        unsigned add = 0;
        unsigned keys_done = 0;

        //Memories
        static word_t _inbuff[SIZE_IN_CHUNK_DATA];
        static ap_uint<32> _outbuff_bit[SIZE_OUT_BIT_DATA];

        // Batching
    batching:
    for (unsigned b = 0; b < key_batch + add; b++)
    {

#pragma HLS loop_tripcount max=1

    go:
        load(_inbuff, in1,
            /* <<--args-->> */
            avg,
            key_length,
            std,
            R,
            L,
            key_batch,
            key_num,
            load_values,
            load_ctrl, b-add);

        compute(_inbuff,
            /* <<--args-->> */
            avg,
            key_length,
            std,
            R,
            L,
            key_batch,
            key_num,
            is_output_ready,
            load_values,
            add,
            _outbuff_bit);

            store(out,
            /* <<--args-->> */
            avg,
            key_length,
            std,
            R,
            L,
            key_batch,
            key_num,
            is_output_ready,
            keys_done,
            store_ctrl, b,
            _outbuff_bit);


        if(keys_done == key_num){

#ifndef __SYNTHESIS__
            std::cout << "TOP : Enough keys were generated " << std::endl;
#endif
            break;

        }
    }
}
