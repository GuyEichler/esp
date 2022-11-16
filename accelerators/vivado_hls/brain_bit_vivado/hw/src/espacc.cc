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
          const unsigned val_num,
          bool &load_values,
          bool &is_keys,
	  dma_info_t &load_ctrl, int batch)
{
load_data:

    const unsigned length = round_up(key_length, VALUES_PER_WORD);
    const unsigned index = length * batch;
    const unsigned val_length = round_up(val_num, VALUES_PER_WORD);

    unsigned dma_length = is_keys ? length / VALUES_PER_WORD : val_length / VALUES_PER_WORD;
    unsigned dma_index = index / VALUES_PER_WORD;

    if(load_values){

        load_ctrl.index = dma_index;
        load_ctrl.length = dma_length;
        load_ctrl.size = SIZE_WORD_T;

#ifndef __SYNTHESIS__

        std::cout << "LOAD : Loading new values " << std::endl;
        std::cout << "LOAD : dma_index "
                  << dma_index
                  << " dma_length "
                  << dma_length
                  << std::endl;

#endif

        for (unsigned i = 0; i < dma_length; i++) {

#pragma HLS loop_tripcount max=512

        load_label0:for(unsigned j = 0; j < VALUES_PER_WORD; j++) {
                _inbuff[i * VALUES_PER_WORD + j] = in1[dma_index + i].word[j];
            }
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
           const unsigned val_num,
           const unsigned tot_iter,
         bool &is_output_ready,
         unsigned &keys_done,
         bool &is_keys,
           dma_info_t &store_ctrl, //int batch,
         ap_uint<32> _outbuff_bit[SIZE_OUT_CHUNK_DATA])
{
store_data:

    unsigned length;
    unsigned val_length;

    // if(is_keys)
        length = round_up(key_length >> DATA_BITWIDTH_LOG, VALUES_PER_WORD);
    // else
        val_length = round_up(val_num, VALUES_PER_WORD);

    // const unsigned length = round_up(key_length >> DATA_BITWIDTH_LOG, VALUES_PER_WORD);
    const unsigned store_offset = round_up(key_length * key_batch * tot_iter, VALUES_PER_WORD);
    const unsigned out_offset = store_offset;

    const unsigned offset = keys_done * length;
    const unsigned index = out_offset + offset;

    unsigned dma_length = is_keys ? (length / VALUES_PER_WORD) : (val_length / VALUES_PER_WORD);
    unsigned dma_index = index / VALUES_PER_WORD;

    if(is_output_ready){

        store_ctrl.index = dma_index;
        store_ctrl.length = dma_length;
        store_ctrl.size = SIZE_WORD_T;

#ifndef __SYNTHESIS__
                std::cout << "STORE : dma_index "
                          << dma_index
                          << " dma_length "
                          << dma_length
                          << std::endl;
#endif

        for (unsigned i = 0; i < dma_length; i++) {

#pragma HLS loop_tripcount max=512

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
         // unsigned &add,
         ap_uint<32> _outbuff_bit[SIZE_OUT_CHUNK_DATA])
{

    const unsigned in_length = round_up(key_length, VALUES_PER_WORD);

    word_t Rs = R * std;
    unsigned result;
    ap_uint<32> result_b;
    static unsigned output_idx = 0;
    static unsigned input_offset = 0;
    unsigned i;

    // static ap_uint<32> bit_val_tot = 0;

#ifndef __SYNTHESIS__
    std::cout << "COMPUTE : Input offset is " << input_offset << std::endl;
#endif

COMPUTE_LOOP:for (i = 0 + input_offset; i < in_length; i++){

#pragma HLS loop_tripcount max=1024

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

            // _outbuff_bit[word][bit] = result_b;
            ap_uint<32> bit_val = _outbuff_bit[word];
            // ap_uint<32> bit_val = bit_val_tot;
            //bit_val[bit] = result_b;
            bit_val = bit_val >> 1;
            bit_val = bit_val | (result_b << (DATA_BITWIDTH - 1));
            _outbuff_bit[word] = bit_val;
            // if(bit == DATA_BITWIDTH - 1)
            //     _outbuff_bit[word] = bit_val;
            // else
            //     bit_val_tot = bit_val;

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
        // add = add + 1;

#ifndef __SYNTHESIS__
        std::cout << "COMPUTE : Load should be blocked " << std::endl;
#endif
// #ifndef __SYNTHESIS__
//         std::cout << "COMPUTE : Adding another batch. Total is " << key_batch + add << std::endl;
// #endif


    }

}

void compute_val(word_t _inbuff[SIZE_IN_CHUNK_DATA],
                 /* <<--compute-params-->> */
                 const word_t avg,
                 const unsigned key_length,
                 const word_t std,
                 const word_t R,
                 const unsigned L,
                 const unsigned key_batch,
                 const unsigned val_num,
                 ap_uint<32> _outbuff_bit[SIZE_OUT_CHUNK_DATA])
{

    // const unsigned in_length = round_up(key_length, VALUES_PER_WORD);

    // word_t Rs = R * std;
    // unsigned result;
    // ap_uint<32> result_b;
    // static unsigned output_idx = 0;
    // static unsigned input_offset = 0;
    // unsigned i;

ASSIGN_LOOP:for(unsigned i = 0; i < val_num; i++){
        word_t word = _inbuff[i];
    BIT_LOOP:for(unsigned j = 0; j < DATA_BITWIDTH; j++){
            bool bit = word[j];
            _outbuff_bit[i](j,j) = word[j];
        }

// #ifndef __SYNTHESIS__
//         std::cout << "COMPUTE VAL : Passed "
//                   << std::bitset<32>(_outbuff_bit[i])
//                   << " Received "
//                   << word.to_string()
//                   << std::endl;
// #endif

    }

}


void top(out_dma_word_t *out, dma_word_t *in1,
            /* <<--params-->> */
            const float conf_info_avg,
            const unsigned conf_info_key_length,
            const float conf_info_std,
            const float conf_info_R,
            const unsigned conf_info_L,
            const unsigned conf_info_key_batch,
            const unsigned conf_info_key_num,
            const unsigned conf_info_val_num,
            const unsigned conf_info_tot_iter,
            dma_info_t &load_ctrl, dma_info_t &store_ctrl)
{

        /* <<--local-params-->> */
        const word_t avg = conf_info_avg;
        const unsigned key_length = conf_info_key_length;
        // const unsigned val_length = conf_info_key_length << DATA_BITWIDTH_LOG;
        const word_t std = conf_info_std;
        const word_t R = conf_info_R;
        const unsigned L = conf_info_L;
        const unsigned key_batch = conf_info_key_batch;
        const unsigned key_num = conf_info_key_num;
        const unsigned val_num = conf_info_val_num;
        const unsigned tot_iter = conf_info_tot_iter;


        //Memories
        // static word_t _inbuff[SIZE_IN_CHUNK_DATA];
        static word_t _inbuff[SIZE_IN_CHUNK_DATA];
        // static ap_uint<32> _outbuff_bit[SIZE_OUT_BIT_DATA];
        static ap_uint<32> _outbuff_bit[SIZE_OUT_CHUNK_DATA];


        for(unsigned i = 0; i < tot_iter; i++){

            bool is_output_ready = false;
            bool load_values = true;

            // unsigned add = 0;
            unsigned keys_done = 0;
            unsigned values_done = 0;

            unsigned b = 0;

            bool is_keys = true;

            // Keys loop
        keys:
            // for (unsigned b = 0; b < key_batch + add; b++)
            while(keys_done != key_num)
            {

//#pragma HLS loop_tripcount max=1

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
                     is_keys,
                     load_ctrl, i * key_batch + b);

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
                        // add,
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
                      val_num,
                      tot_iter,
                      is_output_ready,
                      keys_done,
                      is_keys,
                      store_ctrl, //b,
                      _outbuff_bit);

                if(load_values)
                    b = b + 1;

#ifndef __SYNTHESIS__
                if(keys_done == key_num){
                    std::cout << "TOP : Enough keys were generated " << std::endl;
                }
#endif

            }

            values_done = keys_done;//start from keys_done for offsetting memory
            is_keys = false;

            //values loop
        values:
            // for (unsigned b_v = 0; b_v < val_num; b_v++)
            // while(values_done != val_num)
            if(val_num != 0)
            {

// #pragma HLS loop_tripcount max=1024

    // go_2:

                bool load_values_raw = true;
                bool is_output_raw = true;

                load(_inbuff, in1,
                     /* <<--args-->> */
                     avg,
                     key_length,
                     std,
                     R,
                     L,
                     key_batch,
                     val_num,
                     load_values_raw,
                     is_keys,
                     load_ctrl, i * key_batch + b);

                compute_val(_inbuff,
                            /* <<--args-->> */
                            avg,
                            key_length,
                            std,
                            R,
                            L,
                            key_batch,
                            val_num,
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
                      val_num,
                      tot_iter,
                      is_output_raw,
                      // values_done,
                      keys_done,
                      is_keys,
                      store_ctrl, //b,
                      _outbuff_bit);

                //b = b + 1;

#ifndef __SYNTHESIS__
        // if(values_done == val_num+key_num){
        //     std::cout << "TOP : Enough values were generated " << std::endl;
        // }
            std::cout << "TOP : Values were generated " << std::endl;
#endif

            }
        }
}
