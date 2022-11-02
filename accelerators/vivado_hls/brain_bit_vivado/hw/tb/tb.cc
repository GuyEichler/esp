// Copyright (c) 2011-2022 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#include "../inc/espacc_config.h"
#include "../inc/espacc.h"

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

int main(int argc, char **argv) {

    printf("****start*****\n");

    /* <<--params-->> */
	 const word_t avg = 2.616848987988529;
	 const unsigned key_length = 128;
	 const word_t std = 47.874003721012755;
	 const word_t R = 1.5;
	 const unsigned L = 1500;
	 const unsigned key_batch = 1;

    uint32_t in_words_adj;
    uint32_t out_words_adj;
    uint32_t in_size;
    uint32_t out_size;
    uint32_t dma_in_size;
    uint32_t dma_out_size;
    uint32_t dma_size;


    in_words_adj = round_up(key_length, VALUES_PER_WORD);
    out_words_adj = round_up(key_length, VALUES_PER_WORD);
    in_size = in_words_adj * (key_batch);
    out_size = out_words_adj * (key_batch);

    dma_in_size = in_size / VALUES_PER_WORD;
    dma_out_size = out_size / VALUES_PER_WORD;
    dma_size = dma_in_size + dma_out_size;

    dma_word_t *mem=(dma_word_t*) malloc(dma_size * sizeof(dma_word_t));
    word_t *inbuff=(word_t*) malloc(in_size * sizeof(word_t));
    word_t *outbuff=(word_t*) malloc(out_size * sizeof(word_t));
    word_t *outbuff_gold= (word_t*) malloc(out_size * sizeof(word_t));
    dma_info_t load;
    dma_info_t store;

    word_t input_arr[10] = { 7.976249119563723,
                             -1.8357592988494529,
                             -11.580534496426118,
                             1.9915349907701394,
                             8.664183581327363,
                             2.4353615908873127,
                             -10.50784361660308,
                             -6.624226414376322,
                             -2.2397914673913797,
                             -1.2090025526478039 };

    word_t output_arr[10] = { 1, 1, 1, 1, 1, 0, 0, 1, 1, 0 };


    // Prepare input data
    for(unsigned i = 0; i < key_batch; i++)
        for(unsigned j = 0; j < key_length; j++)
            if(i * in_words_adj + j < 10){
                inbuff[i * in_words_adj + j] = (word_t) input_arr[i * in_words_adj + j];
            }
            else
                inbuff[i * in_words_adj + j] = (word_t) j;

    for(unsigned i = 0; i < dma_in_size; i++)
	for(unsigned k = 0; k < VALUES_PER_WORD; k++)
	    mem[i].word[k] = inbuff[i * VALUES_PER_WORD + k];


    word_t Rs = R * std;
    unsigned result;

    // Set golden output
    for(unsigned i = 0; i < key_batch; i++)
        for(unsigned j = 0; j < key_length; j++)
        {
            word_t val = inbuff[i * in_words_adj + j];
            bool filter = (fabs((float)val - (float)avg) >= (float)Rs);
            if(!filter){
                result = floor((float)(((val - (avg - Rs)) / (2*Rs)) * L));
                result = result % 2;
                outbuff_gold[i * out_words_adj + j] = (word_t) result;
            }
            else
                outbuff_gold[i * out_words_adj + j] = 3;
        }


    // Call the TOP function
    top(mem, mem,
        /* <<--args-->> */
	 	 avg,
	 	 key_length,
	 	 std,
                 R,
	 	 L,
	 	 key_batch,
        load, store);

    // Validate
    uint32_t out_offset = dma_in_size;
    for(unsigned i = 0; i < dma_out_size; i++)
	for(unsigned k = 0; k < VALUES_PER_WORD; k++)
	    outbuff[i * VALUES_PER_WORD + k] = mem[out_offset + i].word[k];

    int errors = 0;
    for(unsigned i = 0; i < key_batch; i++)
        for(unsigned j = 0; j < key_length; j++){
            word_t val = outbuff[i * out_words_adj + j];
            std::cout << "Calculated value " << val << " Golden value " << outbuff_gold[i * out_words_adj + j] << " for index " << i * out_words_adj + j << std::endl;
	    if (outbuff[i * out_words_adj + j] != outbuff_gold[i * out_words_adj + j])
		errors++;
        }

    if (errors)
	std::cout << "Test FAILED with " << errors << " errors." << std::endl;
    else
	std::cout << "Test PASSED." << std::endl;

    // Free memory

    free(mem);
    free(inbuff);
    free(outbuff);
    free(outbuff_gold);

    return 0;
}
