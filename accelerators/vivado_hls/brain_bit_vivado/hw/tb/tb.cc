// Copyright (c) 2011-2022 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0
#include "../inc/espacc_config.h"
#include "../inc/espacc.h"

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <string>
#include <bitset>

int main(int argc, char **argv) {

    printf("****start*****\n");

    /* <<--params-->> */
    const word_t avg = 3.0677295382679177;
    const float avg_f = 3.0677295382679177;
    const unsigned key_length = 128;
    const word_t std = 38.626628825256695;
    const float std_f = 38.626628825256695;
    const word_t R = 1.5;
    const float R_f = 1.5;
    const unsigned L = 1500;
    const unsigned key_batch = 20;
    const unsigned key_num = 15;
    const unsigned val_num = 2;

    uint32_t in_words_adj;
    uint32_t out_words_adj;
    uint32_t in_size;
    uint32_t out_size;
    uint32_t dma_in_size;
    uint32_t dma_out_size;
    uint32_t dma_size;

    // FILE* fp = fopen("/home/geichler/Desktop/esp_accelerators/esp_guy/esp/accelerators/vivado_hls/brain_bit_vivado/hw/tb/raw_values.txt", "r");
    // FILE* fp = fopen("raw_values.txt", "r");
    const int file_length = key_length*(key_batch+1);//78736896;
    // std::ifstream inFile("raw_values.txt");
    float val_arr[file_length];

    std::string inFileName = "/home/geichler/Desktop/esp_accelerators/esp_guy/esp/accelerators/vivado_hls/brain_bit_vivado/hw/tb/raw_values.txt";
    std::ifstream inFile;
    inFile.open(inFileName.c_str());

    // printf("checking file.");

    if (inFile.is_open()){
        for (int i = 0; i < file_length; i++) {
            inFile >> val_arr[i];
            // std::cout.precision(19);
            // if(i == 0)
            //     std::cout << "[ ";
            // std::cout << val_arr[i] << ", ";
            // if(i == file_length - 1)
            //     std::cout << "  ]";
        }
        // std::cout << std::endl;

        inFile.close();
    }
    else { //Error message
        std::cerr << "Can't find input file " << inFileName << std::endl;
    }


    // if (fp == NULL)
    // {
    //     printf("Error! Could not open file\n");
    //     perror("fopen");
    //     exit(-1); // must include stdlib.h
    // }
    // else
    //     printf("file is not null.");


    // for(int i = 0; i < file_length; i++){
    //     fscanf(fp, "%f ", val_arr[i]);
    //     printf("%f ", val_arr[i]);
    // }

    // fclose(fp);

    in_words_adj = round_up(key_length, VALUES_PER_WORD);
    out_words_adj = round_up(key_length, VALUES_PER_WORD);
    in_size = in_words_adj * (key_batch);
    out_size = out_words_adj * (key_batch);

    dma_in_size = in_size / VALUES_PER_WORD;
    dma_out_size = out_size / VALUES_PER_WORD;
    dma_size = dma_in_size + dma_out_size;

    dma_word_t *mem=(dma_word_t*) malloc(dma_size * sizeof(dma_word_t));
    out_dma_word_t *mem_out=(out_dma_word_t*) malloc(dma_size * sizeof(out_dma_word_t));
    word_t *inbuff=(word_t*) malloc(in_size * sizeof(word_t));
    word_t *outbuff=(word_t*) malloc(out_size * sizeof(word_t));
    ap_uint<32> *outbuff_bit=(ap_uint<32>*) malloc(out_size * sizeof(word_t));
    word_t *outbuff_gold= (word_t*) malloc(out_size * sizeof(word_t));
    dma_info_t load;
    dma_info_t store;

    // word_t input_arr[10] = { 7.976249119563723,
    //                          -1.8357592988494529,
    //                          -11.580534496426118,
    //                          1.9915349907701394,
    //                          8.664183581327363,
    //                          2.4353615908873127,
    //                          -10.50784361660308,
    //                          -6.624226414376322,
    //                          -2.2397914673913797,
    //                          -1.2090025526478039 };

    // word_t output_arr[10] = { 1, 1, 1, 1, 1, 0, 0, 1, 1, 0 };

    word_t output_arr[10] = {1, 1, 1, 0, 1, 0, 0, 0, 0, 0};
    unsigned key_offset = 0;

    // Prepare input data
    for(unsigned i = 0; i < key_batch; i++)
        for(unsigned j = 0; j < key_length; j++){
            // if(i * in_words_adj + j < 10){
            // if((i * in_words_adj + j != 10) && (i * in_words_adj + j != 20))
                inbuff[i * in_words_adj + j] = (word_t) val_arr[i * in_words_adj + j + key_offset];
            // else
            //     inbuff[i * in_words_adj + j] = 1000;
            // }
            // else
            //     inbuff[i * in_words_adj + j] = (word_t) j;
        }

    for(unsigned i = 0; i < dma_in_size; i++)
	for(unsigned k = 0; k < VALUES_PER_WORD; k++)
	    mem[i].word[k] = inbuff[i * VALUES_PER_WORD + k];


    float Rs = R_f * std_f;
    unsigned result;
    unsigned added = 0;

    // Set golden output
    for(unsigned i = 0; i < key_batch; i++){
        for(unsigned j = 0; j < key_length; j++)
        {
            float val = val_arr[i * in_words_adj + j + key_offset];
            // if((i * in_words_adj + j == 10) || (i * in_words_adj + j == 20))
            //     val = 1000;
            bool filter = (fabs((float)val - avg_f) >= Rs);
            if(!filter){
                result = floor((float)(((val - (avg_f - Rs)) / (2*Rs)) * L));
                result = result % 2;
                outbuff_gold[i * out_words_adj + j] = (word_t) result;
            }
            else{
                outbuff_gold[i * out_words_adj + j] = 3;
                added += 1;
            }

            // std::cout << "Golden index " << i * out_words_adj + j << " val " << outbuff_gold[i * out_words_adj + j] << std::endl;
        }
    }


    //out_dma_word_t* mem_out = (out_dma_word_t*)&mem[out_offset];

    // Call the TOP function
    top(mem_out, mem,
        /* <<--args-->> */
        avg_f,
        key_length,
        std_f,
        R_f,
        L,
        key_batch,
        key_num,
        val_num,
        load, store);

    // Validate
    uint32_t out_offset = dma_in_size;
    for(unsigned i = 0; i < dma_out_size; i++)
	for(unsigned k = 0; k < VALUES_PER_WORD; k++){
	    // outbuff[i * VALUES_PER_WORD + k] = mem[out_offset + i].word[k];
	    outbuff_bit[i * VALUES_PER_WORD + k] = mem_out[out_offset + i].word[k];
	    //outbuff[i * VALUES_PER_WORD + k] = mem_out[out_offset + i].word[k];
            // std::cout << " mem val is " << std::bitset<32>(outbuff_bit[i * VALUES_PER_WORD + k]) << std::endl;
        }


    std::cout << "\n ^*^*^* VALIDATION ^*^*^* \n" << std::endl;

    int errors = 0;
    int skip = 0;
    int key_counter = 0;

    // for(unsigned i = 0; i < key_batch; i++)
    //     for(unsigned j = 0; j < key_length; j++){
    //         if(key_counter == key_num) break;
    //         unsigned index = i * out_words_adj + j;
    //         word_t val = outbuff[index - skip];
    //         word_t gold_val = outbuff_gold[index];
    //         if(gold_val != 3){
    //             if(!(i == key_batch - (ceil((float)skip/key_length)) && (j > skip - 1) )){
    //                 std::cout << "Calculated value " << val << " Golden value " << outbuff_gold[index] << " for index " << index - skip << std::endl;
    //                 if (outbuff[index - skip] != outbuff_gold[index]){
    //                     errors++;
    //                     std::cout << "ERROR" << std::endl;
    //                 }
    //             }
    //         }
    //         else{
    //             std::cout << "SKIPPING" << std::endl;
    //             skip += 1;
    //         }

    //         if((index - skip + 1) % key_length == 0 && index != 0){
    //             key_counter++;
    //             std::cout << "\n----------KEY " << key_counter << " DONE----------\n" << std::endl;
    //         }
    //     }

    int val_counter = 0;
    word_t* outbuff_val = (word_t*) &outbuff_bit[0];

    for(unsigned i = 0; i < key_batch; i++)
        for(unsigned j = 0; j < key_length; j++){
            if(key_counter != key_num){
            unsigned index = i * out_words_adj + j;
            // word_t val = outbuff[index - skip];
            unsigned bit = (index - skip) % 32;
            unsigned word = (index - skip) >> 5;
            ap_uint<1> val = outbuff_bit[word][bit];
            // std::cout << " word is " << std::bitset<32>(outbuff_bit[word])
            //           << " word " << word
            //           << " bit " << bit
            //           << std::endl;
            word_t gold_val = outbuff_gold[index];
            if(gold_val != 3){
                if(!(i == key_batch - (ceil((float)skip/key_length)) && (j > skip - 1) )){
                    std::cout << "Calculated value " << std::dec << val << " Golden value " << outbuff_gold[index] << " for index " << std::dec << index - skip << std::endl;
                    if (val != gold_val){
                        errors++;
                        std::cout << "ERROR" << std::endl;
                    }
                }
            }
            else{
                std::cout << "SKIPPING" << std::endl;
                skip += 1;
            }

            if((index - skip + 1) % key_length == 0 && index != 0){
                key_counter++;
                std::cout << "\n----------KEY " << std::dec << key_counter << " DONE----------" << std::endl;
                std::cout << "\nKEY IS: [ ";
                    for(int k = key_length / 32 - 1; k >= 0; k--)
                        std::cout << std::hex << outbuff_bit[word-k] << " ";
                    std::cout << "]\n" << std::endl;
            }
            }
            else if(val_counter != val_num){
                unsigned index = i * out_words_adj + j - skip;
                ap_uint<32> val = outbuff_bit[index];
                word_t val_word = 0;
                word_t gold_val = val_arr[index + key_offset];
                for(int b = 0; b < DATA_BITWIDTH; b++){
                    ap_uint<1> val_bit = val[b];
                    val_word[b] = gold_val[b];
                }
                //word_t gold_val = val_arr[index + key_offset];
                if(val_word != gold_val)
                    std::cout << "Calculated value " << std::dec << std::bitset<32>(val_word) << " Golden value " << gold_val << " for index " << std::dec << index << std::endl;
                if((index + 1) % key_length == 0 && index != 0)
                    val_counter++;
            }
        }


    float total = 100 * (float) errors / (key_length*key_batch);

    if (total > 1)
	std::cout << "Test FAILED with " << std::dec << total << "% errors." << std::endl;
    else{
	std::cout << "Keys generated: " << std::dec << key_counter << std::endl;
	std::cout << "Number of bit errors: " << std::dec << errors << std::endl;
	std::cout << "Test PASSED." << std::endl;
    }

    // Free memory

    free(mem);
    free(mem_out);
    free(inbuff);
    free(outbuff);
    free(outbuff_gold);

    return 0;
}
