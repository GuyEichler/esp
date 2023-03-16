// Copyright (c) 2011-2022 Columbia University, System Level Design Group
// SPDX-License-Identifier: Apache-2.0

#include "libesp.h"
#include "cfg.h"
#include <stdio.h>
#include <sys/time.h>
#include "cfg.h"

typedef int32_t token_t;

//#define SYS_RAND 
#define DEV_RAND 
#define KEY_LEN 256
#define KEY_BYTES 32
#define MAX_ITER  100

#define RAND_MODE SYS_RAND

const unsigned num_rand_words = KEY_LEN/ (sizeof(token_t) * 8); 
const unsigned key_len_bytes = KEY_LEN / 8;
unsigned char *keys;
char *key_chars;
struct timeval hw_begin_0, hw_end_0;
int i, k;
long seconds = 0, microseconds = 0;
double elapsed = 0.0, avg_elapsed_us = 0.0, total_time = 0.0;

int run_dev_rand()
{
    FILE *f;
 
    keys = malloc(key_len_bytes);
    if(keys == NULL) {
        printf("Error allocating mem \n");
        return -1;
    }
    
    //average time to read random numbers 
    for(k = 0; k < MAX_ITER; k++) {
        
        //start timer
        gettimeofday(&hw_begin_0, 0);
        
        f = fopen("/dev/random", "r"); 
        if(f == NULL) {
            printf("Error opening /dev/random \n");
            return -1;
        }

        fread((void *) keys, key_len_bytes, 1, f);
        
        fclose(f);
        
        //end timer
        gettimeofday(&hw_end_0, 0);
        
        //print rand number is hex
        for(i = 0; i < num_rand_words; i++)
            printf("%x", (unsigned int) keys[i * sizeof(token_t)]);
        
        printf("\n");
 
        seconds = hw_end_0.tv_sec - hw_begin_0.tv_sec;
        microseconds = hw_end_0.tv_usec - hw_begin_0.tv_usec;
        elapsed += seconds + microseconds * 1e-6;
    }
    
    free(keys);

    total_time = (double) elapsed * 1e6;
    avg_elapsed_us = total_time / MAX_ITER;
     
    printf("Time measured: %u Iterations, avg time %.4f us.\n", MAX_ITER, avg_elapsed_us);
    
    return 0;
}

int run_dev_urand()
{
    FILE *f;
 
    keys = malloc(key_len_bytes);
    if(keys == NULL) {
        printf("Error allocating mem \n");
        return -1;
    }
    
    //average time to read random numbers 
    for(k = 0; k < MAX_ITER; k++) {
        
        //start timer
        gettimeofday(&hw_begin_0, 0);
        
        f = fopen("/dev/urandom", "r"); 
        if(f == NULL) {
            printf("Error opening /dev/urandom \n");
            return -1;
        }

        fread((void *) keys, key_len_bytes, 1, f);
        
        fclose(f);
        
        //end timer
        gettimeofday(&hw_end_0, 0);
        
        //print rand number is hex
        for(i = 0; i < num_rand_words; i++)
            printf("%x", (unsigned int) keys[i * sizeof(token_t)]);
        
        printf("\n");
 
        seconds = hw_end_0.tv_sec - hw_begin_0.tv_sec;
        microseconds = hw_end_0.tv_usec - hw_begin_0.tv_usec;
        elapsed += seconds + microseconds * 1e-6;
    }
    
    free(keys);

    total_time = (double) elapsed * 1e6;
    avg_elapsed_us = total_time / MAX_ITER;
     
    printf("Time measured: %u Iterations, avg time %.4f us.\n", MAX_ITER, avg_elapsed_us);
    
    return 0;
}

int run_rand()
{
   keys = malloc(key_len_bytes);

    for (k = 0; k < MAX_ITER; k++) {
        gettimeofday(&hw_begin_0, 0);

        srand(k);
        for(i = 0; i < num_rand_words; i++){
            keys[i * sizeof(token_t)] = rand();
        }

        gettimeofday(&hw_end_0, 0);

        for(i = 0; i < num_rand_words; i++)
            printf("%x", (unsigned int) keys[i * sizeof(token_t)]);

        printf("\n");
        seconds = hw_end_0.tv_sec - hw_begin_0.tv_sec;
        microseconds = hw_end_0.tv_usec - hw_begin_0.tv_usec;
        elapsed += seconds + microseconds*1e-6;
    }
    
    avg_elapsed_us = (double) (elapsed * 1e6) / MAX_ITER;
    
    printf("Time measured: %0.4fus.\n", avg_elapsed_us);
    
    free(keys);
    return 0;
}

int main(int argc, char *argv[]) 
{   
    if(argc < 2 || *argv[1] == '1') {
        printf("running /dev/random \n");
        run_dev_rand();
    }    
    else if (*argv[1] == '2') {
        printf("running /dev/urandom \n");
        run_dev_urand();
    }
    else if (*argv[1] == '3') {
        printf("running rand \n"); 
        run_rand();
    }

    return 0;
}
