/* Copyright (c) 2011-2022 Columbia University, System Level Design Group */
/* SPDX-License-Identifier: Apache-2.0 */

#include <stdio.h>
#ifndef __riscv
#include <stdlib.h>
#endif

#include <esp_accelerator.h>
#include <esp_probe.h>
#include <fixed_point.h>
//#include <math.h>
#include <stdbool.h>
#include "monitors.h"

#include "prc_utils.h"

//#include "input.h"
#include "input_full.h"
//#include "input_10mil_full.h"
//#include "input_1mil_full.h"
#define DATA_BITWIDTH 32

typedef int32_t token_t;
/* typedef float token_t; */

static unsigned DMA_WORD_PER_BEAT(unsigned _st)
{
        return (sizeof(void *) / _st);
}


#define TH_DEC_PT 14
#define TH_INC_PT 23

#define BB_0 0
#define BB_1 1
#define BB_2 2

#define BB_0_TILE 2
#define BB_1_TILE 4
#define BB_2_TILE 5
#define CPU_TILE 1

#define SLD_BRAIN_BIT_ALT 0x1f4
#define DEV_NAME "sld,brain_bit_alt_vivado"

/* <<--params-->> */
const float avg = 3.0677295382679177;
unsigned* avg_ptr = (unsigned*)&avg;
const int32_t key_length = 256;
const float std = 38.626628825256695;
unsigned* std_ptr = (unsigned*)&std;
const float R = 1.5;
unsigned* R_ptr = (unsigned*)&R;
/* const int32_t L = 1500; */
const int32_t key_batch = 1000;
int32_t key_num = 100;
const int32_t val_num = 0;
const int32_t tot_iter = 1;
const int32_t d = 7;
const int32_t h = 12;
float Rs = 57.939943238; //R * std;

static unsigned in_words_adj;
static unsigned out_words_adj;
static unsigned in_len;
static unsigned out_len;
static unsigned in_size;
static unsigned out_size;
static unsigned out_offset;
static unsigned mem_size;
/* static token_t avg_fixed; */
/* static token_t std_fixed; */
/* static token_t R_fixed; */

/* Size of the contiguous chunks for scatter/gather */
#define CHUNK_SHIFT 20
#define CHUNK_SIZE BIT(CHUNK_SHIFT)
#define NCHUNK(_sz) ((_sz % CHUNK_SIZE == 0) ?		\
			(_sz / CHUNK_SIZE) :		\
			(_sz / CHUNK_SIZE) + 1)

/* User defined registers */
/* <<--regs-->> */
#define BRAIN_BIT_ALT_H_REG 0x64
#define BRAIN_BIT_ALT_D_REG 0x60
#define BRAIN_BIT_ALT_TOT_ITER_REG 0x5c
#define BRAIN_BIT_ALT_VAL_NUM_REG 0x58
#define BRAIN_BIT_ALT_KEY_NUM_REG 0x54
#define BRAIN_BIT_ALT_AVG_REG 0x50
#define BRAIN_BIT_ALT_KEY_LENGTH_REG 0x4c
#define BRAIN_BIT_ALT_STD_REG 0x48
#define BRAIN_BIT_ALT_R_REG 0x44
/* #define BRAIN_BIT_ALT_L_REG 0x44 */
#define BRAIN_BIT_ALT_KEY_BATCH_REG 0x40


/* static int validate_buf(token_t *out, token_t *gold) */
static int validate_buf(token_t *out, token_t *gold)
{
	int i;
	int j;
	int k;
	unsigned errors = 0;
	int skip = 0;
	int key_counter = 0;
	/* int val_counter = 0; */
	int offset = 0;
	bool done = false;

	for (i = 0; i < key_batch; i++)
		for (j = 0; j < key_length; j++){
			if(key_counter != key_num){
				unsigned index = i * out_words_adj + j;
				/* token_t val = out[index - skip]; */
				int bit = (index - skip) % 32;
				int word = (index - skip) >> 5;
				token_t mask = out[word] >> bit;
				/* printf("out val is %d for word %d bit %d\n", out[word], word, bit); */
				token_t val = mask & 1;
				token_t gold_val = gold[index];
				unsigned reduce = (float)(skip/key_length);
				if(gold_val != 3){
					if(!(i == key_batch - reduce && (j > skip - 1) )){
						if (gold_val != val){
							printf("Calculated value %x Golden value %d for index %d \n", val, gold_val, (index-skip));
							errors++;
							printf("ERROR\n");
						}
					}
				}
				else{
					printf("SKIPPING\n");
					skip += 1;
				}

				if((index - skip + 1) % (key_length*(key_counter+1)) == 0 && index != 0){
					key_counter++;
					printf("\n----------KEY %d DONE----------\n", key_counter);
					printf("\nKEY IS: [ ");
					for(int k = key_length / 32 - 1; k >= 0; k--)
						printf("0x%x ", out[word-k]);
					printf("]\n\n");
				}
			}
			else if(!done){
				done = true;
				unsigned max_chunk = 1024;
				offset = i * out_words_adj + j - (skip % key_length);
				if(key_length > max_chunk && skip > key_length)
					offset += key_length - max_chunk;
				else if(key_length > max_chunk && skip < key_length)
					offset -= key_length % max_chunk;
			}
		}

	int index_offset = key_num * key_length / DATA_BITWIDTH;

	for(int i = 0; i < val_num; i++){
		unsigned index = index_offset + i + DATA_BITWIDTH / key_length;
                token_t val = out[index];
                token_t gold_val = (token_t) float_to_fixed32(val_arr[offset+i], 12);

		if(val != gold_val)
			printf("Calculated value %x Golden value %x for index %d \n", val, gold_val, index);
	}


	return errors;
}


static void init_buf (token_t *in, token_t * gold)
{
	int i;
	int j;

	for (i = 0; i < key_batch; i++)
		for (j = 0; j < key_length; j++){
                        float val = val_arr[i * in_words_adj + j];
                        in[i * in_words_adj + j] = (token_t) float_to_fixed32(val, 12);
                        //in[i * in_words_adj + j] = (token_t) val;
                }

        for (i = 0; i < key_batch; i++)
		for (j = 0; j < key_length; j++){
			float val = val_arr[i * in_words_adj + j];
                        bool filter = (float)(val - avg) >= Rs;
			int mul = 2* d;
			int mod = 2* h;
			if(!filter){
				//result can be a negative number here
				/* int result_alt = floor((float)(((val - avg) / (2*Rs)) * Rs * mul)); */
				int result_alt = (float)(((val - avg) / 2) * mul);
				result_alt = result_alt % mod;
				//but result can only be a positive number
				unsigned result = result_alt + mod;
				result = result % mod;
				// result = ((result_alt % mod) + mod) % mod;
				unsigned sum_result = 0;
				for(unsigned k = 0; k < h; k++)
					sum_result = sum_result + ((result >> k) % 2);
				result = sum_result;
				result = result % 2;
				gold[i * out_words_adj + j] = (token_t) result;
			}
			else{
				gold[i * out_words_adj + j] = 3;
			}
		}
}

static inline uint64_t get_counter() {
        uint64_t counter;
        asm volatile (
                "li t0, 0;"
                "csrr t0, mcycle;"
                "mv %0, t0"
                : "=r" ( counter )
                :
                : "t0"
        );
        return counter;
}

int main(int argc, char * argv[])
{
	int i;
	int n;
    int k;
	int ndev;
	struct esp_device *espdevs;
	struct esp_device *dev, *dev_bb_1, *dev_bb_2;
	unsigned done, done_1, done_both;
	unsigned **ptable;
	token_t *mem;
	token_t *gold;
	unsigned errors = 0;
	unsigned coherence;
            
    //Monitor variables bb_0
    esp_monitor_args_t bb_0_mon_args_lsb, bb_0_mon_args_msb;
    const int bb_0_ACC_TILE_IDX = BB_0_TILE;
    bb_0_mon_args_lsb.read_mode = ESP_MON_READ_SINGLE;
    bb_0_mon_args_lsb.tile_index = bb_0_ACC_TILE_IDX;
    bb_0_mon_args_lsb.mon_index = 16;//MON_DVFS_BASE_INDEX + 3;
    bb_0_mon_args_msb.read_mode = ESP_MON_READ_SINGLE;
    bb_0_mon_args_msb.tile_index = bb_0_ACC_TILE_IDX;
    bb_0_mon_args_msb.mon_index = 17;//MON_DVFS_BASE_INDEX + 3;
    unsigned long int bb_0_cycles_start_l, bb_0_cycles_end_l, bb_0_cycles_diff_l;
    unsigned long int bb_0_cycles_start_m, bb_0_cycles_end_m, bb_0_cycles_diff_m;
    unsigned long int bb_0_cycles_diff_l_prev;

    
    //Monitor variables bb_1_
    esp_monitor_args_t bb_1_mon_args_lsb, bb_1_mon_args_msb;
    const int bb_1_ACC_TILE_IDX = BB_1_TILE;
    bb_1_mon_args_lsb.read_mode = ESP_MON_READ_SINGLE;
    bb_1_mon_args_lsb.tile_index = bb_1_ACC_TILE_IDX;
    bb_1_mon_args_lsb.mon_index = 16;//MON_DVFS_BASE_INDEX + 3;
    bb_1_mon_args_msb.read_mode = ESP_MON_READ_SINGLE;
    bb_1_mon_args_msb.tile_index = bb_1_ACC_TILE_IDX;
    bb_1_mon_args_msb.mon_index = 17;//MON_DVFS_BASE_INDEX + 3;
    unsigned long int bb_1_cycles_start_l, bb_1_cycles_end_l, bb_1_cycles_diff_l;
    unsigned long int bb_1_cycles_start_m, bb_1_cycles_end_m, bb_1_cycles_diff_m;
    unsigned long int bb_1_cycles_diff_l_prev;
    
    /*
    bb_2_mon_args_msb.read_mode = ESP_MON_READ_SINGLE;
    bb_2_mon_args_msb.tile_index = CPU_TILE;
    bb_2_mon_args_msb.mon_index = 17;//MON_DVFS_BASE_INDEX + 3;
*/  unsigned long int bb_2_cycles_start_l, bb_2_cycles_end_l, bb_2_cycles_diff_l;
    unsigned long int bb_2_cycles_start_m, bb_2_cycles_end_m, bb_2_cycles_diff_m;
    unsigned long int bb_2_cycles_diff_l_prev;
    
    unsigned long long total_time = 0;

    float max_dpr_threshold = 1.3;
    float min_dpr_threshold = 0.7;

	unsigned long long tmp = 0, tmp_2 = 0;
    unsigned long long acc_run_time = 0, prev_acc_runtime = 0;
    const float new_R = 0.5;
    unsigned* new_R_ptr = (unsigned*) &new_R;

    unsigned long num_accs;
    unsigned long long exec_time_us, tot_exec_time_us = 0;
    float throughput, prev_throughput, tot_throughput = 0;
    unsigned long long total_cycles = 0;
    unsigned scaling_factor = 1e6;
     
    unsigned do_reconfig = 0, do_last_th = 0;

    unsigned key_num_step;
    //data to hold values
    const unsigned N_ITERS = 20;
    unsigned long long data[N_ITERS][7];
    unsigned N_runs = 12;
    uint64_t app_exec_cycles;
    
    //N_runs = N_ITERS;

    if (DMA_WORD_PER_BEAT(sizeof(token_t)) == 0) {
		in_words_adj = key_length;
		out_words_adj = key_length;
	} else {
		in_words_adj = round_up(key_length, DMA_WORD_PER_BEAT(sizeof(token_t)));
		out_words_adj = round_up(key_length, DMA_WORD_PER_BEAT(sizeof(token_t)));
	}
	in_len = in_words_adj * (key_batch);
	out_len = out_words_adj * (key_batch);
	in_size = in_len * sizeof(token_t);
	out_size = out_len * sizeof(token_t);
	out_offset  = in_len;
	mem_size = (out_offset * sizeof(token_t)) + out_size;

        /* avg_fixed = float_to_fixed32(avg, 12); */
        /* std_fixed = float_to_fixed32(std, 12); */
        /* R_fixed = float_to_fixed32(R, 12); */

	// Search for the device
	printf("Scanning device tree... \n");

	ndev = probe(&espdevs, VENDOR_SLD, SLD_BRAIN_BIT_ALT, DEV_NAME);
	if (ndev == 0) {
		printf("brain_bit_alt not found\n");
		return 0;
	}

	//for (n = 2; n < 3; n++) {
        //configure bb_0
        float wave = 1.0;
        n = 0;
        num_accs = 1;

		printf("**************** %s.%d ****************\n", DEV_NAME, n);
        struct esp_device decoupler;

        decoupler.addr = 0x60090580;
        iowrite32(&decoupler, 0X30, 0);

		dev = &espdevs[n];
        //reconfigure_FPGA(&espdevs[0], BB_0);

		// Check DMA capabilities
		if (ioread32(dev, PT_NCHUNK_MAX_REG) == 0) {
			printf("  -> scatter-gather DMA is disabled. Abort.\n");
			return 0;
		}

		if (ioread32(dev, PT_NCHUNK_MAX_REG) < NCHUNK(mem_size)) {
			printf("  -> Not enough TLB entries available. Abort.\n");
			return 0;
		}

		// Allocate memory
		gold = aligned_malloc(out_size);
		mem = aligned_malloc(mem_size);
		printf("  memory buffer base-address = %p\n", mem);

		// Alocate and populate page table
		ptable = aligned_malloc(NCHUNK(mem_size) * sizeof(unsigned *));
		for (i = 0; i < NCHUNK(mem_size); i++)
			ptable[i] = (unsigned *) &mem[i * (CHUNK_SIZE / sizeof(token_t))];

		printf("  ptable = %p\n", ptable);
		printf("  nchunk = %lu\n", NCHUNK(mem_size));

#ifndef __riscv
        for (coherence = ACC_COH_NONE; coherence <= ACC_COH_RECALL; coherence++) {
#else
		{
			/* TODO: Restore full test once ESP caches are integrated */
            coherence = ACC_COH_NONE;
#endif
			printf("  --------------------\n");
			printf("  Generate input...\n");

			/* for(int k = 0; k < in_len+out_len; k++) */
			/* 	mem[k] = 0; */

			init_buf(mem, gold);

			// Pass common configuration parameters

			iowrite32(dev, SELECT_REG, ioread32(dev, DEVID_REG));
			iowrite32(dev, COHERENCE_REG, coherence);

#ifndef __sparc
			iowrite32(dev, PT_ADDRESS_REG, (unsigned long long) ptable);
#else
			iowrite32(dev, PT_ADDRESS_REG, (unsigned) ptable);
#endif
			iowrite32(dev, PT_NCHUNK_REG, NCHUNK(mem_size));
			iowrite32(dev, PT_SHIFT_REG, CHUNK_SHIFT);

			// Use the following if input and output data are not allocated at the default offsets
			iowrite32(dev, SRC_OFFSET_REG, 0x0);
			iowrite32(dev, DST_OFFSET_REG, 0x0);

			// Pass accelerator-specific configuration parameters
			/* <<--regs-config-->> */
			iowrite32(dev, BRAIN_BIT_ALT_AVG_REG, *avg_ptr);
                        /* iowrite32(dev, BRAIN_BIT_ALT_AVG_REG, avg_fixed); */
			iowrite32(dev, BRAIN_BIT_ALT_KEY_LENGTH_REG, key_length);
                        iowrite32(dev, BRAIN_BIT_ALT_STD_REG, *std_ptr);
			/* iowrite32(dev, BRAIN_BIT_ALT_STD_REG, std_fixed); */
                        iowrite32(dev, BRAIN_BIT_ALT_R_REG, *R_ptr);
			/* iowrite32(dev, BRAIN_BIT_ALT_R_REG, R_fixed); */
			/* iowrite32(dev, BRAIN_BIT_ALT_L_REG, L); */
			iowrite32(dev, BRAIN_BIT_ALT_KEY_BATCH_REG, key_batch);
			iowrite32(dev, BRAIN_BIT_ALT_KEY_NUM_REG, key_num);
			iowrite32(dev, BRAIN_BIT_ALT_VAL_NUM_REG, val_num);
			iowrite32(dev, BRAIN_BIT_ALT_TOT_ITER_REG, tot_iter);
			iowrite32(dev, BRAIN_BIT_ALT_D_REG, d);
			iowrite32(dev, BRAIN_BIT_ALT_H_REG, h);

			// Flush (customize coherence model here)
			esp_flush(coherence);

            for(key_num_step = 0; key_num_step < N_ITERS; key_num_step++) {
				uint64_t start_ctr =  get_counter(); 

			for(k = 0; k < N_runs; k++){
				// Start accelerators
				//printf("  Start...\n");

				bb_0_cycles_start_l = esp_monitor(bb_0_mon_args_lsb, NULL);
				bb_0_cycles_start_m = esp_monitor(bb_0_mon_args_msb, NULL);

				iowrite32(dev, CMD_REG, CMD_MASK_START);
                /*
                //check if reconfiguration is requested in the previous run
                if(do_reconfig) {
                    dev_bb_1 = &espdevs[1];
                    reconfigure_FPGA(dev_bb_1, BB_1);
                    
                    //check if the accelerator is done
                    done = ioread32(dev, STATUS_REG);
                    done &= ioread32(dev, STATUS_REG);
                     
                    while (!done) {
                        done = ioread32(dev, STATUS_REG);
                        done &= STATUS_MASK_DONE;
                    }
                    
                    //iowrite32(dev, CMD_REG, 0x0);
                    do_reconfig = 0;
                    do_last_th = 1;
                    //goto do_last_calc; //break;
                }
               */ 

				// Wait for completion
				done = 0;
				while (!done) {
					done = ioread32(dev, STATUS_REG);
					done &= STATUS_MASK_DONE;
				}

				bb_0_cycles_end_l = esp_monitor(bb_0_mon_args_lsb, NULL);
				bb_0_cycles_diff_l = sub_monitor_vals(bb_0_cycles_start_l, bb_0_cycles_end_l);
				/* printf("Monitr time is %u %u %u \n", cycles_start, cycles_end, cycles_diff); */

				/* unsigned nano = cycles_diff * 20; //50MHz - vcu707 */
				/* total_time += nano; */
				/* printf("Accelerator runtime: %u ns \n", nano); */

				//printf("Accelerator 0 runtime: %u cycles LSB\n", bb_0_cycles_diff_l);
				//printf("Accelerator 0 runtime: %u cycles MSB\n", bb_0_cycles_diff_m);
				tmp = (bb_0_cycles_diff_m << 32) + bb_0_cycles_diff_l;
                
                exec_time_us = (13 * tmp) / 1e3;
                throughput = (num_accs * key_num * scaling_factor) / exec_time_us;
				int th_int = (unsigned long) throughput;
                //tmp = cycles_diff_l;
				//printf("Accelerator total runtime: %llu cycles %u(us)  throughput -> %d(keys/s) \n", tmp, exec_time_us,  th_int);
				
                total_time += tmp;
                tot_exec_time_us += exec_time_us;
                tot_throughput += throughput;
                //Check if throughout is reduced
				//if(k > 0 && bb_0_cycles_diff_l > max_dpr_threshold * bb_0_cycles_diff_l_prev) {
				//if(k > 0 && (throughput < min_dpr_threshold * prev_throughput || do_last_th == 1)) {
				/*if(k > 0 && (throughput < min_dpr_threshold * prev_throughput)) {
					printf("Use DPR to increase throughput!!!!\n");
                    if(do_last_th == 1) {
                        do_last_th = 0;
                        iowrite32(dev, CMD_REG, 0x0);
                        break;
                    }
                    else
                        do_reconfig = 1;
                    
                    //dev_bb_1 = &espdevs[1];
                    //reconfigure_FPGA(dev_bb_1, BB_1);
				    //iowrite32(dev, CMD_REG, 0x0);
                    //break;
                }*/
                
                //Check if throughout is increased
				//if(k > 0 && bb_0_cycles_diff_l < min_dpr_threshold * bb_0_cycles_diff_l_prev) 
					//printf("Use DPR to decrease throughput!!!!\n");

				//Saving the previous cycles_diff_l to compare if there's a reduction of throughput
				bb_0_cycles_diff_l_prev = bb_0_cycles_diff_l;
                prev_throughput = throughput;

				iowrite32(dev, CMD_REG, 0x0);

				//printf("  Done run %d\n", k);
				//printf("  validating...\n");
               /* 
                if(k == TH_DEC_PT){
					printf("  Reduce throughput manually \n");
					//const float new_R = 0.5;
					Rs = new_R * std;
					//unsigned* new_R_ptr = (unsigned*)&new_R;
					iowrite32(dev, BRAIN_BIT_ALT_R_REG, *new_R_ptr);
				}
		   */	/*	
                if(k == TH_INC_PT){
					printf("  Increase throughput manually \n");
					const float new_R = 1.5;
					Rs = new_R * std;
					unsigned* new_R_ptr = (unsigned*)&new_R;
					iowrite32(dev, BRAIN_BIT_ALT_R_REG, *new_R_ptr);
				}
			*/
            }

            //bb_2_cycles_end_m = esp_monitor(bb_2_mon_args_msb, NULL);
            //bb_2_cycles_diff_m = sub_monitor_vals(bb_2_cycles_start_m, bb_2_cycles_end_m);
             
            uint64_t end_ctr = get_counter();
            app_exec_cycles = end_ctr - start_ctr;

            data[key_num_step][0] = total_time; //total number of cycles
            data[key_num_step][1] = tot_exec_time_us; 
            data[key_num_step][2] = tot_throughput;
            data[key_num_step][3] = key_num; 
            data[key_num_step][4] = k; //number of iterations
            data[key_num_step][5] = app_exec_cycles;
 
            tot_exec_time_us = 0;
            tot_throughput = 0;
            total_time = 0;
            
            key_num += 50;
            iowrite32(dev, BRAIN_BIT_ALT_KEY_NUM_REG, key_num);
        }
            
            for(int m = 0; m < N_ITERS; m++) {
                float avg_th = data[m][2]/k;
                int avg_th_int = avg_th;
                float avg_cycles = data[m][0]/k;
                int avg_cycles_int = avg_cycles;
                unsigned long avg_time = data[m][1]/k;
                
                float app_exec_time_us = (13 * data[key_num_step][5]) / 1e3;
                uint64_t app_exec_time_int = app_exec_time_us;
                float th_app = (num_accs * data[m][3] * scaling_factor) / app_exec_time_us;
                int th_app_int =  th_app;
 
                int avg_time_int = avg_time;
                printf("keys/invocation, %u, num_iterations, %u, avg cycles, %d, avg_time, %d, avg_througput, %d, app_exec_cycles, %lu, app_exec_time, %llu \n", 
                       data[m][3], data[m][4], avg_cycles_int, avg_time_int, avg_th_int, app_exec_cycles, app_exec_time_int);
            }
/*
            int temp_th = tot_throughput;
			int temp_exec = tot_exec_time_us;

            total_time = total_time / N_runs;
			printf("Average runtime: %u %u %u cycles th %d exec %d\n", key_num_step, k, total_time, temp_th, temp_exec);
*/
			float total = 100 * (float) errors / (key_length*key_batch);
			if (total > 1)
				printf("  ... FAIL\n");
			else
				printf("  ... PASS\n");
}		
        
        //reset previous acc time for bb0
        //bb_0_cycles_diff_l_prev = 0;

//#define RUN_SECOND_ACC

    printf("******************** Restarting execution with 2 accs ************************ \n\n");
#ifdef RUN_SECOND_ACC		
#ifndef __riscv
		for (coherence = ACC_COH_NONE; coherence <= ACC_COH_RECALL; coherence++) {
#else
		{
			/* TODO: Restore full test once ESP caches are integrated */
			coherence = ACC_COH_NONE;
#endif
			// Pass common configuration parameters

			iowrite32(dev_bb_1, SELECT_REG, ioread32(dev_bb_1, DEVID_REG));
			iowrite32(dev_bb_1, COHERENCE_REG, coherence);

#ifndef __sparc
			iowrite32(dev_bb_1, PT_ADDRESS_REG, (unsigned long long) ptable);
#else
			iowrite32(dev_bb_1, PT_ADDRESS_REG, (unsigned) ptable);
#endif
			iowrite32(dev_bb_1, PT_NCHUNK_REG, NCHUNK(mem_size));
			iowrite32(dev_bb_1, PT_SHIFT_REG, CHUNK_SHIFT);

			// Use the following if input and output data are not allocated at the default offsets
			iowrite32(dev_bb_1, SRC_OFFSET_REG, 0x0);
			iowrite32(dev_bb_1, DST_OFFSET_REG, 0x0);

			// Pass accelerator-specific configuration parameters
			/* <<--regs-config-->> */
			iowrite32(dev_bb_1, BRAIN_BIT_ALT_AVG_REG, *avg_ptr);
            /* iowrite32(dev_bb_1, BRAIN_BIT_ALT_AVG_REG, avg_fixed); */
			iowrite32(dev_bb_1, BRAIN_BIT_ALT_KEY_LENGTH_REG, key_length);
            iowrite32(dev_bb_1, BRAIN_BIT_ALT_STD_REG, *std_ptr);
			/* iowrite32(dev_bb_1, BRAIN_BIT_ALT_STD_REG, std_fixed); */
            iowrite32(dev_bb_1, BRAIN_BIT_ALT_R_REG, *new_R_ptr);
			/* iowrite32(dev_bb_1, BRAIN_BIT_ALT_R_REG, R_fixed); */
			/* iowrite32(dev_bb_1, BRAIN_BIT_ALT_L_REG, L); */
			iowrite32(dev_bb_1, BRAIN_BIT_ALT_KEY_BATCH_REG, key_batch);
			iowrite32(dev_bb_1, BRAIN_BIT_ALT_KEY_NUM_REG, key_num);
			iowrite32(dev_bb_1, BRAIN_BIT_ALT_VAL_NUM_REG, val_num);
			iowrite32(dev_bb_1, BRAIN_BIT_ALT_TOT_ITER_REG, tot_iter);
			iowrite32(dev_bb_1, BRAIN_BIT_ALT_D_REG, d);
			iowrite32(dev_bb_1, BRAIN_BIT_ALT_H_REG, h);

			// Flush (customize coherence model here)
			esp_flush(coherence);
            
            k += 1;
			for(; k < N_runs; k++){
				// Start accelerators
				printf("  Start...\n");

				bb_0_cycles_start_l = esp_monitor(bb_0_mon_args_lsb, NULL);
				bb_0_cycles_start_m = esp_monitor(bb_0_mon_args_msb, NULL);
				bb_1_cycles_start_l = esp_monitor(bb_1_mon_args_lsb, NULL);
				bb_1_cycles_start_m = esp_monitor(bb_1_mon_args_msb, NULL);
                
                //start both accelerators
				iowrite32(dev, CMD_REG, CMD_MASK_START);
                iowrite32(dev_bb_1, CMD_REG, CMD_MASK_START);
               
				// Wait for completion
				done = 0;
				done_1 = 0;
                done_both = 0;
                while (!done_both) {
					done = ioread32(dev, STATUS_REG);
					done &= STATUS_MASK_DONE;
					done_1 = ioread32(dev_bb_1, STATUS_REG);
					done_1 &= STATUS_MASK_DONE;
				    done_both = done & done_1;
                }

				bb_0_cycles_end_l = esp_monitor(bb_0_mon_args_lsb, NULL);
				bb_0_cycles_diff_l = sub_monitor_vals(bb_0_cycles_start_l, bb_0_cycles_end_l);
				bb_1_cycles_end_l = esp_monitor(bb_1_mon_args_lsb, NULL);
                bb_1_cycles_diff_l = sub_monitor_vals(bb_1_cycles_start_l, bb_1_cycles_end_l);
				/* printf("Monitr time is %u %u %u \n", cycles_start, cycles_end, cycles_diff); */

				/* unsigned nano = cycles_diff * 20; //50MHz - vcu707 */
				/* total_time += nano; */
				/* printf("Accelerator runtime: %u ns \n", nano); */

				printf("Accelerator 0 runtime: %u cycles LSB\n", bb_0_cycles_diff_l);
				printf("Accelerator 0 runtime: %u cycles MSB\n", bb_0_cycles_diff_m);
				tmp = (bb_0_cycles_diff_m << 32) + bb_0_cycles_diff_l;
				//printf("Accelerator 0 total runtime: %llu cycles\n", tmp);
				//tmp = cycles_diff_l;
				
                printf("Accelerator 1 runtime: %u cycles LSB\n", bb_1_cycles_diff_l);
				printf("Accelerator 1 runtime: %u cycles MSB\n", bb_1_cycles_diff_m);
				tmp_2 = (bb_1_cycles_diff_m << 32) + bb_1_cycles_diff_l;
				//tmp = cycles_diff_l;
				
                //calcuate throughput
                num_accs = 2;
                acc_run_time = (tmp + tmp_2) / 2; 
                exec_time_us = (13 * acc_run_time) / 1e3;
                throughput = (num_accs * key_num * scaling_factor) / exec_time_us;
                  
				int th_int = (unsigned long) throughput;
                printf("Accelerator total runtime: %llu cycles %u(us)  throughput -> %u(keys/s) \n", acc_run_time, exec_time_us, th_int);	
                //tmp = cycles_diff_l;
				
	            total_time += acc_run_time;
                tot_exec_time_us += acc_run_time;
                tot_throughput += throughput;
                //Check if throughout is reduced

                //prev_acc_runtime = (bb_0_cycles_diff_l_prev + bb_1_cycles_diff_l_prev) / 2;
                //Check if throughout is reduced
				/*if(k > TH_DEC_PT && acc_run_time > max_dpr_threshold * prev_acc_runtime)
					printf("Use DPR to increase throughput!!!!\n");
                */
/*                //Check if throughout is increased
				if(k > TH_INC_PT && throughput > max_dpr_threshold * prev_throughput) {
					printf("Use DPR to decrease throughput!!!!\n");
				    do_reconfig = 1;
                    iowrite32(dev_bb_1, CMD_REG, 0x0);
				    iowrite32(dev, CMD_REG, 0x0);
                    break;
                }
i*/
				//Saving the previous cycles_diff_l to compare if there's a reduction of throughput
				bb_0_cycles_diff_l_prev = bb_0_cycles_diff_l;
				bb_1_cycles_diff_l_prev = bb_1_cycles_diff_l;
                prev_throughput = throughput;

				iowrite32(dev, CMD_REG, 0x0);
				iowrite32(dev_bb_1, CMD_REG, 0x0);

				printf("  Done run %d\n", k);
				//printf("  validating...\n");

				//Manually control throughput by tuning R
				/*if(k == 2){
					printf("  Reduce throughput manually \n");
					const float new_R = 0.5;
					Rs = new_R * std;
					unsigned* new_R_ptr = (unsigned*)&new_R;
					iowrite32(dev, BRAIN_BIT_ALT_R_REG, *new_R_ptr);
				}*/
			/*	if(k == TH_INC_PT){
					printf("  Increase throughput manually \n");
					const float new_R_th = 1.5;
					Rs = new_R_th * std;
					unsigned* new_R_ptr_th = (unsigned*)&new_R_th;
					iowrite32(dev, BRAIN_BIT_ALT_R_REG, *new_R_ptr_th);
					iowrite32(dev_bb_1, BRAIN_BIT_ALT_R_REG, *new_R_ptr_th);
				}*/
			}

            data[key_num_step][0] = total_time; //total number of cycles
            data[key_num_step][1] = tot_exec_time_us; 
            data[key_num_step][2] = tot_throughput;
            data[key_num_step][3] = key_num; 
            data[key_num_step][4] = k; //number of iterations
             
            tot_exec_time_us = 0;
            tot_throughput = 0;
            total_time = 0;
            
        }
            key_num += 50;
            iowrite32(dev, BRAIN_BIT_ALT_KEY_NUM_REG, key_num);
            iowrite32(dev, BRAIN_BIT_ALT_R_REG, *R_ptr);    


    }
            
            for(int m = 0; m < N_ITERS; m++) {
                float avg_th = data[m][2]/k;
                int avg_th_int = avg_th;
                float avg_cycles = data[m][0]/k;
                int avg_cycles_int = avg_cycles;
                unsigned long avg_time = data[m][1]/k;
                int avg_time_int = avg_time;
                printf("keys/invocation, %u, num_iterations, %u, avg cycles, %d, avg_time, %d, avg_througput, %d \n",  
                       data[m][3], data[m][4], avg_cycles_int, avg_time_int, avg_th_int);

           
            } 
		/*	
            total_time = total_time / N_runs;
			printf("Average runtime: %u cycles \n", total_time);

			float total = 100 * (float) errors / (key_length*key_batch);
			if (total > 1)
				printf("  ... FAIL\n");
			else
				printf("  ... PASS\n");
		*/
        }
#endif

#ifdef THIRD_LOOP
    printf("******************** Restarting execution with just 1 acc ************************ \n\n");

#ifndef __riscv
        for (coherence = ACC_COH_NONE; coherence <= ACC_COH_RECALL; coherence++) {
#else
		{
			/* TODO: Restore full test once ESP caches are integrated */
            coherence = ACC_COH_NONE;
#endif
			printf("  --------------------\n");
			printf("  Generate input...\n");

			/* for(int k = 0; k < in_len+out_len; k++) */
			/* 	mem[k] = 0; */

			// Flush (customize coherence model here)
			esp_flush(coherence);
            
            num_accs = 1;
			for(; k < N_runs; k++){
				// Start accelerators
				printf("  Start...\n");

				bb_0_cycles_start_l = esp_monitor(bb_0_mon_args_lsb, NULL);
				bb_0_cycles_start_m = esp_monitor(bb_0_mon_args_msb, NULL);

				iowrite32(dev, CMD_REG, CMD_MASK_START);
                if(do_reconfig == 1) {
                    dev_bb_1 = &espdevs[1];
                    //reconfigure_FPGA(dev_bb_1, BB_1); 
                    do_reconfig = 0;
                }

				// Wait for completion
				done = 0;
				while (!done) {
					done = ioread32(dev, STATUS_REG);
					done &= STATUS_MASK_DONE;
				}

				bb_0_cycles_end_l = esp_monitor(bb_0_mon_args_lsb, NULL);
				bb_0_cycles_diff_l = sub_monitor_vals(bb_0_cycles_start_l, bb_0_cycles_end_l);
				/* printf("Monitr time is %u %u %u \n", cycles_start, cycles_end, cycles_diff); */

				/* unsigned nano = cycles_diff * 20; //50MHz - vcu707 */
				/* total_time += nano; */
				/* printf("Accelerator runtime: %u ns \n", nano); */

				printf("Accelerator 0 runtime: %u cycles LSB\n", bb_0_cycles_diff_l);
				printf("Accelerator 0 runtime: %u cycles MSB\n", bb_0_cycles_diff_m);
				tmp = (bb_0_cycles_diff_m << 32) + bb_0_cycles_diff_l;
                
                exec_time_us = (12 * tmp) / 1e3;
                throughput = (num_accs * key_batch * scaling_factor) / exec_time_us;
				//tmp = cycles_diff_l;
				printf("Accelerator total runtime: %llu cycles %u(us)  throughput -> %u(keys/s) \n", tmp, exec_time_us,  throughput);
				total_time += tmp;
                
                //Check if throughout is reduced
				//if(k > 0 && bb_0_cycles_diff_l > max_dpr_threshold * bb_0_cycles_diff_l_prev) {
				if(k > 0 && throughput < min_dpr_threshold * prev_throughput)
					printf("Use DPR to increase throughput!!!!\n");

                    //dev_bb_1 = &espdevs[1];
                    //reconfigure_FPGA(dev_bb_1, BB_1);
				    //iowrite32(dev, CMD_REG, 0x0);
                    //break;
                
                                //Check if throughout is increased
				if(k > 0 && bb_0_cycles_diff_l < min_dpr_threshold * bb_0_cycles_diff_l_prev) 
					printf("Use DPR to decrease throughput!!!!\n");

				//Saving the previous cycles_diff_l to compare if there's a reduction of throughput
				bb_0_cycles_diff_l_prev = bb_0_cycles_diff_l;
                prev_throughput = throughput;

				iowrite32(dev, CMD_REG, 0x0);

				printf("  Done run %d\n", k);
				//printf("  validating...\n");

				//Manually control throughput by tuning R
				if(k == TH_DEC_PT){
					printf("  Reduce throughput manually \n");
					//const float new_R = 0.5;
					Rs = new_R * std;
					//unsigned* new_R_ptr = (unsigned*)&new_R;
					iowrite32(dev, BRAIN_BIT_ALT_R_REG, *new_R_ptr);
				}
				/*
                if(k == 7){
					printf("  Increase throughput manually \n");
					const float new_R = 1.5;
					Rs = new_R * std;
					unsigned* new_R_ptr = (unsigned*)&new_R;
					iowrite32(dev, BRAIN_BIT_ALT_R_REG, *new_R_ptr);
				}*/
			}
            
            data[key_num_step][0] = total_time; //total number of cycles
            data[key_num_step][1] = tot_exec_time_us; 
            data[key_num_step][2] = tot_throughput;
            data[key_num_step][3] = key_num; 
            data[key_num_step][4] = k; //number of iterations
            
            tot_exec_time_us = 0;
            tot_throughput = 0;
            total_time = 0;
            
            key_num += 50;
            iowrite32(dev, BRAIN_BIT_ALT_KEY_NUM_REG, key_num);
            
        }
            
            for(int m = 0; m < N_ITERS; m++) {
                float avg_th = data[m][2]/k;
                int avg_th_int = avg_th;
                float avg_cycles = data[m][0]/k;
                int avg_cycles_int = avg_cycles;
                unsigned long avg_time = data[m][1]/k;
                int avg_time_int = avg_time;
                printf("keys/invocation, %u, num_iterations, %u, avg cycles, %d, avg_time, %d, avg_througput, %d \n",  
                       data[m][3], data[m][4], avg_cycles_int, avg_time_int, avg_th_int);

			total_time = total_time / N_runs;
			printf("Average runtime: %u cycles \n", total_time);

			float total = 100 * (float) errors / (key_length*key_batch);
			if (total > 1)
				printf("  ... FAIL\n");
			else
				printf("  ... PASS\n");
		}
#endif        
        //reset previous acc time for bb0
        bb_0_cycles_diff_l_prev = 0;


        aligned_free(ptable);
		aligned_free(mem);
		aligned_free(gold);
//	}
	return 0;
}
