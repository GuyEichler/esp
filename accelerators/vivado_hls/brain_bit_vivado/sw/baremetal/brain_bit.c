/* Copyright (c) 2011-2022 Columbia University, System Level Design Group */
/* SPDX-License-Identifier: Apache-2.0 */

#include <stdio.h>
#ifndef __riscv
#include <stdlib.h>
#endif

#include <esp_accelerator.h>
#include <esp_probe.h>
#include <fixed_point.h>
#include <math.h>
#include <stdbool.h>
#include "monitors.h"

//#include "input.h"
#include "input_full.h"
//#include "input_10mil_full.h"
#define DATA_BITWIDTH 32

typedef int32_t token_t;
/* typedef float token_t; */

static unsigned DMA_WORD_PER_BEAT(unsigned _st)
{
        return (sizeof(void *) / _st);
}


#define SLD_BRAIN_BIT 0x1e4
#define DEV_NAME "sld,brain_bit_vivado"

/* <<--params-->> */
const float avg = 3.0677295382679177;
unsigned* avg_ptr = (unsigned*)&avg;
const int32_t key_length = 1536;
const float std = 38.626628825256695;
unsigned* std_ptr = (unsigned*)&std;
const float R = 1.5;
unsigned* R_ptr = (unsigned*)&R;
const int32_t L = 1500;
const int32_t key_batch = 5;
const int32_t key_num = 1;
const int32_t val_num = 0;
const int32_t tot_iter = 1;
const float Rs = R * std;

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
#define BRAIN_BIT_TOT_ITER_REG 0x60
#define BRAIN_BIT_VAL_NUM_REG 0x5c
#define BRAIN_BIT_KEY_NUM_REG 0x58
#define BRAIN_BIT_AVG_REG 0x54
#define BRAIN_BIT_KEY_LENGTH_REG 0x50
#define BRAIN_BIT_STD_REG 0x4c
#define BRAIN_BIT_R_REG 0x48
#define BRAIN_BIT_L_REG 0x44
#define BRAIN_BIT_KEY_BATCH_REG 0x40


/* static int validate_buf(token_t *out, token_t *gold) */
static int validate_buf(token_t *out, token_t *gold)
{
	int i;
	int j;
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
				unsigned reduce = (ceil((float)skip/key_length));
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
                        bool filter = (fabs((float)val - avg) >= Rs);
			if(!filter){
				int32_t result = floor((float)(((val - (avg - Rs)) / (2*Rs)) * L));
				result = result % 2;
				gold[i * out_words_adj + j] = (token_t) result;
			}
			else{
				gold[i * out_words_adj + j] = 3;
			}
		}
}


int main(int argc, char * argv[])
{
	int i;
	int n;
	int ndev;
	struct esp_device *espdevs;
	struct esp_device *dev;
	unsigned done;
	unsigned **ptable;
	token_t *mem;
	token_t *gold;
	unsigned errors = 0;
	unsigned coherence;

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

	ndev = probe(&espdevs, VENDOR_SLD, SLD_BRAIN_BIT, DEV_NAME);
	if (ndev == 0) {
		printf("brain_bit not found\n");
		return 0;
	}

	for (n = 0; n < ndev; n++) {

		printf("**************** %s.%d ****************\n", DEV_NAME, n);

		dev = &espdevs[n];

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
			iowrite32(dev, BRAIN_BIT_AVG_REG, *avg_ptr);
                        /* iowrite32(dev, BRAIN_BIT_AVG_REG, avg_fixed); */
			iowrite32(dev, BRAIN_BIT_KEY_LENGTH_REG, key_length);
                        iowrite32(dev, BRAIN_BIT_STD_REG, *std_ptr);
			/* iowrite32(dev, BRAIN_BIT_STD_REG, std_fixed); */
                        iowrite32(dev, BRAIN_BIT_R_REG, *R_ptr);
			/* iowrite32(dev, BRAIN_BIT_R_REG, R_fixed); */
			iowrite32(dev, BRAIN_BIT_L_REG, L);
			iowrite32(dev, BRAIN_BIT_KEY_BATCH_REG, key_batch);
			iowrite32(dev, BRAIN_BIT_KEY_NUM_REG, key_num);
			iowrite32(dev, BRAIN_BIT_VAL_NUM_REG, val_num);
			iowrite32(dev, BRAIN_BIT_TOT_ITER_REG, tot_iter);

			// Flush (customize coherence model here)
			esp_flush(coherence);

			esp_monitor_args_t mon_args_lsb, mon_args_msb;
			const int ACC_TILE_IDX = 6;
			mon_args_lsb.read_mode = ESP_MON_READ_SINGLE;
			mon_args_lsb.tile_index = ACC_TILE_IDX;
			mon_args_lsb.mon_index = 16;//MON_DVFS_BASE_INDEX + 3;
			mon_args_msb.read_mode = ESP_MON_READ_SINGLE;
			mon_args_msb.tile_index = ACC_TILE_IDX;
			mon_args_msb.mon_index = 17;//MON_DVFS_BASE_INDEX + 3;
			unsigned long int cycles_start_l, cycles_end_l, cycles_diff_l;
			unsigned long int cycles_start_m, cycles_end_m, cycles_diff_m;

			unsigned long long total_time = 0;
			unsigned N_runs = 10;

			for(int k = 0; k < N_runs; k++){
				// Start accelerators
				printf("  Start...\n");

				unsigned long long tmp = 0;
				cycles_start_l = esp_monitor(mon_args_lsb, NULL);
				cycles_start_m = esp_monitor(mon_args_msb, NULL);

				iowrite32(dev, CMD_REG, CMD_MASK_START);

				// Wait for completion
				done = 0;
				while (!done) {
					done = ioread32(dev, STATUS_REG);
					done &= STATUS_MASK_DONE;
				}

				cycles_end_l = esp_monitor(mon_args_lsb, NULL);
				cycles_diff_l = sub_monitor_vals(cycles_start_l, cycles_end_l);
				cycles_end_m = esp_monitor(mon_args_msb, NULL);
				cycles_diff_m = sub_monitor_vals(cycles_start_m, cycles_end_m);
				/* printf("Monitr time is %u %u %u \n", cycles_start, cycles_end, cycles_diff); */

				/* unsigned nano = cycles_diff * 20; //50MHz - vcu707 */
				/* total_time += nano; */
				/* printf("Accelerator runtime: %u ns \n", nano); */

				printf("Accelerator runtime: %u cycles LSB\n", cycles_diff_l);
				printf("Accelerator runtime: %u cycles MSB\n", cycles_diff_m);
				tmp = (cycles_diff_m << 32) + cycles_diff_l;
				//tmp = cycles_diff_l;
				printf("Accelerator total runtime: %llu cycles\n", tmp);
				total_time += tmp;

				iowrite32(dev, CMD_REG, 0x0);

				printf("  Done\n");
				printf("  validating...\n");
			}

			/* Validation */
			errors = validate_buf(&mem[out_offset], gold);

			/* // Start accelerators */
			/* printf("  Start...\n"); */
			/* iowrite32(dev, CMD_REG, CMD_MASK_START); */

			/* // Wait for completion */
			/* done = 0; */
			/* while (!done) { */
			/* 	done = ioread32(dev, STATUS_REG); */
			/* 	done &= STATUS_MASK_DONE; */
			/* } */
			/* iowrite32(dev, CMD_REG, 0x0); */

			/* printf("  Done\n"); */
			/* printf("  validating...\n"); */

			/* /\* Validation *\/ */
			/* errors = validate_buf(&mem[out_offset], gold); */

			total_time = total_time / N_runs;
			printf("Average runtime: %llu cycles \n", total_time);

			float total = 100 * (float) errors / (key_length*key_batch);
			if (total > 1)
				printf("  ... FAIL\n");
			else
				printf("  ... PASS\n");
		}
		aligned_free(ptable);
		aligned_free(mem);
		aligned_free(gold);
	}

	return 0;
}
