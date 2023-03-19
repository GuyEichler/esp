#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <time.h>
#include <sys/time.h>
//#include "brain_bit_primality.h"

#include "libesp.h"
#include "cfg_prime.h"
#include "input_1mil_full.h"
//#include "input_full.h"

#define COMPOSITE        0
#define PROBABLE_PRIME   1
#define PRINT 0
#define PRINT_COMPOSITE 0

#define KEY_CONCAT_ITER KEY_LENGTH / sizeof(token_t)

static unsigned in_words_adj;
static unsigned out_words_adj;
static unsigned in_len;
static unsigned out_len;
static unsigned in_size;
static unsigned out_size;
static unsigned out_offset;
static unsigned size;
static int key_counter;

struct timeval hw_begin_0, hw_end_0;
struct timeval prime_begin_0, prime_end_0;

unsigned long long elapsed_bb_ns = 0;
long seconds = 0, microseconds = 0;
double elapsed = 0.0;
double total_elapsed_urand = 0.0;
double total_elapsed_rand = 0.0;
double total_elapsed_prime = 0.0;
double total_elapsed = 0.0;

static unsigned key_len_bytes;// = (KEY_LENGTH / 8 );
static unsigned num_rand_words;// = KEY_LENGTH / DATA_BITWIDTH;

static int bb_idx;

static void start_timer(struct timeval *begin)
{
    gettimeofday(begin, 0);
}

static void end_timer(struct timeval *begin, struct timeval *end)
{
    gettimeofday(end, 0);

    seconds = end->tv_sec - begin->tv_sec;
    microseconds = end->tv_usec - begin->tv_usec;
    elapsed = seconds + microseconds*1e-6;
    total_elapsed  += elapsed * 1e9;
}

/* User-defined code */
static int validate_buffer(token_t *out, token_t *gold)
{
	int i;
	int j;
	unsigned errors = 0;

	int skip = 0;
	key_counter = 0;
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
			/* printf("out val is %x for word %d bit %d\n", out[word], word, bit); */
			token_t val = mask & 1;
			token_t gold_val = gold[index];
			unsigned reduce = (ceil((float)skip/key_length));
			if(gold_val != 3){
				if(!(i == key_batch - reduce && (j > skip - 1) )){
					if (gold_val != val){
						if(PRINT)
							printf("Calculated value %x Golden value %d for index %d \n",val, gold_val, (index-skip));
						errors++;
						if(PRINT)
							printf("ERROR\n");
					}
				}
			}
			else{
				if(PRINT)
					printf("SKIPPING\n");
				skip += 1;
			}

			if((index - skip + 1) % (key_length*(key_counter+1)) == 0 && index != 0){
				key_counter++;
				/* printf("\n----------KEY %d DONE----------\n", key_counter); */
				/* printf("\nKEY IS: [ "); */
				/* for(int k = key_length / 32 - 1; k >= 0; k--) */
				/* 	printf("0x%x ", out[word-k]); */
				/* printf("]\n\n"); */
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

                if(val != gold_val){
			printf("Calculated value %x Golden value %x for index %d \n", val, gold_val, index);
			printf("ERROR\n");
			errors += key_length;
		}
                /* if((index + 1) % key_length == 0 && index != 0) */
                /*         val_counter++; */
         }


	return errors;
}


/* User-defined code */
static void init_buffer(token_t *in, token_t * gold)
{
	int i;
	int j;

	int jump = bb_idx * key_batch * key_length;


	for (i = 0; i < key_batch; i++)
		for (j = 0; j < key_length; j++){
                        float val = val_arr[jump + i * in_words_adj + j];
                        in[i * in_words_adj + j] = (token_t) float_to_fixed32(val, 12);
                        //in[i * in_words_adj + j] = (token_t) val;
                        /* printf("Generated value %f\n", fixed32_to_float(in[i * out_words_adj + j] , 12)); */
                }

        for (i = 0; i < key_batch; i++)
		for (j = 0; j < key_length; j++){
			float val = val_arr[jump + i * in_words_adj + j];
                        bool filter = (fabs((float)val - avg) >= Rs);
			int mul = pow(2, d);
			int mod = pow(2, h);
			if(!filter){
				//result can be a negative number here
				int result_alt = floor((float)(((val - avg) / (2*Rs)) * Rs * mul));
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
                                //printf("Generated golden value %d\n", gold[i * out_words_adj + j]);
			}
			else{
				gold[i * out_words_adj + j] = 3;
			}
		}
}


/* User-defined code */
static void init_parameters()
{
	if (DMA_WORD_PER_BEAT(sizeof(token_t)) == 0) {
		in_words_adj = key_length;
		out_words_adj = key_length;
	} else {
		in_words_adj = round_up(key_length, DMA_WORD_PER_BEAT(sizeof(token_t)));
		out_words_adj = round_up(key_length, DMA_WORD_PER_BEAT(sizeof(token_t)));
	}
	in_len = in_words_adj * (key_batch);
	out_len =  out_words_adj * (key_batch);
	in_size = in_len * sizeof(token_t);
	out_size = out_len * sizeof(token_t);
	out_offset = in_len;
	size = (out_offset * sizeof(token_t)) + out_size;
}


int brain_bit_main(token_t *buf)
{
	int errors;

	token_t *gold;
	cfg_000[0].hw_buf = (token_t *) buf;

	//printf("Inside brain_bit buf %p  passed ptr %p \n", buf, cfg_000[0].hw_buf);

	gold = malloc(out_size);

	init_buffer(buf, gold);

	((struct brain_bit_alt_vivado_access*) cfg_000[0].esp_desc)->key_length = key_length;

	printf("\n====== %s ======\n\n", cfg_000[0].devname);
	/* <<--print-params-->> */
	printf("  .avg = %f\n", avg);
	printf("  .key_length = %d\n", key_length);
	printf("  .std = %f\n", std);
	printf("  .R = %f\n", R);
	printf("  .d = %d\n", D);
	printf("  .h = %d\n", H);
	printf("  .key_batch = %d\n", key_batch);
	printf("  .key_num = %d\n", key_num);
	printf("\n  ** START **\n");

	unsigned avg_u = *avg_ptr;
	unsigned std_u = *std_ptr;
	unsigned R_u = *R_ptr;

	((struct brain_bit_alt_vivado_access*) cfg_000[0].esp_desc)->avg = avg_u;
	((struct brain_bit_alt_vivado_access*) cfg_000[0].esp_desc)->std = std_u;
	((struct brain_bit_alt_vivado_access*) cfg_000[0].esp_desc)->R = R_u;

	token_t* out_location = &buf[out_offset];
	for(int i = 0; i < out_len; i++)
		out_location[i] = 0;

	esp_run_no_print(cfg_000, NACC);

	//get accelerator time
	elapsed_bb_ns += cfg_000[0].hw_ns;

	printf("\n  ** DONE **\n");

	errors = validate_buffer(&buf[out_offset], gold);

	free(gold);
	//esp_free(buf);

        float total = 100 * (float) errors / (key_length*key_batch);
	if (total <= 1){
                printf("- Keys generated: %d\n", key_counter);
		printf("- Number of bit errors: %d\n", errors);
		printf("+ Test PASSED\n");
	}
	else{
                printf("+ Test FAILED with %f%% errors\n", total);
	}

	printf("\n====== %s ======\n\n", cfg_000[0].devname);

	//print accelerator runtime
	if(PRINT)
	    printf("Accelerator runtime %llu \n", elapsed_bb_ns);

	return errors;
}

int miller_rabin_pass(mpz_t a, mpz_t n) {
    int i, s, result;
    mpz_t a_to_power, d, n_minus_one;
    mpz_init(n_minus_one);
    mpz_sub_ui(n_minus_one, n, 1);

    s = 0;
    mpz_init_set(d, n_minus_one);
    while (mpz_even_p(d)) {
        mpz_fdiv_q_2exp(d, d, 1);
        s++;
    }

    mpz_init(a_to_power);
    mpz_powm(a_to_power, a, d, n);
    if (mpz_cmp_ui(a_to_power, 1) == 0)  {
        result=PROBABLE_PRIME; goto exit;
    }
    for(i=0; i < s-1; i++) {
        if (mpz_cmp(a_to_power, n_minus_one) == 0) {
            result=PROBABLE_PRIME; goto exit;
        }
        mpz_powm_ui (a_to_power, a_to_power, 2, n);
    }
    if (mpz_cmp(a_to_power, n_minus_one) == 0) {
        result=PROBABLE_PRIME; goto exit;
    }
    result = COMPOSITE;
exit:
    mpz_clear(a_to_power);
    mpz_clear(d);
    mpz_clear(n_minus_one);
    return result;
}

int miller_rabin(mpz_t n, gmp_randstate_t rand_state) {
    mpz_t a;
    int repeat;
    mpz_init(a);
    for(repeat=0; repeat<20; repeat++) {
        do {
            mpz_urandomm(a, rand_state, n);
        } while (mpz_sgn(a) == 0);
        if (miller_rabin_pass(a, n) == COMPOSITE) {
            return COMPOSITE;
        }
    }
    return PROBABLE_PRIME;
}

int run_rand(char* keys, int k)
{
	if(keys == NULL) {
		printf("Error allocating memory for rand() \n");
		return -1;
	}

        start_timer(&hw_begin_0);

        srand(k);
        for(int i = 0; i < num_rand_words; i++){
            keys[i * sizeof(token_t)] = rand();
        }

        end_timer(&hw_begin_0, &hw_end_0);

        for(int i = 0; i < num_rand_words; i++)
	    if(PRINT)
	        printf("%x", (unsigned int) keys[i * sizeof(token_t)]);

	if(PRINT)
		printf("\n");

        /* seconds = hw_end_0.tv_sec - hw_begin_0.tv_sec; */
        /* microseconds = hw_end_0.tv_usec - hw_begin_0.tv_usec; */
        /* elapsed += seconds + microseconds*1e-6; */


	/* avg_elapsed_us = (double) (elapsed * 1e6) / MAX_ITER; */

	/* printf("Randomness source: rand()\n"); */
	/* printf("Number of random bits: %d\n", KEY_LEN); */
	/* printf("Time measured: %u Iterations, avg time %.4f us.\n", MAX_ITER, avg_elapsed_us); */

	return 0;
}

int run_dev_urand(char *keys)
{
    int i;
    FILE *f;

    if(keys == NULL) {
        printf("Error allocating mem \n");
        return -1;
    }

        //start timer
        start_timer(&hw_begin_0);

        f = fopen("/dev/urandom", "r");
        if(f == NULL) {
            printf("Error opening /dev/urandom \n");
            return -1;
        }

        fread((void *) keys, key_len_bytes, 1, f);
        end_timer(&hw_begin_0, &hw_end_0);

	if(PRINT)
            printf("key_len_bytes %u num_rand_words %u \n", key_len_bytes, num_rand_words);

        fclose(f);

        //print rand number is hex
	if(PRINT)
            for(i = 0; i < num_rand_words; i++)
                printf("0x%x  ",  *(unsigned int*) (keys +  i * sizeof(token_t)));

	if(PRINT)
	    printf("\n");

	return 0;
}


//int main(int argc, char* argv[])
int prime_check(unsigned key_length_out, int* total, int idx)
{
    bb_idx = idx;
    key_len_bytes = (key_length_out / 8 );
    num_rand_words = key_length_out / DATA_BITWIDTH;
    key_length = key_length_out;
    elapsed_bb_ns = 0;
    total_elapsed_urand = 0.0;
    total_elapsed_rand = 0.0;
    total_elapsed_prime = 0.0;

    int i;
    int tot_primes_urand = 0, tot_primes_bb = 0, tot_primes_rand = 0;

    mpz_t n;
    gmp_randstate_t rand_state;
    gmp_randinit_default(rand_state);
    gmp_randseed_ui(rand_state, time(NULL));

    /* char key[200]; */
    /* char tmp[200]; */

    int max_digits = 500;
    int key_size = num_rand_words*max_digits +1;
    char* key = malloc(key_size);
    //memset(key, 0, key_size);

    char *str = &key[0];
    /* char *str_tmp = &tmp[0]; */

    token_t *buf;

    if((key_length == 256 && bb_idx < 3) || bb_idx < 1){
    /* if((key_length == 256 && bb_idx < 35) || (key_length == 512 && bb_idx < 17) || (key_length == 1024 && bb_idx < 8)){ */
    printf("starting brain bit alt primailty test \n");

    init_parameters();
    buf = (token_t *) esp_alloc(size);

    brain_bit_main(buf);

    printf("**** Brain Bit alt DONE ****\n");
    token_t* out_location = &buf[out_offset];

    printf("\nChecking primality on brain_bit_alt keys.....\n");

    i = 0;
    for(i = 0; i < KEY_NUM; i++){

	memset(key, 0, key_size);

	int pos = 0;
	for(int j = 0; j < num_rand_words; j++){
		pos += sprintf(&key[pos], "%u", out_location[i * num_rand_words + j]);
	}

	if(PRINT_COMPOSITE){
	    printf("\nbrain_bit alt key is = %s\n", str);
	    printf("words were = \n");
	    for(int j = 0; j < num_rand_words; j++){
		printf("%u ", out_location[i*num_rand_words+j]);
	    }
	    printf("\n");
	}

        //start timer
        start_timer(&prime_begin_0);

        mpz_init_set_str(n, str, 10);
        if(miller_rabin(n, rand_state) == PROBABLE_PRIME) {
	    end_timer(&prime_begin_0, &prime_end_0);
            printf("PRIME Brain_bit alt %s\n", str);
            tot_primes_bb++;
        }
        else{
	    end_timer(&prime_begin_0, &prime_end_0);
	    if(PRINT_COMPOSITE)
	        printf("COMPOSITE Brain_bit_alt\n");
	}


        //Reset timer
	total_elapsed_prime += total_elapsed;
	total_elapsed = 0;

        /* mpz_init_set_str(n, str_tmp, 10); */
        /* if(miller_rabin(n, rand_state) == PROBABLE_PRIME) { */
        /*     printf("PRIME Brain_bit\n"); */
        /*     tot_primes_bb++; */
        /* } */
        /* else */
	/*     if(PRINT) */
	/*         printf("COMPOSITE Brain_bit\n"); */

    }

    total[0] += tot_primes_bb;
    esp_free(buf);
    }

/*    
    printf("********************************************************* \n");    
    printf("Report: %u Prime Numbers out of %u Keys \n", tot_primes_bb, KEY_NUM);
    printf("********************************************************* \n");    
*/ 

    printf("\nChecking primality on /dev/urandom keys.....\n");

    char key_urand[key_len_bytes];

    //key_urand = malloc(key_len_bytes + 100);

    str = &key[0];

    i = 0;
    for(i = 0; i < KEY_NUM ; i++) {

        run_dev_urand(&key_urand[0]);
	memset(key, 0, key_size);

	int pos = 0;
        for (int j = 0; j < num_rand_words; j++)
	    pos += sprintf(&key[pos], "%u", *(unsigned int*) (key_urand + j * sizeof(token_t)));


	if(PRINT_COMPOSITE){
	    printf("\nurandom key is = %s\n", str);
	    printf("words were = \n");
	    for(int j = 0; j < num_rand_words; j++){
		printf("%u ", *(unsigned int*) (key_urand + j * sizeof(token_t)));
	    }
	    printf("\n");
	}

        //Reset timer
	total_elapsed_urand += total_elapsed;
	total_elapsed = 0;

        //start timer
        start_timer(&prime_begin_0);

        mpz_init_set_str(n, str, 10);
        if(miller_rabin(n, rand_state) == PROBABLE_PRIME) {
	    end_timer(&prime_begin_0, &prime_end_0);
            printf("PRIME /dev/urandom %s\n", str);
            tot_primes_urand++;
        }
        else{
	    end_timer(&prime_begin_0, &hw_end_0);
	    if(PRINT_COMPOSITE)
		    printf("COMPOSITE /dev/urandom\n");
	}


        //Reset timer
	total_elapsed_prime += total_elapsed;
	total_elapsed = 0;
    }


    total[1] += tot_primes_urand;

    printf("\nChecking primality on rand() keys.....\n");

    char key_rand[key_len_bytes];


    str = &key[0];

    i = 0;
    for(i = 0; i < KEY_NUM ; i++) {

        run_rand(&key_rand[0], i);
	memset(key, 0, key_size);

	int pos = 0;
        for (int j = 0; j < num_rand_words; j++)
	    pos += sprintf(&key[pos], "%u", *(unsigned int*) (key_rand + j * sizeof(token_t)));


	if(PRINT_COMPOSITE){
	    printf("\nrand() key is = %s\n", str);
	    printf("words were = \n");
	    for(int j = 0; j < num_rand_words; j++){
		printf("%u ", *(unsigned int*) (key_rand + j * sizeof(token_t)));
	    }
	    printf("\n");
	}

        //Reset timer
	total_elapsed_rand += total_elapsed;
	total_elapsed = 0;

        //start timer
        start_timer(&prime_begin_0);

        mpz_init_set_str(n, str, 10);
        if(miller_rabin(n, rand_state) == PROBABLE_PRIME) {
	    end_timer(&prime_begin_0, &prime_end_0);
            printf("PRIME rand() %s\n", str);
            tot_primes_rand++;

	    if(PRINT_COMPOSITE){
		    printf("\nrand() prime is = %s\n", str);
		    printf("words were = \n");
		    for(int j = 0; j < num_rand_words; j++){
			    printf("%u ", *(unsigned int*) (key_rand + j * sizeof(token_t)));
		    }
		    printf("\n");
	    }

        }
        else{
	    end_timer(&prime_begin_0, &prime_end_0);
	    if(PRINT_COMPOSITE)
		    printf("COMPOSITE /dev/urandom\n");
	}

        //Reset timer
	total_elapsed_prime += total_elapsed;
	total_elapsed = 0;
    }

    total[2] += tot_primes_rand;

    printf("*************************************************************************************************** \n");
    printf("Report: %u Prime Numbers from brain_bit_alt out of %u total keys total keys gen time(us) %f\n", tot_primes_bb, KEY_NUM, (float) (elapsed_bb_ns / 1e3));
    printf("****************************************************************************************************** \n");

    printf("************************************************************************************************** \n");
    printf("Report: %u Prime Numbers from /dev/urandom out of %u total keys gen(us) %f\n", tot_primes_urand, KEY_NUM, (float) (total_elapsed_urand / 1e3));
    printf("*************************************************************************************************** \n");

    printf("************************************************************************************************** \n");
    printf("Report: %u Prime Numbers from rand() out of %u total keys gen(us) %f\n", tot_primes_rand, KEY_NUM, (float) (total_elapsed_rand / 1e3));
    printf("*************************************************************************************************** \n");

    printf("*************************************************************************************************** \n");
    printf("Report: Average primality check for %d key size and %d keys time(us) %f\n", key_length, KEY_NUM*3, (float) (total_elapsed_prime / 1e3) / (KEY_NUM*3));
    printf("****************************************************************************************************** \n");

    free(key);
    /* esp_free(buf); */

    return 0;
}

int main(int argc, char* argv[]){

	unsigned key = 256;
	int num = 0;
	int iter = 100;
	int total_prime_bb_256 = 0, total_prime_urand_256 = 0, total_prime_rand_256 = 0;
	int total_prime_bb_512 = 0, total_prime_urand_512 = 0, total_prime_rand_512 = 0;
	int total_prime_bb_1024 = 0, total_prime_urand_1024 = 0, total_prime_rand_1024 = 0;
	int total[3];
	total[0] = 0;
	total[1] = 0;
	total[2] = 0;

	for(int i = 0; i < iter; i++){
		key = 256;
		num = prime_check(key, &total[0], i);

		total_prime_bb_256 += total[0];
		total_prime_urand_256 += total[1];
		total_prime_rand_256 += total[2];
		total[0] = 0;
		total[1] = 0;
		total[2] = 0;

		key = 512;
		num = prime_check(key, &total[0], i);

		total_prime_bb_512 += total[0];
		total_prime_urand_512 += total[1];
		total_prime_rand_512 += total[2];
		total[0] = 0;
		total[1] = 0;
		total[2] = 0;

		key = 1024;
		num = prime_check(key, &total[0], i);

		total_prime_bb_1024 += total[0];
		total_prime_urand_1024 += total[1];
		total_prime_rand_1024 += total[2];
		total[0] = 0;
		total[1] = 0;
		total[2] = 0;
	}

	printf("256-bit numbers:\n");
	printf("Total prime numbers from brain bit alt: %d, average per %d iterations: %f\n", total_prime_bb_256, KEY_NUM, (float) total_prime_bb_256/iter);
	printf("Total prime numbers from /dev/urandom: %d, average per %d iterations: %f\n", total_prime_urand_256, KEY_NUM, (float) total_prime_urand_256/iter);
	printf("Total prime numbers from rand(): %d, average per %d iterations: %f\n", total_prime_rand_256, KEY_NUM, (float) total_prime_rand_256/iter);

	printf("512-bit numbers:\n");
	printf("Total prime numbers from brain bit alt: %d, average per %d iterations: %f\n", total_prime_bb_512, KEY_NUM, (float) total_prime_bb_512/iter);
	printf("Total prime numbers from /dev/urandom: %d, average per %d iterations: %f\n", total_prime_urand_512, KEY_NUM, (float) total_prime_urand_512/iter);
	printf("Total prime numbers from rand(): %d, average per %d iterations: %f\n", total_prime_rand_512, KEY_NUM, (float) total_prime_rand_512/iter);

	printf("1024-bit numbers:\n");
	printf("Total prime numbers from brain bit alt: %d, average per %d iterations: %f\n", total_prime_bb_1024, KEY_NUM, (float) total_prime_bb_1024/iter);
	printf("Total prime numbers from /dev/urandom: %d, average per %d iterations: %f\n", total_prime_urand_1024, KEY_NUM, (float) total_prime_urand_1024/iter);
	printf("Total prime numbers from rand(): %d, average per %d iterations: %f\n", total_prime_rand_1024, KEY_NUM, (float) total_prime_rand_1024/iter);

return num;

}
