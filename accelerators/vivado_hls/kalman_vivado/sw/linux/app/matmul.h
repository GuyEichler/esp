/*
 * matmul.h
 *
 *  Created on: Aug 15, 2023
 *      Author: geichler
 */

#ifndef MATMUL_H_
#define MATMUL_H_

#include "KalmanFilter.h"

void matmul(float a[SIZE][SIZE], float b[SIZE][SIZE], float c[SIZE][SIZE])
{
	int m = N; //A rows
	int n = N; //A columns, B rows
	int q = N; //B columns

	//int a[SIZE][SIZE], b[SIZE][SIZE], c[SIZE][SIZE];

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < q; j++) {
            c[i][j] = 0;
            for (int k = 0; k < n; k++) {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }

}

template<int DIM>
void matmul(float a[DIM][DIM], float b[DIM][DIM], float c[DIM][DIM])
{
	int m = DIM; //A rows
	int n = DIM; //A columns, B rows
	int q = DIM; //B columns

	//int a[SIZE][SIZE], b[SIZE][SIZE], c[SIZE][SIZE];

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < q; j++) {
            c[i][j] = 0;
            for (int k = 0; k < n; k++) {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }

}

template<int ROWA, int COLA, int COLB>
void matmul(float a[ROWA][COLA], float b[COLA][COLB], float c[ROWA][COLB])
{
	int m = ROWA; //A rows
	int n = COLA; //A columns, B rows
	int q = COLB; //B columns

	//int a[SIZE][SIZE], b[SIZE][SIZE], c[SIZE][SIZE];

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < q; j++) {
            c[i][j] = 0;
            for (int k = 0; k < n; k++) {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }

}



template<int ROWA, int COLA, int ROWB>
void matmul_trans(float a[ROWA][COLA], float b[ROWB][COLA], float c[ROWA][ROWB])
{
	//transpose matrix B
	int m = ROWA; //A rows
	int n = COLA; //A columns, B columns
	int q = ROWB; //B rows

	//int a[SIZE][SIZE], b[SIZE][SIZE], c[SIZE][SIZE];

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < q; j++) {
            c[i][j] = 0;
            for (int k = 0; k < n; k++) {
                c[i][j] += a[i][k] * b[j][k];
            }
        }
    }

}

template<int ROW, int COL>
void matadd(float a[ROW][COL], float b[ROW][COL], float c[ROW][COL])
{
	//transpose matrix B
	int m = ROW; //A rows
	int n = COL; //A columns


	//int a[SIZE][SIZE], b[SIZE][SIZE], c[SIZE][SIZE];

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
        	c[i][j] = a[i][j] + b[i][j];
        }
    }

}

template<int ROW, int COL>
void matsub(float a[ROW][COL], float b[ROW][COL], float c[ROW][COL])
{
	//transpose matrix B
	int m = ROW; //A rows
	int n = COL; //A columns


	//int a[SIZE][SIZE], b[SIZE][SIZE], c[SIZE][SIZE];

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
        	c[i][j] = a[i][j] - b[i][j];
        }
    }

}



#endif /* MATMUL_H_ */
