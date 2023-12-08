#include <stdlib.h>
#include <math.h>
#include "vector.h"
#include "config.h"
#include <stdio.h>

//compute: Updates the positions and locations of the objects in the system based on gravity.
//Parameters: None
//Returns: None
//Side Effect: Modifies the hPos and hVel arrays with the new positions and accelerations after 1 INTERVAL
__global__ void compute(double *d_mass, vector3 *d_hPos, vector3 *d_hVel, vector3* d_accels, int d_numObjects, vector3* d_accels_sum){
	//make an acceleration matrix which is NUMENTITIES squared in size;
	int blockRow = blockIdx.y;
	int blockCol = blockIdx.x;
	int row = threadIdx.y;
	int col = threadIdx.x;
	int k = threadIdx.z;

	// __shared__ vector3 accelsSub[BLOCK_SIZE * BLOCK_SIZE];
	// __shared__ vector3 hPosSub[2*BLOCK_SIZE];
	// __shared__ double massSub[BLOCK_SIZE];
	// __shared__ vector3 accel_sum[BLOCK_SIZE];
	// hPosSub[row][k] = d_hPos[blockRow + row][k];
	// hPosSub[col + BLOCK_SIZE][k] = d_hPos[blockCol + col][k];
	// massSub[col] = d_mass[blockCol + col];	
	// __syncthreads();

	if (blockRow + row==blockCol + col) {
		d_accels[(blockRow + row)*d_numObjects + blockCol + col][k] = 0;
	}
	else{
		vector3 distance;
		distance[k]=d_hPos[blockRow + row][k]-d_hPos[blockCol + col][k];
		double magnitude_sq=distance[0]*distance[0]+distance[1]*distance[1]+distance[2]*distance[2];
		double magnitude=sqrt(magnitude_sq);
		double accelmag=-1*GRAV_CONSTANT*d_mass[blockCol + col]/magnitude_sq;
		d_accels[(blockRow + row)*d_numObjects + blockCol + col][k] = accelmag*distance[k]/magnitude;
	}
	d_accels_sum[blockRow + row][k] = 0;
	__syncthreads();

	d_accels_sum[blockRow + row][k]+=d_accels[(blockRow + row)*d_numObjects + blockCol + col][k];
	__syncthreads();
	d_hVel[blockRow + row][k]+=d_accels_sum[blockRow + row][k]*INTERVAL;
	d_hPos[blockRow + row][k]+=d_hVel[blockRow + row][k]*INTERVAL;
}

// __global__ void sumAccels(vector3 *d_hPos, vector3 *d_hVel, vector3 **d_accels, vector3 *d_accels_sum){
// 	int blockRow = blockIdx.y;
// 	int blockCol = blockIdx.x;
// 	int row = threadIdx.y;
// 	int col = threadIdx.x;
// 	int k = threadIdx.z;

// 	d_accels_sum[blockRow + row][k] = 0;
// 	__syncthreads();

// 	d_accels_sum[blockRow + row][k]+=d_accels[blockRow+row][blockCol + col][k];
// 	__syncthreads();
// 	d_hVel[blockRow + row][k]+=d_accels_sum[blockRow + row][k]*INTERVAL;
// 	d_hPos[blockRow + row][k]+=d_hVel[blockRow + row][k]*INTERVAL;
// }
