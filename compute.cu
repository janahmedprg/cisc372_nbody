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
	int i = blockIdx.y*blockDim.y + threadIdx.y;
	int j = blockIdx.x*blockDim.x + threadIdx.x;
	int k = threadIdx.z;
	int x = threadIdx.x;
	int y = threadIdx.y;

	if(i>d_numObjects){
		return;
	}
	if(j>d_numObjects){
		return;
	}

	// __shared__ vector3 accelsSub[BLOCK_SIZE * BLOCK_SIZE];
	// __shared__ vector3 hPosSub[2*BLOCK_SIZE];
	// __shared__ double massSub[BLOCK_SIZE];
	// __shared__ vector3 accel_sum[BLOCK_SIZE];
	// hPosSub[row][k] = d_hPos[blockRow + row][k];
	// hPosSub[col + BLOCK_SIZE][k] = d_hPos[blockCol + col][k];
	// massSub[col] = d_mass[blockCol + col];	
	// __syncthreads();

	__shared__ vector3 distance[16][16];
	distance[x][y][k]=d_hPos[i][k]-d_hPos[j][k];
	__syncthreads();
	
	if (i==j) {
		d_accels[i*d_numObjects + j][k] = 0;
	}
	else{
		double magnitude_sq=distance[x][y][0]*distance[x][y][0]+distance[x][y][1]*distance[x][y][1]+distance[x][y][2]*distance[x][y][2];
		double magnitude=sqrt(magnitude_sq);
		double accelmag=-1*GRAV_CONSTANT*d_mass[j]/magnitude_sq;
		d_accels[i*d_numObjects + j][k] = accelmag*distance[x][y][k]/magnitude;
	}
	d_accels_sum[i][k] = 0;
	__syncthreads();

	d_accels_sum[i][k]+=d_accels[i*d_numObjects + j][k];
	__syncthreads();
	d_hVel[i][k]+=d_accels_sum[i][k]*INTERVAL;
	d_hPos[i][k]+=d_hVel[i][k]*INTERVAL;
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
