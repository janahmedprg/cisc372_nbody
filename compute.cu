#include <stdlib.h>
#include <math.h>
#include "vector.h"
#include "config.h"
#include <stdio.h>

//compute: Updates the positions and locations of the objects in the system based on gravity.
//Parameters: None
//Returns: None
//Side Effect: Modifies the hPos and hVel arrays with the new positions and accelerations after 1 INTERVAL
__global__ void compute(double *d_mass, vector3 *d_hPos, vector3 *d_hVel, vector3* d_accels, vector3* d_accels_sum){
	//make an acceleration matrix which is NUMENTITIES squared in size;
	int i = blockIdx.y*blockDim.y + threadIdx.y;
	int j = blockIdx.x*blockDim.x + threadIdx.x;
	int k = threadIdx.z;
	int x = threadIdx.y;
	int y = threadIdx.x;

	if(i>NUMENTITIES){
		return;
	}
	if(j>NUMENTITIES){
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

	__shared__ vector3 distance[BLOCK_SIZE][BLOCK_SIZE];
	// __shared__ vector3 accels_sum_sub[BLOCK_SIZE];
	distance[x][y][k]=d_hPos[i][k]-d_hPos[j][k];
	__syncthreads();
	
	double accelmag, magnitude, magnitude_sq;
	if (i==j) {
		d_accels[i*NUMENTITIES + j][k] = 0;
	}
	else{
		magnitude_sq=distance[x][y][0]*distance[x][y][0]+distance[x][y][1]*distance[x][y][1]+distance[x][y][2]*distance[x][y][2];
		magnitude=sqrt(magnitude_sq);
		accelmag=-1*GRAV_CONSTANT*d_mass[j]/magnitude_sq;
	}
	__syncthreads();
	if(i != j){
		d_accels[i*NUMENTITIES + j][k] = accelmag*distance[x][y][k]/magnitude;
	}
	d_accels_sum[i][k] = 0;
}

__global__ void sumAccels(vector3 *d_accels, vector3 *d_accels_sum){
	int i = blockIdx.y*blockDim.y + threadIdx.y;
	int k = threadIdx.z;
	int x = threadIdx.x;

	int size = ((NUMENTITIES+blockDim.x-1)/blockDim.x) * blockDim.x;
	__shared__ vector3 subRow[blockDim.x];
	int idx;
	for (idx = 0; idx<size; idx += blockDim.x){
		subRow[x][k] += (x + idx)<NUMENTITIES ? d_accels[i * NUMENTITIES + x + idx][k] : 0;
	}
	__syncthreads();

	for(int stride = 1; stride<blockDim.x; stride *= 2){
		int arrIdx = x*stride*2;	
		if(x+ stride<blockDim.x){
			subRow[arrIdx][k] += subRow[arrIdx + stride][k];
		}
		__syncthreads();
	}
	if (x == 0){
		d_accels_sum[i][k]=subRow[0][k];
	}
}

__global__ void updateVelPos(vector3 *d_hPos, vector3 *d_hVel, vector3* d_accels_sum){
	int i = blockIdx.y*blockDim.y + threadIdx.y;
	int k = threadIdx.z;
	d_hVel[i][k]+=d_accels_sum[i][k]*INTERVAL;
	d_hPos[i][k]+=d_hVel[i][k]*INTERVAL;
}