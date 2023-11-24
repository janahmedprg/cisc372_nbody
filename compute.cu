#include <stdlib.h>
#include <math.h>
#include "vector.h"
#include "config.h"

//compute: Updates the positions and locations of the objects in the system based on gravity.
//Parameters: None
//Returns: None
//Side Effect: Modifies the hPos and hVel arrays with the new positions and accelerations after 1 INTERVAL
__global__ void compute(double *d_mass, vector3 *d_hPos, vector3 *d_hVel){
	//make an acceleration matrix which is NUMENTITIES squared in size;
	int blockRow = blockIdx.y;
	int blockCol = blockIdx.x; 
	int row = threadIdx.y;
	int col = threadIdx.x;
	int k;
	// vector3** accels=(vector3**)malloc(sizeof(vector3*)*NUMENTITIES);
	// accels[row]=&values[row*NUMENTITIES];
	//first compute the pairwise accelerations.  Effect is on the first argument.
	__shared__ vector3 accelsSub[BLOCK_SIZE * BLOCK_SIZE];
	__shared__ vector3 hPosSub[2*BLOCK_SIZE];
	__shared__ double massSub[BLOCK_SIZE];
	__shared__ vector3 accel_sum[BLOCK_SIZE];
	hPosSub[row][0] = d_hPos[blockRow + row][0];
	hPosSub[row][1] = d_hPos[blockRow + row][1];
	hPosSub[row][2] = d_hPos[blockRow + row][2];
	hPosSub[col + BLOCK_SIZE][0] = d_hPos[blockCol + col][0];
	hPosSub[col + BLOCK_SIZE][1] = d_hPos[blockCol + col][1];
	hPosSub[col + BLOCK_SIZE][2] = d_hPos[blockCol + col][2];
	massSub[col] = d_mass[blockCol + col];	
	__syncthreads();

	if (blockRow + row==blockCol + col) {
		FILL_VECTOR(accelsSub[row * BLOCK_SIZE + col],0,0,0);
	}
	else{
		vector3 distance;
		for (k=0;k<3;k++) distance[k]=hPosSub[row][k]-hPosSub[col + BLOCK_SIZE][k];
		double magnitude_sq=distance[0]*distance[0]+distance[1]*distance[1]+distance[2]*distance[2];
		double magnitude=sqrt(magnitude_sq);
		double accelmag=-1*GRAV_CONSTANT*massSub[col]/magnitude_sq;
		FILL_VECTOR(accelsSub[row * BLOCK_SIZE + col], accelmag*distance[0]/magnitude,accelmag*distance[1]/magnitude,accelmag*distance[2]/magnitude);
	}
	__syncthreads();
	//sum up the rows of our matrix to get effect on each entity, then update velocity and position.
	for (k=0;k<3;k++)
		accel_sum[row][k]+=accelsSub[row * BLOCK_SIZE + col][k];
		//compute the new velocity based on the acceleration and time interval
		//compute the new position based on the velocity and time interval
	for (k=0;k<3;k++){
		d_hVel[row][k]+=accel_sum[row][k]*INTERVAL;
		d_hPos[row][k]+=d_hVel[row][k]*INTERVAL;
	}
	// free(accels);
	// free(values);
}