#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ReadFile.h"
#include "ReadFile.c"
#define LINE 1024
#define NAME 1024
#include "time.h"
#include "omp.h"
#include "cuda_runtime.h"

#define GRID_BLOCK 1
#define BLOCK_SIZE 1
#define GTHREAD_N ( GRID_BLOCK * BLOCK_SIZE )
#define GIGA 1073741824
#define MP_THREAD 10

int defect_gain_correct(char *fin, char *gain, char *fout, MrcHeader *head,int threads);
int dispatcher_gpu(float *coord_l, float *gain_l , long size_x, long size_y, long slice_n, int type);

__global__ void mutiplier_kernel_type_c_exp(float *coord, float *gain, long long, long long, long long);
__global__ void mutiplier_kernel_type_c_exp_mini_2(float *coord, float *gain, long long, long long, long long);
void mutiplier_kernel_type_c_exp_test(float *coord, float *gain, long long total_s, long long single_s, long long unit_n);


int main(int argc, char *argv[])
{

	//printf("raw_name %s, out_name %s, gain_name %s\n",raw_name,out_name,gain_name);
	char *raw_name = NULL, *gain_name = NULL, *out_name = NULL;
	MrcHeader *inhead=NULL;
    defect_gain_correct(raw_name,gain_name,out_name,inhead, MP_THREAD );
	/*
	free(raw_name);
	free(out_name);
	free(gain_name);
	*/
	return 0;
}

int defect_gain_correct(char *fin, char *gain, char *fout, MrcHeader *head,int threads){
	
	//revise input image header then written for output file
	int size_x=5000;
	int size_y=2000;
	int slice_n=2;
    int indicator=0;
	int result=0;
	int sum=0;
	float *coor_xy = NULL;
	float *gain_xy = NULL;

	coor_xy = (float *)malloc( sizeof(int)*size_x*size_y*slice_n );
	gain_xy = (float *)malloc( sizeof(int)*size_x*size_y );

	for (int j=0; j < size_x*size_y; j++ ){
		gain_xy[j] = 2;
	}
	for (int k=0; k < size_x*size_y*slice_n; k++ ){
		if ( k < size_x*size_y  ) {
			coor_xy[k] = 1;
		}
		else {
			coor_xy[k] = 2;
		}
	}

	indicator=dispatcher_gpu( coor_xy, gain_xy, size_x, size_y, slice_n, 0 );

	for(int i=0; i<size_x*size_y*slice_n; i++){
		if ( i < size_x*size_y && coor_xy[i] != 2 ) {
			result+=1;
		}
		else if( i > size_x*size_y && coor_xy[i] != 4 ) {
			result+=1;
		}
		else;
		sum += coor_xy[i];
	}

	printf("result: %ld, %d\n", result, sum );
	return 0;
}

/**********************/

int dispatcher_gpu(float *coord_l, float *gain_l , long size_x, long size_y, long slice_n, int type){
	//set up cuda 
	cudaSetDevice(0);	
	void *device_coord=NULL; //potential: change all to double
	void *device_gain=NULL;
	//long *d_x=NULL, *d_y=NULL,*d_n=NULL;
	//void *temp_coord=NULL;
	//void *temp_gain=NULL;
	cudaError_t f1, f2;

	long long total_s=0, single_s=0, unit_n=0;
	single_s = size_x * size_y;
	total_s = single_s * slice_n;

	printf("pointer check(host): %p -- %p \n", coord_l, gain_l);
	printf("pointer check(device): %p -- %p \n", device_coord, device_gain);

	printf("Initialize allocation.\n");
	printf("Required total(GB): %ld \n", (sizeof(float)*total_s + sizeof(float)*single_s)/GIGA);
	printf("Required each: %ld -- %ld \n", sizeof(float)*total_s, sizeof(float)*single_s);
	f1 = cudaMalloc( (void **)&device_coord,  sizeof(float)*total_s  );
	f2 = cudaMalloc( (void **)&device_gain,  sizeof(float)*single_s );
	//f3 = cudaMalloc( (void **)&device_src, sizeof(char*)*src_size );
	//GRAM allocation check
	if ( f1 != cudaSuccess || f2 != cudaSuccess ){
		printf("cuda memory allocation failed: %s -- %s \n", f1, f2);
		exit(4);
	}
	
	/*
	cudaMalloc( (void **)&device_coord,  sizeof(float)*size_x*size_y*slice_n  );
	cudaMalloc( (void **)&device_gain,  sizeof(float)*size_x*size_y );
	cudaMalloc( (void **)&device_src, sizeof(char*)*src_size );
	*/
	printf("Allocation done.\n");
	printf("status: %s -- %s \n", f1, f2);

	printf("pointer check(host): %p -- %p \n", coord_l, gain_l);
	printf("pointer check(device): %p -- %p \n", device_coord, device_gain);
	/*
	cudaMalloc( (void **)&d_x, sizeof(long) );
	cudaMalloc( (void **)&d_y, sizeof(long) );
	cudaMalloc( (void **)&d_n, sizeof(long) );
	*/
	//data transfer to device
	printf("Initialize Memcpy: host->device\n");
	
	///* Kernel validation
	f1 = cudaMemcpy( device_coord, coord_l, sizeof(float)*total_s, cudaMemcpyHostToDevice );
	f2 = cudaMemcpy( device_gain, gain_l, sizeof(float)*single_s, cudaMemcpyHostToDevice );
	//GRAM allocation check
	if ( f1 != cudaSuccess || f2 != cudaSuccess ){
		printf("cudaMemcpy failed(H->D): %s -- %s \n", f1, f2);
		exit(5);
	}
	//*/
	printf("Memcpy Done: host->device\n");

	printf("pointer check(host): %p -- %p \n", coord_l, gain_l);
	printf("pointer check(device): %p -- %p \n", device_coord, device_gain);

	/*
	cudaMemcpy( d_x, &size_x, sizeof(long), cudaMemcpyHostToDevice);
	cudaMemcpy( d_y, &size_y, sizeof(long), cudaMemcpyHostToDevice);
	cudaMemcpy( d_n, &slice_n, sizeof(long), cudaMemcpyHostToDevice);
	*/
	//activate kernel
	unit_n = (total_s/( GTHREAD_N )) + 1; 
	switch(type){
		case 0://c
			printf("waypoint c\n");
			//mutiplier_kernel_type_c_exp<<< GRID_BLOCK, BLOCK_SIZE >>>( (float*)device_coord, (float*)device_gain, total_s, single_s, unit_n );
			mutiplier_kernel_type_c_exp_mini_2<<< GRID_BLOCK, BLOCK_SIZE >>>( (float*)device_coord, (float*)device_gain, total_s, single_s, unit_n );
			/*
			unit_n = total_s;
			mutiplier_kernel_type_c_exp<<< 1, 1 >>>( (float*)device_coord, (float*)device_gain, total_s, single_s, unit_n );
			*/
			////Kernel validation
			/*
			unit_n = total_s;
			mutiplier_kernel_type_c_exp_test( coord_l, gain_l, total_s, single_s, unit_n);
			*/
			printf("C out\n");
			break;
		case 1://s
			//mutiplier_kernel_type_su<<< GRID_BLOCK, BLOCK_SIZE >>>( (float*)device_coord, (float*)device_gain, device_src, size_x, size_y, slice_n );			
			break;
		case 2://f
			//mutiplier_kernel_type_f<<< GRID_BLOCK, BLOCK_SIZE >>>( (float*)device_coord, (float*)device_gain, device_src, size_x, size_y, slice_n );
			break;
		case 6://u
			//mutiplier_kernel_type_su<<< GRID_BLOCK, BLOCK_SIZE >>>( (float*)device_coord, (float*)device_gain, device_src, size_x, size_y, slice_n );			
			break;
		default:
			printf("Well...Houston, we've got some problem...\n");
			exit(1);
	}
	//thread synchronization 
	f1 = cudaThreadSynchronize();
	if ( f1 != cudaSuccess ){
		printf("cuda sync failed: %s \n", f1 );
		exit(9);
	}
	printf("Sync Done\n");

	printf("pointer check(host): %p -- %p \n", coord_l, gain_l);
	printf("pointer check(device): %p -- %p \n", device_coord, device_gain);

	//data transfer from device
	///* Kernel validation
	printf("Initialize Memcpy: device->host\n");
	
	f1 = cudaMemcpy( coord_l, device_coord, sizeof(float)*total_s, cudaMemcpyDeviceToHost );

	//GRAM allocation check
	if ( f1 != cudaSuccess ){
		printf("cudaMemcpy failed(D->H): %s \n", f1 );
		exit(6);
	}
	printf("Memcpy Done: device->host\n");

	printf("pointer check(host): %p -- %p \n", coord_l, gain_l);
	printf("pointer check(device): %p -- %p \n", device_coord, device_gain);
	//*/
	//cudaMemcpy( coord_l, device_coord, size_x*size_y*slice_n, cudaMemcpyDeviceToHost );
	//cudaMemcpy( gain_l, device_gain, size_x*size_y,cudaMemcpyDeviceToHost );
	
	//clean up
	printf("Initialize cudaFree\n");
	f1 = cudaFree( device_coord );
	f2 = cudaFree( device_gain );
	if ( f1 != cudaSuccess || f2 != cudaSuccess ){
		printf("cudaFree failed: %s -- %s \n", f1, f2);
		exit(7);
	}
	printf("cudaFree Done\n");
	printf("pointer check(host): %p -- %p \n", coord_l, gain_l);
	printf("pointer check(device): %p -- %p \n", device_coord, device_gain);

	//cudaFree( device_src );
	/*
	cudaFree( d_x );
	cudaFree( d_y );
	cudaFree( d_n );
	*/
	printf("Dispatcher out\n");
	return 117;
}

__global__ void mutiplier_kernel_type_c_exp(float *coord, float *gain, long long total_s, long long single_s, long long unit_n){
	//execute all those mutiplications of type_c data
	long long index=0;
	long long i = 0; 
	//long long unit_n = (size_x*size_y*slice_n/( GTHREAD_N ))+1;
	//long long unit_n = unit_s;
	//long unit_n = UNIT_N;
	//long unit_n = (size_x*size_y*slice_n);
	index = blockDim.x*blockIdx.x + threadIdx.x;  //such fun
	long long slice_c=0;

	printf("index check: %d\n ", index);
	//printf("Size check: %d \n", size_x*size_y*slice_n);
	printf("Size check: %d \n", total_s);
	printf("unit_n check: %d \n", unit_n);
	printf("input check: %p -- %p \n", coord, gain);

	for (i=0; i<unit_n; i++){
		//printf("head-> %d\n", i);
		if (index*unit_n + i >= total_s ){ //boundary check
			printf("X: %d\n", index*unit_n + i); //POTENTIAL: return a counter through pointer
			return;
		}
		//printf("progress: %d \n", i);
		//slice_c = (index*unit_n+i)%(size_x*size_y);
		slice_c = (index*unit_n+i)%(single_s);
		coord[index*unit_n+i] = (coord[index*unit_n+i] * gain[slice_c]);
		//printf("tail-> %d\n", i);
	}
	//__threadfence()
	printf("Kernel out: %d\n", i);
	return;
}

__global__ void mutiplier_kernel_type_c_exp_mini_2(float *coord, float *gain, long long total_s, long long single_s, long long unit_n){
	//execute all those mutiplications of type_c data
	long long index=0;
	long long i = 0; 
	index = blockDim.x*blockIdx.x + threadIdx.x;  //such fun
	long long slice_c=0;
	
	long long index_1=0;

	printf("index check: %d\n ", index);
	printf("Size check: %d \n", total_s);
	printf("unit_n check: %d \n", unit_n);
	printf("input check: %p -- %p \n", coord, gain);

	index_1 = index*unit_n;
	printf("Local iteration start\n");
	for (i=0; i<unit_n; i++){
		
		if (index_1 >= total_s ){ //boundary check
			printf("Y: %d\n", index_1); //POTENTIAL: return a counter through pointer
			return;
		}
		
		slice_c = (index_1)%(single_s);
		//coord[index_1] = (coord[index_1] * gain[slice_c]);
		index_1+=1;
	}
	
	printf("Kernel_mini out: %d\n", i);
	return;
}

void mutiplier_kernel_type_c_exp_test(float *coord, float *gain, long long total_s, long long single_s, long long unit_n){
	//execute all those mutiplications of type_c data
	long long index=0; 
	//long long unit_n = (size_x*size_y*slice_n/( GTHREAD_N ))+1;
	//long long unit_n = unit_s;
	//long unit_n = UNIT_N;
	//long unit_n = (size_x*size_y*slice_n);
	//index = blockDim.x*blockIdx.x + threadIdx.x;  //such fun
	long slice_c=0;

	printf("index check: %d\n ", index);
	//printf("Size check: %d \n", size_x*size_y*slice_n);
	printf("Size check: %d \n", total_s);
	printf("unit_n check: %d \n", unit_n);

	for (int i=0; i<unit_n; i++){
		//printf("head-> %d\n", i);
		if (index*unit_n + i >= total_s ){ //boundary check
			printf("X: %d", index*unit_n + i);
			return;
		}
		//printf("progress: %d \n", i);
		//slice_c = (index*unit_n+i)%(size_x*size_y);
		slice_c = (index*unit_n+i)%(single_s);
		coord[index*unit_n+i] = (coord[index*unit_n+i] * gain[slice_c]);
		//printf("tail-> %d\n", i);
	}
	//__threadfence()
	printf("Kernel out\n");
	return;
}

