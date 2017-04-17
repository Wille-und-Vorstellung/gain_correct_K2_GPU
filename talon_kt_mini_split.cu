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

__global__ void mutiplier_kernel_type_c_exp_split(float *coord, float *gain, long long, long long, long long, long long);
void mutiplier_kernel_type_c_exp_test(float *coord, float *gain, long long total_s, long long single_s, long long unit_n);

__global__ void mutiplier_kernel_type_c(float *coord, float *gain, void* src, long size_x, long size_y, long slice_n);
__global__ void mutiplier_kernel_type_su(float *coord, float *gain, void* src,long size_x, long size_y, long slice_n);
__global__ void mutiplier_kernel_type_f(float *coord, float *gain, void* src,long size_x, long size_y, long slice_n);
void mutiplier_kernel_test(float *coord, float *gain, void* src, long size_x, long size_y, long slice_n);

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
	int size_y=200;
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
	void *device_coord_1=NULL, *device_coord_2=NULL;
	void *device_gain_1=NULL, *device_gain_2=NULL;
	int left_over = 0;
	cudaError_t f1, f2, f3, f4;
	long long total_s=0, single_s=0, unit_n=0;
	long long offset = 0; 

	single_s = size_x * size_y;
	total_s = single_s * slice_n;
	if ( total_s%2 != 0 ){
		left_over = 1;
	}
	offset = total_s/2;
	/*
	printf("pointer check(host): %p -- %p \n", coord_l, gain_l);
	printf("pointer check(device): %p -- %p \n", device_coord, device_gain);
	*/
	printf("Initialize allocation.\n");
	printf("Required total(GB): %ld \n", (sizeof(float)*total_s + sizeof(float)*single_s)/GIGA);
	printf("Required each: %ld -- %ld \n", sizeof(float)*total_s, sizeof(float)*single_s);
	f1 = cudaMalloc( (void **)&device_coord_1,  sizeof(float)*total_s/2  );
	f2 = cudaMalloc( (void **)&device_gain_1,  sizeof(float)*single_s );
	f3 = cudaMalloc( (void **)&device_coord_2,  sizeof(float)*(total_s/2 + left_over)  );
	//f4 = cudaMalloc( (void **)&device_gain_2,  sizeof(float)*single_s );
	//GRAM allocation check
	if ( f1 != cudaSuccess || f2 != cudaSuccess || f3 != cudaSuccess ){
		printf("cuda memory allocation failed: %s -- %s -- %s \n", f1, f2, f3 );
		exit(4);
	}
	printf("Allocation done.\n");
	printf("status: %s -- %s -- %s \n", f1, f2, f3);
	/*
	printf("pointer check(host): %p -- %p \n", coord_l, gain_l);
	printf("pointer check(device): %p -- %p \n", device_coord, device_gain);
	*/
	//data transfer to device
	printf("Initialize Memcpy: host->device\n");
	
	///* Kernel validation
	f1 = cudaMemcpy( device_coord_1, coord_l, sizeof(float)*total_s/2, cudaMemcpyHostToDevice );
	f2 = cudaMemcpy( device_gain_1, gain_l, sizeof(float)*single_s, cudaMemcpyHostToDevice );
	f3 = cudaMemcpy( device_coord_2, (coord_l + total_s/2), sizeof(float)*(total_s/2 + left_over), cudaMemcpyHostToDevice );
	//GRAM allocation check
	if ( f1 != cudaSuccess || f2 != cudaSuccess ){
		printf("cudaMemcpy failed(H->D): %s -- %s \n", f1, f2);
		exit(5);
	}
	//*/
	printf("Memcpy Done: host->device\n");
	/*
	printf("pointer check(host): %p -- %p \n", coord_l, gain_l);
	printf("pointer check(device): %p -- %p \n", device_coord, device_gain);
	*/
	
	//activate kernel
	unit_n = (total_s/(2*( GTHREAD_N ))) + 1; 
	


			printf("waypoint C\n");

			mutiplier_kernel_type_c_exp_split<<< GRID_BLOCK, BLOCK_SIZE >>>( (float*)device_coord_1, (float*)device_gain_1, total_s/2, single_s, unit_n, 0 ); 
			//HAZARD: slice_c might be a problem... + a offset, yes offset will do.BUT NOT FULLY VERIFIED YET
			//can switch cpu here and remove the sync below
			///*
			f1 = cudaThreadSynchronize();
			if ( f1 != cudaSuccess ){
				printf("cuda sync(mid-way) failed: %s \n", f1 );
				exit(12);
			}
			//*/
			mutiplier_kernel_type_c_exp_split<<< GRID_BLOCK, BLOCK_SIZE >>>( (float*)device_coord_2, (float*)device_gain_1, total_s/2+left_over, single_s, unit_n, offset );


			printf("C out\n");



	//thread synchronization 
	f1 = cudaThreadSynchronize();
	if ( f1 != cudaSuccess ){
		printf("cuda sync failed: %s \n", f1 );
		exit(9);
	}
	printf("Sync Done\n");
	/*
	printf("pointer check(host): %p -- %p \n", coord_l, gain_l);
	printf("pointer check(device): %p -- %p \n", device_coord, device_gain);
	*/
	//data transfer from device
	///* Kernel validation
	printf("Initialize Memcpy: device->host\n");
	
	f1 = cudaMemcpy( coord_l, device_coord_1, sizeof(float)*total_s/2, cudaMemcpyDeviceToHost );
	printf("Memcpy: device->host, half-way through\n");
	////Narrator: replace cudaMemcpy manually,if cudaMemcpy doesn't work for too long addresses
	f2 = cudaMemcpy( (coord_l + total_s/2), device_coord_2, sizeof(float)*(total_s/2 + left_over), cudaMemcpyDeviceToHost );
	/*
	for (long long i=0; i < sizeof(float)*total_s/2 + left_over; i++){
		coord_l[i+total_s/2] = device_coord_2[i];
	}
	*/
	//Narrator: if manual copy doesn't work either, then we declare another agent and do a manual cascade.

	//GRAM allocation check
	if ( f1 != cudaSuccess || f2 != cudaSuccess ){
		printf("cudaMemcpy failed(D->H): %s -- %s \n", f1, f2 );
		exit(6);
	}
	printf("Memcpy Done: device->host\n");
	/*
	printf("pointer check(host): %p -- %p \n", coord_l, gain_l);
	printf("pointer check(device): %p -- %p \n", device_coord, device_gain);
	*/
	
	//clean up
	printf("Initialize cudaFree\n");
	f1 = cudaFree( device_coord_1 );
	f2 = cudaFree( device_gain_1 );
	f3 = cudaFree( device_coord_2 );
	if ( f1 != cudaSuccess || f2 != cudaSuccess ){
		printf("cudaFree failed: %s -- %s \n", f1, f2);
		exit(7);
	}
	printf("cudaFree Done\n");
	/*
	printf("pointer check(host): %p -- %p \n", coord_l, gain_l);
	printf("pointer check(device): %p -- %p \n", device_coord, device_gain);
	*/
	
	printf("Dispatcher out\n");
	return 117;
}

__global__ void mutiplier_kernel_type_c_exp_split(float *coord, float *gain, long long total_s, long long single_s, long long unit_n, long long offset){
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
		slice_c = (index*unit_n+i+offset)%(single_s);
		coord[index*unit_n+i] = (coord[index*unit_n+i] * gain[slice_c]);
		//printf("tail-> %d\n", i);
	}
	//__threadfence()
	printf("Kernel out: %d\n", i);
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

__global__ void mutiplier_kernel_type_c(float *coord, float *gain, void* src, long size_x, long size_y, long slice_n){
	//execute all those mutiplications of type_c data
	long index=0; 
	//long unit_n = (size_x*size_y*slice_n/(32*64))+1;
	//long unit_n = UNIT_N;
	long unit_n = (size_x*size_y*slice_n);
	index = blockDim.x*blockIdx.x + threadIdx.x;  //such fun
	long slice_c=0;

	printf("index %d\n ", index);
	printf("Size check: %d \n", size_x*size_y*slice_n);
	printf("unit_n check: %d \n", unit_n);

	for (int i=0; i<unit_n; i++){
		printf("head %d\n", i);
		if (index*unit_n + i >= size_x*size_y*slice_n ){ //boundary check
			printf("X: %d", index*unit_n + i);
			return;
		}
		printf("progress: %d \n", i);
		//slice_c = (index*unit_n+i)%(size_x*size_y);
		coord[index*unit_n+i] = (float)(*(unsigned char*)(((unsigned char*)src)+index*unit_n+i)) * gain[slice_c];
		printf("tail %d\n", i);
	}
	//__threadfence()
	printf("Kernel out\n");
	return;
}
__global__ void mutiplier_kernel_type_su(float *coord, float *gain, void* src, long size_x, long size_y, long slice_n){//execute all those mutiplications of type_s data
	//execute all those mutiplications of type_c data
	long index=0; 
	long unit_n = size_x*size_y*slice_n/GRID_BLOCK;
	//long unit_n = UNIT_N;
	index = blockDim.x*blockIdx.x + threadIdx.x;  //such fun
	long slice_c=0;

	for (int i=0; i<unit_n; i++){
		if (index*unit_n + i >= size_x*size_y*slice_n ){ //boundary check
			return;
		}
		slice_c = (index*unit_n+i)%(size_x*size_y);
		coord[index*unit_n+i] = (*(((short*)src)+index*unit_n+i))*gain[slice_c];
	}
	//__threadfence()
	return;
}
__global__ void mutiplier_kernel_type_f(float *coord, float *gain, void* src, long size_x, long size_y, long slice_n){//execute all those mutiplications of type_f data
	//execute all those mutiplications of type_c data
	long index=0; 
	long unit_n = size_x*size_y*slice_n/GRID_BLOCK;
	//long unit_n = UNIT_N;
	index = blockDim.x*blockIdx.x + threadIdx.x;  //such fun
	long slice_c=0;

	for (int i=0; i<unit_n; i++){
		if (index*unit_n + i >= size_x*size_y*slice_n ){ //boundary check
			return;
		}
		slice_c = (index*unit_n+i)%(size_x*size_y);
		coord[index*unit_n+i] = (*(((float*)src)+index*unit_n+i))*gain[slice_c];
	}
	//__threadfence()
	return;
}

void mutiplier_kernel_test(float *coord, float *gain, void* src, long size_x, long size_y, long slice_n){
	///////////////////////////////
	printf("Activating kernel test\n");
	//execute all those mutiplications of type_c data
	long index=0; 
	//long unit_n = size_x*size_y*slice_n/GRID_BLOCK;
	//long unit_n = UNIT_N;
	long unit_n = size_x*size_y*slice_n;
	//index = blockDim.x*blockIdx.x + threadIdx.x;  //such fun
	long slice_c=0;

	for (int i=0; i<unit_n; i++){
		
		if (index*unit_n + i >= size_x*size_y*slice_n ){ //boundary check
			return;
		}
		
		slice_c = (index*unit_n+i)%(size_x*size_y);
		//slice_c = (i)%(size_x*size_y);
		coord[index*unit_n+i] = (*(((unsigned char*)src)+index*unit_n+i))*gain[slice_c];
		//coord[i] = (*(((unsigned char*)src)+i))*gain[slice_c];
	}
	//__threadfence()
	printf("Kernel test done\n");
	return;
}
