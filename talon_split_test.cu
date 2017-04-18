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

#define GRID_BLOCK 16
#define BLOCK_SIZE 32
#define GTHREAD_N ( GRID_BLOCK * BLOCK_SIZE )
#define GIGA 1073741824
#define MP_THREAD 10

int defect_gain_correct(char *fin, char *gain, char *fout, MrcHeader *head,int threads);
int dispatcher_gpu(float *coord_l, float *gain_l , long size_x, long size_y, long slice_n, int type);

__global__ void mutiplier_kernel_type_c_exp(float *coord, float *gain, long long, long long, long long);
void mutiplier_kernel_type_c_exp_test(float *coord, float *gain, long long total_s, long long single_s, long long unit_n);

__global__ void mutiplier_kernel_type_c(float *coord, float *gain, void* src, long size_x, long size_y, long slice_n);
__global__ void mutiplier_kernel_type_su(float *coord, float *gain, void* src,long size_x, long size_y, long slice_n);
__global__ void mutiplier_kernel_type_f(float *coord, float *gain, void* src,long size_x, long size_y, long slice_n);
void mutiplier_kernel_test(float *coord, float *gain, void* src, long size_x, long size_y, long slice_n);

/////talon_split
__global__ void mutiplier_kernel_type_c_exp_split(float *coord, float *gain, long long, long long, long long, long long);
int dispatcher_gpu_split(float *coord_l, float *gain_l , long size_x, long size_y, long slice_n, int type);
/////

int main(int argc, char *argv[])
{
    char *raw_name,*out_name,*gain_name;
    FILE *file,*gain;
    //int threads;
    MrcHeader *inhead;

    raw_name=(char *)malloc(NAME*sizeof(char));
    out_name=(char *)malloc(NAME*sizeof(char));
    gain_name=(char *)malloc(NAME*sizeof(char));

    inhead=(MrcHeader *)malloc(sizeof(MrcHeader));

    if(argc!=4)
    {
        printf("Please input: raw_image gain_name out_name threads\n");
        return;
    }
    raw_name=argv[1];
    gain_name=argv[2];
    out_name=argv[3];
    //threads=atoi(argv[4]);

    file=fopen(raw_name,"rb");	
    mrc_read_head(file,inhead);
    fclose(file);

//    printf("raw_name %s, out_name %s, gain_name %s\n",raw_name,out_name,gain_name);
    defect_gain_correct(raw_name,gain_name,out_name,inhead, MP_THREAD );
	/*
	free(raw_name);
	free(out_name);
	free(gain_name);
	*/
	return 0;
}

int defect_gain_correct(char *fin, char *gain, char *fout, MrcHeader *head,int threads)
	{
	unsigned char lbuf[LINE];
	float fbuf[LINE];
	float *tmp_array,*gain_xy=NULL;
	//unsigned char *coor_xy;
	float *coor_xy=NULL;
	int i,j,x,y,p,k,w,w_row,w_length,m,range,point_s,point_e,n_file,size_bit,*pxy_array;
	long int n;
	FILE *input,*output,*gain_f, *temp_head;
	char *input_byte_c=NULL;short *input_byte_s=NULL;float *input_byte_f=NULL;short *input_byte_u=NULL;
	srand((unsigned) time(NULL));
	//void *source =NULL;

	//open input output and gain file
	printf("Start defect and gain correction|\ninput %s output %s gain %s\n",fin,fout,gain);
	input=fopen(fin,"rb");
	gain_f=fopen(gain,"rb");

	//revise input image header then written for output file
	int size_x=head->nx;
	int size_y=head->ny;
	int slice_n=head->nz;
        int file_type=head->mode;
	head->mode=2;
        switch(file_type)
                {
                case 0:size_bit=1;break;
                case 1:size_bit=2;break;
                case 2:size_bit=4;break;
                case 6:size_bit=2;break;
                default:printf("File type error!");
                }

	//calculate the size of the input and gain file by byte.
	fseek(input,0,SEEK_END);
	long int input_size=ftell(input);
	rewind(input);

	fseek(gain_f,0,SEEK_END);
	long int gain_f_size=ftell(gain_f);
	rewind(gain_f);

	//skip header
	fseek(input,1024,0);
	fseek(gain_f,1024,0);
	//printf("input_size %d file_type %d size_bit %d\n",input_size,file_type,size_bit);

	//malloc memory for the file pointer
	float *gain_byte=(float*)malloc(sizeof(char*)*(gain_f_size-1024));

	//malloc memory for pxy_array
        pxy_array=(int*)malloc(sizeof(int*)*4000);

	//malloc memory for variate 
	tmp_array=(float*)malloc(sizeof(float)*size_x);

	printf("pointer check(pre-allocate): %p -- %p \n", coor_xy, gain_xy);
	gain_xy=(float*)malloc(sizeof(float)*size_y*size_x);
	//coor_xy=(unsigned char*)malloc(sizeof(unsigned char)*size_y*size_x*slice_n);
	coor_xy=(float*)malloc(sizeof(float)*size_y*size_x*slice_n);
	printf("pointer check(post-allocate): %p -- %p \n", coor_xy, gain_xy);

	//read the input and gain file into memory
	//note: the byte size(second parameter) is basen on the type of the file

	switch(file_type)
                {
                case 0:size_bit=1;input_byte_c=(char*)malloc(sizeof(char)*(input_size-1024));fread(input_byte_c,size_bit,(input_size-1024)/size_bit,input);
				//source = input_byte_c;
				break;
                case 1:size_bit=2;input_byte_s=(short*)malloc(sizeof(char)*(input_size-1024));fread(input_byte_s,size_bit,(input_size-1024)/size_bit,input);
				//source = input_byte_s;
				break;
                case 2:size_bit=4;input_byte_f=(float*)malloc(sizeof(char)*(input_size-1024));fread(input_byte_f,size_bit,(input_size-1024)/size_bit,input);
				//source = input_byte_f;
				break;
                case 6:size_bit=2;input_byte_u=(short*)malloc(sizeof(char)*(input_size-1024));fread(input_byte_u,size_bit,(input_size-1024)/size_bit,input);
				//source = input_byte_u;
				break;
                default:printf("File type error!");
                }

	fread(gain_byte,sizeof(float),(gain_f_size-1024)/4,gain_f);

	//read gain 
	float badcut=6.0;
	long int num_gain=0;
	int defect_num=0;
	clock_t start,finish, mid;

	long int total_num=0,gain_num=0,gain_m=0;
	
	start=clock();
	for(n=0; n<size_x*size_y; n++)
		{
		gain_xy[n]=(*(float*)(gain_byte+n))*10.0;
		}
	finish=clock();
	printf("read gain_time %d \n",finish-start);
	start=clock();

	#pragma omp parallel for num_threads(threads) private(gain_num)
	//#pragma omp parallel for schedule(dynamic)
	for(n=0; n<size_x*size_y*slice_n; n++){
			gain_num=n%(size_x*size_y);
            switch(file_type){
                case 0:coor_xy[n]=(*(unsigned char*)(input_byte_c+n));
				break;
                case 1:coor_xy[n]=(*(short*)(input_byte_s+n));
				break;
                case 2:coor_xy[n]=(*(float*)(input_byte_f+n));
				break;
                case 6:coor_xy[n]=(*(short*)(input_byte_u+n));
				break;
                default:printf("File type error!\n");
            }
	}

	mid=clock();
	int indicator=0;
	printf("pointer check(pre-invoke): %p -- %p \n", coor_xy, gain_xy);
	indicator=dispatcher_gpu_split( coor_xy, gain_xy, size_x, size_y, slice_n, file_type );

	//defect correction for points
	finish=clock();
	printf("gain time: omp %ds, GPU %ds, total %ds  \n",(mid-start)/CLOCKS_PER_SEC, (finish-mid)/CLOCKS_PER_SEC, (finish-start)/CLOCKS_PER_SEC);

	fclose(input);
	fclose(gain_f);
	free(gain_byte);
	switch(file_type)
               {
               case 0:free(input_byte_c);break;
               case 1:free(input_byte_s);break;
               case 2:free(input_byte_f);break;
               case 6:free(input_byte_u);break;
               default:printf("File type error!\n");
               }
	free(tmp_array);


	//write output file
	output=fopen(fout,"wb");
	//write head for output file
	fwrite(head,(sizeof(MrcHeader)+head->next),1,output);
	start=clock();
	//for(n=0;n<size_x*size_y*slice_n;n++)
		{
		//fwrite(&coor_xy[n],sizeof(unsigned char),1,output);
		//fwrite(&coor_xy[n],sizeof(float),1,output);
		fwrite(coor_xy,sizeof(float),size_x*size_y*slice_n,output);
		}
		printf("entrice written: %ld\n", size_x*size_y*slice_n);
		printf("head length: %d\n", (sizeof(MrcHeader)+head->next));
	finish=clock();
	printf("write_time %d \n",finish-start);

	start=clock();
	free(coor_xy);
	fclose(output);
	finish=clock();
	printf("free time %d \n",finish-start);
	printf("Defect and gain correction finished!\n");

	printf("En Taro Tassadar!!!\n");
	return 0;
}


////talon_split
/**********************/
int dispatcher_gpu_split(float *coord_l, float *gain_l , long size_x, long size_y, long slice_n, int type){
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
	


			printf("waypoint Serpent\n");

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


			printf("Serpent out\n");



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

////talon_split
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
			mutiplier_kernel_type_c_exp<<< GRID_BLOCK, BLOCK_SIZE >>>( (float*)device_coord, (float*)device_gain, total_s, single_s, unit_n );
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
			printf("X: %d", index*unit_n + i); //POTENTIAL: return a counter through pointer
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
