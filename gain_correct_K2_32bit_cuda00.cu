#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ReadFile.h"
#include "ReadFile.c"
#define LINE 1024
#define NAME 1024
#include "time.h"
//#include "omp.h"
#include "cuda_runtime.h"

#define GRID_BLOCK 32 
#define BLOCK_SIZE 64
#define UNIT_N 1024*1024

int defect_gain_correct(char *fin, char *gain, char *fout, MrcHeader *head,int threads);
int dispatcher_gpu(float *coord_l, float *gain_l , long size_x, long size_y, long slice_n, int type, void* source, long int src_size);

__global__ void mutiplier_kernel_type_c(float *coord, float *gain, void* src, long size_x, long size_y, long slice_n);
__global__ void mutiplier_kernel_type_su(float *coord, float *gain, void* src,long size_x, long size_y, long slice_n);
__global__ void mutiplier_kernel_type_f(float *coord, float *gain, void* src,long size_x, long size_y, long slice_n);
void mutiplier_kernel_test(float *coord, float *gain, void* src, long size_x, long size_y, long slice_n);

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
    defect_gain_correct(raw_name,gain_name,out_name,inhead,0);
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
	float *tmp_array,*gain_xy;
	//unsigned char *coor_xy;
	float *coor_xy;
	int i,j,x,y,p,k,w,w_row,w_length,m,range,point_s,point_e,n_file,size_bit,*pxy_array;
	long int n;
	FILE *input,*output,*gain_f, *temp_head;
	char *input_byte_c;short *input_byte_s;float *input_byte_f;short *input_byte_u;
	srand((unsigned) time(NULL));
	void *source =NULL;

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

	gain_xy=(float*)malloc(sizeof(float)*size_y*size_x);
	//coor_xy=(unsigned char*)malloc(sizeof(unsigned char)*size_y*size_x*slice_n);
	coor_xy=(float*)malloc(sizeof(float)*size_y*size_x*slice_n);


	//read the input and gain file into memory
	//note: the byte size(second parameter) is basen on the type of the file

	switch(file_type)
                {
                case 0:size_bit=1;input_byte_c=(char*)malloc(sizeof(char*)*(input_size-1024));fread(input_byte_c,size_bit,(input_size-1024)/size_bit,input);
				source = input_byte_c;
				break;
                case 1:size_bit=2;input_byte_s=(short*)malloc(sizeof(char*)*(input_size-1024));fread(input_byte_s,size_bit,(input_size-1024)/size_bit,input);
				source = input_byte_s;
				break;
                case 2:size_bit=4;input_byte_f=(float*)malloc(sizeof(char*)*(input_size-1024));fread(input_byte_f,size_bit,(input_size-1024)/size_bit,input);
				source = input_byte_f;
				break;
                case 6:size_bit=2;input_byte_u=(short*)malloc(sizeof(char*)*(input_size-1024));fread(input_byte_u,size_bit,(input_size-1024)/size_bit,input);
				source = input_byte_u;
				break;
                default:printf("File type error!");
                }

	fread(gain_byte,sizeof(float),(gain_f_size-1024)/4,gain_f);

	//read gain 
	float badcut=6.0;
	long int num_gain=0;
	int defect_num=0;
	clock_t start,finish;

	long int total_num=0,gain_num=0,gain_m=0;
	
	start=clock();
	for(n=0; n<size_x*size_y; n++)
		{
		gain_xy[n]=(*(float*)(gain_byte+n))*10.0;
		}
	finish=clock();
	printf("read gain_time %d \n",finish-start);
	start=clock();
	//#pragma omp parallel for num_threads(threads) private(gain_num)
	//#pragma omp parallel for schedule(dynamic)
	/*
		for(n=0; n<size_x*size_y*slice_n; n++)
			{
			gain_num=n%(size_x*size_y);
                        switch(file_type)
                                {
                                case 0:coor_xy[n]=(*(unsigned char*)(input_byte_c+n));coor_xy[n]=(coor_xy[n])*(gain_xy[gain_num]);break;
                                case 1:coor_xy[n]=(*(short*)(input_byte_s+n));coor_xy[n]=(coor_xy[n])*(gain_xy[gain_num]);break;
                                case 2:coor_xy[n]=(*(float*)(input_byte_f+n));coor_xy[n]=(coor_xy[n])*(gain_xy[gain_num]);break;
                                case 6:coor_xy[n]=(*(short*)(input_byte_u+n));coor_xy[n]=(coor_xy[n])*(gain_xy[gain_num]);break;
                                default:printf("File type error!\n");
                                }
			
			}
	*/
	int indicator=0;
	indicator=dispatcher_gpu( coor_xy, gain_xy, size_x, size_y, slice_n, file_type, source, input_size-1024 );

	//defect correction for points
	finish=clock();
	printf("gain time %d s \n",(finish-start)/CLOCKS_PER_SEC);

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

	printf("En Taro Tassaddar!!!\n");
	return 0;
}

/**********************/

int dispatcher_gpu(float *coord_l, float *gain_l , long size_x, long size_y, long slice_n, int type, void* source, long int src_size){
	//set up cuda 
	cudaSetDevice(1);	
	void *device_coord=NULL;
	void *device_gain=NULL;
	void *device_src=NULL;
	//void *temp_coord=NULL;
	//void *temp_gain=NULL;
	
	cudaMalloc( (void **)&device_coord,  sizeof(float)*size_x*size_y*slice_n  );
	cudaMalloc( (void **)&device_gain,  sizeof(float)*size_x*size_y );
	cudaMalloc( (void **)&device_src, sizeof(char*)*src_size );
	//data transfer to device
	cudaMemcpy( device_src, source, sizeof(char*)*src_size, cudaMemcpyHostToDevice );
	cudaMemcpy( device_gain, gain_l, sizeof(float)*size_x*size_y, cudaMemcpyHostToDevice );

	//activate kerneld
	switch(type){
		case 0://c
			mutiplier_kernel_type_c<<< GRID_BLOCK, BLOCK_SIZE >>>( (float*)device_coord, (float*)device_gain, device_src, size_x, size_y, slice_n );
			//mutiplier_kernel_test( coord_l, gain_l, source, size_x, size_y, slice_n );
			break;
		case 1://s
			mutiplier_kernel_type_su<<< GRID_BLOCK, BLOCK_SIZE >>>( (float*)device_coord, (float*)device_gain, device_src, size_x, size_y, slice_n );			
			break;
		case 2://f
			mutiplier_kernel_type_f<<< GRID_BLOCK, BLOCK_SIZE >>>( (float*)device_coord, (float*)device_gain, device_src, size_x, size_y, slice_n );
			break;
		case 6://u
			mutiplier_kernel_type_su<<< GRID_BLOCK, BLOCK_SIZE >>>( (float*)device_coord, (float*)device_gain, device_src, size_x, size_y, slice_n );			
			break;
		default:
			printf("Well well, Houston, we have some problem....\n");
			exit(1);
	}
	//data transfer from device
	cudaMemcpy( coord_l, device_coord, sizeof(float)*size_x*size_y*slice_n, cudaMemcpyDeviceToHost );
	//cudaMemcpy( coord_l, device_coord, size_x*size_y*slice_n, cudaMemcpyDeviceToHost );
	//cudaMemcpy( gain_l, device_gain, size_x*size_y,cudaMemcpyDeviceToHost );
	
	//clean up
	cudaFree( device_coord );
	cudaFree( device_gain );

	return 117;
}

__global__ void mutiplier_kernel_type_c(float *coord, float *gain, void* src, long size_x, long size_y, long slice_n){
	//execute all those mutiplications of type_c data
	long index=0; 
	//long unit_n = size_x*size_y*slice_n/GRID_BLOCK;
	long unit_n = UNIT_N;
	index = blockDim.x*blockIdx.x + threadIdx.x;  //such fun
	long slice_c=0;

	for (int i=0; i<unit_n; i++){
		if (index*unit_n + i >= size_x*size_y*slice_n ){ //boundary check
			return;
		}
		slice_c = (index*unit_n+i)%(size_x*size_y);
		coord[index*unit_n+i] = (*(((unsigned char*)src)+index*unit_n+i))*gain[slice_c];
	}
	//__threadfence()
	return;
}
__global__ void mutiplier_kernel_type_su(float *coord, float *gain, void* src, long size_x, long size_y, long slice_n){//execute all those mutiplications of type_s data
	//execute all those mutiplications of type_c data
	long index=0; 
	//long unit_n = size_x*size_y*slice_n/GRID_BLOCK;
	long unit_n = UNIT_N;
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
	//long unit_n = size_x*size_y*slice_n/GRID_BLOCK;
	long unit_n = UNIT_N;
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
