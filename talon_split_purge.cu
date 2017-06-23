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

#define GRID_BLOCK 16 //16
#define BLOCK_SIZE 32 //32
#define GTHREAD_N ( GRID_BLOCK * BLOCK_SIZE )
#define GIGA 1073741824
#define MP_THREAD 10
#define TEST_RUN false

#define SIZE 2000 //debug only, alpha_test()
#define B_SIZE 28125 // debug only, beta_test()
#define B_GAIN_S 5625

int defect_gain_correct(char *fin, char *gain, char *fout, MrcHeader *head,int threads);

__global__ void mutiplier_kernel_split_purge(float *coord, float *gain, long long, long long, long long, long long);

__global__ void transfer_kernel( float *coord, char *src, long long total_s, long long unit_n, int step );

float dispatcher_gpu_split_purge(float *coord_l, void *src, float *gain_l , long size_x, long size_y, long slice_n, int type, long );

bool cudaErrorCheck( cudaError_t, int );

bool alpha_test(); //type 2 input

bool beta_test( int type = 0 ); //type 0, 1, 6

int main(int argc, char *argv[]){
	bool k = true, h = true;
	if ( TEST_RUN ){
		//k = alpha_test();
		h = beta_test();
		if ( k == false ){
			printf("Alpha Failed \n");
			return 1;
		}
		if ( h == false ){
			printf( "Beta Failed \n" );
			return 2;
		}
		printf( "Passed \n" );
		return 0;
	}

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

	//printf("raw_name %s, out_name %s, gain_name %s\n",raw_name,out_name,gain_name);
    defect_gain_correct(raw_name,gain_name,out_name,inhead, MP_THREAD );
	/*
	free(raw_name);
	free(out_name);
	free(gain_name);
	*/
	return 0;
}

int defect_gain_correct(char *fin, char *gain, char *fout, MrcHeader *head,int threads){
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

	printf("pointer check(pre-allocate): %p -- %p \n", coor_xy, gain_xy);
	gain_xy=(float*)malloc(sizeof(float)*size_y*size_x);
	//coor_xy=(unsigned char*)malloc(sizeof(unsigned char)*size_y*size_x*slice_n);
	coor_xy=(float*)malloc(sizeof(float)*size_y*size_x*slice_n);
	printf("pointer check(post-allocate): %p -- %p \n", coor_xy, gain_xy);

	//read the input and gain file into memory
	//note: the byte size(second parameter) is basen on the type of the file
	printf( "--> Input Mrc type: %d \n", file_type );
	printf("Before dispatch, statues: src_size - %d, size_x - %d, size_y - %d, slice_n %d \n", input_size-1024, size_x, size_y, slice_n );

	switch(file_type)
                {
                case 0:size_bit=1;input_byte_c=(char*)malloc(sizeof(char)*(input_size-1024));fread(input_byte_c,size_bit,(input_size-1024)/size_bit,input);
				source = input_byte_c;
				break;
                case 1:size_bit=2;input_byte_s=(short*)malloc(sizeof(char)*(input_size-1024));fread(input_byte_s,size_bit,(input_size-1024)/size_bit,input);
				source = input_byte_s;
				break;
                case 2:size_bit=4;input_byte_f=(float*)malloc(sizeof(char)*(input_size-1024));fread(input_byte_f,size_bit,(input_size-1024)/size_bit,input);
				source = input_byte_f;
				break;
                case 6:size_bit=2;input_byte_u=(short*)malloc(sizeof(char)*(input_size-1024));fread(input_byte_u,size_bit,(input_size-1024)/size_bit,input);
				source = input_byte_u;
				break;
                default:printf("File type error!");
                }

	fread(gain_byte,sizeof(float),(gain_f_size-1024)/4,gain_f);

	//read gain 
	float badcut=6.0;
	long int num_gain=0;
	int defect_num=0;
	clock_t start = 0 ,finish = 0, mid = 0;

	long int total_num=0,gain_num=0,gain_m=0;
	
	start=clock();
	for(n=0; n<size_x*size_y; n++)
		{
		gain_xy[n]=(*(float*)(gain_byte+n))*10.0;
		}
	finish=clock();
	printf("read gain_time %ds \n",(finish-start)/CLOCKS_PER_SEC);
	mid=clock();
	
	float indicator=0;
	//printf("pointer check(pre-invoke): %p -- %p \n", coor_xy, gain_xy);
	indicator=dispatcher_gpu_split_purge( coor_xy, source, gain_xy, size_x, size_y, slice_n, file_type, input_size-1024 );

	//defect correction for points
	finish=clock();
	printf( "gain time: memcpy %ds, calc %ds\n",(indicator)/CLOCKS_PER_SEC, (finish-mid)/CLOCKS_PER_SEC );

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
		{
		
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
	printf("free time %ds \n",(finish-start)/CLOCKS_PER_SEC);
	printf("Defect and gain correction finished!\n");

	printf("En Taro Tassadar!!!\n");
	return 0;
}


////
/**********************/
float dispatcher_gpu_split_purge( float *coord_l, void *src,float *gain_l , long size_x, long size_y, long slice_n, int type, long src_size ){
	////NOTICE: src_size = sizeof( 'src_entry' )*size_x*size_y*slice_n
	//type check	
	if ( type != 0 && type != 1 && type != 2 && type != 6 ) {
		printf( "dispatcher_gpu_split_purge: wrong input Mrc type: %d \n", type );
		return -1;
	}

	//set up cuda 
	cudaSetDevice(0);	
	void *device_coord_1=NULL, *device_coord_2=NULL;
	void *device_gain_1=NULL, *device_gain_2=NULL;
	void *device_src=NULL;
	int left_over = 0;
	cudaError_t f1, f2, f3, f4;
	long long total_s=0, single_s=0, unit_n=0;
	long long offset = 0; 
	clock_t t1=0, t2=0;
	float ret = 0;

	single_s = size_x * size_y;
	total_s = single_s * slice_n;
	if ( total_s%2 != 0 ){
		left_over = 1;
	}
	offset = total_s/2;
	
	printf("Initialize allocation.\n");
	printf("Required total(GB): %ld \n", (sizeof(float)*total_s + sizeof(float)*single_s)/GIGA);
	printf("Required each: %ld -- %ld \n", sizeof(float)*total_s, sizeof(float)*single_s);
	f1 = cudaMalloc( (void **)&device_coord_1,  sizeof(float)*total_s/2  );
	f2 = cudaMalloc( (void **)&device_src, src_size );
	f4 = cudaMemcpy( device_src, src, src_size, cudaMemcpyHostToDevice );

	f3 = cudaMalloc( (void **)&device_coord_2,  sizeof(float)*(total_s/2 + left_over)  );
	//GRAM allocation check
	if ( f1 != cudaSuccess || f2 != cudaSuccess || f3 != cudaSuccess || f4 != cudaSuccess ){
		printf("cuda memory allocation failed: %s -- %s -- %s \n", f1, f2, f3 );
		exit(4);
	}
	printf("Allocation done.\n");
	
	printf("Before any Kernel, statues: leftover - %d, total_s - %d, offset - %d, src_size - %d, size_x - %d, size_y - %d, slice_n %d \n", left_over, total_s, offset, src_size, size_x, size_y, slice_n );
	t1 = clock();
	if (type != 2){ /////////this branch is yet to be reviewed
			printf(" type: %d \n", type);
			long unit_temp = ceil(total_s/(2.0*( GTHREAD_N )));
			int step=0;
			long offset_temp=0;
			switch( type ){
				case 0:
					step = 1;
					break;
				case 1:
					step = 2;
					break;
				case 6:
					step = 2;
					break;
				default:
					printf("What have you done?!...");
					exit(123);
			}

			offset_temp = (total_s/2)*step;
			printf("offest_temp: %d, unit_temp: %d, step: %d \n", offset_temp, unit_temp, step );
			//printf("--offest_temp: %f, unit_temp: %f, step: %f \n", (float)offset_temp, (float)unit_temp, (float)step );
			transfer_kernel<<< GRID_BLOCK, BLOCK_SIZE >>>( (float*)device_coord_1, (char *)device_src, total_s/2, unit_temp, step);

			transfer_kernel<<< GRID_BLOCK, BLOCK_SIZE >>>( (float*)device_coord_2,((char*)device_src + offset_temp), total_s/2 + left_over, unit_temp, step);

			f1 = cudaThreadSynchronize();
			cudaErrorCheck(f1, 321);
	}else{
		printf("type: 2\n");
		f1 = cudaMemcpy( device_coord_1, src, sizeof(float)*total_s/2, cudaMemcpyHostToDevice );
		f3 = cudaMemcpy( device_coord_2, (((float*)src) + total_s/2), sizeof(float)*(total_s/2 + left_over), cudaMemcpyHostToDevice );/////////////////////////////pH
		if ( f1 != cudaSuccess || f3 != cudaSuccess ){
			printf("purge finish-off failed: %s \n", f1 );
			exit(121);
		}
	}
	t2 = clock();
	ret = (t2-t1);
	f2 = cudaFree(device_src);
	if ( f2 != cudaSuccess ){
		printf("purge finish-off failed: %s \n", f1 );
		exit(121);
	}

	//end{mem transfer}

	//data transfer to device
	printf("Initialize Memcpy: host->device\n");
	f4 = cudaMalloc( (void **)&device_gain_1,  sizeof(float)*single_s );
	//f1 = cudaMemcpy( device_coord_1, coord_l, sizeof(float)*total_s/2, cudaMemcpyHostToDevice );
	f2 = cudaMemcpy( device_gain_1, gain_l, sizeof(float)*single_s, cudaMemcpyHostToDevice );
	//f3 = cudaMemcpy( device_coord_2, (coord_l + total_s/2), sizeof(float)*(total_s/2 + left_over), cudaMemcpyHostToDevice );
	//GRAM allocation check
	if (  f2 != cudaSuccess || f4 != cudaSuccess ){
		printf("cudaMemcpy failed(H->D): %s -- %s \n", f2, f4);
		exit(5);
	}
	//*/
	printf("Memcpy Done: host->device\n");
	
	//activate kernel
	unit_n = ceil(total_s/(2.0*( GTHREAD_N ))); 
	printf("Insider: the unit_n is %d \n", unit_n);	


			printf("waypoint Serpent\n");

			mutiplier_kernel_split_purge<<< GRID_BLOCK, BLOCK_SIZE >>>( (float*)device_coord_1, (float*)device_gain_1, total_s/2, single_s, unit_n, 0 ); 
			//HAZARD: slice_c might be a problem... + a offset, yes offset will do.BUT NOT FULLY VERIFIED YET
			//can switch cpu here and remove the sync below
			///*
			f1 = cudaThreadSynchronize();
			if ( f1 != cudaSuccess ){
				printf("cuda sync(mid-way) failed: %s \n", f1 );
				exit(12);
			}
			//*/
			mutiplier_kernel_split_purge<<< GRID_BLOCK, BLOCK_SIZE >>>( (float*)device_coord_2, (float*)device_gain_1, total_s/2+left_over, single_s, unit_n, offset );


			printf("Serpent out\n");

	//thread synchronization 
	f1 = cudaThreadSynchronize();
	if ( f1 != cudaSuccess ){
		printf("cuda sync failed: %s \n", f1 );
		exit(9);
	}
	printf("Sync Done\n");
	
	//data transfer from device
	///* Kernel validation
	printf("Initialize Memcpy: device->host\n");
	
	f1 = cudaMemcpy( coord_l, device_coord_1, sizeof(float)*total_s/2, cudaMemcpyDeviceToHost );
	printf("Memcpy: device->host, half-way through\n");
	////Narrator: replace cudaMemcpy manually,if cudaMemcpy doesn't work for too long addresses
	f2 = cudaMemcpy( (coord_l + total_s/2), device_coord_2, sizeof(float)*(total_s/2 + left_over), cudaMemcpyDeviceToHost );
	//Narrator: if manual copy doesn't work either, then we declare another agent and do a manual cascade.

	//GRAM allocation check
	if ( f1 != cudaSuccess || f2 != cudaSuccess ){
		printf("cudaMemcpy failed(D->H): %s -- %s \n", f1, f2 );
		exit(6);
	}
	printf("Memcpy Done: device->host\n");
	
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
	
	printf("Dispatcher out\n");
	return ret;
}

__global__ void mutiplier_kernel_split_purge(float *coord, float *gain, long long total_s, long long single_s, long long unit_n, long long offset){
	//execute all those mutiplications of type_c data
	long long index=0;
	long long i = 0; 
	index = blockDim.x*blockIdx.x + threadIdx.x;  //such fun
	long long slice_c=0;

	for (i=0; i<unit_n; i++){
		if (index*unit_n + i >= total_s ){
			return;
		}
		slice_c = (index*unit_n+i+offset)%(single_s);
		//printf("In M kernel: %d, gain %d at %d \n", coord[index*unit_n+i], gain[slice_c], index*unit_n+i );
		coord[index*unit_n+i] = (coord[index*unit_n+i] * gain[slice_c]);
		//printf("In M kernel: %d \n", coord[index*unit_n+i] );
	}
	
	return;
}

__global__ void transfer_kernel( float *coord, char *src, long long total_s, long long unit_n, int step ){
	long long index = 0;
	long long i=0;
	index = blockDim.x*blockIdx.x + threadIdx.x;

	for ( i=0; i<unit_n; i++ ){
		if ( index*unit_n+i >= total_s ){
			return;
		}

		coord[index*unit_n+i] = *(src + (index*unit_n+i)*step);
		//printf("In transfer kernel: %d at \n", coord[index*unit_n+i], index*unit_n+i );

	}

	return;
}

bool cudaErrorCheck( cudaError_t x, int error_code ){
	if ( x != cudaSuccess ){
		exit(error_code);
	}
	return true;
}

bool alpha_test(){
	char *log_name = "tsp_alpha_test.log";
	FILE *log_f = NULL, *ori_f = NULL, *out_f = NULL;
	
	float source[SIZE] = {0};
	float gain[SIZE/2] = {0};
	float des[SIZE] = {0};
	float sum_check = 0;
	int x_size = 10;
	long i=0, j=0;
	bool flag = true;
	int misfit_fh = 0, misfit_lh = 0;

	//initialization
	for ( j=0; j<SIZE; j++ ) {
		if ( j < SIZE/2 ){
			source[j] = 2;
			gain[j] = 2;
		}
		else if ( j < SIZE ){
			source[j] = 4;
		}
	}

	/*
	float dispatcher_gpu_split_purge(float *coord_l, void *src, float *gain_l , long size_x, long size_y, long slice_n, int type, long );
	*/
	dispatcher_gpu_split_purge( des, source, gain, SIZE/(2*x_size), x_size, 2, 2, SIZE*4 );
	
	//validation
	for ( i=0; i<SIZE; i++ ) {
		if ( i <SIZE/2 && des[i] != 4){
			flag = false;
			misfit_fh += 1;
		}
		else if( i >= SIZE/2 && i < SIZE && des[i] != 8 ){
			flag = false;
			misfit_lh += 1;
		}
		sum_check += des[i];
		continue;
	}
	
	printf( "sum_check: %f \n", sum_check );
	printf("Misfit, first-half: %d; last-half: %d \n", misfit_fh, misfit_lh );
	log_f = fopen( log_name, "w");
	ori_f = fopen( "tsp_alpha_ori.log", "w" );
	out_f = fopen( "tsp_alpha_out.log", "w" );

	for ( j=0; j<SIZE ; j++){
		fprintf( log_f, "%.0f ", des[j]-2*source[j] );
		fprintf( ori_f, "%.0f ", source[j] );
		fprintf( out_f, "%.0f ", des[j] );
	}
	fclose( log_f );
	fclose( ori_f );
	fclose( out_f );
	printf( "--" );
	return flag;
}
bool beta_test( int type ){
	if ( type != 0 && type != 1 && type != 6 ){
		printf("Wrong type: %d \n", type);
		return false;
	}

	char *log_name = "tsp_beta_test.log";
	FILE *log_f = NULL, *ori_f = NULL, *out_f = NULL;
	
	char source[B_SIZE] = {0};
	float gain[B_GAIN_S] = {0};
	float des[B_SIZE] = {0};
	float sum_check = 0;

	int x_size = 10;
	int slice_N = 2;

	long i=0, j=0;
	bool flag = true;
	int misfit_fh = 0, misfit_lh = 0;

	//initialization
	for (j = 0; j < B_GAIN_S; j++ ){
		gain[j] = 2;
	}

	for ( j=0; j<B_SIZE; j++ ) {
		if ( j < B_SIZE/2 ){
			
			source[j] = '!';
		}
		else{
			source[j] = '#';
		}
		des[j] = -1;
	}

	/*
	float dispatcher_gpu_split_purge(float *coord_l, void *src, float *gain_l , long size_x, long size_y, long slice_n, int type, long );
	*/
	//dispatcher_gpu_split_purge( des, source, gain, B_SIZE/(slice_N*x_size), x_size, slice_N, type, B_SIZE*sizeof(char) );
	dispatcher_gpu_split_purge( des, source, gain, 75, 75, 5, type, 75*75*5*sizeof(char) );
	
	//validation
	for ( i=0; i<B_SIZE; i++ ) {
		if (  i< B_SIZE/2 && (int)des[i] != 66){
			flag = false;
			misfit_fh += 1;
		}
		else if ( i >= B_SIZE/2 && (int)des[i] != 70 ){
			flag = false;
			misfit_lh += 1;
		}
	}
	
	printf("Misfit, first-half: %d, last-half: %d \n", misfit_fh, misfit_lh );
	log_f = fopen( log_name, "w");
	ori_f = fopen( "tsp_beta_ori.log", "w" );
	out_f = fopen( "tsp_beta_out.log", "w" );

	for ( j=0; j<B_SIZE ; j++){
		fprintf( log_f, "%.0f ", des[j]-2*(int)source[j] );
		fprintf( ori_f, "%c ", source[j] );
		fprintf( out_f, "%d ", (int)des[j] );
	}
	fclose( log_f );
	fclose( ori_f );
	fclose( out_f );
	printf( "--" );
	return flag;
}

////
/**********************/