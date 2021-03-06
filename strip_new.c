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

/************************NOTICE**************/
//support mode 0 mrcs only, you've been warned.

int defect_gain_correct(char *fin, char *gain, char *fout, MrcHeader *head,int threads);

int main(int argc, char *argv[])
{
    char *raw_name,*out_name,*gain_name;
    FILE *file,*gain;
    int threads;
    MrcHeader *inhead;

    raw_name=(char *)malloc(NAME*sizeof(char));
    out_name=(char *)malloc(NAME*sizeof(char));
    gain_name=(char *)malloc(NAME*sizeof(char));

    inhead=(MrcHeader *)malloc(sizeof(MrcHeader));

    if(argc!=5)
    {
        printf("Please input: raw_image gain_name out_name threads\n");
        return 1;
    }
    raw_name=argv[1];
    gain_name=argv[2];
    out_name=argv[3];
    threads=atoi(argv[4]);

    file=fopen(raw_name,"rb");	
    mrc_read_head(file,inhead);
    fclose(file);

//    printf("raw_name %s, out_name %s, gain_name %s\n",raw_name,out_name,gain_name);
    defect_gain_correct(raw_name,gain_name,out_name,inhead,threads);
	return 1;
}

int defect_gain_correct(char *fin, char *gain, char *fout, MrcHeader *head,int threads)
	{
	unsigned char lbuf[LINE];
	float fbuf[LINE];
	float *tmp_array,*gain_xy;
	//unsigned char *coor_xy;
	char *coor_xy;
	int i,j,x,y,p,k,w,w_row,w_length,m,range,point_s,point_e,n_file,size_bit,*pxy_array;
	long int n;
	FILE *input,*output,*gain_f, *temp_head;
	char *input_byte_c;short *input_byte_s;float *input_byte_f;short *input_byte_u;
	srand((unsigned) time(NULL));
	

	//open input output and gain file
	printf("Start defect and gain correction|\ninput %s output %s gain %s\n",fin,fout,gain);
	input=fopen(fin,"rb");
	gain_f=fopen(gain,"rb");

	//revise input image header then written for output file
	int size_x=head->nx;
	int size_y=head->ny;
	int slice_n=head->nz;
        int file_type=head->mode;
	head->mode=0;
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
	//coor_xy=(float*)malloc(sizeof(float)*size_y*size_x*slice_n);
	coor_xy=(char*)malloc(sizeof(char)*size_y*size_x);

	//read the input and gain file into memory
	//note: the byte size(second parameter) is basen on the type of the file

	switch(file_type)
                {
                case 0:size_bit=1;input_byte_c=(char*)malloc(sizeof(char)*(input_size-1024));fread(input_byte_c,size_bit,(input_size-1024)/size_bit,input);break;
                case 1:size_bit=2;input_byte_s=(short*)malloc(sizeof(char)*(input_size-1024));fread(input_byte_s,size_bit,(input_size-1024)/size_bit,input);break;
                case 2:size_bit=4;input_byte_f=(float*)malloc(sizeof(char)*(input_size-1024));fread(input_byte_f,size_bit,(input_size-1024)/size_bit,input);break;
                case 6:size_bit=2;input_byte_u=(short*)malloc(sizeof(char)*(input_size-1024));fread(input_byte_u,size_bit,(input_size-1024)/size_bit,input);break;
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
	#pragma omp parallel for num_threads(threads) private(gain_num)
	//#pragma omp parallel for schedule(dynamic)
		for(n=0; n<size_x*size_y/**slice_n*/; n++)
			{
			gain_num=n%(size_x*size_y);
                        switch(file_type)
                                {
                                case 0:coor_xy[n]=(*(unsigned char*)(input_byte_c+n));//coor_xy[n]=(coor_xy[n])*(gain_xy[gain_num]);
								break;
                                case 1:coor_xy[n]=(*(short*)(input_byte_s+n));
								//coor_xy[n]=(coor_xy[n])*(gain_xy[gain_num]);
								break;
                                case 2:coor_xy[n]=(*(float*)(input_byte_f+n));
								//coor_xy[n]=(coor_xy[n])*(gain_xy[gain_num]);
								break;
                                case 6:coor_xy[n]=(*(short*)(input_byte_u+n));
								//coor_xy[n]=(coor_xy[n])*(gain_xy[gain_num]);
								break;
                                default:printf("File type error!\n");
                                }
			
			}
	//defect correction for points
	finish=clock();
	printf("gain time %d \n",finish-start);

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
	head->nz = 1;
	printf("header slice_n: %d\n", head->nz);
	printf("header mode: %d\n", head->mode);
	fwrite(head,(sizeof(MrcHeader)+head->next),1,output);
	start=clock();
	//for(n=0;n<size_x*size_y*slice_n;n++)
		{
		//fwrite(&coor_xy[n],sizeof(unsigned char),1,output);
		//fwrite(&coor_xy[n],sizeof(float),1,output);
		fwrite(coor_xy,sizeof(char),size_x*size_y*1/*extract 1 slice only*/,output);
		}
	finish=clock();
	printf("write_time %d \n",finish-start);

	start=clock();
	free(coor_xy);
	fclose(output);
	finish=clock();
	printf("free time %d \n",finish-start);
	printf("Defect and gain correction finished!\n");
	}
