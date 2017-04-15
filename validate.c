#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ReadFile.h"
#include "ReadFile.c"

#define HEAD_L 1024
#define ARG_N 3
#define NAME_N 1024

const float EPSINON = 0.000001;

long compare(void*, void*, int ,int ,int, int type );

int main(int argc, char *argv[]){
    char *temp_name,*obj_name;
    FILE *temp_f, *obj_f;
    MrcHeader *inhead;

    temp_name=(char *)malloc(NAME_N*sizeof(char));
    obj_name=(char *)malloc(NAME_N*sizeof(char));
    inhead=(MrcHeader *)malloc(sizeof(MrcHeader));

    if(argc!=ARG_N){
        printf("Please input: mics1 mics2 \n");
        return 4;
    }
    temp_name=argv[1];
    obj_name=argv[2];
    
    temp_f=fopen(temp_name,"rb");	
    mrc_read_head(temp_f,inhead);
    obj_f = fopen(obj_name, "rb");
    rewind( temp_f );
    rewind( obj_f );
    //data loading and etc
    int size_x = inhead->nx;
    int size_y = inhead->ny;
    int slice_n = inhead->nz;
    int file_type = inhead->mode;
    inhead->mode = 2;
    
    long temp_length = 0;
    long obj_length = 0;
    fseek(temp_f, 0, SEEK_END);
	temp_length = ftell(temp_f);
	rewind(temp_f);
	fseek(obj_f, 0, SEEK_END);
	obj_length = ftell(obj_f);
	rewind(obj_f);

    //check point
    if ( temp_length != obj_length ){
        printf( "Eh...Houston, we've got a problem here...\n" );
        exit(2);
    }

    int bit_size=0;
    void *temp_content=NULL;
    void *obj_content=NULL;
    temp_content = malloc( sizeof(char*)*(temp_length-HEAD_L) );
    obj_content = malloc( sizeof(char*)*(temp_length-HEAD_L) );
    fseek( temp_f, HEAD_L, 0 );
    fseek( obj_f, HEAD_L, 0 );
    switch(file_type){
        case 0:
            bit_size=1;
            break;
        case 1:
            bit_size=2;
            break;
        case 2:
            bit_size=4;
            break;
        case 6:
            bit_size=2;
            break;
        default:
            printf("File type error!");
    }
    fread( temp_content, bit_size, (temp_length-HEAD_L)/bit_size, temp_f );
    fread( obj_content, bit_size, (obj_length-HEAD_L)/bit_size, obj_f );
    
    fclose(temp_f);
    fclose(obj_f);
    //pixel-wise compare 
    long diff_count=0;
    diff_count = compare( temp_content, obj_content, size_x, size_y, slice_n, file_type );

    printf("Different pixel numbers: %ld \n", diff_count);

    return 0;
}

/////////////////////////////////////////
long compare(void *temp, void *obj, int size_x, int size_y, int slice_n, int type ){
    long diff_n=0;
    //compare every entry
    float temp_1=0;
    float temp_2=0;
    long k=0;
    for ( k=0; k<size_x*size_y*slice_n; k++ ){
        switch(type){
            case 0:
                temp_1 = (float)*(((unsigned char*)temp)+k);
                temp_2 = (float)*(((unsigned char*)obj)+k);
                break;
            case 1:
                temp_1 = (float)*(((short*)temp)+k);
                temp_2 = (float)*(((short*)obj)+k);
                break;
            case 2:
                temp_1 = (float)*(((float*)temp)+k);
                temp_2 = (float)*(((float*)obj)+k);
                break;
            case 6:
                temp_1 = (float)*(((short*)temp)+k);
                temp_2 = (float)*(((short*)obj)+k);
                break;
            default:
                printf("what have you done....\n");
                exit(3);
        }
        if ( ( temp_1 - temp_2 > EPSINON ) || ( temp_1 - temp_2 < -EPSINON ) ){
            diff_n += 1;
        }
    }

    return diff_n;
}



