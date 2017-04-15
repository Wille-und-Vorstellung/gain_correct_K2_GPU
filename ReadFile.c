#include "ReadFile.h"

char geo_init_fn[] = "geo_init.par";
char geo_ref_fn[] = "geo.par";
char xyshift_fn[] = "xyshift.par";
char xyshift2_fn[] = "xyshift2.par";
char ang_init_fn[] = "ang_init.par";

long long get_file_size(FILE *fin) {
	fseek(fin, 0, SEEK_END);

	return ftell(fin);

}

/*******************************************************************************************/
int mrc_read_head(FILE *fin, MrcHeader *head) {

	if (ftell(fin) != 0)
		rewind(fin);
	//if(ftello64(fin)!=0)rewind(fin);
	fread(head, (sizeof(MrcHeader) + head->next), 1, fin); //HAZARD? without initializing head->next?

	if (!(head->cmap[0] == 'M' && head->cmap[1] == 'A' && head->cmap[2] == 'P')) {
		printf(
				"Error with function 'mrc_read_head()'! Warning: Not MRC format! \n");
		return -1;
		//exit(1);
	}

	return 0;

}

/*******************************************************************************************/
int mrc_write_head(FILE *fout, MrcHeader *head) {

	if (ftell(fout) != 0)
		rewind(fout);

	if (!head
			|| !(head->cmap[0] == 'M' && head->cmap[1] == 'A'
					&& head->cmap[2] == 'P')) {
		printf(
				"Error with function 'mrc_write_head()'! Can not write the MrcHeader!\n");
		exit(1);
	}
	fwrite(head, (sizeof(MrcHeader) + head->next), 1, fout);

	return 0;
}

/*******************************************************************************************/
int mrc_init_head(MrcHeader *head) {

	head->nx = 1;
	head->ny = 1;
	head->nz = 1;

	head->mode = MRC_MODE_FLOAT;

	head->nxstart = 0;
	head->nystart = 0;
	head->nzstart = 0;

	head->mx = 1;
	head->my = 1;
	head->mz = 1;

	head->xlen = 1;
	head->ylen = 1;
	head->zlen = 1;

	head->alpha = 90;
	head->beta = 90;
	head->gamma = 90;

	head->mapc = 1;
	head->mapr = 2;
	head->maps = 3;

	head->amin = 0.0;
	head->amax = 255.0;
	head->amean = 128.0;

	head->ispg = 1;
	head->nsymbt = 0;

	head->next = 0;

	head->creatid = 1000;
	head->cmap[0] = 'M';
	head->cmap[1] = 'A';
	head->cmap[2] = 'P';

	head->stamp[0] = 'D';

	return 0;
}

/*******************************************************************************************/
MrcHeader *mrc_new_head() {

	MrcHeader *head;
	head = (MrcHeader *) malloc(sizeof(MrcHeader));
	mrc_init_head(head);

	return head;
}

/*******************************************************************************************/
int mrc_replace_head(char *outf, MrcHeader *head) {

	FILE *fout;
	if ((fout = fopen(outf, "r+")) == NULL) {
		printf("Cannot open file strike any key exit!\n");
		exit(1);
	}
	mrc_write_head(fout, head);
	fclose(fout);

	return TRUE;

}

/*******************************************************************************************/
int mrc_check_head(char *inf) {
	MrcHeader *head;
	head = (MrcHeader *) malloc(sizeof(MrcHeader));

	FILE *fin;

	if ((fin = fopen(inf, "r+")) == NULL) {
		printf("Cannot open file strike any key exit!\n");
		exit(1);
	}

	mrc_read_head(fin, head);
	fclose(fin);

	if (!head
			|| !(head->cmap[0] == 'M' && head->cmap[1] == 'A'
					&& head->cmap[2] == 'P')) {
		printf(
				"Error with function 'mrc_write_head()'! Can not write the MrcHeader!");
		exit(1);
	}

	return TRUE;

}

/*******************************************************************************************/
//only MRC_MODE_BYTE MRC_MODE_SHORT MRC_MODE_FLOAT will be considered in this function
int mrc_update_head(char *inoutf) {

	MrcHeader *head;
	head = (MrcHeader *) malloc(sizeof(MrcHeader));

	FILE *finout;

	if ((finout = fopen(inoutf, "r+")) == NULL) {
		printf("Cannot open file strike any key exit!\n");
		exit(1);

	}
	mrc_init_head(head);
	mrc_read_head(finout, head);
	double sum, sum_xy, amin, amax, amean;
	int k, pNum, pLen;
	//pNum is the number of pixels in one XoY  slice,
	//pLen is the
	unsigned long site;
	unsigned char *p_uchar;
	short *p_short;
	float *p_float;

	fseek(finout, (sizeof(MrcHeader) + head->next), SEEK_SET);

	pNum = head->nx * head->ny;
	switch (head->mode) {
	//switch start

	/**********case MRC_MODE_BYTE ***********/
	case MRC_MODE_BYTE:

		pLen = pNum * sizeof(unsigned char);

		if ((p_uchar = (unsigned char *) malloc(pLen)) == NULL) {
			printf("Function 'malloc' erro, while updating head!\n");
			exit(1);

		}
		printf("updating head!\n");
		fflush(stdout);
		fread(p_uchar, pLen, 1, finout);

		sum = sum_xy = 0.0;
		amin = amax = p_uchar[0];

		for (site = 0; site < pNum; site++) {
			if (p_uchar[site] > amax)
				amax = p_uchar[site];
			if (p_uchar[site] < amin)
				amin = p_uchar[site];
			sum_xy += p_uchar[site];

		}

		sum += sum_xy;

		for (k = 1; k < head->nz; k++) {

			sum_xy = 0.0;

			fread(p_uchar, pLen, 1, finout);

			for (site = 0; site < pNum; site++) {
				if (p_uchar[site] > amax)
					amax = p_uchar[site];
				if (p_uchar[site] < amin)
					amin = p_uchar[site];
				sum_xy += p_uchar[site];

			}

			sum += sum_xy;
		}
		amean = sum / ((long) head->nx * head->ny * head->nz);
		free(p_uchar);

		break;
		/**********case MRC_MODE_SHORT ***********/
	case MRC_MODE_SHORT:
		pLen = pNum * sizeof(short);
		if ((p_short = (short *) malloc(pLen)) == NULL) {
			printf("Function 'malloc' erro, while updating head!\n");
			exit(1);

		}
		printf("updating head!\n");
		fflush(stdout);
		fread(p_short, pLen, 1, finout);

		sum = sum_xy = 0.0;
		amin = amax = p_short[0];

		for (site = 0; site < pNum; site++) {
			if (p_short[site] > amax)
				amax = p_short[site];
			if (p_short[site] < amin)
				amin = p_short[site];
			sum_xy += p_short[site];

		}

		sum += sum_xy;

		for (k = 1; k < head->nz; k++) {

			sum_xy = 0.0;

			fread(p_short, pLen, 1, finout);

			for (site = 0; site < pNum; site++) {
				if (p_short[site] > amax)
					amax = p_short[site];
				if (p_short[site] < amin)
					amin = p_short[site];
				sum_xy += p_short[site];

			}

			sum += sum_xy;

		}

		amean = sum / ((long) head->nx * head->ny * head->nz);

		free(p_short);

		break;

		/**********case MRC_MODE_FLOAT ***********/
	case MRC_MODE_FLOAT:
		pLen = pNum * sizeof(float);
		if ((p_float = (float *) malloc(pLen)) == NULL) {
			printf("Function 'malloc' erro, while updating head!\n");
			fflush(stdout);
			exit(1);
		}
		printf("updating head!\n");
		fflush(stdout);
		fread(p_float, pLen, 1, finout);

		sum = sum_xy = 0.0;
		amin = amax = p_float[0];

		for (site = 0; site < pNum; site++) {
			if (p_float[site] > amax)
				amax = p_float[site];
			if (p_float[site] < amin)
				amin = p_float[site];
			sum_xy += p_float[site];

		}
		sum += sum_xy;

		for (k = 1; k < head->nz; k++) {
			sum_xy = 0.0;

			fread(p_float, pLen, 1, finout);

			for (site = 0; site < pNum; site++) {
				if (p_float[site] > amax)
					amax = p_float[site];
				if (p_float[site] < amin)
					amin = p_float[site];
				/*if(p_float[site]>10000000){
					printf("%d %d \n",k,site);
				}*/
				sum_xy += p_float[site];
			}
			sum += sum_xy;
		}
		amean = sum / ((long) head->nx * head->ny * head->nz);

		free(p_float);

		break;

	} //switch end
	head->amin = amin;
	head->amax = amax;
	head->amean = amean;

	fclose(finout);
	mrc_replace_head(inoutf, head);
	free(head);
	printf("updating head finished!\n");
	fflush(stdout);
	return 0;
}

/****************************************************************************************/
void mrc_new_file(char *newf, MrcHeader *outhead) {
	FILE *fnew;
	if ((fnew = fopen(newf, "w+")) == NULL) {
		printf("\nCannot open file strike any key exit!");
	}
	mrc_write_head(fnew, outhead);
	fclose(fnew);
}

/*******************************************************************************************/
int mrc_read_pixel(float *pix_gray, FILE *fin, MrcHeader *head, int x, int y,
		int z) {

	switch (head->mode) {
	case MRC_MODE_BYTE:

		fseek(fin,
				(sizeof(MrcHeader) + head->next)
						+ (z * head->nx * head->ny + y * head->nx + x)
								* sizeof(char), SEEK_SET);

		if ((fread(pix_gray, sizeof(char), 1, fin) == 0)) {
			printf(
					"Error with Function 'mrc_read_pixel()'! Reading file failed!");
			return -1;
		}
		break;

	case MRC_MODE_SHORT:

		fseek(fin,
				(sizeof(MrcHeader) + head->next)
						+ (z * head->nx * head->ny + y * head->nx + x)
								* sizeof(short), SEEK_SET);

		if ((fread(pix_gray, sizeof(short), 1, fin) == 0)) {
			printf(
					"Error with Function 'mrc_read_pixel()'! Reading file failed!");
			return -1;
		}
		break;

	case MRC_MODE_FLOAT:

		fseek(fin,
				(sizeof(MrcHeader) + head->next)
						+ (z * head->nx * head->ny + y * head->nx + x)
								* sizeof(float), SEEK_SET);

		if ((fread(pix_gray, sizeof(float), 1, fin) == 0)) {
			printf(
					"Error with Function 'mrc_read_pixel()'! Reading file failed!");
			return -1;
		}
		break;

	}

	return 0;

}

/*******************************************************************************************/
/*******slcN couts from 0 to N-1, so if you want to read the first slice slcN shoud be 0****/

int mrc_read_slice(FILE *fin, MrcHeader *head, int slcN, char axis,
		float *slcdata) {
//check the mrc file to make sure the size is exact in register with the head
	switch (head->mode) {
	case MRC_MODE_BYTE:
		if (get_file_size(fin) - (sizeof(MrcHeader) + head->next)
				!= head->nx * head->ny * head->nz * sizeof(char)) {
			printf("Error with Function 'mrc_read_slic()'! File size erro!\n");
		}
		break;

	case MRC_MODE_SHORT:
	case MRC_MODE_USHORT:
		if (get_file_size(fin) - (sizeof(MrcHeader) + head->next)
				!= head->nx * head->ny * head->nz * sizeof(short)) {
			printf("Error with Function 'mrc_read_slice()'! File size erro!\n");
		}
		break;

	case MRC_MODE_FLOAT:
		if (get_file_size(fin) - (sizeof(MrcHeader) + head->next)
				!= (long int)head->nx * head->ny * head->nz * sizeof(float)) {
			long int total_size = get_file_size(fin);
			int headsize = sizeof(MrcHeader) + head->next;
			printf("remainder:%lld shouldbe:%lld\n",get_file_size(fin) - (sizeof(MrcHeader) + head->next), (long long int)head->nx * head->ny * head->nz * sizeof(float));
			printf("%d %d %d %d \n", head->nx, head->ny, head->nz, head->next);
			fflush(stdout);
			printf("Error with Function 'mrc_read_slice()'! File size erro!\n");
		}
		break;

	default:
		printf("Error with Function 'mrc_read_slice()'! File type unknown!\n");

		break;
	}

	size_t psize;
	//int psize;
	short buf_short;
	unsigned short buf_ushort;
	unsigned char buf_byte;
	float buf_float;
	int i, k;

	switch (head->mode) {
	case MRC_MODE_BYTE:
		psize = sizeof(unsigned char);

		break;

	case MRC_MODE_SHORT:
	case MRC_MODE_USHORT:
		psize = sizeof(short);

		break;

	case MRC_MODE_FLOAT:
		psize = sizeof(float);

		break;
	}

	switch (axis) {

	/***********************************X************************************/
	case 'x':
	case 'X':

		fseek(fin, (sizeof(MrcHeader) + head->next) + slcN * psize, SEEK_SET);

		switch (head->mode) {
		case MRC_MODE_BYTE:
			for (i = 0; i < head->ny * head->nz; i++) {
				fread(&buf_byte, psize, 1, fin);
				slcdata[i] = (float) buf_byte;
				fseek(fin, (head->nx - 1) * psize, SEEK_CUR);
			}

			break;

		case MRC_MODE_SHORT:
			for (i = 0; i < head->ny * head->nz; i++) {
				fread(&buf_short, psize, 1, fin);
				slcdata[i] = (float) (buf_short);
				fseek(fin, (head->nx - 1) * psize, SEEK_CUR);
			}

			break;

		case MRC_MODE_USHORT:
			for (i = 0; i < head->ny * head->nz; i++) {
				fread(&buf_ushort, psize, 1, fin);
				slcdata[i] = (float) (buf_ushort);
				fseek(fin, (head->nx - 1) * psize, SEEK_CUR);
			}

			break;

		case MRC_MODE_FLOAT:
			for (i = 0; i < head->ny * head->nz; i++) {
				fread(&buf_float, psize, 1, fin);
				slcdata[i] = buf_float;
				fseek(fin, (head->nx - 1) * psize, SEEK_CUR);
			}
			break;

		}

		break;

		/***********************************Y************************************/
	case 'y':
	case 'Y':

		for (k = 0; k < head->nz; k++) {
			fseek(fin,
					(sizeof(MrcHeader) + head->next)
							+ (long int)psize
									* (k * (long int)head->nx * head->ny + head->nx * slcN),
					SEEK_SET);

			switch (head->mode) {
			case MRC_MODE_BYTE:
				for (i = 0; i < head->nx; i++) {
					fread(&buf_byte, psize, 1, fin);
					slcdata[k * head->nx + i] = (float) buf_byte;
				}

				break;

			case MRC_MODE_SHORT:
				for (i = 0; i < head->nx; i++) {
					fread(&buf_short, psize, 1, fin);
					slcdata[k * head->nx + i] = (float) (buf_short);
				}

				break;

			case MRC_MODE_USHORT:
				for (i = 0; i < head->nx; i++) {
					fread(&buf_ushort, psize, 1, fin);
					slcdata[k * head->nx + i] = (float) (buf_ushort);
				}

				break;

			case MRC_MODE_FLOAT:
				for (i = 0; i < head->nx; i++) {
					fread(&buf_float, psize, 1, fin);
					slcdata[k * head->nx + i] = buf_float;
				}

				break;

			}

		}
		break;

		/***********************************Z************************************/
	case 'z':
	case 'Z':
		fseek(fin,
				(sizeof(MrcHeader) + head->next)
						+ psize * slcN * head->nx * head->ny, SEEK_SET);

		if (head->mode == MRC_MODE_FLOAT)
			fread(slcdata, psize * head->nx * head->ny, 1, fin);

		else if (head->mode == MRC_MODE_BYTE) {
			for (i = 0; i < head->nx * head->ny; i++) {
				fread(&buf_byte, psize, 1, fin);
				slcdata[i] = (float) buf_byte;
			}
		}

		else if (head->mode == MRC_MODE_SHORT) {
			for (i = 0; i < head->nx * head->ny; i++) {
				fread(&buf_short, psize, 1, fin);
				slcdata[i] = (float) buf_short;
			}
		}

		else if (head->mode == MRC_MODE_USHORT) {
			for (i = 0; i < head->nx * head->ny; i++) {
				fread(&buf_ushort, psize, 1, fin);
				slcdata[i] = (float) buf_ushort;
			}
		}

		break;

	}

	return 0;
}

/*******************************************************************************************/
/*******slcN couts from 0 to N-1, so if you want to read the first slice slcN shoud be 0****/

int mrc_read_slice_double(FILE *fin, MrcHeader *head, int slcN, char axis,
		double *slcdata) {

//check the mrc file to make sure the size is exact in register with the head
	switch (head->mode) {
	case MRC_MODE_BYTE:
		if (get_file_size(fin) - (sizeof(MrcHeader) + head->next)
				!= head->nx * head->ny * head->nz * sizeof(char)) {
			printf(
					"Error with Function 'mrc_read_slice_double()'! File size erro!");
		}
		break;

	case MRC_MODE_SHORT:
	case MRC_MODE_USHORT:
		if (get_file_size(fin) - (sizeof(MrcHeader) + head->next)
				!= head->nx * head->ny * head->nz * sizeof(short)) {
			printf(
					"Error with Function 'mrc_read_slice_double()'! File size erro!");
		}
		break;

	case MRC_MODE_FLOAT:
		if (get_file_size(fin) - (sizeof(MrcHeader) + head->next)
				!= head->nx * head->ny * head->nz * sizeof(float)) {
			printf(
					"Error with Function 'mrc_read_slice_double()'! File size erro!");
		}
		break;

	default:
		printf(
				"Error with Function 'mrc_read_slice_double()'! File type unknown!");
		return -1;
		break;
	}

	int psize;
	short buf_short;
	unsigned char buf_byte;
	float buf_float;
	int i, k;

	switch (head->mode) {
	case MRC_MODE_BYTE:
		psize = sizeof(unsigned char);
		break;

	case MRC_MODE_SHORT:
	case MRC_MODE_USHORT:
		psize = sizeof(short);
		break;

	case MRC_MODE_FLOAT:
		psize = sizeof(float);
		break;
	}

	switch (axis) {

	/***********************************X************************************/
	case 'x':
	case 'X':

		fseek(fin, (sizeof(MrcHeader) + head->next) + slcN * psize, SEEK_SET);

		switch (head->mode) {
		case MRC_MODE_BYTE:
			for (i = 0; i < head->ny * head->nz; i++) {
				fread(&buf_byte, psize, 1, fin);
				slcdata[i] = (double) buf_byte;
				fseek(fin, (head->nx - 1) * psize, SEEK_CUR);
			}

			break;

		case MRC_MODE_SHORT:
			for (i = 0; i < head->ny * head->nz; i++) {
				fread(&buf_short, psize, 1, fin);
				slcdata[i] = (double) (buf_short);
				fseek(fin, (head->nx - 1) * psize, SEEK_CUR);
			}

			break;

		case MRC_MODE_FLOAT:
			for (i = 0; i < head->ny * head->nz; i++) {
				fread(&buf_float, psize, 1, fin);
				slcdata[i] = (double) buf_float;
				fseek(fin, (head->nx - 1) * psize, SEEK_CUR);
			}
			break;

		}

		break;

		/***********************************Y************************************/
	case 'y':
	case 'Y':

		for (k = 0; k < head->nz; k++) {
			fseek(fin,
					(sizeof(MrcHeader) + head->next)
							+ (k * head->nx * head->ny + head->nx * slcN)
									* psize, SEEK_SET);

			switch (head->mode) {
			case MRC_MODE_BYTE:
				for (i = 0; i < head->nx; i++) {
					fread(&buf_byte, psize, 1, fin);
					slcdata[k * head->nx + i] = (double) buf_byte;
				}

				break;

			case MRC_MODE_SHORT:
				for (i = 0; i < head->nx; i++) {
					fread(&buf_short, psize, 1, fin);
					slcdata[k * head->nx + i] = (double) (buf_short);
				}

				break;

			case MRC_MODE_FLOAT:
				for (i = 0; i < head->nx; i++) {
					fread(&buf_float, psize, 1, fin);
					slcdata[k * head->nx + i] = (double) buf_float;
				}

				break;

			}

		}
		break;

		/***********************************Z************************************/
	case 'z':
	case 'Z':
		fseek(fin,
				(sizeof(MrcHeader) + head->next)
						+ slcN * head->nx * head->ny * psize, SEEK_SET);

		if (head->mode == MRC_MODE_FLOAT) {
			for (i = 0; i < head->nx * head->ny; i++) {
				fread(&buf_float, psize, 1, fin);
				slcdata[i] = (double) buf_float;
			}
		}

		else if (head->mode == MRC_MODE_BYTE) {
			for (i = 0; i < head->nx * head->ny; i++) {
				fread(&buf_byte, psize, 1, fin);
				slcdata[i] = (double) buf_byte;
			}
		}

		else if (head->mode == MRC_MODE_SHORT) {
			for (i = 0; i < head->nx * head->ny; i++) {
				fread(&buf_short, psize, 1, fin);
				slcdata[i] = (double) buf_short;
			}
		}
		break;

	}

	return 0;
}

int mrc_read_slice_byte(FILE *fin, MrcHeader *head, int slcN, char axis,
		unsigned char *slcdata) {

//check the mrc file to make sure the size is exact in register with the head
	switch (head->mode) {
	case MRC_MODE_BYTE:
		if (get_file_size(fin) - (sizeof(MrcHeader) + head->next)
				!= head->nx * head->ny * head->nz * sizeof(char)) {
			printf(
					"Error with Function 'mrc_read_slice_double()'! File size erro!");
		}
		break;

	case MRC_MODE_SHORT:
	case MRC_MODE_USHORT:
		if (get_file_size(fin) - (sizeof(MrcHeader) + head->next)
				!= head->nx * head->ny * head->nz * sizeof(short)) {
			printf(
					"Error with Function 'mrc_read_slice_double()'! File size erro!");
		}
		break;

	case MRC_MODE_FLOAT:
		if (get_file_size(fin) - (sizeof(MrcHeader) + head->next)
				!= head->nx * head->ny * head->nz * sizeof(float)) {
			printf(
					"Error with Function 'mrc_read_slice_double()'! File size erro!");
		}
		break;

	default:
		printf(
				"Error with Function 'mrc_read_slice_double()'! File type unknown!");
		return -1;
		break;
	}

	int psize;
	short buf_short;
	unsigned char buf_byte;
	float buf_float;
	int i, k;

	float min = head->amin;
	float max = head->amax;

	float norm_coef = 255.0 / (max - min);
	switch (head->mode) {
	case MRC_MODE_BYTE:
		psize = sizeof(unsigned char);
		break;

	case MRC_MODE_SHORT:
	case MRC_MODE_USHORT:
		psize = sizeof(short);
		break;

	case MRC_MODE_FLOAT:
		psize = sizeof(float);
		break;
	}

	switch (axis) {

	/***********************************X************************************/
	case 'x':
	case 'X':

		fseek(fin, (sizeof(MrcHeader) + head->next) + slcN * psize, SEEK_SET);

		switch (head->mode) {
		case MRC_MODE_BYTE:
			for (i = 0; i < head->ny * head->nz; i++) {
				fread(&buf_byte, psize, 1, fin);
				slcdata[i] = (unsigned char) (norm_coef * (max - buf_byte));
				fseek(fin, (head->nx - 1) * psize, SEEK_CUR);
			}

			break;

		case MRC_MODE_SHORT:
			for (i = 0; i < head->ny * head->nz; i++) {
				fread(&buf_short, psize, 1, fin);
				slcdata[i] = (unsigned char) (norm_coef * (max - buf_short));
				fseek(fin, (head->nx - 1) * psize, SEEK_CUR);
			}

			break;

		case MRC_MODE_FLOAT:
			for (i = 0; i < head->ny * head->nz; i++) {
				fread(&buf_float, psize, 1, fin);
				slcdata[i] = (unsigned char) (norm_coef * (max - buf_float));
				fseek(fin, (head->nx - 1) * psize, SEEK_CUR);
			}
			break;

		}

		break;

		/***********************************Y************************************/
	case 'y':
	case 'Y':

		for (k = 0; k < head->nz; k++) {
			fseek(fin,
					(sizeof(MrcHeader) + head->next)
							+ (k * head->nx * head->ny + head->nx * slcN)
									* psize, SEEK_SET);

			switch (head->mode) {
			case MRC_MODE_BYTE:
				for (i = 0; i < head->nx; i++) {
					fread(&buf_byte, psize, 1, fin);
					slcdata[k * head->nx + i] = (unsigned char) (norm_coef
							* (max - buf_byte));
					buf_byte;
				}

				break;

			case MRC_MODE_SHORT:
				for (i = 0; i < head->nx; i++) {
					fread(&buf_short, psize, 1, fin);
					slcdata[k * head->nx + i] = (unsigned char) (norm_coef
							* (max - buf_short));
				}

				break;

			case MRC_MODE_FLOAT:
				for (i = 0; i < head->nx; i++) {
					fread(&buf_float, psize, 1, fin);
					slcdata[k * head->nx + i] = (unsigned char) (norm_coef
							* (max - buf_float));
				}

				break;

			}

		}
		break;

		/***********************************Z************************************/
	case 'z':
	case 'Z':
		fseek(fin,
				(sizeof(MrcHeader) + head->next)
						+ slcN * head->nx * head->ny * psize, SEEK_SET);

		if (head->mode == MRC_MODE_FLOAT) {
			for (i = 0; i < head->nx * head->ny; i++) {
				fread(&buf_float, psize, 1, fin);
				slcdata[i] = (unsigned char) (norm_coef * (max - buf_float));
			}
		}

		else if (head->mode == MRC_MODE_BYTE) {
			for (i = 0; i < head->nx * head->ny; i++) {
				fread(&buf_byte, psize, 1, fin);
				slcdata[i] = (unsigned char) (norm_coef * (max - buf_byte));
			}
		}

		else if (head->mode == MRC_MODE_SHORT) {
			for (i = 0; i < head->nx * head->ny; i++) {
				fread(&buf_short, psize, 1, fin);
				slcdata[i] = (unsigned char) (norm_coef * (max - buf_short));
			}
		}
		break;

	}

	return 0;
}

/*****************************************************************************************************/

int mrc_read_all(FILE *fin, MrcHeader *head, float *mrc_data_all) {

//check the mrc file to make sure the size is exact in register with the head
	switch (head->mode) {
	case MRC_MODE_BYTE:
		if (get_file_size(fin) - (sizeof(MrcHeader) + head->next)
				!= head->nx * head->ny * head->nz * sizeof(char)) {
			printf("Error with Function 'mrc_read_all()'! File size erro!");
		}
		break;

	case MRC_MODE_SHORT:
	case MRC_MODE_USHORT:
		if (get_file_size(fin) - (sizeof(MrcHeader) + head->next)
				!= head->nx * head->ny * head->nz * sizeof(short)) {
			printf("Error with Function 'mrc_read_all()'! File size erro!");
		}
		break;

	case MRC_MODE_FLOAT:
		if (get_file_size(fin) - (sizeof(MrcHeader) + head->next)
				!= head->nx * head->ny * head->nz * sizeof(float)) {
			printf("Error with Function 'mrc_read_all()'! File size erro!");
		}
		break;

	default:
		printf("Error with Function 'mrc_read_all()'! File type unknown!");

		break;
	}

	long i;
	unsigned char buf_byte;
	short buf_short;
	short buf_ushort;

	fseek(fin, (sizeof(MrcHeader) + head->next), SEEK_SET);

	switch (head->mode) {
	case MRC_MODE_BYTE:

		for (i = 0; i < head->nx * head->ny * head->nz; i++) {
			fread(&buf_byte, sizeof(char), 1, fin);
			mrc_data_all[i] = (float) buf_byte;
		}
		break;

	case MRC_MODE_SHORT:

		for (i = 0; i < head->nx * head->ny * head->nz; i++) {
			fread(&buf_short, sizeof(short), 1, fin);
			mrc_data_all[i] = (float) buf_short;
		}
		break;

	case MRC_MODE_USHORT:

		for (i = 0; i < head->nx * head->ny * head->nz; i++) {
			fread(&buf_ushort, sizeof(short), 1, fin);
			mrc_data_all[i] = (float) buf_ushort;
		}
		break;

	case MRC_MODE_FLOAT:

		if ((fread(mrc_data_all, head->nx * head->ny * head->nz * sizeof(float),
				1, fin) == 0)) {
			printf(
					"Error with Function 'mrc_read_all()'! Reading file failed!");
			return -1;
		}
		break;

	default:
		printf("Error with Function 'mrc_read_all()'! File type unknown!");
		break;

	}

	return 0;
}

/*******************************************************************************************/
int mrc_add_sliceN(FILE *fout, MrcHeader *headout, float *slcdata, int slcN) {

	fseek(fout,
			(sizeof(MrcHeader) + headout->next)
					+ (long int )sizeof(float) * headout->nx * headout->ny * slcN,
			SEEK_SET);
	fwrite(slcdata, sizeof(float) * headout->nx * headout->ny, 1, fout);
	return 0;
}

/*******************************************************************************************/
int mrc_add_slice(FILE *fout, MrcHeader *headout, float *slcdata) {

	fseek(fout, 0, SEEK_END);
	fwrite(slcdata, headout->nx * headout->ny * sizeof(float), 1, fout);
	return 0;
}

/*******************************************************************************************/
int mrc_add_slice_byte(FILE *fout, MrcHeader *headout, unsigned char *slcdata) {

	fseek(fout, 0, SEEK_END);
	fwrite(slcdata, headout->nx * headout->ny * sizeof(float), 1, fout);
	return 0;
}

/****************************************************************************************/
void mrc_update_slcX(FILE *fout, int slcN, float *data) {
	int i;
	MrcHeader *inhead;
	inhead = (MrcHeader *) malloc(sizeof(MrcHeader));
	mrc_read_head(fout, inhead);

	fseek(fout, (sizeof(MrcHeader) + inhead->next) + slcN * sizeof(float),
			SEEK_SET);
	for (i = 0; i < inhead->ny * inhead->nz; i++) {
		fwrite(data + i, sizeof(float), 1, fout);
		fseek(fout, sizeof(float) * (inhead->nx - 1), SEEK_CUR);
	}
	free(inhead);
}

/*******************************************************************************************/
int mrc_flipyz(char *inf, char *outf) {
	printf("\nBegin flipping:\n");
	fflush(stdout);
	MrcHeader *inhead, *outhead;
	inhead = (MrcHeader *) malloc(sizeof(MrcHeader));
	outhead = (MrcHeader *) malloc(sizeof(MrcHeader));

	FILE *fin, *fout;
	if ((fin = fopen(inf, "r")) == NULL) {
		printf("\nCannot open file strike any key exit!\n");
	}

	if ((fout = fopen(outf, "w+")) == NULL) {
		printf("\nCannot open file strike any key exit!\n");
	}

	fflush(stdout);
	mrc_init_head(inhead);
	mrc_read_head(fin, inhead);
	//mrc_read_head(fin,outhead);
	mrc_init_head(outhead);
	memcpy(outhead, inhead, sizeof(MrcHeader));
	//memcpy(outhead, inhead, sizeof(MrcHeader));
	outhead->nx = inhead->nx;
	outhead->ny = inhead->nz;
	outhead->nz = inhead->ny;
	outhead->mode = MRC_MODE_FLOAT;
	mrc_write_head(fout, outhead);

	float *buf;
	buf = (float *) malloc(sizeof(float) * inhead->nx * inhead->nz);
	int j;

	for (j = 0; j < inhead->ny; j++) {
		printf("%d \n", j);
		fflush(stdout);
		mrc_read_slice(fin, inhead, j, 'y', buf);
		mrc_add_slice(fout, outhead, buf);
	}

	free(inhead);
	free(outhead);
	free(buf);

	fclose(fin);
	fclose(fout);
	printf("\nflipping finished!\n");
	fflush(stdout);
	return 0;
}
int mrc_write_slice(FILE *fout, MrcHeader  *head, int slcN,char axis,float *slcdata)
{
  int psize;
  if (head->mode==MRC_MODE_FLOAT)
  psize=sizeof(float);
  else {
        printf ("outfile headmode is error!\n");
        return 1;
        }

  int i,k;
  off_t offset;

  switch(axis)
  {

/***********************************X************************************/
    case 'x':
    case 'X':

      fseeko(fout,(sizeof(MrcHeader) + head->next)+slcN*psize,SEEK_SET);

          for(i=0;i<head->ny*head->nz;i++)
            {
            fwrite(slcdata+i,psize,1,fout);
            fseeko(fout,(head->nx-1)*psize,SEEK_CUR);
            }
    break;

/***********************************Y************************************/
    case 'y':
    case 'Y':
      fseeko(fout,(sizeof(MrcHeader) + head->next)+slcN*head->nx*psize,SEEK_SET);

          for(k=0;k<head->nz;k++)
            {
             fwrite(slcdata+k*head->nx,psize,head->nx,fout);
             fseeko(fout,head->nx*(head->ny-1)*psize,SEEK_CUR);
            }
    break;

/***********************************Z************************************/
    case 'z':
    case 'Z'://problem
      //fseeko(fout,1024+slcN*head->nx*head->ny*psize,SEEK_SET);
      offset=head->nx*head->ny;
      offset*=(slcN*psize);
      offset+=(sizeof(MrcHeader) + head->next);
      fseeko(fout, offset, SEEK_SET);
      fwrite(slcdata,psize,head->nx*head->ny,fout);

    break;

    }
return 0;
}

