#line 228 "/home/sajid2/benchmarks/splash2/codes/null_macros/c.m4.null"

#line 1 "main.C"
/*************************************************************************/
/*                                                                       */
/*  Copyright (c) 1994 Stanford University                               */
/*                                                                       */
/*  All rights reserved.                                                 */
/*                                                                       */
/*  Permission is given to use, copy, and modify this software for any   */
/*  non-commercial purpose as long as this copyright notice is not       */
/*  removed.  All other uses, including redistribution in whole or in    */
/*  part, are forbidden without prior written permission.                */
/*                                                                       */
/*  This software is provided with absolutely no warranty and no         */
/*  support.                                                             */
/*                                                                       */
/*************************************************************************/

/*************************************************************************
*                                                                        *
*     main.c:  Starting point for rendering system.                      *
*                                                                        *
      Usage:  VOLREND num_processes input_file [-a]

      where input_file is head for the head data set. i.e. the filename
          without a suffix.
      and the -a option enables adaptive sampling of pixels.

*************************************************************************/

#include "incl.h"
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <limits.h>
#include "tiffio.h"

#define SH_MEM_AMT 60000000


#line 38
#include <pthread.h>
#line 38
#include <sys/time.h>
#line 38
#include <unistd.h>
#line 38
#include <stdlib.h>
#line 38
#define MAX_THREADS 32
#line 38
pthread_t PThreadTable[MAX_THREADS];
#line 38


#include "anl.h"

struct GlobalMemory *Global;

long image_section[NI];
long voxel_section[NM];

long num_nodes,frame;
long num_blocks,num_xblocks,num_yblocks;
PIXEL *image_address;
MPIXEL *mask_image_address;
PIXEL *image_block,*mask_image_block;
PIXEL *shd_address;
BOOLEAN *sbit_address;
long shd_length;
long image_len[NI], mask_image_len[NI];
int  image_length;
long mask_image_length;
char filename[FILENAME_STRING_SIZE];

void mclock(long stoptime, long starttime, long *exectime)
{
  if (stoptime < starttime)
    *exectime = ((ULONG_MAX - starttime) + stoptime)/1000;
  else
    *exectime = (stoptime - starttime)/1000;
}


int main(int argc, char *argv[])
{
  if ((argc < 3) || (strncmp(argv[1],"-h",strlen("-h")) == 0) || (strncmp(argv[1],"-h",strlen("-H")) == 0)){
    printf("usage:  VOLREND num_processes input_file\n");
    exit(-1);
  }

  {;};

  num_nodes = atol(argv[1]);

  strcpy(filename,argv[2]);

  if (argc == 4) {
    if (strncmp(argv[3],"-a",strlen("-a")) == 0)
      adaptive = YES;
    else {
      printf("usage:  VOLREND num_processes input_file [-a] \n");
      exit(-1);
    }
  }

  Frame();

/*  if (num_nodes > 1)
    {
#line 94
	unsigned long	i, Error;
#line 94
	for (i = 0; i < (num_nodes-1) - 1; i++) {
#line 94
		Error = pthread_join(PThreadTable[i], NULL);
#line 94
		if (Error != 0) {
#line 94
			printf("Error in pthread_join().\n");
#line 94
			exit(-1);
#line 94
		}
#line 94
	}
#line 94
};*/
  if (num_nodes > 1)
    {
#line 96
	unsigned long	i, Error;
#line 96
	for (i = 0; i < (num_nodes) - 1; i++) {
#line 96
		Error = pthread_join(PThreadTable[i], NULL);
#line 96
		if (Error != 0) {
#line 96
			printf("Error in pthread_join().\n");
#line 96
			exit(-1);
#line 96
		}
#line 96
	}
#line 96
};
  {exit(0);};
}


void Frame()
{
  long starttime,stoptime,exectime,i;

  Init_Options();

  printf("*****Entering init_decomposition with num_nodes = %ld\n",num_nodes);
  fflush(stdout);

  Init_Decomposition();

  printf("*****Exited init_decomposition with num_nodes = %ld\n",num_nodes);
  fflush(stdout);



  Global = (struct GlobalMemory *)malloc(sizeof(struct GlobalMemory));;
  {
#line 118
	unsigned long	Error;
#line 118

#line 118
	Error = pthread_mutex_init(&(Global->SlaveBarrier).mutex, NULL);
#line 118
	if (Error != 0) {
#line 118
		printf("Error while initializing barrier.\n");
#line 118
		exit(-1);
#line 118
	}
#line 118

#line 118
	Error = pthread_cond_init(&(Global->SlaveBarrier).cv, NULL);
#line 118
	if (Error != 0) {
#line 118
		printf("Error while initializing barrier.\n");
#line 118
		pthread_mutex_destroy(&(Global->SlaveBarrier).mutex);
#line 118
		exit(-1);
#line 118
	}
#line 118

#line 118
	(Global->SlaveBarrier).counter = 0;
#line 118
	(Global->SlaveBarrier).cycle = 0;
#line 118
};
  {
#line 119
	unsigned long	Error;
#line 119

#line 119
	Error = pthread_mutex_init(&(Global->TimeBarrier).mutex, NULL);
#line 119
	if (Error != 0) {
#line 119
		printf("Error while initializing barrier.\n");
#line 119
		exit(-1);
#line 119
	}
#line 119

#line 119
	Error = pthread_cond_init(&(Global->TimeBarrier).cv, NULL);
#line 119
	if (Error != 0) {
#line 119
		printf("Error while initializing barrier.\n");
#line 119
		pthread_mutex_destroy(&(Global->TimeBarrier).mutex);
#line 119
		exit(-1);
#line 119
	}
#line 119

#line 119
	(Global->TimeBarrier).counter = 0;
#line 119
	(Global->TimeBarrier).cycle = 0;
#line 119
};
  {pthread_mutex_init(&(Global->IndexLock), NULL);};
  {pthread_mutex_init(&(Global->CountLock), NULL);};
  {
#line 122
	unsigned long	i, Error;
#line 122

#line 122
	for (i = 0; i < MAX_NUMPROC+1; i++) {
#line 122
		Error = pthread_mutex_init(&Global->QLock[i], NULL);
#line 122
		if (Error != 0) {
#line 122
			printf("Error while initializing array of locks.\n");
#line 122
			exit(-1);
#line 122
		}
#line 122
	}
#line 122
};

  /* load dataset from file to each node */
#ifndef RENDER_ONLY
  {
#line 126
	struct timeval	FullTime;
#line 126

#line 126
	gettimeofday(&FullTime, NULL);
#line 126
	(starttime) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 126
};
  Load_Map(filename);
  {
#line 128
	struct timeval	FullTime;
#line 128

#line 128
	gettimeofday(&FullTime, NULL);
#line 128
	(stoptime) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 128
};
  mclock(stoptime,starttime,&exectime);
  printf("wall clock execution time to load map:  %lu ms\n", exectime);
#endif

  {
#line 133
	struct timeval	FullTime;
#line 133

#line 133
	gettimeofday(&FullTime, NULL);
#line 133
	(starttime) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 133
};
#ifndef RENDER_ONLY
  Compute_Normal();
#ifdef PREPROCESS
  Store_Normal(filename);
#endif
#else
  Load_Normal(filename);
#endif
  {
#line 142
	struct timeval	FullTime;
#line 142

#line 142
	gettimeofday(&FullTime, NULL);
#line 142
	(stoptime) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 142
};
  mclock(stoptime,starttime,&exectime);
  printf("wall clock execution time to compute normal:  %lu ms\n", exectime);

  {
#line 146
	struct timeval	FullTime;
#line 146

#line 146
	gettimeofday(&FullTime, NULL);
#line 146
	(starttime) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 146
};
#ifndef RENDER_ONLY
  Compute_Opacity();
#ifdef PREPROCESS
  Store_Opacity(filename);
#endif
#else
  Load_Opacity(filename);
#endif
  {
#line 155
	struct timeval	FullTime;
#line 155

#line 155
	gettimeofday(&FullTime, NULL);
#line 155
	(stoptime) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 155
};
  mclock(stoptime,starttime,&exectime);
  printf("wall clock execution time to compute opacity:  %lu ms\n", exectime);

  Compute_Pre_View();
  shd_length = LOOKUP_SIZE;
  Allocate_Shading_Table(&shd_address,shd_length);
  /* allocate space for image */
  image_len[X] = frust_len;
  image_len[Y] = frust_len;
  image_length = image_len[X] * image_len[Y];
  Allocate_Image(&image_address,image_length);

  if (num_nodes == 1) {
    block_xlen = image_len[X];
    block_ylen = image_len[Y];
    num_blocks = 1;
    num_xblocks = 1;
    num_yblocks = 1;
    image_block = image_address;
  }
  else {
    num_xblocks = ROUNDUP((float)image_len[X]/(float)block_xlen);
    num_yblocks = ROUNDUP((float)image_len[Y]/(float)block_ylen);
    num_blocks = num_xblocks * num_yblocks;
    Lallocate_Image(&image_block,block_xlen*block_ylen);
  }

  {
#line 183
	struct timeval	FullTime;
#line 183

#line 183
	gettimeofday(&FullTime, NULL);
#line 183
	(starttime) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 183
};
#ifndef RENDER_ONLY
  Compute_Octree();
#ifdef PREPROCESS
  Store_Octree(filename);
#endif
#else
  Load_Octree(filename);
#endif
  {
#line 192
	struct timeval	FullTime;
#line 192

#line 192
	gettimeofday(&FullTime, NULL);
#line 192
	(stoptime) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 192
};
  mclock(stoptime,starttime,&exectime);
  printf("wall clock execution time to compute octree:  %lu ms\n", exectime);

#ifdef PREPROCESS
  return;
#endif

  if (adaptive) {
    printf("1.\n");
    for (i=0; i<NI; i++) {
      mask_image_len[i] = image_len[i];
    }
    mask_image_length = image_length;
    Allocate_MImage(&mask_image_address, mask_image_length);
    if (num_nodes == 1)
      mask_image_block = (PIXEL *)mask_image_address;
    else
      Lallocate_Image(&mask_image_block, block_xlen*block_ylen);
    printf("2.\n");
  }

#ifndef RENDER_ONLY
  Deallocate_Map(&map_address);
#endif

  Global->Index = NODE0;

  printf("\nRendering...\n");
  printf("node\tframe\ttime\titime\trays\thrays\tsamples trilirped\n");

  {
#line 223
	long	i, Error;
#line 223

#line 223
	for (i = 0; i < (num_nodes) - 1; i++) {
#line 223
		Error = pthread_create(&PThreadTable[i], NULL, (void * (*)(void *))(Render_Loop), NULL);
#line 223
		if (Error != 0) {
#line 223
			printf("Error in pthread_create().\n");
#line 223
			exit(-1);
#line 223
		}
#line 223
	}
#line 223

#line 223
	Render_Loop();
#line 223
};
}


void Render_Loop()
{
  long step,i;
  PIXEL *local_image_address;
  MPIXEL *local_mask_image_address;
  char outfile[FILENAME_STRING_SIZE];
  long image_partition,mask_image_partition;
  float inv_num_nodes;
  long my_node;

  {pthread_mutex_lock(&(Global->IndexLock));};
  my_node = Global->Index++;
  {pthread_mutex_unlock(&(Global->IndexLock));};
  my_node = my_node%num_nodes;

  {;};
  {;};

/*  POSSIBLE ENHANCEMENT:  Here's where one might bind the process to a
    processor, if one wanted to.
*/

  inv_num_nodes = 1.0/(float)num_nodes;
  image_partition = ROUNDUP(image_length*inv_num_nodes);
  mask_image_partition = ROUNDUP(mask_image_length*inv_num_nodes);

#ifdef DIM
  for (dim=0; dim<NM; dim++) {
#endif

    for (step=0; step<ROTATE_STEPS; step++) { /* do rotation sequence */

/*  POSSIBLE ENHANCEMENT:  Here is where one might reset statistics, if
    		one wanted to.
*/

      frame = step;
      /* initialize images here */
      local_image_address = image_address + image_partition * my_node;
      local_mask_image_address = mask_image_address +
	mask_image_partition * my_node;

      {
#line 269
	unsigned long	Error, Cycle;
#line 269
	long		Cancel, Temp;
#line 269

#line 269
	Error = pthread_mutex_lock(&(Global->SlaveBarrier).mutex);
#line 269
	if (Error != 0) {
#line 269
		printf("Error while trying to get lock in barrier.\n");
#line 269
		exit(-1);
#line 269
	}
#line 269

#line 269
	Cycle = (Global->SlaveBarrier).cycle;
#line 269
	if (++(Global->SlaveBarrier).counter != (num_nodes)) {
#line 269
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 269
		while (Cycle == (Global->SlaveBarrier).cycle) {
#line 269
			Error = pthread_cond_wait(&(Global->SlaveBarrier).cv, &(Global->SlaveBarrier).mutex);
#line 269
			if (Error != 0) {
#line 269
				break;
#line 269
			}
#line 269
		}
#line 269
		pthread_setcancelstate(Cancel, &Temp);
#line 269
	} else {
#line 269
		(Global->SlaveBarrier).cycle = !(Global->SlaveBarrier).cycle;
#line 269
		(Global->SlaveBarrier).counter = 0;
#line 269
		Error = pthread_cond_broadcast(&(Global->SlaveBarrier).cv);
#line 269
	}
#line 269
	pthread_mutex_unlock(&(Global->SlaveBarrier).mutex);
#line 269
};

      if (my_node == num_nodes-1) {
	for (i=image_partition*my_node; i<image_length; i++)
	  *local_image_address++ = background;
	if (adaptive)
	  for (i=mask_image_partition*my_node; i<mask_image_length; i++)
	    *local_mask_image_address++ = NULL_PIXEL;
      }
      else {
	for (i=0; i<image_partition; i++)
	  *local_image_address++ = background;
	if (adaptive)
	  for (i=0; i<mask_image_partition; i++)
	    *local_mask_image_address++ = NULL_PIXEL;
      }

      if (my_node == ROOT) {
#ifdef DIM
	Select_View((float)STEP_SIZE, dim);
#else
        Select_View((float)STEP_SIZE, Y);
#endif
}

      {
#line 294
	unsigned long	Error, Cycle;
#line 294
	long		Cancel, Temp;
#line 294

#line 294
	Error = pthread_mutex_lock(&(Global->SlaveBarrier).mutex);
#line 294
	if (Error != 0) {
#line 294
		printf("Error while trying to get lock in barrier.\n");
#line 294
		exit(-1);
#line 294
	}
#line 294

#line 294
	Cycle = (Global->SlaveBarrier).cycle;
#line 294
	if (++(Global->SlaveBarrier).counter != (num_nodes)) {
#line 294
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 294
		while (Cycle == (Global->SlaveBarrier).cycle) {
#line 294
			Error = pthread_cond_wait(&(Global->SlaveBarrier).cv, &(Global->SlaveBarrier).mutex);
#line 294
			if (Error != 0) {
#line 294
				break;
#line 294
			}
#line 294
		}
#line 294
		pthread_setcancelstate(Cancel, &Temp);
#line 294
	} else {
#line 294
		(Global->SlaveBarrier).cycle = !(Global->SlaveBarrier).cycle;
#line 294
		(Global->SlaveBarrier).counter = 0;
#line 294
		Error = pthread_cond_broadcast(&(Global->SlaveBarrier).cv);
#line 294
	}
#line 294
	pthread_mutex_unlock(&(Global->SlaveBarrier).mutex);
#line 294
};

      Global->Counter = num_nodes;
      Global->Queue[num_nodes][0] = num_nodes;
      Global->Queue[my_node][0] = 0;

      Render(my_node);

      if (my_node == ROOT) {
	if (ROTATE_STEPS > 1) {
#ifdef DIM
	  sprintf(outfile, "%s_%ld",filename, 1000+dim*ROTATE_STEPS+step);
#else
	  sprintf(outfile, "%s_%ld.tiff",filename, 1000+step);
#endif
/*	  Store_Image(outfile);
          p = image_address;
          for (zz = 0;zz < image_length;zz++) {
            tiff_image[zz] = (long) ((*p)*256*256*256 + (*p)*256*256 +
				    (*p)*256 + (*p));
	    p++;
          }
tiff_save_rgba(outfile,tiff_image,image_len[X],image_len[Y]);  */
WriteGrayscaleTIFF(outfile, image_len[X],image_len[Y],image_len[X], image_address);
	} else {
/*	  Store_Image(filename);
	  p = image_address;
          for (zz = 0;zz < image_length;zz++) {
            tiff_image[zz] = (long) ((*p)*256*256*256 + (*p)*256*256 +
				    (*p)*256 + (*p));
	    p++;
          }
tiff_save_rgba(filename,tiff_image,image_len[X],image_len[Y]);    */
          strcat(filename,".tiff");
WriteGrayscaleTIFF(filename, image_len[X],image_len[Y],image_len[X], image_address);
	}
      }
    }
#ifdef DIM
  }
#endif
}

#if 0
void Error(char string[], char *arg1, char *arg2, char *arg3, char *arg4, char *arg5, char *arg6, char *arg7, char *arg8)
{
  fprintf(stderr,string,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8);
  exit(1);
}
#else
void Error(char string[], ...)
{
	va_list	ap;
	char	*arg1 = NULL, *arg2 = NULL, *arg3 = NULL, *arg4 = NULL, *arg5 = NULL, *arg6 = NULL, *arg7 = NULL, *arg8 = NULL;

	va_start(ap, string);
	arg1 = va_arg(ap, char *);
	arg1 = va_arg(ap, char *);
	arg1 = va_arg(ap, char *);
	arg1 = va_arg(ap, char *);
	arg1 = va_arg(ap, char *);
	arg1 = va_arg(ap, char *);
	arg1 = va_arg(ap, char *);
	arg1 = va_arg(ap, char *);
	va_end(ap);
	fprintf(stderr,string,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8);
	exit(1);
}
#endif

void Allocate_Image(PIXEL **address, long length)
{
  long i;

  printf("    Allocating image of %ld bytes...\n", length*sizeof(PIXEL));

  *address = (PIXEL *)malloc(length*sizeof(PIXEL));;

  if (*address == NULL)
	  Error("    No space available for image.\n", NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

  for (i=0; i<length; i++) *(*address+i) = 0;

}


void Allocate_MImage(MPIXEL **address, long length)
{
  long i;

  printf("    Allocating image of %ld bytes...\n", length*sizeof(MPIXEL));

  *address = (MPIXEL *)malloc(length*sizeof(MPIXEL));;

  if (*address == NULL)
    Error("    No space available for image.\n");

/*  POSSIBLE ENHANCEMENT:  Here's where one might distribute the
    opacity map among physical memories if one wanted to.
*/

  for (i=0; i<length; i++) *(*address+i) = 0;

}


void Lallocate_Image(PIXEL **address, long length)
{
  printf("    Allocating image of %ld bytes...\n", length*sizeof(PIXEL));
  *address = (PIXEL *)calloc(length,sizeof(PIXEL));
  if (*address == NULL)
    Error("    No space available for image.\n");

}


void Store_Image(char filename[])
{
  char local_filename[FILENAME_STRING_SIZE];
  long fd;
  short pix_version;
  short local_image_len[NI+1]; /* dimension larger than NI for backwards     */
                               /* compatibility of .pix file with no color   */

  local_image_len[X] = image_len[X];
  local_image_len[Y] = image_len[Y];
  local_image_len[2] = 1;

  pix_version = PIX_CUR_VERSION;
  strcpy(local_filename,filename);
  strcat(local_filename,".pix");
  fd = Create_File(local_filename);
  Write_Shorts(fd,(unsigned char *)&pix_version,(long)sizeof(pix_version));

  Write_Shorts(fd,(unsigned char *)local_image_len,(long)sizeof(local_image_len));
  Write_Longs(fd,(unsigned char *)&image_length,(long)sizeof(image_length));

  Write_Bytes(fd,image_address,(long)(image_length*sizeof(PIXEL)));
  Close_File(fd);
}


void Allocate_Shading_Table(PIXEL **address1, long length)
{
  long i;

  printf("    Allocating shade lookup table of %ld bytes...\n",
	 length*sizeof(PIXEL));

/*  POSSIBLE ENHANCEMENT:  If you want to replicate the shade table,
    replace the macro with a simple malloc in the line below */

  *address1 = (PIXEL *)malloc(length);;

  if (*address1 == NULL)
    Error("    No space available for table.\n");

/*  POSSIBLE ENHANCEMENT:  Here's where one might distribute the
    shading table among physical memories if one wanted to.
*/

  for (i=0; i<length; i++) *(*address1+i) = 0;

}


void Init_Decomposition()
{
  long factors[MAX_NUMPROC];
  double processors,newfactor;
  long i,sq,cu,maxcu,count;

  /* figure out what to divide dimensions of image and volume by to */
  /* partition data and computation to processors                   */
  if (num_nodes == 1) {
    image_section[X] = 1;
    image_section[Y] = 1;
    voxel_section[X] = 1;
    voxel_section[Y] = 1;
    voxel_section[Z] = 1;
  }
  else {
    count = 1;
    processors = (double)num_nodes;
    sq = (long)sqrt(processors);
    cu = (long)pow(processors,1.0/3.0);
    factors[0] = 1;

    for (i=2; i<sq; i++) {
      if (FRACT(processors/(double)i) == 0.0) {
        factors[count++] = i;
        if (i <= cu) {
          maxcu = i;
          newfactor = (double)(num_nodes/i);
        }
      }
    }
    count--;
    image_section[X] = factors[count];
    image_section[Y] = num_nodes/factors[count];

    sq = (long)sqrt(newfactor);
    count = 1;

    for (i=2; i<sq; i++) {
      if (FRACT(newfactor/(double)i) == 0.0)
        factors[count++] = i;
    }
    count--;
    voxel_section[X] = MIN(factors[count],maxcu);
    voxel_section[Y] = MAX(factors[count],maxcu);
    voxel_section[Z] = (long)newfactor/factors[count];
  }
}

/*
 * WriteGrayscaleTIFF
 *
 * Create a grayscale TIFF image.  The input is a sequence of bytes in an
 * array (each byte representing one pixel).  This is converted into a
 * compressed TIFF image file.
 *
 * Return value is 1 for success, 0 for failure.
 */

long WriteGrayscaleTIFF(char *filename, long width, long height, long scanbytes, unsigned char *data)
{
    long y;
    double factor;
    long c;
    unsigned long cmap[256];           /* output color map */
    TIFF *outimage;                     /* TIFF image handle */

    /* create a grayscale ramp for the output color map */
    factor = (double)((1 << 16) - 1) / (double)((1 << 8) - 1);
    for (c = 0; c < 256; c++)
        cmap[c] = (long)(c * factor);

    /* open and initialize output file */
    if ((outimage = TIFFOpen(filename, "w")) == NULL)
        return(0);
    TIFFSetField(outimage, TIFFTAG_IMAGEWIDTH, width);
    TIFFSetField(outimage, TIFFTAG_IMAGELENGTH, height);
    TIFFSetField(outimage, TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(outimage, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(outimage, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(outimage, TIFFTAG_COMPRESSION, COMPRESSION_LZW);
    TIFFSetField(outimage, TIFFTAG_ORIENTATION, ORIENTATION_BOTLEFT);
    TIFFSetField(outimage, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_PALETTE);
    TIFFSetField(outimage, TIFFTAG_COLORMAP, cmap, cmap, cmap);

    /* write the image data */
    for (y = 0; y < height; y++) {
        if (!TIFFWriteScanline(outimage, data, y, 0)) {
            TIFFClose(outimage);
            return(0);
        }
        data += scanbytes;
    }

    /* close the file */
    TIFFClose(outimage);
    return(1);
}
