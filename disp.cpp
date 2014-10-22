/*   CS580 HW   */
#include    "stdafx.h"  
#include	"Gz.h"
#include	"disp.h"


int GzNewFrameBuffer(char** framebuffer, int width, int height)
{
	if(!framebuffer || width <= 0 || width > MAXXRES || height <= 0 || height > MAXYRES)
		return GZ_FAILURE;

	//allocate memory for framebuffer : 3 * (sizeof)char * width * height
	char* fb =(char*) malloc(NUM_COMPONENTS*sizeof(char)*width*height);
	if(!fb)
		return GZ_FAILURE;
	
	//pass back pointer
	*framebuffer = fb;

	return GZ_SUCCESS;
}

int GzNewDisplay(GzDisplay	**display, int xRes, int yRes)
{

	if(!display || xRes <= 0 || xRes > MAXXRES || yRes <= 0 || yRes > MAXYRES)
		return GZ_FAILURE;

	//allocate memory for indicated resolution
	GzDisplay *disp;
	disp =(GzDisplay*) malloc(sizeof(GzDisplay));
	if(!disp)
		return GZ_FAILURE;

	disp->fbuf =(GzPixel*) malloc(sizeof(GzPixel)*xRes*yRes);
	if(!disp->fbuf){
		free(disp);
		return GZ_FAILURE;
	}

	disp->xres = xRes;
	disp->yres = yRes;

	//pass back pointer to GzDisplay object in display
	*display = disp;
	return GZ_SUCCESS;
}

//clean up, free memory
int GzFreeDisplay(GzDisplay	*display)
{
	if(display){
		if(display->fbuf)
			free(display->fbuf);
		free(display);
	}

	return GZ_SUCCESS;
}

//pass back values for a display
int GzGetDisplayParams(GzDisplay *display, int *xRes, int *yRes)
{
	if(!display || !xRes || !yRes)
		return GZ_FAILURE;
	
	*xRes = display->xres;
	*yRes = display->yres;

	return GZ_SUCCESS;
}

//Definitions for background color
#define SILVER_R 3083
#define SILVER_G 3083
#define SILVER_B 3083

//set everything to some default values - start a new frame	
int GzInitDisplay(GzDisplay	*display)
{
	if(!display)
		return GZ_FAILURE;
	const unsigned int max_i = display->xres*display->yres;
	for(unsigned int i = 0; i < max_i; ++i){
		display->fbuf[i].z = MAXINT;
		display->fbuf[i].alpha = MAX_INTENSITY;
		display->fbuf[i].red = SILVER_R;
		display->fbuf[i].green = SILVER_G;
		display->fbuf[i].blue = SILVER_B;
	}
	
	return GZ_SUCCESS;
}

//write pixel values into the display
int GzPutDisplay(GzDisplay *display, int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z)
{
	if(!display)
		return GZ_FAILURE;
	//Clip unnecessary pixels
	if(i < 0 || i >= display->xres || j < 0 || j >= display->yres)
		return GZ_SUCCESS;

	display->fbuf[ARRAY(i,j)].z = z;
	display->fbuf[ARRAY(i,j)].alpha = max(MIN_INTENSITY, min(MAX_INTENSITY, a));
	display->fbuf[ARRAY(i,j)].red = max(MIN_INTENSITY, min(MAX_INTENSITY, r));
	display->fbuf[ARRAY(i,j)].green = max(MIN_INTENSITY, min(MAX_INTENSITY, g));
	display->fbuf[ARRAY(i,j)].blue = max(MIN_INTENSITY, min(MAX_INTENSITY, b));

	return GZ_SUCCESS;
}

//pass back pixel value in the display
int GzGetDisplay(GzDisplay *display, int i, int j, GzIntensity *r, GzIntensity *g, GzIntensity *b, GzIntensity *a, GzDepth *z)
{
	if(!display || i < 0 || i >= display->xres || j < 0 || j >= display->yres || !r || !g || !b || !a || !z)
		return GZ_FAILURE;
	*r = display->fbuf[ARRAY(i,j)].red;
	*g = display->fbuf[ARRAY(i,j)].green;
	*b = display->fbuf[ARRAY(i,j)].blue;
	*a = display->fbuf[ARRAY(i,j)].alpha;
	*z = display->fbuf[ARRAY(i,j)].z;

	return GZ_SUCCESS;
}

//write pixels to ppm file 
int GzFlushDisplay2File(FILE* outfile, GzDisplay *display)
{
	if(!outfile || !display)
		return GZ_FAILURE;

	//Header: "P6 W H 255\n"
	fprintf(outfile, "P6 %d %d 255\n", display->xres, display->yres);
	//Raster
	const unsigned int max_i = display->xres*display->yres;
	for(unsigned int i = 0; i < max_i; ++i){	
		fprintf(outfile, "%c%c%c", (char) (display->fbuf[i].red >> 4), (char) (display->fbuf[i].green >> 4), (char) (display->fbuf[i].blue >> 4));
	}

	return GZ_SUCCESS;
}

//write pixels to framebuffer
int GzFlushDisplay2FrameBuffer(char* framebuffer, GzDisplay *display)
{
	if(!framebuffer || !display)
		return GZ_FAILURE;

	//store the pixel to the frame buffer as the order of blue, green, and red 
	const unsigned int max_i = display->xres*display->yres;
	for(unsigned int i = 0; i < max_i; ++i){
		framebuffer[i*3] =(char) (display->fbuf[i].blue >> 4);
		framebuffer[i*3+1] =(char) (display->fbuf[i].green >> 4);
		framebuffer[i*3+2] =(char) (display->fbuf[i].red >> 4);
	}

	return GZ_SUCCESS;
}