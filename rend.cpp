/* CS580 Homework 3 */

#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"
#include	"float.h"

#define M_PI 3.14159265358979323846f 

#ifndef DDA_FS
#define DDA_FS
typedef struct{
	GzCoord start, end, current;
	float slopeX, slopeZ;
} GzDDA_FS;
#endif

#ifndef DDA_GS
#define DDA_GS
typedef struct{
	GzCoord start, end, current;
	GzColor startColor, endColor, currentColor;

	float slopeX, slopeZ;
	GzColor slopeColor;
} GzDDA_GS;
#endif

#ifndef DDA_PS
#define DDA_PS
typedef struct{
	GzCoord start, end, current;
	GzCoord startNormal, endNormal, currentNormal;

	float slopeX, slopeZ;
	GzCoord slopeNormal;
} GzDDA_PS;
#endif

#ifndef SPAN_FS
#define SPAN_FS
typedef struct{
	GzCoord start, end, current;
	float slopeZ;
} GzSpan_FS;
#endif

#ifndef SPAN_GS
#define SPAN_GS
typedef struct{
	GzCoord start, end, current;
	GzColor startColor, endColor, currentColor;

	float slopeZ;
	GzColor slopeColor;
} GzSpan_GS;
#endif

#ifndef SPAN_PS
#define SPAN_PS
typedef struct{
	GzCoord start, end, current;
	GzCoord startNormal, endNormal, currentNormal;

	float slopeZ;
	GzCoord slopeNormal;
} GzSpan_PS;
#endif

#define NUM_DIMENSIONS 3
#define NUM_VERTICES 3
#define NUM_EDGES 3

//Functions prototypes external to API
//Rasterizer's functions
int GzFrustumCull(GzRender *render, GzCoord vertexList[NUM_VERTICES], bool &culled);
int GzRasterizerFlatShading(GzRender *render, GzCoord vertexList[NUM_VERTICES]);
int GzRasterizerGouraudShading(GzRender *render, GzCoord vertexList[NUM_VERTICES], GzCoord normalList[NUM_VERTICES]);
int GzRasterizerPhongShading(GzRender *render, GzCoord vertexList[NUM_VERTICES], GzCoord normalList[NUM_VERTICES]);
void GzShading(GzRender *render, GzCoord normal, GzColor &color);
void GzTriangleNormal(GzCoord vertexList[NUM_VERTICES], GzCoord normal);
int GzPixelToDisplay(GzCoord *point, GzDisplay *display, GzColor *flatcolor);

//Transforms' functions
int GzSetupXsp(GzRender *render);
int GzSetupXpi(GzRender *render);
int GzSetupXiw(GzRender *render);

//Algebra functions
float GzVectorNorm(GzCoord *v);
float GzDotProduct(GzCoord *a, GzCoord *b);
void GzMatrixMult(GzMatrix A, GzMatrix B, GzMatrix &res);
void GzMatVecMult(GzMatrix A, float v[4], float (&res)[4]);

//Usefull functions
short ctoi(float color);
int GzInitDefaultCamera(GzRender *render);
int GzXformTriangle(GzRender *render, GzCoord *vertexList);
int GzXformNormal(GzRender *render, GzCoord *normalList);
void GzPreprocessMatrix(GzMatrix matrix, GzMatrix &res);

//Identity matrix
GzMatrix identityMat = {
	1.0f,	0.0f,	0.0f,	0.0f,
	0.0f,	1.0f,	0.0f,	0.0f,
	0.0f,	0.0f,	1.0f,	0.0f,
	0.0f,	0.0f,	0.0f,	1.0f
};

int GzRotXMat(float degree, GzMatrix mat)
{
// Create rotate matrix : rotate along x axis
// Pass back the matrix using mat value
	float theta = degree * M_PI / 180.0f;
	mat[0][0] = 1.0f; mat[0][1] = 0.0f; mat[0][2] = 0.0f; mat[0][3] = 0.0f;
	mat[1][0] = 0.0f; mat[1][1] = cos(theta); mat[1][2] = - sin(theta); mat[1][3] = 0.0f;
	mat[2][0] = 0.0f; mat[2][1] = sin(theta); mat[2][2] = cos(theta); mat[2][3] = 0.0f;
	mat[3][0] = 0.0f; mat[3][1] = 0.0f; mat[3][2] = 0.0f; mat[3][3] = 1.0f;

	return GZ_SUCCESS;
}


int GzRotYMat(float degree, GzMatrix mat)
{
// Create rotate matrix : rotate along y axis
// Pass back the matrix using mat value
	float theta = degree * M_PI / 180.0f;
	mat[0][0] = cos(theta); mat[0][1] = 0.0f; mat[0][2] = sin(theta); mat[0][3] = 0.0f;
	mat[1][0] = 0.0f; mat[1][1] = 1.0f; mat[1][2] = 0.0f; mat[1][3] = 0.0f;
	mat[2][0] = -sin(theta); mat[2][1] = 0.0f; mat[2][2] = cos(theta); mat[2][3] = 0.0f;
	mat[3][0] = 0.0f; mat[3][1] = 0.0f; mat[3][2] = 0.0f; mat[3][3] = 1.0f;

	return GZ_SUCCESS;
}


int GzRotZMat(float degree, GzMatrix mat)
{
// Create rotate matrix : rotate along z axis
// Pass back the matrix using mat value
	float theta = degree * M_PI / 180.0f;
	mat[0][0] = cos(theta); mat[0][1] = -sin(theta); mat[0][2] = 0.0f; mat[0][3] = 0.0f;
	mat[1][0] = sin(theta); mat[1][1] = cos(theta); mat[1][2] = 0.0f; mat[1][3] = 0.0f;
	mat[2][0] = 0.0f; mat[2][1] = 0.0f; mat[2][2] = 1.0f; mat[2][3] = 0.0f;
	mat[3][0] = 0.0f; mat[3][1] = 0.0f; mat[3][2] = 0.0f; mat[3][3] = 1.0f;

	return GZ_SUCCESS;
}


int GzTrxMat(GzCoord translate, GzMatrix mat)
{
// Create translation matrix
// Pass back the matrix using mat value
	mat[0][0] = 1.0f; mat[0][1] = 0.0f; mat[0][2] = 0.0f; mat[0][3] = translate[0];
	mat[1][0] = 0.0f; mat[1][1] = 1.0f; mat[1][2] = 0.0f; mat[1][3] = translate[1];
	mat[2][0] = 0.0f; mat[2][1] = 0.0f; mat[2][2] = 1.0f; mat[2][3] = translate[2];
	mat[3][0] = 0.0f; mat[3][1] = 0.0f; mat[3][2] = 0.0f; mat[3][3] = 1.0f;
	
	return GZ_SUCCESS;
}


int GzScaleMat(GzCoord scale, GzMatrix mat)
{
// Create scaling matrix
// Pass back the matrix using mat value
	mat[0][0] = scale[0]; mat[0][1] = 0.0f; mat[0][2] = 0.0f; mat[0][3] = 0.0f;
	mat[1][0] = 0.0f; mat[1][1] = scale[1]; mat[1][2] = 0.0f; mat[1][3] = 0.0f;
	mat[2][0] = 0.0f; mat[2][1] = 0.0f; mat[2][2] = scale[2]; mat[2][3] = 0.0f;
	mat[3][0] = 0.0f; mat[3][1] = 0.0f; mat[3][2] = 0.0f; mat[3][3] = 1.0f;
	
	return GZ_SUCCESS;
}


//----------------------------------------------------------
// Begin main functions

int GzNewRender(GzRender **render, GzDisplay	*display)
{
	if(!render || !display)
		return GZ_FAILURE;
	
	//malloc a renderer struct
	GzRender* pRender =(GzRender*) malloc(sizeof(GzRender));
	if(!pRender)
		return GZ_FAILURE;

	//span interpolator needs pointer to display for pixel writes
	pRender->display = display;

	//setup Xsp
	if(GzSetupXsp(pRender))
		return GZ_FAILURE;

	//init default camera
	if(GzInitDefaultCamera(pRender))
		return GZ_FAILURE;

	pRender->interp_mode = GZ_FLAT;
	pRender->numlights = 0;

	*render = pRender;
	return GZ_SUCCESS;
}


int GzFreeRender(GzRender *render)
{
	//free all renderer resources
	if(render) free(render);
	return GZ_SUCCESS;
}


int GzBeginRender(GzRender *render)
{

	//setup for start of each frame - init frame buffer color,alpha,z
	if(GzInitDisplay(render->display))
		return GZ_FAILURE;  /* init for new frame */

	//compute Xiw and projection xform Xpi from camera definition
	if(GzSetupXpi(render) || GzSetupXiw(render))
		return GZ_FAILURE;

	//init Ximage - put Xsp at base of stack, push on Xpi and Xiw 
	render->matlevel = -1;
	if(GzPushMatrix(render, render->Xsp) || GzPushMatrix(render, render->camera.Xpi) || GzPushMatrix(render, render->camera.Xiw))
		return GZ_FAILURE;

	//now stack contains Xsw and app can push model Xforms when needed 
	return GZ_SUCCESS;
}

int GzPutCamera(GzRender *render, GzCamera *camera)
{
	if(!render || !camera)
		return GZ_FAILURE;
/*
- overwrite renderer camera structure with new camera definition
*/
	memcpy(render->camera.position, camera->position, sizeof(GzCoord));
	memcpy(render->camera.lookat, camera->lookat, sizeof(GzCoord));
	memcpy(render->camera.worldup, camera->worldup, sizeof(GzCoord));
	render->camera.FOV = camera->FOV;

	return GZ_SUCCESS;	
}

//Push a matrix into the Ximage stack
int GzPushMatrix(GzRender *render, GzMatrix	matrix)
{

	//Check for errors and stack overflow
	if(!render || render->matlevel >= MATLEVELS)
		return GZ_FAILURE;
	
	//1 - Vertices' Transforms
	if(render->matlevel == -1)
	{
		render->matlevel++;
		memcpy(render->Ximage[render->matlevel], matrix, sizeof(GzMatrix));
	} else
	{
		//Matrix multiplication
		GzMatrix tmpMat;
		GzMatrixMult(render->Ximage[render->matlevel], matrix, tmpMat);
		//Increase Top of Stack and copy tmpMat into top position
		render->matlevel++;
		memcpy(render->Ximage[render->matlevel], tmpMat, sizeof(GzMatrix));
	}

	//2 - Normals' Transforms
	if(render->matlevel == 0 || render->matlevel == 1){
		//If renderer is pushing Xsp or Xpi, copy identity matrix to TOS
		memcpy(render->Xnorm[render->matlevel], identityMat, sizeof(GzMatrix));
	}else if(render->matlevel == 2){
		//If renderer is pushing Xiw, copy Xiw to TOS
		
		//Preprocessing input matrix
		// - Removing translation component
		GzMatrix processedMat;
		GzPreprocessMatrix(matrix, processedMat);
		memcpy(render->Xnorm[render->matlevel], processedMat, sizeof(GzMatrix));
	}else{
		//If renderer is pushing Xwm, push the matrix in the stack
		//Matrix multiplication
		render->matlevel--;
		GzMatrix tmpMat;
		//Preprocessing input matrix
		GzMatrix processedMat;
		GzPreprocessMatrix(matrix, processedMat);
		GzMatrixMult(render->Xnorm[render->matlevel], processedMat, tmpMat);

		//Increase Top of Stack and copy tmpMat into top position
		render->matlevel++;
		memcpy(render->Xnorm[render->matlevel], tmpMat, sizeof(GzMatrix));
	}

	return GZ_SUCCESS;
}

//Preprocess a matrix to transform the normals
void GzPreprocessMatrix(GzMatrix matrix, GzMatrix &res){
	memcpy(res, matrix, sizeof(GzMatrix));
	
	// - Removing translation component
	res[0][3] = 0.0f; res[1][3] = 0.0f; res[2][3] = 0.0f;

	float K = 1 / (sqrt(pow(res[0][0],2) + pow(res[0][1],2) + pow(res[0][2],2)));
	res[0][0] *= K; res[0][1] *= K; res[0][2] *= K;
	res[1][0] *= K; res[1][1] *= K; res[1][2] *= K;
	res[2][0] *= K; res[2][1] *= K; res[2][2] *= K;
}


//pop a matrix off the Ximage stack
int GzPopMatrix(GzRender *render)
{
	//Check for errors and stack underflow
	if(!render || render->matlevel < 0)
		return GZ_FAILURE;
	
	//Decrease Top of Stack
	render->matlevel--;
	return GZ_SUCCESS;
}


int GzPutAttribute(GzRender	*render, int numAttributes, GzToken	*nameList, 
	GzPointer	*valueList) /* void** valuelist */
{
	if(!render || numAttributes <= 0 || !nameList || !valueList)
		return GZ_FAILURE;

	int status = 0;
	for(unsigned int i = 0; i < numAttributes; ++i){
		switch(nameList[i]){
		case GZ_NULL_TOKEN: break;
		case GZ_RGB_COLOR:
			memcpy(render->flatcolor, valueList[i], sizeof(GzColor));
			break;
		case GZ_INTERPOLATE:
			render->interp_mode = *((int*)(valueList[i]));
			break;
		case GZ_DIRECTIONAL_LIGHT:
			//Checking for lights overflow
			if(render->numlights == MAX_LIGHTS){
				status = GZ_FAILURE;
				break;
			}
			memcpy(&(render->lights[render->numlights]), valueList[i], sizeof(GzLight));
			++(render->numlights);
			break;
		case GZ_AMBIENT_LIGHT:
			memcpy(&(render->ambientlight), valueList[i], sizeof(GzLight));
			break;
		case GZ_AMBIENT_COEFFICIENT:
			memcpy(render->Ka, valueList[i], sizeof(GzColor));
			break;
		case GZ_DIFFUSE_COEFFICIENT:
			memcpy(render->Kd, valueList[i], sizeof(GzColor));
			break;
		case GZ_SPECULAR_COEFFICIENT:
			memcpy(render->Ks, valueList[i], sizeof(GzColor));
			break;
		case GZ_DISTRIBUTION_COEFFICIENT:
			render->spec = *((float*)(valueList[i]));
			break;
		default: status = GZ_FAILURE;
		}
		if(status) return GZ_FAILURE;
	}

	return GZ_SUCCESS;
}

void GzTriangleNormal(GzCoord vertexList[NUM_VERTICES], GzCoord normal){
	GzCoord vY, vZ;
	float size;
	vY[X] = vertexList[1][X] - vertexList[0][X];
	vY[Y] = vertexList[1][Y] - vertexList[0][Y];
	vY[Z] = vertexList[1][Z] - vertexList[0][Z];
	//size = 1.0f / sqrt(pow(vY[X],2)+pow(vY[Y],2)+pow(vY[Z],2));
	//vY[X] *= size; vY[Y] *= size; vY[Z] *= size;
	
	vZ[X] = vertexList[2][X] - vertexList[1][X];
	vZ[Y] = vertexList[2][Y] - vertexList[1][Y];
	vZ[Z] = vertexList[2][Z] - vertexList[1][Z];
	//size = 1.0f / sqrt(pow(vZ[X],2)+pow(vZ[Y],2)+pow(vZ[Z],2));
	//vZ[X] *= size; vZ[Y] *= size; vZ[Z] *= size;

	//Computing normal
	normal[X] = vY[Y] * vZ[Z] - vY[Z] * vZ[Y];
	normal[Y] = vY[Z] * vZ[X] - vY[X] * vZ[Z];
	normal[Z] = vY[X] * vZ[Y] - vY[Y] * vZ[X];

	size = 1.0f / sqrt(normal[X]*normal[X] + normal[Y]*normal[Y] + normal[Z]*normal[Z]);
	normal[X] *= size; normal[Y] *= size; normal[Z] *= size;
}

int GzPutTriangle(GzRender	*render, int numParts, GzToken *nameList, GzPointer	*valueList)
{
/* numParts : how many names and values */

/* 
pass in a triangle description with tokens and values corresponding to
      GZ_NULL_TOKEN:		do nothing - no values
      GZ_POSITION:		3 vert positions 
	  GZ_NORMAL:		normal direction
*/
	if(!render || numParts <= 0 || !nameList || !valueList)
		return GZ_FAILURE;

	int status = GZ_SUCCESS;
	bool hasPosition = false, hasNormal = false;
	GzCoord vertexList[NUM_VERTICES];
	GzCoord normalList[NUM_VERTICES];
	GzCoord triNorm;
	//Loop on every part in the lists
	for(unsigned int p = 0; p < numParts; ++p){
		switch(nameList[p]){
		case GZ_NULL_TOKEN: break;
		case GZ_POSITION:
			//Recording vertices positions
			memcpy(vertexList, valueList[p], NUM_VERTICES*sizeof(GzCoord));
			//calculating triangle norm (for Flat Shading)
			GzTriangleNormal(vertexList, triNorm);
			for(unsigned int v = 0; v < NUM_VERTICES; ++v){
				//Xform positions of verts using matrix on top of stack 
				if(GzXformTriangle(render, &(vertexList[v])))
					return GZ_FAILURE;
				// Clip - just discard any triangle with any vert(s) behind view plane 
				if(vertexList[v][Z] < 0.0f)
					return GZ_SUCCESS;
			}
			//Test for triangles with all three verts off-screen (trivial frustum cull)
			bool culled;
			if(GzFrustumCull(render, vertexList, culled))
				return GZ_FAILURE;
			if(culled)
				return GZ_SUCCESS;
			
			hasPosition = true;
			break;
		case GZ_NORMAL:
			//Recording normals positions
			memcpy(normalList, valueList[p], NUM_VERTICES*sizeof(GzCoord));
			for(unsigned int v = 0; v < NUM_VERTICES; ++v){
				//Xform positions of verts using matrix on top of stack 
				if(GzXformNormal(render, &(normalList[v])))
					return GZ_FAILURE;
			}
			hasNormal = true;
			break;
		default: return GZ_FAILURE;
		}
	}

	//Switching on shader
	switch(render->interp_mode){
	case GZ_FLAT:
		//Flat Shading
		if(hasPosition){
			GzXformNormal(render, &(triNorm));
			GzShading(render, triNorm, render->flatcolor);
			return GzRasterizerFlatShading(render, vertexList);
		}
		else return GZ_FAILURE;
		break;
	case GZ_COLOR:
		//Gouraud Shading
		if(hasPosition && hasNormal) return GzRasterizerGouraudShading(render, vertexList, normalList);
		else return GZ_FAILURE;
		break;
	case GZ_NORMALS:
		//Phong Shading
		if(hasPosition && hasNormal) return GzRasterizerPhongShading(render, vertexList, normalList);
		else return GZ_FAILURE;
		break;
	default:
		return GZ_FAILURE;
		break;
	};

	return GZ_SUCCESS;
}

/* NOT part of API - just for general assistance */
short	ctoi(float color)		/* convert float color to GzIntensity short */
{
  return(short)((int)(color * ((1 << 12) - 1)));
}


int GzXformTriangle(GzRender *render, GzCoord *vertex)
{
	if(!render || !vertex)
		return GZ_FAILURE;

	//Getting Xsm
	GzMatrix Xsm;
	memcpy(Xsm, render->Ximage[render->matlevel], sizeof(GzMatrix));

	//Building 4D vector
	float modVert[4];
	memcpy(modVert, *vertex, sizeof(GzCoord));
	modVert[3] = 1.0f;

	//Applying transform
	float scrVert[4];
	GzMatVecMult(Xsm, modVert, scrVert);

	//Converting 4D vector to 3D
	for(unsigned int i = 0; i < 3; ++i){
		(*vertex)[i] = scrVert[i] / scrVert[3];
	}

	return GZ_SUCCESS;
}

int GzXformNormal(GzRender *render, GzCoord *normal)
{
	if(!render || !normal)
		return GZ_FAILURE;

	//Getting Xn
	GzMatrix Xn;
	memcpy(Xn, render->Xnorm[render->matlevel], sizeof(GzMatrix));

	//Building 4D vector
	float modNorm[4];
	memcpy(modNorm, *normal, sizeof(GzCoord));
	modNorm[3] = 1.0f;

	//Applying transform
	float scrNorm[4];
	GzMatVecMult(Xn, modNorm, scrNorm); 

	memcpy(*normal, scrNorm, sizeof(GzCoord));
	return GZ_SUCCESS;
}

/* Shader: Calculate color based on lights, camera position, and vertex normal */
void GzShading(GzRender *render, GzCoord normal, GzColor &color){
	GzCoord E = {0.0f, 0.0f, -1.0f};

	//Ambient light
	color[RED] = render->ambientlight.color[RED] * render->Ka[RED];
	color[GREEN] = render->ambientlight.color[GREEN] * render->Ka[GREEN];
	color[BLUE] = render->ambientlight.color[BLUE] * render->Ka[BLUE];
	
	//Lights
	for(unsigned int l = 0; l < render->numlights; ++l){
		GzCoord N; memcpy(N, normal, sizeof(GzCoord));
		GzCoord L; memcpy(L, render->lights[l].direction, sizeof(GzCoord));
		GzColor Ie; memcpy(Ie, render->lights[l].color, sizeof(GzColor));

		float NdotL = GzDotProduct(&N, &L);
		float NdotE = GzDotProduct(&N, &E);
		//Special cases
		//check if light and eye are on opposite sides of the surface
		if((NdotL > 0.0f && NdotE < 0.0f) || (NdotL < 0.0f && NdotE > 0.0f))
			continue;
		//check if we need to invert the normal
		if(NdotL < 0.0f && NdotE < 0.0f){
			N[X] *= -1.0f; N[Y] *= -1.0f; N[Z] *= -1.0f;
			NdotL *= -1;
		}

		//Calculating vector R
		GzCoord R;
		R[X] = 2 * NdotL * N[X] - L[X];
		R[Y] = 2 * NdotL * N[Y] - L[Y];
		R[Z] = 2 * NdotL * N[Z] - L[Z];

		//Specular component
		float RdotE = GzDotProduct(&R, &E);
		RdotE = min(max(RdotE, 0.0f), 1.0f);
		color[RED] += render->Ks[RED] * Ie[RED] * pow(RdotE, render->spec);
		color[GREEN] += render->Ks[GREEN] * Ie[GREEN] * pow(RdotE, render->spec);
		color[BLUE] += render->Ks[BLUE] * Ie[BLUE] * pow(RdotE, render->spec);

		//Diffuse component
		color[RED] += render->Kd[RED] * Ie[RED] * NdotL;
		color[GREEN] += render->Kd[GREEN] * Ie[GREEN] * NdotL;
		color[BLUE] += render->Kd[BLUE] * Ie[BLUE] * NdotL;

	}
	
	//Checking for overflow
	color[RED] = min(color[RED], 1.0f);
	color[GREEN] = min(color[GREEN], 1.0f);
	color[BLUE] = min(color[BLUE], 1.0f);
}

/////////
// RASTERIZER'S FUNCTIONS
/////////
/* Flat-Shading Scan-Line Rasterizer algorithm */
int GzRasterizerFlatShading(GzRender *render, GzCoord vertexList[NUM_VERTICES]){
	if(!render)
		return GZ_FAILURE;

	// 1 - Sorting vertices
	GzCoord tmpVertex;
	for(unsigned int i = 0; i < NUM_VERTICES-1; ++i){
		for(unsigned int j = i+1; j < NUM_VERTICES; ++j){
			if(vertexList[i][Y] > vertexList[j][Y]){
				memcpy(tmpVertex, vertexList[i], sizeof(GzCoord));
				memcpy(vertexList[i], vertexList[j], sizeof(GzCoord));
				memcpy(vertexList[j], tmpVertex, sizeof(GzCoord));
			}
		}
	}
	float Ymax = vertexList[NUM_VERTICES-1][Y];
	float Ymin = vertexList[0][Y];

	// 2 - Setup edge DDAs
	GzDDA_FS DDAList[NUM_EDGES];
	unsigned int num_h = 0, num_v = 0;
	GzDDA_FS* DDA_H = NULL; //Horizontal DDA
	bool DDA_H_isTop; //Horizontal DDA is a top edge
	GzDDA_FS* DDA_NC[NUM_EDGES]; //Unclassified DDA
	unsigned int DDA_NC_size = 0;
	{
		//Edge 0
		memcpy(DDAList[0].start, vertexList[0], sizeof(GzCoord));
		memcpy(DDAList[0].current, vertexList[0], sizeof(GzCoord));
		memcpy(DDAList[0].end, vertexList[1], sizeof(GzCoord));
		if(DDAList[0].end[Y] != DDAList[0].start[Y]){
			DDAList[0].slopeX = (DDAList[0].end[X] - DDAList[0].start[X])/(DDAList[0].end[Y] - DDAList[0].start[Y]);
			DDAList[0].slopeZ = (DDAList[0].end[Z] - DDAList[0].start[Z])/(DDAList[0].end[Y] - DDAList[0].start[Y]);
			DDA_NC[DDA_NC_size++] = &(DDAList[0]);
		}else{
			DDAList[0].slopeX = DDAList[0].slopeZ = 0;
			DDA_H = &DDAList[0];
			++num_h;
			DDA_H_isTop = true;
		}
		if(DDAList[0].end[X] == DDAList[0].start[X]) ++num_v;

		//Edge 1
		memcpy(DDAList[1].start, vertexList[0], sizeof(GzCoord));
		memcpy(DDAList[1].current, vertexList[0], sizeof(GzCoord));
		memcpy(DDAList[1].end, vertexList[2], sizeof(GzCoord));
		if(DDAList[1].end[Y] != DDAList[1].start[Y]){
			DDAList[1].slopeX = (DDAList[1].end[X] - DDAList[1].start[X])/(DDAList[1].end[Y] - DDAList[1].start[Y]);
			DDAList[1].slopeZ = (DDAList[1].end[Z] - DDAList[1].start[Z])/(DDAList[1].end[Y] - DDAList[1].start[Y]);
			DDA_NC[DDA_NC_size++] = &(DDAList[1]);
		}else{
			DDAList[1].slopeX = DDAList[1].slopeZ = 0;
			DDA_H = &DDAList[1];
			++num_h;
		}
		if(DDAList[1].end[X] == DDAList[1].start[X]) ++num_v;

		//Edge 2
		memcpy(DDAList[2].start, vertexList[1], sizeof(GzCoord));
		memcpy(DDAList[2].current, vertexList[1], sizeof(GzCoord));
		memcpy(DDAList[2].end, vertexList[2], sizeof(GzCoord));
		if(DDAList[2].end[Y] != DDAList[2].start[Y]){
			DDAList[2].slopeX = (DDAList[2].end[X] - DDAList[2].start[X])/(DDAList[2].end[Y] - DDAList[2].start[Y]);
			DDAList[2].slopeZ = (DDAList[2].end[Z] - DDAList[2].start[Z])/(DDAList[2].end[Y] - DDAList[2].start[Y]);
			DDA_NC[DDA_NC_size++] = &(DDAList[2]);
		}else{
			DDAList[2].slopeX = DDAList[2].slopeZ = 0;
			DDA_H = &DDAList[2];
			++num_h;
			DDA_H_isTop = false;
		}
		if(DDAList[2].end[X] == DDAList[2].start[X]) ++num_v;
	}
	//Check for ill-defined triangles
	if(num_v > 1 || num_h > 1)
		return GZ_SUCCESS;

	
	// 3 - Sort edges by L or R
	GzDDA_FS* DDA_L[NUM_EDGES-1]; //Left DDA
	GzDDA_FS* DDA_R[NUM_EDGES-1]; //Right DDA
	if(num_h == 0){
		if(DDA_NC[0]->slopeX < DDA_NC[1]->slopeX){
			DDA_L[0] = DDA_NC[0];
			DDA_R[0] = DDA_NC[1];
		}else{
			DDA_L[0] = DDA_NC[1];
			DDA_R[0] = DDA_NC[0];
		}

		if(DDA_L[0]->end[Y] < DDA_R[0]->end[Y]){
			DDA_L[1] = DDA_NC[2];
			DDA_R[1] = NULL;
		}else{
			DDA_R[1] = DDA_NC[2];
			DDA_L[1] = NULL;
		}
	}else{
		//Horizontal edge
		DDA_L[1] = NULL;
		DDA_R[1] = NULL;
		if(DDA_H_isTop){
			if(DDA_NC[0]->slopeX > DDA_NC[1]->slopeX){
				DDA_L[0] = DDA_NC[0];
				DDA_R[0] = DDA_NC[1];
			}else{
				DDA_L[0] = DDA_NC[1];
				DDA_R[0] = DDA_NC[0];
			}
		}else{
			if(DDA_NC[0]->slopeX < DDA_NC[1]->slopeX){
				DDA_L[0] = DDA_NC[0];
				DDA_R[0] = DDA_NC[1];
			}else{
				DDA_L[0] = DDA_NC[1];
				DDA_R[0] = DDA_NC[0];
			}
		}
	}

	//Rasterizing loop
	unsigned int pos_L = 0, pos_R = 0;
	float Ypos = ceil(Ymin);
	while(Ypos < Ymax){
		//Avoiding useless iterations
		if(Ypos < 0.0f)
			Ypos = 0.0f;
		if(Ypos > render->display->yres)
			break;

		//Jumping to the next edge when required
		if(Ypos > DDA_L[pos_L]->end[Y]){
			++pos_L;
			continue;
		}
		if(Ypos > DDA_R[pos_R]->end[Y]){
			++pos_R;
			continue;
		}

		// 5 - Advance DDAs current positions to y-scan line
		float deltaY;
		if(DDA_L[pos_L]->current[Y] != Ypos){
			deltaY = Ypos - DDA_L[pos_L]->current[Y];
			DDA_L[pos_L]->current[X] += deltaY * DDA_L[pos_L]->slopeX;
			DDA_L[pos_L]->current[Y] = Ypos;
			DDA_L[pos_L]->current[Z] += deltaY * DDA_L[pos_L]->slopeZ;
		}
		if(DDA_R[pos_R]->current[Y] != Ypos){
			deltaY = Ypos - DDA_R[pos_R]->current[Y];
			DDA_R[pos_R]->current[X] += deltaY * DDA_R[pos_R]->slopeX;
			DDA_R[pos_R]->current[Y] = Ypos;
			DDA_R[pos_R]->current[Z] += deltaY * DDA_R[pos_R]->slopeZ;
		}

		// 7 - Set span DDA current and end positions to right and left edge values
		GzSpan_FS span;
		memcpy(span.start, DDA_L[pos_L]->current, sizeof(GzCoord));
		memcpy(span.end, DDA_R[pos_R]->current, sizeof(GzCoord));
		memcpy(span.current, DDA_L[pos_L]->current, sizeof(GzCoord));
		span.slopeZ = (DDA_R[pos_R]->current[Z] - DDA_L[pos_L]->current[Z])/(DDA_R[pos_R]->current[X] - DDA_L[pos_L]->current[X]);

		// 8 - Advance span current position to left-most covered pixel (ceiling)
		span.current[X] = (ceil(span.start[X]) < 0.0f)?0.0f:ceil(span.start[X]);
		float deltaX = span.current[X] - span.start[X];
		span.current[Z] += deltaX * span.slopeZ;

		//NOT DRAWING pixels covered by bottom horizontal edges
		if(!(num_h == 1 && !DDA_H_isTop && Ypos == Ymax)){
			// 9 - Interpolate span position and parameters (Z) until current position > end
			//NOT DRAWING pixels covered by right edges
			while(span.current[X] < span.end[X] && span.current[X] <= render->display->xres){
				if(GzPixelToDisplay(&(span.current), render->display, &(render->flatcolor)))
					return GZ_FAILURE;
				span.current[X] += 1.0f;
				span.current[Z] += span.slopeZ;		
			}
		}
		//Advance to the next scan line
		Ypos += 1.0f;
	} //END (Ypos < Ymax)

	return GZ_SUCCESS;
}

/* Gouraud-Shading Scan-Line Rasterizer algorithm */
int GzRasterizerGouraudShading(GzRender *render, GzCoord vertexList[NUM_VERTICES], GzCoord normalList[NUM_VERTICES]){

	if(!render)
		return GZ_FAILURE;
	
	// Gouraud-Shading Calculating color at vertices
	GzColor colorList[NUM_VERTICES];
	for(unsigned int i = 0; i < NUM_VERTICES; ++i){
		GzShading(render, normalList[i], colorList[i]);
	}

	// 1 - Sorting vertices
	GzCoord tmpVertex;
	GzColor tmpColor;
	for(unsigned int i = 0; i < NUM_VERTICES-1; ++i){
		for(unsigned int j = i+1; j < NUM_VERTICES; ++j){
			if(vertexList[i][Y] > vertexList[j][Y]){
				memcpy(tmpVertex, vertexList[i], sizeof(GzCoord));
				memcpy(vertexList[i], vertexList[j], sizeof(GzCoord));
				memcpy(vertexList[j], tmpVertex, sizeof(GzCoord));

				memcpy(tmpColor, colorList[i], sizeof(GzCoord));
				memcpy(colorList[i], colorList[j], sizeof(GzCoord));
				memcpy(colorList[j], tmpColor, sizeof(GzCoord));
			}
		}
	}
	float Ymax = vertexList[NUM_VERTICES-1][Y];
	float Ymin = vertexList[0][Y];

	// 2 - Setup edge DDAs
	GzDDA_GS DDAList[NUM_EDGES];
	unsigned int num_h = 0, num_v = 0;
	GzDDA_GS* DDA_H = NULL; //Horizontal DDA
	bool DDA_H_isTop; //Horizontal DDA is a top edge
	GzDDA_GS* DDA_NC[NUM_EDGES]; //Unclassified DDA
	unsigned int DDA_NC_size = 0;
	{
		//Edge 0
		memcpy(DDAList[0].start, vertexList[0], sizeof(GzCoord));
		memcpy(DDAList[0].current, vertexList[0], sizeof(GzCoord));
		memcpy(DDAList[0].end, vertexList[1], sizeof(GzCoord));
		memcpy(DDAList[0].startColor, colorList[0], sizeof(GzColor));
		memcpy(DDAList[0].currentColor, colorList[0], sizeof(GzColor));
		memcpy(DDAList[0].endColor, colorList[1], sizeof(GzColor));
		if(DDAList[0].end[Y] != DDAList[0].start[Y]){
			DDAList[0].slopeX = (DDAList[0].end[X] - DDAList[0].start[X])/(DDAList[0].end[Y] - DDAList[0].start[Y]);
			DDAList[0].slopeZ = (DDAList[0].end[Z] - DDAList[0].start[Z])/(DDAList[0].end[Y] - DDAList[0].start[Y]);
			DDAList[0].slopeColor[RED] = (DDAList[0].endColor[RED] - DDAList[0].startColor[RED])/(DDAList[0].end[Y] - DDAList[0].start[Y]);
			DDAList[0].slopeColor[GREEN] = (DDAList[0].endColor[GREEN] - DDAList[0].startColor[GREEN])/(DDAList[0].end[Y] - DDAList[0].start[Y]);
			DDAList[0].slopeColor[BLUE] = (DDAList[0].endColor[BLUE] - DDAList[0].startColor[BLUE])/(DDAList[0].end[Y] - DDAList[0].start[Y]);
			DDA_NC[DDA_NC_size++] = &(DDAList[0]);
		}else{
			DDAList[0].slopeX = DDAList[0].slopeZ = DDAList[0].slopeColor[RED] = DDAList[0].slopeColor[GREEN] = DDAList[0].slopeColor[BLUE] = 0;
			DDA_H = &DDAList[0];
			++num_h;
			DDA_H_isTop = true;
		}
		if(DDAList[0].end[X] == DDAList[0].start[X]) ++num_v;

		//Edge 1
		memcpy(DDAList[1].start, vertexList[0], sizeof(GzCoord));
		memcpy(DDAList[1].current, vertexList[0], sizeof(GzCoord));
		memcpy(DDAList[1].end, vertexList[2], sizeof(GzCoord));
		memcpy(DDAList[1].startColor, colorList[0], sizeof(GzColor));
		memcpy(DDAList[1].currentColor, colorList[0], sizeof(GzColor));
		memcpy(DDAList[1].endColor, colorList[2], sizeof(GzColor));
		if(DDAList[1].end[Y] != DDAList[1].start[Y]){
			DDAList[1].slopeX = (DDAList[1].end[X] - DDAList[1].start[X])/(DDAList[1].end[Y] - DDAList[1].start[Y]);
			DDAList[1].slopeZ = (DDAList[1].end[Z] - DDAList[1].start[Z])/(DDAList[1].end[Y] - DDAList[1].start[Y]);
			DDAList[1].slopeColor[RED] = (DDAList[1].endColor[RED] - DDAList[1].startColor[RED])/(DDAList[1].end[Y] - DDAList[1].start[Y]);
			DDAList[1].slopeColor[GREEN] = (DDAList[1].endColor[GREEN] - DDAList[1].startColor[GREEN])/(DDAList[1].end[Y] - DDAList[1].start[Y]);
			DDAList[1].slopeColor[BLUE] = (DDAList[1].endColor[BLUE] - DDAList[1].startColor[BLUE])/(DDAList[1].end[Y] - DDAList[1].start[Y]);
			DDA_NC[DDA_NC_size++] = &(DDAList[1]);
		}else{
			DDAList[1].slopeX = DDAList[1].slopeZ = DDAList[1].slopeColor[RED] = DDAList[1].slopeColor[GREEN] = DDAList[1].slopeColor[BLUE] = 0;
			DDA_H = &DDAList[1];
			++num_h;
		}
		if(DDAList[1].end[X] == DDAList[1].start[X]) ++num_v;

		//Edge 2
		memcpy(DDAList[2].start, vertexList[1], sizeof(GzCoord));
		memcpy(DDAList[2].current, vertexList[1], sizeof(GzCoord));
		memcpy(DDAList[2].end, vertexList[2], sizeof(GzCoord));
		memcpy(DDAList[2].startColor, colorList[1], sizeof(GzColor));
		memcpy(DDAList[2].currentColor, colorList[1], sizeof(GzColor));
		memcpy(DDAList[2].endColor, colorList[2], sizeof(GzColor));
		if(DDAList[2].end[Y] != DDAList[2].start[Y]){
			DDAList[2].slopeX = (DDAList[2].end[X] - DDAList[2].start[X])/(DDAList[2].end[Y] - DDAList[2].start[Y]);
			DDAList[2].slopeZ = (DDAList[2].end[Z] - DDAList[2].start[Z])/(DDAList[2].end[Y] - DDAList[2].start[Y]);
			DDAList[2].slopeColor[RED] = (DDAList[2].endColor[RED] - DDAList[2].startColor[RED])/(DDAList[2].end[Y] - DDAList[2].start[Y]);
			DDAList[2].slopeColor[GREEN] = (DDAList[2].endColor[GREEN] - DDAList[2].startColor[GREEN])/(DDAList[2].end[Y] - DDAList[2].start[Y]);
			DDAList[2].slopeColor[BLUE] = (DDAList[2].endColor[BLUE] - DDAList[2].startColor[BLUE])/(DDAList[2].end[Y] - DDAList[2].start[Y]);
			DDA_NC[DDA_NC_size++] = &(DDAList[2]);
		}else{
			DDAList[2].slopeX = DDAList[2].slopeZ  = DDAList[2].slopeColor[RED] = DDAList[2].slopeColor[GREEN] = DDAList[2].slopeColor[BLUE] = 0;
			DDA_H = &DDAList[2];
			++num_h;
			DDA_H_isTop = false;
		}
		if(DDAList[2].end[X] == DDAList[2].start[X]) ++num_v;
	}
	//Check for ill-defined triangles
	if(num_v > 1 || num_h > 1)
		return GZ_SUCCESS;

	
	// 3 - Sort edges by L or R
	GzDDA_GS* DDA_L[NUM_EDGES-1]; //Left DDA
	GzDDA_GS* DDA_R[NUM_EDGES-1]; //Right DDA
	if(num_h == 0){
		if(DDA_NC[0]->slopeX < DDA_NC[1]->slopeX){
			DDA_L[0] = DDA_NC[0];
			DDA_R[0] = DDA_NC[1];
		}else{
			DDA_L[0] = DDA_NC[1];
			DDA_R[0] = DDA_NC[0];
		}

		if(DDA_L[0]->end[Y] < DDA_R[0]->end[Y]){
			DDA_L[1] = DDA_NC[2];
			DDA_R[1] = NULL;
		}else{
			DDA_R[1] = DDA_NC[2];
			DDA_L[1] = NULL;
		}
	}else{
		//Horizontal edge
		DDA_L[1] = NULL;
		DDA_R[1] = NULL;
		if(DDA_H_isTop){
			if(DDA_NC[0]->slopeX > DDA_NC[1]->slopeX){
				DDA_L[0] = DDA_NC[0];
				DDA_R[0] = DDA_NC[1];
			}else{
				DDA_L[0] = DDA_NC[1];
				DDA_R[0] = DDA_NC[0];
			}
		}else{
			if(DDA_NC[0]->slopeX < DDA_NC[1]->slopeX){
				DDA_L[0] = DDA_NC[0];
				DDA_R[0] = DDA_NC[1];
			}else{
				DDA_L[0] = DDA_NC[1];
				DDA_R[0] = DDA_NC[0];
			}
		}
	}

	//Rasterizing loop
	unsigned int pos_L = 0, pos_R = 0;
	float Ypos = ceil(Ymin);
	while(Ypos < Ymax){
		//Avoiding useless iterations
		if(Ypos < 0.0f)
			Ypos = 0.0f;
		if(Ypos > render->display->yres)
			break;

		//Jumping to the next edge when required
		if(Ypos > DDA_L[pos_L]->end[Y]){
			++pos_L;
			continue;
		}
		if(Ypos > DDA_R[pos_R]->end[Y]){
			++pos_R;
			continue;
		}

		// 5 - Advance DDAs current positions to y-scan line
		float deltaY;
		if(DDA_L[pos_L]->current[Y] != Ypos){
			deltaY = Ypos - DDA_L[pos_L]->current[Y];
			DDA_L[pos_L]->current[X] += deltaY * DDA_L[pos_L]->slopeX;
			DDA_L[pos_L]->current[Y] = Ypos;
			DDA_L[pos_L]->current[Z] += deltaY * DDA_L[pos_L]->slopeZ;
			DDA_L[pos_L]->currentColor[RED] += deltaY * DDA_L[pos_L]->slopeColor[RED];
			DDA_L[pos_L]->currentColor[GREEN] += deltaY * DDA_L[pos_L]->slopeColor[GREEN];
			DDA_L[pos_L]->currentColor[BLUE] += deltaY * DDA_L[pos_L]->slopeColor[BLUE];
		}
		if(DDA_R[pos_R]->current[Y] != Ypos){
			deltaY = Ypos - DDA_R[pos_R]->current[Y];
			DDA_R[pos_R]->current[X] += deltaY * DDA_R[pos_R]->slopeX;
			DDA_R[pos_R]->current[Y] = Ypos;
			DDA_R[pos_R]->current[Z] += deltaY * DDA_R[pos_R]->slopeZ;
			DDA_R[pos_R]->currentColor[RED] += deltaY * DDA_R[pos_R]->slopeColor[RED];
			DDA_R[pos_R]->currentColor[GREEN] += deltaY * DDA_R[pos_R]->slopeColor[GREEN];
			DDA_R[pos_R]->currentColor[BLUE] += deltaY * DDA_R[pos_R]->slopeColor[BLUE];
		}

		// 7 - Set span DDA current and end positions to right and left edge values
		GzSpan_GS span;
		memcpy(span.start, DDA_L[pos_L]->current, sizeof(GzCoord));
		memcpy(span.end, DDA_R[pos_R]->current, sizeof(GzCoord));
		memcpy(span.current, DDA_L[pos_L]->current, sizeof(GzCoord));
		memcpy(span.startColor, DDA_L[pos_L]->currentColor, sizeof(GzCoord));
		memcpy(span.endColor, DDA_R[pos_R]->currentColor, sizeof(GzCoord));
		memcpy(span.currentColor, DDA_L[pos_L]->currentColor, sizeof(GzCoord));
		span.slopeZ = (DDA_R[pos_R]->current[Z] - DDA_L[pos_L]->current[Z])/(DDA_R[pos_R]->current[X] - DDA_L[pos_L]->current[X]);
		span.slopeColor[RED] = (DDA_R[pos_R]->currentColor[RED] - DDA_L[pos_L]->currentColor[RED])/(DDA_R[pos_R]->current[X] - DDA_L[pos_L]->current[X]);
		span.slopeColor[GREEN] = (DDA_R[pos_R]->currentColor[GREEN] - DDA_L[pos_L]->currentColor[GREEN])/(DDA_R[pos_R]->current[X] - DDA_L[pos_L]->current[X]);
		span.slopeColor[BLUE] = (DDA_R[pos_R]->currentColor[BLUE] - DDA_L[pos_L]->currentColor[BLUE])/(DDA_R[pos_R]->current[X] - DDA_L[pos_L]->current[X]);

		// 8 - Advance span current position to left-most covered pixel (ceiling)
		span.current[X] = (ceil(span.start[X]) < 0.0f)?0.0f:ceil(span.start[X]);
		float deltaX = span.current[X] - span.start[X];
		span.current[Z] += deltaX * span.slopeZ;
		span.currentColor[RED] += deltaX * span.slopeColor[RED];
		span.currentColor[GREEN] += deltaX * span.slopeColor[GREEN];
		span.currentColor[BLUE] += deltaX * span.slopeColor[BLUE];

		//NOT DRAWING pixels covered by bottom horizontal edges
		if(!(num_h == 1 && !DDA_H_isTop && Ypos == Ymax)){
			// 9 - Interpolate span position and parameters (Z) until current position > end
			//NOT DRAWING pixels covered by right edges
			while(span.current[X] < span.end[X] && span.current[X] <= render->display->xres){
				if(GzPixelToDisplay(&(span.current), render->display, &(span.currentColor)))
					return GZ_FAILURE;
				span.current[X] += 1.0f;
				span.current[Z] += span.slopeZ;
				span.currentColor[RED] += span.slopeColor[RED];
				span.currentColor[GREEN] += span.slopeColor[GREEN];
				span.currentColor[BLUE] += span.slopeColor[BLUE];
			}
		}
		//Advance to the next scan line
		Ypos += 1.0f;
	} //END (Ypos < Ymax)

	return GZ_SUCCESS;
}

/* Phong-Shading Scan-Line Rasterizer algorithm */
int GzRasterizerPhongShading(GzRender *render, GzCoord vertexList[NUM_VERTICES], GzCoord normalList[NUM_VERTICES]){
	if(!render)
		return GZ_FAILURE;

	// 1 - Sorting vertices
	GzCoord tmpVertex;
	GzColor tmpNormal;
	for(unsigned int i = 0; i < NUM_VERTICES-1; ++i){
		for(unsigned int j = i+1; j < NUM_VERTICES; ++j){
			if(vertexList[i][Y] > vertexList[j][Y]){
				memcpy(tmpVertex, vertexList[i], sizeof(GzCoord));
				memcpy(vertexList[i], vertexList[j], sizeof(GzCoord));
				memcpy(vertexList[j], tmpVertex, sizeof(GzCoord));

				memcpy(tmpNormal, normalList[i], sizeof(GzCoord));
				memcpy(normalList[i], normalList[j], sizeof(GzCoord));
				memcpy(normalList[j], tmpNormal, sizeof(GzCoord));
			}
		}
	}
	float Ymax = vertexList[NUM_VERTICES-1][Y];
	float Ymin = vertexList[0][Y];

	// 2 - Setup edge DDAs
	GzDDA_PS DDAList[NUM_EDGES];
	unsigned int num_h = 0, num_v = 0;
	GzDDA_PS* DDA_H = NULL; //Horizontal DDA
	bool DDA_H_isTop; //Horizontal DDA is a top edge
	GzDDA_PS* DDA_NC[NUM_EDGES]; //Unclassified DDA
	unsigned int DDA_NC_size = 0;
	{
		//Edge 0
		memcpy(DDAList[0].start, vertexList[0], sizeof(GzCoord));
		memcpy(DDAList[0].current, vertexList[0], sizeof(GzCoord));
		memcpy(DDAList[0].end, vertexList[1], sizeof(GzCoord));
		memcpy(DDAList[0].startNormal, normalList[0], sizeof(GzColor));
		memcpy(DDAList[0].currentNormal, normalList[0], sizeof(GzColor));
		memcpy(DDAList[0].endNormal, normalList[1], sizeof(GzColor));
		if(DDAList[0].end[Y] != DDAList[0].start[Y]){
			DDAList[0].slopeX = (DDAList[0].end[X] - DDAList[0].start[X])/(DDAList[0].end[Y] - DDAList[0].start[Y]);
			DDAList[0].slopeZ = (DDAList[0].end[Z] - DDAList[0].start[Z])/(DDAList[0].end[Y] - DDAList[0].start[Y]);
			DDAList[0].slopeNormal[X] = (DDAList[0].endNormal[X] - DDAList[0].startNormal[X])/(DDAList[0].end[Y] - DDAList[0].start[Y]);
			DDAList[0].slopeNormal[Y] = (DDAList[0].endNormal[Y] - DDAList[0].startNormal[Y])/(DDAList[0].end[Y] - DDAList[0].start[Y]);
			DDAList[0].slopeNormal[Z] = (DDAList[0].endNormal[Z] - DDAList[0].startNormal[Z])/(DDAList[0].end[Y] - DDAList[0].start[Y]);
			DDA_NC[DDA_NC_size++] = &(DDAList[0]);
		}else{
			DDAList[0].slopeX = DDAList[0].slopeZ = DDAList[0].slopeNormal[X] = DDAList[0].slopeNormal[Y] = DDAList[0].slopeNormal[Z] = 0;
			DDA_H = &DDAList[0];
			++num_h;
			DDA_H_isTop = true;
		}
		if(DDAList[0].end[X] == DDAList[0].start[X]) ++num_v;

		//Edge 1
		memcpy(DDAList[1].start, vertexList[0], sizeof(GzCoord));
		memcpy(DDAList[1].current, vertexList[0], sizeof(GzCoord));
		memcpy(DDAList[1].end, vertexList[2], sizeof(GzCoord));
		memcpy(DDAList[1].startNormal, normalList[0], sizeof(GzColor));
		memcpy(DDAList[1].currentNormal, normalList[0], sizeof(GzColor));
		memcpy(DDAList[1].endNormal, normalList[2], sizeof(GzColor));
		if(DDAList[1].end[Y] != DDAList[1].start[Y]){
			DDAList[1].slopeX = (DDAList[1].end[X] - DDAList[1].start[X])/(DDAList[1].end[Y] - DDAList[1].start[Y]);
			DDAList[1].slopeZ = (DDAList[1].end[Z] - DDAList[1].start[Z])/(DDAList[1].end[Y] - DDAList[1].start[Y]);
			DDAList[1].slopeNormal[X] = (DDAList[1].endNormal[X] - DDAList[1].startNormal[X])/(DDAList[1].end[Y] - DDAList[1].start[Y]);
			DDAList[1].slopeNormal[Y] = (DDAList[1].endNormal[Y] - DDAList[1].startNormal[Y])/(DDAList[1].end[Y] - DDAList[1].start[Y]);
			DDAList[1].slopeNormal[Z] = (DDAList[1].endNormal[Z] - DDAList[1].startNormal[Z])/(DDAList[1].end[Y] - DDAList[1].start[Y]);
			DDA_NC[DDA_NC_size++] = &(DDAList[1]);
		}else{
			DDAList[1].slopeX = DDAList[1].slopeZ = DDAList[1].slopeNormal[X] = DDAList[1].slopeNormal[Y] = DDAList[1].slopeNormal[Z] = 0;
			DDA_H = &DDAList[1];
			++num_h;
		}
		if(DDAList[1].end[X] == DDAList[1].start[X]) ++num_v;

		//Edge 2
		memcpy(DDAList[2].start, vertexList[1], sizeof(GzCoord));
		memcpy(DDAList[2].current, vertexList[1], sizeof(GzCoord));
		memcpy(DDAList[2].end, vertexList[2], sizeof(GzCoord));
		memcpy(DDAList[2].startNormal, normalList[1], sizeof(GzColor));
		memcpy(DDAList[2].currentNormal, normalList[1], sizeof(GzColor));
		memcpy(DDAList[2].endNormal, normalList[2], sizeof(GzColor));
		if(DDAList[2].end[Y] != DDAList[2].start[Y]){
			DDAList[2].slopeX = (DDAList[2].end[X] - DDAList[2].start[X])/(DDAList[2].end[Y] - DDAList[2].start[Y]);
			DDAList[2].slopeZ = (DDAList[2].end[Z] - DDAList[2].start[Z])/(DDAList[2].end[Y] - DDAList[2].start[Y]);
			DDAList[2].slopeNormal[X] = (DDAList[2].endNormal[X] - DDAList[2].startNormal[X])/(DDAList[2].end[Y] - DDAList[2].start[Y]);
			DDAList[2].slopeNormal[Y] = (DDAList[2].endNormal[Y] - DDAList[2].startNormal[Y])/(DDAList[2].end[Y] - DDAList[2].start[Y]);
			DDAList[2].slopeNormal[Z] = (DDAList[2].endNormal[Z] - DDAList[2].startNormal[Z])/(DDAList[2].end[Y] - DDAList[2].start[Y]);
			DDA_NC[DDA_NC_size++] = &(DDAList[2]);
		}else{
			DDAList[2].slopeX = DDAList[2].slopeZ  = DDAList[2].slopeNormal[X] = DDAList[2].slopeNormal[Y] = DDAList[2].slopeNormal[Z] = 0;
			DDA_H = &DDAList[2];
			++num_h;
			DDA_H_isTop = false;
		}
		if(DDAList[2].end[X] == DDAList[2].start[X]) ++num_v;
	}
	//Check for ill-defined triangles
	if(num_v > 1 || num_h > 1)
		return GZ_SUCCESS;

	
	// 3 - Sort edges by L or R
	GzDDA_PS* DDA_L[NUM_EDGES-1]; //Left DDA
	GzDDA_PS* DDA_R[NUM_EDGES-1]; //Right DDA
	if(num_h == 0){
		if(DDA_NC[0]->slopeX < DDA_NC[1]->slopeX){
			DDA_L[0] = DDA_NC[0];
			DDA_R[0] = DDA_NC[1];
		}else{
			DDA_L[0] = DDA_NC[1];
			DDA_R[0] = DDA_NC[0];
		}

		if(DDA_L[0]->end[Y] < DDA_R[0]->end[Y]){
			DDA_L[1] = DDA_NC[2];
			DDA_R[1] = NULL;
		}else{
			DDA_R[1] = DDA_NC[2];
			DDA_L[1] = NULL;
		}
	}else{
		//Horizontal edge
		DDA_L[1] = NULL;
		DDA_R[1] = NULL;
		if(DDA_H_isTop){
			if(DDA_NC[0]->slopeX > DDA_NC[1]->slopeX){
				DDA_L[0] = DDA_NC[0];
				DDA_R[0] = DDA_NC[1];
			}else{
				DDA_L[0] = DDA_NC[1];
				DDA_R[0] = DDA_NC[0];
			}
		}else{
			if(DDA_NC[0]->slopeX < DDA_NC[1]->slopeX){
				DDA_L[0] = DDA_NC[0];
				DDA_R[0] = DDA_NC[1];
			}else{
				DDA_L[0] = DDA_NC[1];
				DDA_R[0] = DDA_NC[0];
			}
		}
	}

	//Rasterizing loop
	unsigned int pos_L = 0, pos_R = 0;
	float Ypos = ceil(Ymin);
	while(Ypos < Ymax){
		//Avoiding useless iterations
		if(Ypos < 0.0f)
			Ypos = 0.0f;
		if(Ypos > render->display->yres)
			break;

		//Jumping to the next edge when required
		if(Ypos > DDA_L[pos_L]->end[Y]){
			++pos_L;
			continue;
		}
		if(Ypos > DDA_R[pos_R]->end[Y]){
			++pos_R;
			continue;
		}

		// 5 - Advance DDAs current positions to y-scan line
		float deltaY;
		if(DDA_L[pos_L]->current[Y] != Ypos){
			deltaY = Ypos - DDA_L[pos_L]->current[Y];
			DDA_L[pos_L]->current[X] += deltaY * DDA_L[pos_L]->slopeX;
			DDA_L[pos_L]->current[Y] = Ypos;
			DDA_L[pos_L]->current[Z] += deltaY * DDA_L[pos_L]->slopeZ;
			DDA_L[pos_L]->currentNormal[X] += deltaY * DDA_L[pos_L]->slopeNormal[X];
			DDA_L[pos_L]->currentNormal[Y] += deltaY * DDA_L[pos_L]->slopeNormal[Y];
			DDA_L[pos_L]->currentNormal[Z] += deltaY * DDA_L[pos_L]->slopeNormal[Z];
		}
		if(DDA_R[pos_R]->current[Y] != Ypos){
			deltaY = Ypos - DDA_R[pos_R]->current[Y];
			DDA_R[pos_R]->current[X] += deltaY * DDA_R[pos_R]->slopeX;
			DDA_R[pos_R]->current[Y] = Ypos;
			DDA_R[pos_R]->current[Z] += deltaY * DDA_R[pos_R]->slopeZ;
			DDA_R[pos_R]->currentNormal[X] += deltaY * DDA_R[pos_R]->slopeNormal[X];
			DDA_R[pos_R]->currentNormal[Y] += deltaY * DDA_R[pos_R]->slopeNormal[Y];
			DDA_R[pos_R]->currentNormal[Z] += deltaY * DDA_R[pos_R]->slopeNormal[Z];
		}

		// 7 - Set span DDA current and end positions to right and left edge values
		GzSpan_PS span;
		memcpy(span.start, DDA_L[pos_L]->current, sizeof(GzCoord));
		memcpy(span.end, DDA_R[pos_R]->current, sizeof(GzCoord));
		memcpy(span.current, DDA_L[pos_L]->current, sizeof(GzCoord));
		memcpy(span.startNormal, DDA_L[pos_L]->currentNormal, sizeof(GzCoord));
		memcpy(span.endNormal, DDA_R[pos_R]->currentNormal, sizeof(GzCoord));
		memcpy(span.currentNormal, DDA_L[pos_L]->currentNormal, sizeof(GzCoord));
		span.slopeZ = (DDA_R[pos_R]->current[Z] - DDA_L[pos_L]->current[Z])/(DDA_R[pos_R]->current[X] - DDA_L[pos_L]->current[X]);
		span.slopeNormal[X] = (DDA_R[pos_R]->currentNormal[X] - DDA_L[pos_L]->currentNormal[X])/(DDA_R[pos_R]->current[X] - DDA_L[pos_L]->current[X]);
		span.slopeNormal[Y] = (DDA_R[pos_R]->currentNormal[Y] - DDA_L[pos_L]->currentNormal[Y])/(DDA_R[pos_R]->current[X] - DDA_L[pos_L]->current[X]);
		span.slopeNormal[Z] = (DDA_R[pos_R]->currentNormal[Z] - DDA_L[pos_L]->currentNormal[Z])/(DDA_R[pos_R]->current[X] - DDA_L[pos_L]->current[X]);

		// 8 - Advance span current position to left-most covered pixel (ceiling)
		span.current[X] = (ceil(span.start[X]) < 0.0f)?0.0f:ceil(span.start[X]);
		float deltaX = span.current[X] - span.start[X];
		span.current[Z] += deltaX * span.slopeZ;
		span.currentNormal[X] += deltaX * span.slopeNormal[X];
		span.currentNormal[Y] += deltaX * span.slopeNormal[Y];
		span.currentNormal[Z] += deltaX * span.slopeNormal[Z];

		//NOT DRAWING pixels covered by bottom horizontal edges
		if(!(num_h == 1 && !DDA_H_isTop && Ypos == Ymax)){
			// 9 - Interpolate span position and parameters (Z) until current position > end
			//NOT DRAWING pixels covered by right edges
			while(span.current[X] < span.end[X] && span.current[X] <= render->display->xres){
				GzColor currentColor;
				GzCoord currentNormal;
				memcpy(currentNormal, span.currentNormal, sizeof(GzCoord));
				float K = 1 / GzVectorNorm(&currentNormal);
				currentNormal[X] *= K; currentNormal[Y] *= K; currentNormal[Z] *= K;
				GzShading(render, currentNormal, currentColor);
				if(GzPixelToDisplay(&(span.current), render->display, &(currentColor)))
					return GZ_FAILURE;
				span.current[X] += 1.0f;
				span.current[Z] += span.slopeZ;
				span.currentNormal[X] += span.slopeNormal[X];
				span.currentNormal[Y] += span.slopeNormal[Y];
				span.currentNormal[Z] += span.slopeNormal[Z];
			}
		}
		//Advance to the next scan line
		Ypos += 1.0f;
	} //END (Ypos < Ymax)

	return GZ_SUCCESS;
}


int GzPixelToDisplay(GzCoord *point, GzDisplay *display, GzColor *flatcolor){
	if(!point || !display || !flatcolor)
		return GZ_FAILURE;

	if( (*point)[Z] <= 0 || 
		(*point)[X] < 0 || (*point)[X] >= display->xres ||
		(*point)[Y] < 0 || (*point)[Y] >= display->yres)
		return GZ_SUCCESS; //clipping pixel

	int pos = ARRAY((int) (*point)[X], (int) (*point)[Y]);

	//Test interpolated-Z against FB-Z for each pixel - low Z wins
	if((*point)[Z] < display->fbuf[pos].z)
		//Write color value into FB pixel  (default or computed color)
		GzPutDisplay(display, (int) (*point)[X], (int) (*point)[Y], 
			ctoi((*flatcolor)[RED]), ctoi((*flatcolor)[GREEN]), ctoi((*flatcolor)[BLUE]), 
			MAX_INTENSITY, (int) (*point)[Z]);

	return GZ_SUCCESS;
}

int GzFrustumCull(GzRender *render, GzCoord vertexList[NUM_VERTICES], bool &culled){
	if(!render)
		return GZ_FAILURE;
	
	culled = false;
	
	if(vertexList[0][X] < 0.0f && vertexList[1][X] < 0.0f && vertexList[2][X] < 0.0F){
		culled = true;
		return GZ_SUCCESS;
	}

	if(vertexList[0][X] > render->display->xres && vertexList[1][X] > render->display->xres && vertexList[2][X] > render->display->xres){
		culled = true;
		return GZ_SUCCESS;
	}

	if(vertexList[0][Y] < 0.0f && vertexList[1][Y] < 0.0f && vertexList[2][Y] < 0.0F){
		culled = true;
		return GZ_SUCCESS;
	}

	if(vertexList[0][Y] > render->display->yres && vertexList[1][Y] > render->display->yres && vertexList[2][Y] > render->display->yres){
		culled = true;
		return GZ_SUCCESS;
	}
	
	return GZ_SUCCESS;
}

/////////
//TRANSFORMS' FUNCTIONS
/////////

//Setup Xsp transform
int GzSetupXsp(GzRender *render)
{
	if(!render || !render->display)
		return GZ_FAILURE;
	
	//Populating by row
	render->Xsp[0][0] = render->display->xres/2.0f;
	render->Xsp[0][1] = 0.0f;
	render->Xsp[0][2] = 0.0f;
	render->Xsp[0][3] = render->display->xres/2.0f;
	
	render->Xsp[1][0] = 0.0f;
	render->Xsp[1][1] = -render->display->yres/2.0f;
	render->Xsp[1][2] = 0.0f;
	render->Xsp[1][3] = render->display->yres/2.0f;

	render->Xsp[2][0] = 0.0f;
	render->Xsp[2][1] = 0.0f;
	render->Xsp[2][2] =(float) MAXINT;
	render->Xsp[2][3] = 0.0f;
	
	render->Xsp[3][0] = 0.0f;
	render->Xsp[3][1] = 0.0f;
	render->Xsp[3][2] = 0.0f;
	render->Xsp[3][3] = 1.0f;

	return GZ_SUCCESS;
}

int GzSetupXpi(GzRender *render)
{
	if(!render)
		return GZ_FAILURE;
	GzCamera *pCam = &(render->camera);
	float dInv = tan((pCam->FOV/2.0f) * M_PI / 180.0f);
	
	//Populating the matrix by row
	pCam->Xpi[0][0] = 1.0f;
	pCam->Xpi[0][1] = 0.0f;
	pCam->Xpi[0][2] = 0.0f;
	pCam->Xpi[0][3] = 0.0f;

	pCam->Xpi[1][0] = 0.0f;
	pCam->Xpi[1][1] = 1.0f;
	pCam->Xpi[1][2] = 0.0f;
	pCam->Xpi[1][3] = 0.0f;
	
	pCam->Xpi[2][0] = 0.0f;
	pCam->Xpi[2][1] = 0.0f;
	pCam->Xpi[2][2] = dInv;
	pCam->Xpi[2][3] = 0.0f;

	pCam->Xpi[3][0] = 0.0f;
	pCam->Xpi[3][1] = 0.0f;
	pCam->Xpi[3][2] = dInv;
	pCam->Xpi[3][3] = 1.0f;

	return GZ_SUCCESS;
}

int GzSetupXiw(GzRender *render)
{
	if(!render)
		return GZ_FAILURE;

	GzCamera *pCam = &(render->camera);
	GzCoord vX, vY, vZ;

	//Defining Z-axis
	for(unsigned int i = 0; i < NUM_DIMENSIONS; ++i){
		vZ[i] = pCam->lookat[i] - pCam->position[i];
	}
	float tmpFloat = GzVectorNorm(&vZ);
	if(tmpFloat == 0.0f)
		return GZ_FAILURE;
	for(unsigned int i = 0; i < NUM_DIMENSIONS; ++i){
		vZ[i] = vZ[i] / tmpFloat; 
	}

	//Defining Y-axis
	tmpFloat = GzDotProduct(&(pCam->worldup), &vZ);
	for(unsigned int i = 0; i < NUM_DIMENSIONS; ++i){
		vY[i] = pCam->worldup[i] - tmpFloat * vZ[i];
	}
	tmpFloat = GzVectorNorm(&vY);
	if(tmpFloat == 0.0f)
		return GZ_FAILURE;
	for(unsigned int i = 0; i < NUM_DIMENSIONS; ++i){
		vY[i] = vY[i] / tmpFloat; 
	}

	//Defining X-axis
	vX[X] = vY[Y] * vZ[Z] - vY[Z] * vZ[Y];
	vX[Y] = vY[Z] * vZ[X] - vY[X] * vZ[Z];
	vX[Z] = vY[X] * vZ[Y] - vY[Y] * vZ[X];

	//Populating Xiw by row
	memcpy(pCam->Xiw[0], vX, sizeof(GzCoord));
	pCam->Xiw[0][3] = - GzDotProduct(&vX, &(pCam->position));
	
	memcpy(pCam->Xiw[1], vY, sizeof(GzCoord));
	pCam->Xiw[1][3] = - GzDotProduct(&vY, &(pCam->position));
	
	memcpy(pCam->Xiw[2], vZ, sizeof(GzCoord));
	pCam->Xiw[2][3] = - GzDotProduct(&vZ, &(pCam->position));

	pCam->Xiw[3][0] = 0.0f;
	pCam->Xiw[3][1] = 0.0f;
	pCam->Xiw[3][2] = 0.0f;
	pCam->Xiw[3][3] = 1.0f;
	
	return GZ_SUCCESS;
}

int GzInitDefaultCamera(GzRender *render)
{
	if(!render)
		return GZ_FAILURE;

	render->camera.position[X] = DEFAULT_IM_X;      
  	render->camera.position[Y] = DEFAULT_IM_Y;
  	render->camera.position[Z] = DEFAULT_IM_Z;

  	render->camera.lookat[X] = 0.0f;
  	render->camera.lookat[Y] = 0.0f;
  	render->camera.lookat[Z] = 0.0f;

  	render->camera.worldup[X] = 0.0f;
  	render->camera.worldup[Y] = 1.0f;
  	render->camera.worldup[Z] = 0.0f;

	render->camera.FOV = DEFAULT_FOV;              /* degrees */

	return GZ_SUCCESS;
}

/////////
//ALGEBRA FUNCTIONS
/////////

float GzVectorNorm(GzCoord *v){
	return sqrt(pow((*v)[X],2) + pow((*v)[Y],2) + pow((*v)[Z],2));
}
float GzDotProduct(GzCoord *a, GzCoord *b){
	return (*a)[X]*(*b)[X] + (*a)[Y]*(*b)[Y] + (*a)[Z]*(*b)[Z];
}
void GzMatrixMult(GzMatrix A, GzMatrix B, GzMatrix &res){
	for(unsigned int i = 0; i < 4; ++i){
			for(unsigned int j = 0; j < 4; ++j){
				res[i][j] = 0.0f;
				for(unsigned int pos = 0; pos < 4; ++pos){
					res[i][j] += A[i][pos] * B[pos][j];
				}
			}
	}
}
void GzMatVecMult(GzMatrix A, float v[4], float (&res)[4]){
	for(unsigned int i = 0; i < 4; ++i){
		res[i] = 0.0f;
		for(unsigned int j = 0; j < 4; ++j){
			res[i] += A[i][j] * v[j];
		}
	}

}