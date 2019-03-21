#include <GL/glew.h>
#include <GL/freeglut.h>

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <math.h> 
#include <cstdlib>.

#include <windows.h> 
#include <wingdi.h>

#include <iostream>

using namespace std;

#define FileString "MRbrain_growing_result.raw"
#define LENGTH 150
#define WIDTH 123
#define HEIGHT 184
#define FileType unsigned char

#define Line_Side_length 10
#define Cell_Side_length 1
#define IsoValue 55

float rotate_x = 0;
float rotate_y = 0;

float old_rot_x = 0;                              //剛按下滑鼠時的視窗座標 
float old_rot_y = 0;

float rot_x = 0;                                 //拖曳後的相對座標，用這決定要旋轉幾度 
float rot_y = 0;

float record_x = 0;                              //紀錄上一次旋轉的角度 
float record_y = 0;

int TotalTriangle = 0;

/*-----Create image space for textures -----*/
unsigned char  Cube_W_L[((WIDTH / 10) + 1)][((LENGTH / 10) + 1)][4];
unsigned char  Cube_H_L[((HEIGHT / 10) + 1)][((LENGTH / 10) + 1)][4];
unsigned char  Cube_H_W[((HEIGHT / 10) + 1)][((WIDTH / 10) + 1)][4];
unsigned int   textName[3];


FileType buf[HEIGHT*WIDTH*LENGTH];
double gGradient[HEIGHT*WIDTH*LENGTH][3];

typedef struct {
	float x;
	float y;
	float z;
} XYZ;

typedef struct tri{
	XYZ p[3];
	XYZ pGradient[3];
	struct tri *next;
} TRIANGLE;

TRIANGLE  *headTrianglePtr = NULL;
TRIANGLE  *nowReadPtr = NULL;

typedef struct {
	XYZ p[8];
	unsigned char val[8];
} GRIDCELL;

GRIDCELL gbuf[(HEIGHT - 1)*(WIDTH - 1)*(LENGTH - 1) / (Cell_Side_length*Cell_Side_length*Cell_Side_length)];

void fileOpen() {
	cout << "Openong file begin." << endl;

	FILE* pInput = fopen(FileString, "rb");
	//FILE* pOutput = NULL;
	if (!pInput) {
		cout << "檔案開啟失敗...\n" << endl;
		exit(1);
	}
	
	fread(buf, sizeof(FileType), LENGTH*WIDTH*HEIGHT, pInput);

	fclose(pInput);

	cout << "Opening file done." << endl;

}

void makeCell(int celllengh) {
	cout << "Making cell begin." << endl;

	for (int k = 0, z = 0; z < (HEIGHT - 1)/ celllengh; k += celllengh, z++) {
		for (int j = 0, y = 0; y < (WIDTH - 1)/celllengh; j += celllengh, y++) {
			for (int i = 0, x = 0; x < (LENGTH - 1)/celllengh; i += celllengh, x++) {
				gbuf[x + (LENGTH-1)/ Cell_Side_length*y + (LENGTH - 1) / Cell_Side_length*(WIDTH-1) / Cell_Side_length*z].p[0].x = i;
				gbuf[x + (LENGTH-1)/ Cell_Side_length*y + (LENGTH - 1) / Cell_Side_length*(WIDTH-1) / Cell_Side_length*z].p[0].y = j;
				gbuf[x + (LENGTH-1)/ Cell_Side_length*y + (LENGTH - 1) / Cell_Side_length*(WIDTH-1) / Cell_Side_length*z].p[0].z = k;
				gbuf[x + (LENGTH-1)/ Cell_Side_length*y + (LENGTH - 1) / Cell_Side_length*(WIDTH-1) / Cell_Side_length*z].val[0] = buf[(i)+LENGTH*(j)+LENGTH*WIDTH*(k)];
				
				gbuf[x + (LENGTH-1)/ Cell_Side_length*y + (LENGTH - 1) / Cell_Side_length*(WIDTH-1) / Cell_Side_length*z].p[1].x = i + celllengh;
				gbuf[x + (LENGTH-1)/ Cell_Side_length*y + (LENGTH - 1) / Cell_Side_length*(WIDTH-1) / Cell_Side_length*z].p[1].y = j;
				gbuf[x + (LENGTH-1)/ Cell_Side_length*y + (LENGTH - 1) / Cell_Side_length*(WIDTH-1) / Cell_Side_length*z].p[1].z = k;
				gbuf[x + (LENGTH-1)/ Cell_Side_length*y + (LENGTH - 1) / Cell_Side_length*(WIDTH-1) / Cell_Side_length*z].val[1] = buf[(i + celllengh) + LENGTH*(j)+LENGTH*WIDTH*(k)];
				
				gbuf[x + (LENGTH-1)/ Cell_Side_length*y + (LENGTH - 1) / Cell_Side_length*(WIDTH-1) / Cell_Side_length*z].p[2].x = i + celllengh;
				gbuf[x + (LENGTH-1)/ Cell_Side_length*y + (LENGTH - 1) / Cell_Side_length*(WIDTH-1) / Cell_Side_length*z].p[2].y = j + celllengh;
				gbuf[x + (LENGTH-1)/ Cell_Side_length*y + (LENGTH - 1) / Cell_Side_length*(WIDTH-1) / Cell_Side_length*z].p[2].z = k;
				gbuf[x + (LENGTH-1)/ Cell_Side_length*y + (LENGTH - 1) / Cell_Side_length*(WIDTH-1) / Cell_Side_length*z].val[2] = buf[(i + celllengh) + LENGTH*(j + celllengh) + LENGTH*WIDTH*(k)];
				
				gbuf[x + (LENGTH-1)/ Cell_Side_length*y + (LENGTH - 1) / Cell_Side_length*(WIDTH-1) / Cell_Side_length*z].p[3].x = i;
				gbuf[x + (LENGTH-1)/ Cell_Side_length*y + (LENGTH - 1) / Cell_Side_length*(WIDTH-1) / Cell_Side_length*z].p[3].y = j + celllengh;
				gbuf[x + (LENGTH-1)/ Cell_Side_length*y + (LENGTH - 1) / Cell_Side_length*(WIDTH-1) / Cell_Side_length*z].p[3].z = k;
				gbuf[x + (LENGTH-1)/ Cell_Side_length*y + (LENGTH - 1) / Cell_Side_length*(WIDTH-1) / Cell_Side_length*z].val[3] = buf[(i)+LENGTH*(j + celllengh) + LENGTH*WIDTH*(k)];
			
				gbuf[x + (LENGTH-1)/ Cell_Side_length*y + (LENGTH - 1) / Cell_Side_length*(WIDTH-1) / Cell_Side_length*z].p[4].x = i;
				gbuf[x + (LENGTH-1)/ Cell_Side_length*y + (LENGTH - 1) / Cell_Side_length*(WIDTH-1) / Cell_Side_length*z].p[4].y = j;
				gbuf[x + (LENGTH-1)/ Cell_Side_length*y + (LENGTH - 1) / Cell_Side_length*(WIDTH-1) / Cell_Side_length*z].p[4].z = k + celllengh;
				gbuf[x + (LENGTH-1)/ Cell_Side_length*y + (LENGTH - 1) / Cell_Side_length*(WIDTH-1) / Cell_Side_length*z].val[4] = buf[(i)+LENGTH*(j)+LENGTH*WIDTH*(k + celllengh)];
		
				gbuf[x + (LENGTH-1)/ Cell_Side_length*y + (LENGTH - 1) / Cell_Side_length*(WIDTH-1) / Cell_Side_length*z].p[5].x = i + celllengh;
				gbuf[x + (LENGTH-1)/ Cell_Side_length*y + (LENGTH - 1) / Cell_Side_length*(WIDTH-1) / Cell_Side_length*z].p[5].y = j;
				gbuf[x + (LENGTH-1)/ Cell_Side_length*y + (LENGTH - 1) / Cell_Side_length*(WIDTH-1) / Cell_Side_length*z].p[5].z = k + celllengh;
				gbuf[x + (LENGTH-1)/ Cell_Side_length*y + (LENGTH - 1) / Cell_Side_length*(WIDTH-1) / Cell_Side_length*z].val[5] = buf[(i + celllengh) + LENGTH*(j)+LENGTH*WIDTH*(k + celllengh)];
			
				gbuf[x + (LENGTH-1)/ Cell_Side_length*y + (LENGTH - 1) / Cell_Side_length*(WIDTH-1) / Cell_Side_length*z].p[6].x = i + celllengh;
				gbuf[x + (LENGTH-1)/ Cell_Side_length*y + (LENGTH - 1) / Cell_Side_length*(WIDTH-1) / Cell_Side_length*z].p[6].y = j + celllengh;
				gbuf[x + (LENGTH-1)/ Cell_Side_length*y + (LENGTH - 1) / Cell_Side_length*(WIDTH-1) / Cell_Side_length*z].p[6].z = k + celllengh;
				gbuf[x + (LENGTH-1)/ Cell_Side_length*y + (LENGTH - 1) / Cell_Side_length*(WIDTH-1) / Cell_Side_length*z].val[6] = buf[(i + celllengh) + LENGTH*(j + celllengh) + LENGTH*WIDTH*(k + celllengh)];

				gbuf[x + (LENGTH-1)/ Cell_Side_length*y + (LENGTH - 1) / Cell_Side_length*(WIDTH-1) / Cell_Side_length*z].p[7].x = i;
				gbuf[x + (LENGTH-1)/ Cell_Side_length*y + (LENGTH - 1) / Cell_Side_length*(WIDTH-1) / Cell_Side_length*z].p[7].y = j + celllengh;
				gbuf[x + (LENGTH-1)/ Cell_Side_length*y + (LENGTH - 1) / Cell_Side_length*(WIDTH-1) / Cell_Side_length*z].p[7].z = k + celllengh;
				gbuf[x + (LENGTH-1)/ Cell_Side_length*y + (LENGTH - 1) / Cell_Side_length*(WIDTH-1) / Cell_Side_length*z].val[7] = buf[(i)+LENGTH*(j + celllengh) + LENGTH*WIDTH*(k + celllengh)];
			}
		}
	}

	cout << "Making cell done." << endl;
}

void pointGradient() {
	cout << "Making gradient begin." << endl;
	for (int h = 0; h < HEIGHT; h++) {
		for (int w = 0; w < WIDTH; w++) {
			for (int l = 0; l < LENGTH; l++) {
				if (l == 0) {
					gGradient[l + w*LENGTH + h*LENGTH*WIDTH][0] = (buf[(l + 1) + w*LENGTH + h*LENGTH*WIDTH] - buf[l + w*LENGTH + h*LENGTH*WIDTH]);
				}
				else if(l == LENGTH - 1) {
					gGradient[l + w*LENGTH + h*LENGTH*WIDTH][0] = (buf[l + w*LENGTH + h*LENGTH*WIDTH] - buf[(l - 1) + w*LENGTH + h*LENGTH*WIDTH]);
				}
				else{
					gGradient[l + w*LENGTH + h*LENGTH*WIDTH][0] = (buf[(l + 1) + w*LENGTH + h*LENGTH*WIDTH] - buf[(l - 1) + w*LENGTH + h*LENGTH*WIDTH]) / 2;
				}
				if (w == 0) {
					gGradient[l + w*LENGTH + h*LENGTH*WIDTH][1] = (buf[l + (w + 1)*LENGTH + h*LENGTH*WIDTH] - buf[l + w*LENGTH + h*LENGTH*WIDTH]);
				}
				else if (w == WIDTH - 1) {
					gGradient[l + w*LENGTH + h*LENGTH*WIDTH][1] = (buf[l + w*LENGTH + h*LENGTH*WIDTH] - buf[l + (w - 1)*LENGTH + h*LENGTH*WIDTH]);
				}
				else {
					gGradient[l + w*LENGTH + h*LENGTH*WIDTH][1] = (buf[l + (w + 1)*LENGTH + h*LENGTH*WIDTH] - buf[l + (w - 1)*LENGTH + h*LENGTH*WIDTH]) / 2;
				}
				if (h == 0) {
					gGradient[l + w*LENGTH + h*LENGTH*WIDTH][2] = (buf[l + w*LENGTH + (h + 1)*LENGTH*WIDTH] - buf[l + w*LENGTH + h*LENGTH*WIDTH]);
				}
				else if (h == HEIGHT - 1) {
					gGradient[l + w*LENGTH + h*LENGTH*WIDTH][2] = (buf[l + w*LENGTH + h*LENGTH*WIDTH] - buf[l + w*LENGTH + (h - 1)*LENGTH*WIDTH]);
				}
				else {
					gGradient[l + w*LENGTH + h*LENGTH*WIDTH][2] = (buf[l + w*LENGTH + (h + 1)*LENGTH*WIDTH] - buf[l + w*LENGTH + (h - 1)*LENGTH*WIDTH]) / 2;
				}
				//Normalize
				double d = sqrt(pow(gGradient[l + w*LENGTH + h*LENGTH*WIDTH][0], 2) + pow(gGradient[l + w*LENGTH + h*LENGTH*WIDTH][1], 2) + pow(gGradient[l + w*LENGTH + h*LENGTH*WIDTH][2], 2));
				if (d != 0) {
					gGradient[l + w*LENGTH + h*LENGTH*WIDTH][0] = (gGradient[l + w*LENGTH + h*LENGTH*WIDTH][0] / d);
					gGradient[l + w*LENGTH + h*LENGTH*WIDTH][1] = (gGradient[l + w*LENGTH + h*LENGTH*WIDTH][1] / d);
					gGradient[l + w*LENGTH + h*LENGTH*WIDTH][2] = (gGradient[l + w*LENGTH + h*LENGTH*WIDTH][2] / d);
				}
			}
		}
	}
	cout << "Making gradient end." << endl;
}

XYZ VertexInterp(double isolevel, XYZ p1, XYZ p2, double valp1, double valp2) {
	double mu;
	XYZ p;

	if (fabs(isolevel - valp1) < 0.00001)
		return(p1);
	if (fabs(isolevel - valp2) < 0.00001)
		return(p2);
	if (fabs(valp1 - valp2) < 0.00001)
		return(p1);
	mu = (isolevel - valp1) / (valp2 - valp1);
	p.x = p1.x + mu * (p2.x - p1.x);
	p.y = p1.y + mu * (p2.y - p1.y);
	p.z = p1.z + mu * (p2.z - p1.z);

	return(p);
}

XYZ triangleGradient(XYZ p1,XYZ p2,XYZ tp) {
	XYZ p;
	double dx = p1.x - p2.x;
	double dy = p1.y - p2.y;
	double dz = p1.z - p2.z;
	if (dx) {
		p.x = fabs((p1.x - tp.x) / dx)*(double)gGradient[(int)(p1.x + p1.y*LENGTH + p1.z*LENGTH*WIDTH)][0] + fabs((p2.x - tp.x) / dx)*(double)gGradient[(int)(p2.x + p2.y*LENGTH + p2.z*LENGTH*WIDTH)][0];
	}
	else {
		p.x = (double)gGradient[(int)(p2.x + p2.y*LENGTH + p2.z*LENGTH*WIDTH)][0];
	}
	if (dy) {
		p.y = fabs((p1.y - tp.y) / dy)*(double)gGradient[(int)(p1.x + p1.y*LENGTH + p1.z*LENGTH*WIDTH)][1] + fabs((p2.y - tp.y) / dy)*(double)gGradient[(int)(p2.x + p2.y*LENGTH + p2.z*LENGTH*WIDTH)][1];
	}
	else {
		p.y = (double)gGradient[(int)(p1.x + p1.y*LENGTH + p1.z*LENGTH*WIDTH)][1];
	}
	if (dz) {
		p.z = fabs((p1.z - tp.z) / dz)*(double)gGradient[(int)(p1.x + p1.y*LENGTH + p1.z*LENGTH*WIDTH)][2] + fabs((p2.z - tp.z) / dz)*(double)gGradient[(int)(p2.x + p2.y*LENGTH + p2.z*LENGTH*WIDTH)][2];
	}
	else {
		p.z = (double)gGradient[(int)(p1.x + p1.y*LENGTH + p1.z*LENGTH*WIDTH)][2];
	}
		
	return p;
}

int Polygonise(GRIDCELL grid, double isolevel, TRIANGLE *triangles)
{
	int i, ntriang;
	int cubeindex;
	XYZ vertlist[12];
	XYZ pG[12];

	int edgeTable[256] = {
		0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
		0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
		0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
		0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
		0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
		0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
		0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
		0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
		0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
		0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
		0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
		0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
		0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
		0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
		0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
		0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
		0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
		0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
		0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
		0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
		0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
		0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
		0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
		0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
		0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
		0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
		0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
		0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
		0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
		0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
		0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
		0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0 };
	int triTable[256][16] =
	{ { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1 },
	{ 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1 },
	{ 3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1 },
	{ 3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1 },
	{ 2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1 },
	{ 8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1 },
	{ 4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1 },
	{ 3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1 },
	{ 4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1 },
	{ 4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1 },
	{ 5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1 },
	{ 2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1 },
	{ 9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1 },
	{ 2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1 },
	{ 10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1 },
	{ 4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1 },
	{ 5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1 },
	{ 5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1 },
	{ 10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1 },
	{ 8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1 },
	{ 2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1 },
	{ 7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1 },
	{ 2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1 },
	{ 11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1 },
	{ 5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1 },
	{ 11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1 },
	{ 11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1 },
	{ 5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1 },
	{ 2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1 },
	{ 5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1 },
	{ 6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1 },
	{ 3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1 },
	{ 6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1 },
	{ 5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1 },
	{ 10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1 },
	{ 6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1 },
	{ 8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1 },
	{ 7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1 },
	{ 3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1 },
	{ 5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1 },
	{ 0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1 },
	{ 9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1 },
	{ 8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1 },
	{ 5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1 },
	{ 0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1 },
	{ 6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1 },
	{ 10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1 },
	{ 10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1 },
	{ 8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1 },
	{ 1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1 },
	{ 3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1 },
	{ 0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1 },
	{ 10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1 },
	{ 3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1 },
	{ 6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1 },
	{ 9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1 },
	{ 8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1 },
	{ 3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1 },
	{ 6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1 },
	{ 10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1 },
	{ 10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1 },
	{ 2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1 },
	{ 7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1 },
	{ 7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1 },
	{ 2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1 },
	{ 1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1 },
	{ 11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1 },
	{ 8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1 },
	{ 0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1 },
	{ 7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1 },
	{ 10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1 },
	{ 2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1 },
	{ 6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1 },
	{ 7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1 },
	{ 2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1 },
	{ 10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1 },
	{ 10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1 },
	{ 0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1 },
	{ 7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1 },
	{ 6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1 },
	{ 8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1 },
	{ 6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1 },
	{ 4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1 },
	{ 10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1 },
	{ 8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1 },
	{ 1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1 },
	{ 8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1 },
	{ 10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1 },
	{ 4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1 },
	{ 10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1 },
	{ 5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1 },
	{ 11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1 },
	{ 9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1 },
	{ 6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1 },
	{ 7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1 },
	{ 3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1 },
	{ 7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1 },
	{ 3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1 },
	{ 6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1 },
	{ 9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1 },
	{ 1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1 },
	{ 4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1 },
	{ 7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1 },
	{ 6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1 },
	{ 3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1 },
	{ 0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1 },
	{ 6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1 },
	{ 0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1 },
	{ 11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1 },
	{ 6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1 },
	{ 5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1 },
	{ 9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1 },
	{ 1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1 },
	{ 10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1 },
	{ 0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1 },
	{ 5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1 },
	{ 10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1 },
	{ 11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1 },
	{ 9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1 },
	{ 7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1 },
	{ 2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1 },
	{ 8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1 },
	{ 9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1 },
	{ 9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1 },
	{ 1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1 },
	{ 5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1 },
	{ 0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1 },
	{ 10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1 },
	{ 2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1 },
	{ 0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1 },
	{ 0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1 },
	{ 9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1 },
	{ 5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1 },
	{ 3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1 },
	{ 5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1 },
	{ 8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1 },
	{ 9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1 },
	{ 1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1 },
	{ 3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1 },
	{ 4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1 },
	{ 9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1 },
	{ 11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1 },
	{ 11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1 },
	{ 2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1 },
	{ 9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1 },
	{ 3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1 },
	{ 1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1 },
	{ 4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1 },
	{ 4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1 },
	{ 3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1 },
	{ 3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1 },
	{ 0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1 },
	{ 9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1 },
	{ 1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ 0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
	{ -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 } };

	/*
	Determine the index into the edge table which
	tells us which vertices are inside of the surface
	*/
	cubeindex = 0;
	if (grid.val[0] < isolevel) cubeindex |= 1;
	if (grid.val[1] < isolevel) cubeindex |= 2;
	if (grid.val[2] < isolevel) cubeindex |= 4;
	if (grid.val[3] < isolevel) cubeindex |= 8;
	if (grid.val[4] < isolevel) cubeindex |= 16;
	if (grid.val[5] < isolevel) cubeindex |= 32;
	if (grid.val[6] < isolevel) cubeindex |= 64;
	if (grid.val[7] < isolevel) cubeindex |= 128;

	/* Cube is entirely in/out of the surface */
	if (edgeTable[cubeindex] == 0)
		return(0);
	/* Find the vertices where the surface intersects the cube */
	if (edgeTable[cubeindex] & 1) {
		vertlist[0] = VertexInterp(isolevel, grid.p[0], grid.p[1], grid.val[0], grid.val[1]);
		pG[0] = triangleGradient(grid.p[0], grid.p[1], vertlist[0]);
	}
	if (edgeTable[cubeindex] & 2) {
		vertlist[1] = VertexInterp(isolevel, grid.p[1], grid.p[2], grid.val[1], grid.val[2]);
		pG[1] = triangleGradient(grid.p[1], grid.p[2], vertlist[1]);
	}
	if (edgeTable[cubeindex] & 4) {
		vertlist[2] = VertexInterp(isolevel, grid.p[2], grid.p[3], grid.val[2], grid.val[3]);
		pG[2] = triangleGradient(grid.p[2], grid.p[3], vertlist[2]);
	}
	if (edgeTable[cubeindex] & 8) {
		vertlist[3] = VertexInterp(isolevel, grid.p[3], grid.p[0], grid.val[3], grid.val[0]);
		pG[3] = triangleGradient(grid.p[3], grid.p[0], vertlist[3]);
	}
	if (edgeTable[cubeindex] & 16) {
		vertlist[4] = VertexInterp(isolevel, grid.p[4], grid.p[5], grid.val[4], grid.val[5]);
		pG[4] = triangleGradient(grid.p[4], grid.p[5], vertlist[4]);
	}
	if (edgeTable[cubeindex] & 32) {
		vertlist[5] = VertexInterp(isolevel, grid.p[5], grid.p[6], grid.val[5], grid.val[6]);
		pG[5] = triangleGradient(grid.p[5], grid.p[6], vertlist[5]);
	}
	if (edgeTable[cubeindex] & 64) {
		vertlist[6] = VertexInterp(isolevel, grid.p[6], grid.p[7], grid.val[6], grid.val[7]);
		pG[6] = triangleGradient(grid.p[6], grid.p[7], vertlist[6]);
	}
	if (edgeTable[cubeindex] & 128) {
		vertlist[7] = VertexInterp(isolevel, grid.p[7], grid.p[4], grid.val[7], grid.val[4]);
		pG[7] = triangleGradient(grid.p[7], grid.p[4], vertlist[7]);
	}
	if (edgeTable[cubeindex] & 256) {
		vertlist[8] = VertexInterp(isolevel, grid.p[0], grid.p[4], grid.val[0], grid.val[4]);
		pG[8] = triangleGradient(grid.p[0], grid.p[4], vertlist[8]);
	}
	if (edgeTable[cubeindex] & 512) {
		vertlist[9] = VertexInterp(isolevel, grid.p[1], grid.p[5], grid.val[1], grid.val[5]);
		pG[9] = triangleGradient(grid.p[1], grid.p[5], vertlist[9]);
	}
	if (edgeTable[cubeindex] & 1024) {
		vertlist[10] = VertexInterp(isolevel, grid.p[2], grid.p[6], grid.val[2], grid.val[6]);
		pG[10] = triangleGradient(grid.p[2], grid.p[6], vertlist[10]);
	}
	if (edgeTable[cubeindex] & 2048) {
		vertlist[11] = VertexInterp(isolevel, grid.p[3], grid.p[7], grid.val[3], grid.val[7]);
		pG[11] = triangleGradient(grid.p[3], grid.p[7], vertlist[11]);
	}

	/* Create the triangle */
	ntriang = 0;
	for (i = 0; triTable[cubeindex][i] != -1; i += 3) {
		TRIANGLE * newTriangle = (TRIANGLE *)malloc(sizeof(TRIANGLE));
		newTriangle->p[0] = vertlist[triTable[cubeindex][i]];
		newTriangle->pGradient[0] = pG[triTable[cubeindex][i]];
		newTriangle->p[1] = vertlist[triTable[cubeindex][i + 1]];
		newTriangle->pGradient[1] = pG[triTable[cubeindex][i + 1]];
		newTriangle->p[2] = vertlist[triTable[cubeindex][i + 2]];
		newTriangle->pGradient[2] = pG[triTable[cubeindex][i + 2]];
		newTriangle->next = NULL;
		if (TotalTriangle == 0) {
			headTrianglePtr = newTriangle;
			nowReadPtr = newTriangle;
		}
		else {
			nowReadPtr->next = newTriangle;
			nowReadPtr = newTriangle;
		}
		ntriang++;
		TotalTriangle++;
	}
	
	return(ntriang);
}

void drawBackCube() {

	cout << "Draw Cube begin." << endl;
	
	//Cube Material
	float cube_material_amb[4] = { 0.2,0.2,0.2,1 };
	float cube_material_dif[4] = { 0.4,0.4,0.4,1 };
	float cube_material_spe[4] = { 0.4,0.4,0.4,1 };

	glMaterialfv(GL_FRONT, GL_AMBIENT, cube_material_amb);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, cube_material_dif);
	glMaterialfv(GL_FRONT, GL_SPECULAR, cube_material_spe);
	glMaterialf(GL_FRONT, GL_SHININESS, 64.0);

	glEnable(GL_CULL_FACE);
	glCullFace(GL_FRONT);

	glBegin(GL_POLYGON);
	glTexCoord2f(0.0, 0.0); glNormal3f(0.0, 0.0, 1.0); glVertex3f(0.0, 0.0, 0.0);
	glTexCoord2f(1.0, 0.0); glNormal3f(0.0, 0.0, 1.0); glVertex3f(1.0, 0.0, 0.0);
	glTexCoord2f(1.0, 1.0); glNormal3f(0.0, 0.0, 1.0); glVertex3f(1.0, 1.0, 0.0);
	glTexCoord2f(0.0, 1.0); glNormal3f(0.0, 0.0, 1.0); glVertex3f(0.0, 1.0, 0.0);
	glEnd();

	glBegin(GL_POLYGON);
	glTexCoord2f(0.0, 0.0); glNormal3f(0.0, 0.0, -1.0); glVertex3f(1.0, 0.0, -1.0);
	glTexCoord2f(1.0, 0.0); glNormal3f(0.0, 0.0, -1.0); glVertex3f(0.0, 0.0, -1.0);
	glTexCoord2f(1.0, 1.0); glNormal3f(0.0, 0.0, -1.0); glVertex3f(0.0, 1.0, -1.0);
	glTexCoord2f(0.0, 1.0); glNormal3f(0.0, 0.0, -1.0); glVertex3f(1.0, 1.0, -1.0);
	glEnd();

	glBegin(GL_POLYGON);
	glTexCoord2f(0.0, 0.0); glNormal3f(1.0, 0.0, 0.0); glVertex3f(1.0, 0.0, 0.0);
	glTexCoord2f(1.0, 0.0); glNormal3f(1.0, 0.0, 0.0); glVertex3f(1.0, 0.0, -1.0);
	glTexCoord2f(1.0, 1.0); glNormal3f(1.0, 0.0, 0.0); glVertex3f(1.0, 1.0, -1.0);
	glTexCoord2f(0.0, 1.0); glNormal3f(1.0, 0.0, 0.0); glVertex3f(1.0, 1.0, 0.0);
	glEnd();

	glBegin(GL_POLYGON);
	glTexCoord2f(0.0, 0.0); glNormal3f(0.0, 1.0, 0.0); glVertex3f(0.0, 1.0, 0.0);
	glTexCoord2f(1.0, 0.0); glNormal3f(0.0, 1.0, 0.0); glVertex3f(1.0, 1.0, 0.0);
	glTexCoord2f(1.0, 1.0); glNormal3f(0.0, 1.0, 0.0); glVertex3f(1.0, 1.0, -1.0);
	glTexCoord2f(0.0, 1.0); glNormal3f(0.0, 1.0, 0.0); glVertex3f(0.0, 1.0, -1.0);
	glEnd();

	glBegin(GL_POLYGON);
	glTexCoord2f(0.0, 0.0); glNormal3f(-1.0, 0.0, 0.0); glVertex3f(0.0, 0.0, -1.0);
	glTexCoord2f(1.0, 0.0); glNormal3f(-1.0, 0.0, 0.0); glVertex3f(0.0, 0.0, 0.0);
	glTexCoord2f(1.0, 1.0); glNormal3f(-1.0, 0.0, 0.0); glVertex3f(0.0, 1.0, 0.0);
	glTexCoord2f(0.0, 1.0); glNormal3f(-1.0, 0.0, 0.0); glVertex3f(0.0, 1.0, -1.0);
	glEnd();

	glBegin(GL_POLYGON);
	glTexCoord2f(0.0, 0.0); glNormal3f(0.0, -1.0, 0.0); glVertex3f(0.0, 0.0, -1.0);
	glTexCoord2f(1.0, 0.0); glNormal3f(0.0, -1.0, 0.0); glVertex3f(1.0, 0.0, -1.0);
	glTexCoord2f(1.0, 1.0); glNormal3f(0.0, -1.0, 0.0); glVertex3f(1.0, 0.0, 0.0);
	glTexCoord2f(0.0, 1.0); glNormal3f(0.0, -1.0, 0.0); glVertex3f(0.0, 0.0, 0.0);
	glEnd();

	glDisable(GL_CULL_FACE);

	cout << "Draw Cube end." << endl;
}

void Caculate_Triangle() {

	cout << "Caculate triangle begin." << endl;
	
	float triangle_material_amb[4] = { 0.0,0.0,0.0,1 };
	float triangle_material_dif[4] = { 0.5,0.5,0.5,1.0 };
	float triangle_material_spe[4] = { 0.8,0.8,0.8,1.0 };

	glMaterialfv(GL_FRONT, GL_AMBIENT, triangle_material_amb);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, triangle_material_dif);
	glMaterialfv(GL_FRONT, GL_SPECULAR, triangle_material_spe);
	glMaterialf(GL_FRONT, GL_SHININESS, 64.0);
	
	for (int i = 0; i < ((HEIGHT - 1)*(WIDTH - 1)*(LENGTH - 1) / (Cell_Side_length*Cell_Side_length*Cell_Side_length)); i++) {
		Polygonise(gbuf[i], IsoValue, nowReadPtr);
	}

	cout << "TotalTriangle:" << TotalTriangle << endl << "Caculate triangle done." << endl;
}

void Light() {


	//左方燈
	float light1_amb[4] = { 0.3,0.3,0.3,1 };
	float light1_dif[4] = { 1.0,1.0,1.0,1 };
	float light1_spe[4] = { 0.3,0.3,0.3,1 };
	float light1_dir[4] = { 0, 0, -600, 0.0 };

	glEnable(GL_LIGHT1);
	glLightfv(GL_LIGHT1, GL_AMBIENT, light1_amb);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_dif);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light1_spe);
	glLightfv(GL_LIGHT1, GL_POSITION, light1_dir);
	glDisable(GL_LIGHT1);


}

void Make_Cube_Texture(int h, int w, int l,int side) {
	cout << "Making Cube Texture begin." << endl;

	//make W_L_Cube
	for (int y = 0; y < w; y++) {
		for (int x = 0; x < l; x++) {
			Cube_W_L[y][x][0] = 1.0 * 255;
			Cube_W_L[y][x][1] = 1.0 * 255;
			Cube_W_L[y][x][2] = 1.0 * 255;
			Cube_W_L[y][x][3] = 255;
		}
	}

	//make H_L_Cube
	for (int y = 0; y < h; y++) {
		for (int x = 0; x < l; x++) {

		}
	}

	//make H_W_Cube
	for (int y = 0; y < h; y++) {
		for (int x = 0; x < w; x++) {

		}
	}




	cout << "Making Cube Texture end." << endl;
}

void Texture()
{

	glPixelStorei(GL_UNPACK_ALIGNMENT, 4);
	glGenTextures(3, textName);

	Make_Cube_Texture(((HEIGHT / 10) + 1), ((WIDTH / 10) + 1), ((LENGTH / 10) + 1), Line_Side_length);
	glBindTexture(GL_TEXTURE_2D, textName[0]);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, ((WIDTH / 10) + 1), ((LENGTH / 10) + 1), 0, GL_RGBA, GL_UNSIGNED_BYTE, Cube_W_L);

	glBindTexture(GL_TEXTURE_2D, textName[1]);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, ((HEIGHT / 10) + 1), ((LENGTH / 10) + 1), 0, GL_RGBA, GL_UNSIGNED_BYTE, Cube_H_L);

	glBindTexture(GL_TEXTURE_2D, textName[2]);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, ((HEIGHT / 10) + 1), ((WIDTH / 10) + 1), 0, GL_RGBA, GL_UNSIGNED_BYTE, Cube_H_W);

}

void Init()
{
	fileOpen();
	pointGradient();
	makeCell(Cell_Side_length);
	Caculate_Triangle();
	Light();
	Texture();
}

void drawTriangle() {
	cout << "Draw Triangle begin."<<endl;
	glBegin(GL_TRIANGLES);
	TRIANGLE * nowPtr = headTrianglePtr;
	while (nowPtr != NULL) {
		glNormal3f(nowPtr->pGradient[0].x, nowPtr->pGradient[0].y, nowPtr->pGradient[0].z);
		glVertex3f(nowPtr->p[0].x, nowPtr->p[0].y, nowPtr->p[0].z);
		glNormal3f(nowPtr->pGradient[1].x, nowPtr->pGradient[1].y, nowPtr->pGradient[1].z);
		glVertex3f(nowPtr->p[1].x, nowPtr->p[1].y, nowPtr->p[1].z);
		glNormal3f(nowPtr->pGradient[2].x, nowPtr->pGradient[2].y, nowPtr->pGradient[2].z);
		glVertex3f(nowPtr->p[2].x, nowPtr->p[2].y, nowPtr->p[2].z);
		if (nowPtr->next != NULL) {
			nowPtr = nowPtr->next;
		}
		else {
			nowPtr = NULL;
		}
	}
	glEnd();
	cout << "Draw Triangle done." << endl;
	return;
}

void Display()
{
	glEnable(GL_DEPTH_TEST);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glShadeModel(GL_SMOOTH);

	glClearColor(0.0, 0.0, 0.0, 1);

	//projection
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45, 1, 50, 5000);

	//Camera
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	gluLookAt(0, 0, 600, 0, 0, 0, 0, 1, 0);

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT1);

	glPushMatrix();
	glTranslated(-LENGTH/2, -WIDTH/2, -HEIGHT/2);
	glScalef(LENGTH * 2, WIDTH * 2, HEIGHT * 2);
	//drawBackCube();
	glPopMatrix();

	glPushMatrix();
	glRotated(rotate_x, 0.0, 1.0, 0.0);
	glRotated(rotate_y, 1.0, 0.0, 0.0);

	glTranslated(-LENGTH / 2, -WIDTH / 2, -HEIGHT / 2);
	//畫座標軸
	//畫等值三角面
	drawTriangle();
	glPopMatrix();
	glDisable(GL_BLEND);

	glDisable(GL_LIGHTING);

	glutSwapBuffers();
}

void Timer(int c)
{
	glutPostRedisplay();
	glutTimerFunc(50, Timer, 0);
}

void KeyboardDown(unsigned char c, int x, int y)
{
	glFinish();
}

void Mouse(int button, int state, int x, int y)
{
	if (state)
	{
		record_x += x - old_rot_x;
		record_y += y - old_rot_y;

		rot_x = 0;
		rot_y = 0;
	}
	else
	{
		old_rot_x = x;
		old_rot_y = y;
	}
	glFinish();
}

void Motion(int x, int y)
{
	rot_x = x - old_rot_x;
	rot_y = y - old_rot_y;
	glFinish();
	rotate_x = (rot_x + record_x);
	rotate_y = (rot_y + record_y);

	glutPostRedisplay();
}

int main(int argc, char** argv) {

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
	glutInitWindowSize(800, 800);
	glutCreateWindow("Iso-Surface Test");
	glutDisplayFunc(Display);
	
	glutKeyboardFunc(KeyboardDown);
	glutMotionFunc(Motion);
	glutMouseFunc(Mouse);
	
	glutTimerFunc(50, Timer, 0);
	glewInit();
	Init();
	glutMainLoop();
	
}
