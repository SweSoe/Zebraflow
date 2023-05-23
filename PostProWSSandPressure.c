#include <getopt.h>
#include <sys/time.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#define PI 3.14159265358979
#define	KM	723//	  //number of x-nodes
#define	JM	547 //number of y-nodes
#define	IM	249//422  //number of z nodes

#define Z0 -27.// Xmax = -27. + 0.5*(723-1) = 334
#define Y0 29.// Ymax = 29. + 0.5*(547-1) = 302
#define X0 -1.// Zmax = -1 + 0.5*(249-1) = 123 
#define DX 0.5

#define iterstart 0
#define itermax 720000//set to the maximum iteration level you have simulated
#define outputfreq 10000
#define NodeM 1107042//298184//1107243
#define CellM 2212824//553209//2213222
#define Pgauge 120
//float WallLIST[3][1100000], pressure[1100000];
float Xm[NodeM], Ym[NodeM], Zm[NodeM], PRESSW[NodeM], WSS[3][NodeM];
int NeighCount1[NodeM], NeighCount2[NodeM];
int POLYTYPE[CellM], POLYPID[3][CellM];
int ValidData[NodeM];

float X[IM+1][JM+1][KM+1], Y[IM+1][JM+1][KM+1], Z[IM+1][JM+1][KM+1], U[IM+1][JM+1][KM+1], V[IM+1][JM+1][KM+1], W[IM+1][JM+1][KM+1], PRESS[IM+1][JM+1][KM+1];
int B_index[IM+1][JM+1][KM+1];
float SRT[9][IM+1][JM+1][KM+1];//shearratetensor 9 members in a 3D physical space
float Norm[3][NodeM],WALLCELLNORMALS[3][CellM];

int main(void){
	FILE *fOUT, *fIN;
	char cin;
	int i,j,k,ii,jj,kk,l,m,iii,jjj,kkk, linecount;
	float UU,VV,WW,DD,Bx,By,Bz;
	char number[15], num[20], name[50], bulk[200], junkchar;
	float xp, yp, zp;
	int i1, i2, j1, j2, k1, k2;
	float distance, mindis;
	int bvalue;
	char zero8[] = "00000000", zero7[] = "0000000", zero6[] = "000000", zero5[]="00000", zero4[]="0000", zero3[]="000", zero2[]="00", zero1[]="0";
	float pp;
	int iter;
	int junkint;
	//WallList[0][N] = x, WallList[1][N] = x, WallList[2][N] = x, WallList[3][N] = pressure
	int check, counter;
	float weight, vtan[3], vnormmag;
	
	///Read in nodes from vtk file
	
	strcpy(name, "InputData/WallMeshNodeCoord.dat");	
	printf("Opening %s \n", name);
	fIN = fopen(name, "r");
		fgets(bulk, 180, fIN);
		for(l=0;l<NodeM;l++){
			fscanf(fIN,"%f",&Xm[l]);
			fscanf(fIN,"%f",&Ym[l]);
			fscanf(fIN,"%f",&Zm[l]);
		}
		
		
	fclose(fIN);
	
	strcpy(name, "InputData/WallMeshNodeNormals.dat");	
	printf("Opening %s \n", name);
	fIN = fopen(name, "r");
		fgets(bulk, 180, fIN);
		fgets(bulk, 180, fIN);
		for(l=0;l<NodeM;l++){
			fscanf(fIN,"%f",&Norm[0][l]);
			fscanf(fIN,"%f",&Norm[1][l]);
			fscanf(fIN,"%f",&Norm[2][l]);
		}
			
	fclose(fIN);
	
	
	strcpy(name, "InputData/WallMeshCells.dat");	
	printf("Opening %s \n", name);
	fIN = fopen(name, "r");
		fgets(bulk, 180, fIN);
		for(l=0;l<CellM;l++){
			fscanf(fIN,"%d",&junkint);
			fscanf(fIN,"%d",&i);
			fscanf(fIN,"%d",&j);
			fscanf(fIN,"%d",&k);
			POLYPID[0][l] = i+1;
			POLYPID[1][l] = j+1;
			POLYPID[2][l] = k+1;
		}
	fclose(fIN);
	
	/// The wall mesh normals file was too big for GitHub so it was split into 2 files///
	linecount = 0;
	strcpy(name, "InputData/WallMeshCellNormals.aa");
		// Open the file 
		fIN = fopen(name, "r"); 
	  
		// Check if file exists 
		if(fIN == NULL) 
		{ 
			printf("Could not open file %s", name); 
			return 0; 
		} 
	  
		// Extract characters from file and store in character c 
		for (cin = getc(fIN); cin != EOF; cin = getc(fIN)){
			if (cin == '\n') // Increment count if this character is newline 
				linecount = linecount + 1; 
		}
		// Close the file 
	fclose(fIN); 
	printf("The file %s has %d lines\n ", name, linecount); 
	linecount = linecount-2;
	
	counter = linecount*3;//3 triangle element data in each row, so total triangles in each file is linecount*3
	strcpy(name, "InputData/WallMeshCellNormals.aa");	
	printf("Opening %s \n", name);
	fIN = fopen(name, "r");
		fgets(bulk, 180, fIN);
		fgets(bulk, 180, fIN);
		for(m=1;m<=counter;m++){
			fscanf(fIN,"%f",&WALLCELLNORMALS[0][m-1]);
			fscanf(fIN,"%f",&WALLCELLNORMALS[1][m-1]);
			fscanf(fIN,"%f",&WALLCELLNORMALS[2][m-1]);
		}	
	fclose(fIN);
	
	strcpy(name, "InputData/WallMeshCellNormals.ab");	
	printf("Opening %s \n", name);
	fIN = fopen(name, "r");
		for(m=counter+1;m<=CellM;m++){
			fscanf(fIN,"%f",&WALLCELLNORMALS[0][m-1]);
			fscanf(fIN,"%f",&WALLCELLNORMALS[1][m-1]);
			fscanf(fIN,"%f",&WALLCELLNORMALS[2][m-1]);
		}	
	fclose(fIN);
	/// //////////////
	
	
	
	for(i=1;i<=IM;i++){
		for(j=1;j<=JM;j++){
			for(k=1;k<=KM;k++){
				B_index[i][j][k] = -1;
				X[i][j][k] = X0 + (i-1)*DX;// x-coordinates from logical i position
				Y[i][j][k] = Y0 + (j-1)*DX;//y-coord from logical j position
				Z[i][j][k] = Z0 + (k-1)*DX;
			}
		}
	}
	
	
	iter = iterstart;
	do{
		sprintf(number, "%d", iter);
			
		if(iter < 10){
			strcpy(num, zero8);
			strcat(num, number);          
		}
		if(iter >= 10 && iter < 100){
			strcpy(num, zero7);
			strcat(num, number);          
		}
		if(iter >= 100 && iter < 1000){
			strcpy(num, zero6);
			strcat(num, number);          
		}
		if(iter >= 1000 && iter < 10000){
			strcpy(num, zero5);	
			strcat(num, number);
		}
		if(iter >= 10000 && iter < 100000){
			strcpy(num, zero4);	
			strcat(num, number);
		}
		if(iter >= 100000 && iter < 1000000){
			strcpy(num, zero3);	
			strcat(num, number);
		}
		if(iter >= 1000000 && iter < 10000000){
			strcpy(num, zero2);	
			strcat(num, number);
		}
		if(iter >= 10000000 && iter < 100000000){
			strcpy(num, zero1);	
			strcat(num, number);
		}
		if(iter >= 100000000)	strcpy(num, number);
		
		
		
		linecount = 0;
		strcpy(name, "SimulationOutputData/FLUIDdomain-");
		strcat(name, num);
		strcat(name, ".csv");
		// Open the file 
		fIN = fopen(name, "r"); 
		  
			// Check if file exists 
			if(fIN == NULL) 
			{ 
				printf("Could not open file %s", name); 
				return 0; 
			} 
		  
			// Extract characters from file and store in character c 
			for (cin = getc(fIN); cin != EOF; cin = getc(fIN)){
				if (cin == '\n') // Increment count if this character is newline 
					linecount = linecount + 1; 
			}
			// Close the file 
		fclose(fIN); 
		printf("The file %s has %d lines\n ", name, linecount); 
		
		strcpy(name, "SimulationOutputData/FLUIDdomain-");
		strcat(name, num);
		strcat(name, ".csv");
		printf("Opening %s \n", name);
		fIN = fopen(name, "r");
			for(m=1;m<=linecount;m++){
				if(m==1) fgets(bulk, 180, fIN);
				else{
					fscanf(fIN,"%f",&xp);
					fscanf(fIN,"%c",&junkchar);
					fscanf(fIN,"%f",&yp);
					fscanf(fIN,"%c",&junkchar);
					fscanf(fIN,"%f",&zp);
					fscanf(fIN,"%c",&junkchar);
					fscanf(fIN,"%d",&bvalue);
					fscanf(fIN,"%c",&junkchar);
					fscanf(fIN,"%f",&UU);
					fscanf(fIN,"%c",&junkchar);
					fscanf(fIN,"%f",&VV);
					fscanf(fIN,"%c",&junkchar);
					fscanf(fIN,"%f",&WW);
					fscanf(fIN,"%c",&junkchar);
					fscanf(fIN,"%f",&DD);					
					fscanf(fIN,"%c",&junkchar);
					fscanf(fIN,"%f",&Bx);
					fscanf(fIN,"%c",&junkchar);
					fscanf(fIN,"%f",&By);
					fscanf(fIN,"%c",&junkchar);
					fscanf(fIN,"%f",&Bz);
					fscanf(fIN,"%c",&junkchar);
					fscanf(fIN,"%d",&junkint);
					fscanf(fIN,"%c",&junkchar);
					fscanf(fIN,"%d",&junkint);
					
					if(bvalue >= 0){
						i = (int)((xp-X0)/DX + 1.01);
						j = (int)((yp-Y0)/DX + 1.01);
						k = (int)((zp-Z0)/DX + 1.01);
						
						B_index[i][j][k] = bvalue;
						
						pp = (DD-1)/3.*10000.;
						
						PRESS[i][j][k] = pp;
						U[i][j][k] = UU;
						V[i][j][k] = VV;
						W[i][j][k] = WW;
					}
				}			
			}
		fclose(fIN);
		
		for(m=0;m<NodeM;m++){
			ii = (int)((Xm[m]-X0)/DX + 1.01);
			jj = (int)((Ym[m]-Y0)/DX + 1.01);
			kk = (int)((Zm[m]-Z0)/DX + 1.01);
			
			i1 = ii - 3; i2 = ii + 3;
			if(i2 > IM) i2 = IM; 
			else if(i1 < 3) i1 = 3;
			j1 = jj - 3; j2 = jj + 3;
			if(j2 > JM) j2 = JM;
			else if(j1 < 1) j1 = 1;
			k1 = kk - 3; k2 = kk + 3;
			if(k2 > KM) k2 = KM;
			else if(k1 < 1) k1 = 1;
			
			mindis = 1.42*DX;
			///find nearest B_index=0 pt. to define the wall pressure on surface mesh
			for(i=i1;i<=i2;i++){
				for(j=j1;j<=j2;j++){
					for(k=k1;k<=k2;k++){
						//if(B_index[i][j][k] > -0.00001 && B_index[i][j][k] < 0.0001){
						if(B_index[i][j][k] > -0.0001){
							distance = sqrtf((Xm[m]-X[i][j][k])*(Xm[m]-X[i][j][k]) + (Ym[m]-Y[i][j][k])*(Ym[m]-Y[i][j][k]) + (Zm[m]-Z[i][j][k])*(Zm[m]-Z[i][j][k]));
							
							if(distance < mindis){
								mindis = distance;
								iii = i; jjj = j; kkk = k;
							}
						}
					}
				}
			}
			
			if(mindis > 1.419999*DX){
				//printf("Error in matching wall mesh to cartesian mesh\n");
				ValidData[m] = 0;
				PRESSW[m] = 0;
			}
			else ValidData[m] = 1;
			
			if(ValidData[m] == 1) PRESSW[m] = PRESS[iii][jjj][kkk]-Pgauge;
			
			if(PRESSW[m] < 1) ValidData[m] = 0;
			
			///WSS calculation
			
			
			xp = Xm[m] + DX*Norm[0][m]; yp = Ym[m] + DX*Norm[1][m]; zp = Zm[m] + DX*Norm[2][m];
			
			///interpolate velocity to xp, yp, zp
			ii = (int)((xp-X0)/DX + 1.01);
			jj = (int)((yp-Y0)/DX + 1.01);
			kk = (int)((zp-Z0)/DX + 1.01);
			
			i1 = ii - 3; i2 = ii + 3;
			if(i2 > IM) i2 = IM;
			else if(i1 < 3) i1 = 3;
			j1 = jj - 3; j2 = jj + 3;
			if(j2 > JM) j2 = JM;
			else if(j1 < 1) j1 = 1;
			k1 = kk - 3; k2 = kk + 3;
			if(k2 > KM) k2 = KM;
			else if(k1 < 1) k1 = 1;
			
			weight = 0;
			check = 1;
			UU = 0; VV = 0; WW = 0;
			for(i=i1;i<=i2;i++){
				for(j=j1;j<=j2;j++){
					for(k=k1;k<=k2;k++){
						distance = sqrtf((xp-X[i][j][k])*(xp-X[i][j][k]) + (yp-Y[i][j][k])*(yp-Y[i][j][k]) + (zp-Z[i][j][k])*(zp-Z[i][j][k]));
						if(distance <= 2*DX){
							weight = (1.+cosf(PI*fabsf(xp-X[i][j][k])/2/DX))*(1.+cosf(PI*fabsf(yp-Y[i][j][k])/2/DX))*(1.+cosf(PI*fabsf(zp-Z[i][j][k])/2/DX))/64.;
							if(B_index[i][j][k] > 0.11 && B_index[i][j][k] < 4997){
								UU += weight*U[i][j][k];
								VV += weight*V[i][j][k];
								WW += weight*W[i][j][k];
							}
						}
					}
				}
			}
			
			
			if(check == 1){
				vnormmag = UU*Norm[0][m] + VV*Norm[1][m] + WW*Norm[2][m];
				
				vtan[0] = UU-vnormmag*Norm[0][m];
				vtan[1] = VV-vnormmag*Norm[1][m];
				vtan[2] = WW-vnormmag*Norm[2][m];
				
				WSS[0][m] = -0.0012*vtan[0]/DX/1.e-06;
				WSS[1][m] = -0.0012*vtan[1]/DX/1.e-06;
				WSS[2][m] = -0.0012*vtan[2]/DX/1.e-06;
			}
			if(check == 0) printf("WSS not calculated on mesh node %d\n",m);
		}
		
		
		
		strcpy(name, "PostProcessedData/Wallpressure");
		strcat(name,"-");
		strcat(name,num);
		strcat(name,".vtk");
		printf("Ouputting surface pressures at iter#%d to %s...\n", iter, name);
		fOUT = fopen(name, "w");
			fprintf(fOUT, "# vtk DataFile Version 4.0\n");
			fprintf(fOUT, "vtk output\n");
			fprintf(fOUT, "ASCII\n");
			fprintf(fOUT, "DATASET POLYDATA\n");
			fprintf(fOUT, "POINTS %d float\n",NodeM);
			
			for(m=0;m<NodeM;m++){
				fprintf(fOUT, "%f %f %f ", Xm[m], Ym[m], Zm[m]);
				if((m+1)%3 == 0) fprintf(fOUT, "\n");	
			}
			
			fprintf(fOUT, "\n");
			
			fprintf(fOUT, "POLYGONS %d %d\n", CellM, CellM*4);
			for(m=0;m<CellM;m++){
				
					fprintf(fOUT, "%d %d %d %d\n", 3, POLYPID[0][m]-1, POLYPID[1][m]-1, POLYPID[2][m]-1);
				//}
			}
			
			fprintf(fOUT, "POINT_DATA %d\n",NodeM);
			fprintf(fOUT, "SCALARS PRESSURE float\n");
			fprintf(fOUT, "LOOKUP_TABLE default\n");
			for(m=0;m<NodeM;m++){
				fprintf(fOUT, "%f ", PRESSW[m]);
				if((m+1)%6 == 0) fprintf(fOUT, "\n");	
			}
			
			fprintf(fOUT, "FIELD FieldData 4\n");
			fprintf(fOUT, "WSSx 1 %d float\n",NodeM);
			for(m=0;m<NodeM;m++){
				fprintf(fOUT, "%f ", WSS[0][m]);
				if((m+1)%6 == 0) fprintf(fOUT, "\n");	
			}
			
			fprintf(fOUT, "WSSy 1 %d float\n",NodeM);
			for(m=0;m<NodeM;m++){
				fprintf(fOUT, "%f ", WSS[1][m]);
				if((m+1)%6 == 0) fprintf(fOUT, "\n");	
			}
			
			fprintf(fOUT, "WSSz 1 %d float\n",NodeM);
			for(m=0;m<NodeM;m++){
				fprintf(fOUT, "%f ", WSS[2][m]);
				if((m+1)%6 == 0) fprintf(fOUT, "\n");	
			}
			
			fprintf(fOUT, "Valid 1 %d int\n",NodeM);
			for(m=0;m<NodeM;m++){
				fprintf(fOUT, "%d ", ValidData[m]);
				if((m+1)%6 == 0) fprintf(fOUT, "\n");	
			}
		fclose(fOUT);
		
		
		
		iter += outputfreq;
	}while(iter <= itermax);
	
	return 0;
}
