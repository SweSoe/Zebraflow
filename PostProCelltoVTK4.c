#include <getopt.h>
#include <sys/time.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#define CapsM 3000
#define NodeM 504
#define EdgeM 1506
#define TriM 1004

#define iterstart 0
#define itermax 720000//set to the maximum iteration level you have simulated

#define outputfreq 10000


float NODECOORD[3][CapsM+1][NodeM+1], NODE_VEL[3][CapsM+1][NodeM+1];
int TRINODE1[TriM+1],TRINODE2[TriM+1],TRINODE3[TriM+1];
int caps, node, tri;
int CellStatus[CapsM+1],iter,CellCount,CellCounter;
int NODE_WALLNEIGH[CapsM+1][NodeM+1];
float F_Total[3][CapsM+1][NodeM+1], NODE_WALLREPUL[3][CapsM+1][NodeM+1], NODE_ARATIO[CapsM+1][NodeM+1];
int REGRESSFLAG[CapsM+1];
int MembIn_Neigh[CapsM+1][NodeM+1], MembIn_NeighFlag[CapsM+1][NodeM+1]/*this flag will be -1 if the node has penetrated into the neighbor surface, otherwise 1 if there is a previously identified neighbor, 0 if there is a previously no identified neighbor*/;
int MembOut_Neigh[CapsM+1][NodeM+1], MembOut_NeighFlag[CapsM+1][NodeM+1]/*this flag will be -1 if the node has penetrated into the neighbor surface, otherwise 1 if there is a previously identified neighbor, 0 if there is a previously no identified neighbor*/;
float F_VOLUME[3][CapsM+1][NodeM+1];

int main(void){
	FILE *fOUT, *fIN;
    char cin;
	int m,nodenumber;
	char number[20], num[40], name[250], bulk[200], junkchar;
	float xp, yp, zp, vx, vy, vz, nwrepulx, nwrepuly, nwrepulz, astrain, fx, fy, fz, fvx, fvy, fvz;
	int nwtri,rflag,MNin,MNinflag,MNout,MNoutflag;
	int cellindex,statval;
	char zero8[] = "00000000", zero7[] = "0000000", zero6[] = "000000", zero5[]="00000", zero4[]="0000", zero3[]="000", zero2[]="00", zero1[]="0";
	int iter;
	
	///Read in nodes from vtk file
	printf("loading ZFISHRBCTRIANGLES.dat file\n");
	fIN = fopen("InputData/ZFISHRBCTRIANGLES.dat", "r");
		fgets(bulk, 30, fIN);//first line states TriM, ignore the line, put into junk char array bulk[]
	
			for(tri=1;tri<=TriM;tri++){
				fscanf(fIN,"%d", &TRINODE1[tri]);
				fscanf(fIN,"%d", &TRINODE2[tri]);
				fscanf(fIN,"%d", &TRINODE3[tri]);
			}
		
	fclose(fIN);
	
	
	iter = iterstart;

	do{
		CellCount = 0;
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
		
		for(caps=1;caps<=CapsM;caps++){
			CellStatus[caps] = -5000;//not in simulation domain
		}
		
		
		strcpy(name, "SimulationOutputData/AllCellMemb-");
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
				CellCount = CellCount + 1; 
		}
		// Close the file 
		fclose(fIN); 
		printf("The file %s has %d lines\n ", name, CellCount); 
		
		CellCount = CellCount - 1;
		CellCount = CellCount/504;
		
		printf("The file %s has %d number of cells\n ", name, CellCount); 
		
		strcpy(name, "SimulationOutputData/AllCellMemb-");
		strcat(name, num);
		strcat(name, ".csv");
		
		printf("Opening %s \n", name);
		fIN = fopen(name, "r");
		fgets(bulk, 180, fIN);
			for(m=1;m<=NodeM*CellCount/*1931*//*969*/;m++){
				if(m==1) nodenumber = 0;
				nodenumber += 1;
				
				fscanf(fIN,"%f",&xp);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&yp);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&zp);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%d",&cellindex);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%d",&statval);//cellstat
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&vx);//vx
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&vy);//vy
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&vz);//vz					
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%d",&nwtri);//vz
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&nwrepulx);//vz
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&nwrepuly);//vz
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&nwrepulz);//vz
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%d",&rflag);//vz
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&astrain);//vz
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&fx);//vx
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&fy);//vy
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&fz);//vz
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&fvx);//vz
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&fvy);//vz
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&fvz);//vz
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%d",&MNoutflag);//vz
				
				CellStatus[cellindex] = statval;
				REGRESSFLAG[cellindex] = rflag;
				NODECOORD[0][cellindex][nodenumber] = xp;
				NODECOORD[1][cellindex][nodenumber] = yp;
				NODECOORD[2][cellindex][nodenumber] = zp;
				
				NODE_VEL[0][cellindex][nodenumber] = vx;
				NODE_VEL[1][cellindex][nodenumber] = vy;
				NODE_VEL[2][cellindex][nodenumber] = vz;
				
				F_Total[0][cellindex][nodenumber] = fx;
				F_Total[1][cellindex][nodenumber] = fy;
				F_Total[2][cellindex][nodenumber] = fz;
				
				F_VOLUME[0][cellindex][nodenumber] = fvx;
				F_VOLUME[1][cellindex][nodenumber] = fvy;
				F_VOLUME[2][cellindex][nodenumber] = fvz;
				
				NODE_WALLNEIGH[cellindex][nodenumber] = nwtri;
				NODE_WALLREPUL[0][cellindex][nodenumber] = nwrepulx;
				NODE_WALLREPUL[1][cellindex][nodenumber] = nwrepuly;
				NODE_WALLREPUL[2][cellindex][nodenumber] = nwrepulz;
				NODE_ARATIO[cellindex][nodenumber] = astrain;
				
				MembIn_Neigh[cellindex][nodenumber] = MNin;
				MembIn_NeighFlag[cellindex][nodenumber] = MNinflag;
				
				MembOut_Neigh[cellindex][nodenumber] = MNout;
				MembOut_NeighFlag[cellindex][nodenumber] = MNoutflag;
				
				if(statval == 4){
					if(zp < 0){//DAP parent
						if(zp > -56) NODECOORD[2][cellindex][nodenumber] += -60.;
					}
					else{//PCV parent
						if(zp < 362.5) NODECOORD[2][cellindex][nodenumber] += 60.;
					}
				}
				
				if(m != 1 && (m%NodeM == 0)) nodenumber = 0;
			}
		fclose(fIN);
		
		
		strcpy(name, "PostProcessedData/AllCellMemb");
		strcat(name,"-");
		strcat(name,num);
		strcat(name,".vtk");
		printf("Ouputting fluid at iter#%d to %s...\n", iter, name);
		fOUT = fopen(name, "w");
			fprintf(fOUT, "# vtk DataFile Version 4.0\n");
			fprintf(fOUT, "vtk output\n");
			fprintf(fOUT, "ASCII\n");
			fprintf(fOUT, "DATASET UNSTRUCTURED_GRID\n");
			fprintf(fOUT, "POINTS %d double\n",CellCount*504);
			for(caps=1;caps<=CapsM;caps++){
				if(CellStatus[caps] != -5000){
					for(node=1;node<=NodeM;node++){
						fprintf(fOUT, "%.7e %.7e %.7e\n",NODECOORD[0][caps][node],NODECOORD[1][caps][node],NODECOORD[2][caps][node]);
					}
				}
			}
					
			fprintf(fOUT, "\nCELLS %d %d\n", CellCount*1004, CellCount*4016);
			CellCounter = 0;
			for(caps=1;caps<=CapsM;caps++){
				if(CellStatus[caps] != -5000){		
					for(tri=1;tri<=TriM;tri++){
						fprintf(fOUT, "3 %d %d %d\n", CellCounter*504+TRINODE1[tri]-1, CellCounter*504+TRINODE2[tri]-1, CellCounter*504+TRINODE3[tri]-1);
					}
					CellCounter += 1;
				}
				
			}
					
			fprintf(fOUT, "\nCELL_TYPES %d\n", CellCount*1004);
			for(caps=1;caps<=CapsM;caps++){
				if(CellStatus[caps] != -5000){	
					for(tri=1;tri<=TriM;tri++){
						fprintf(fOUT, "5\n");
					}
				}
			}
			
			
			fprintf(fOUT, "\nPOINT_DATA %d\n",CellCount*504);
			fprintf(fOUT, "SCALARS RBC int\n");
			fprintf(fOUT, "LOOKUP_TABLE default\n");
			for(caps=1;caps<=CapsM;caps++){
				if(CellStatus[caps] != -5000){	
					for(node=1;node<=NodeM;node++){
						fprintf(fOUT, "%d\n", caps);
					}
				}
			}
			
			fprintf(fOUT, "\nFIELD FieldData 11\n");
			fprintf(fOUT, "vx 1 %d double\n",CellCount*504);
			for(caps=1;caps<=CapsM;caps++){
				if(CellStatus[caps] != -5000){	
					for(node=1;node<=NodeM;node++){
						fprintf(fOUT, "%.7e\n", NODE_VEL[0][caps][node]);
					}
				}
			}
			
			fprintf(fOUT, "\nvy 1 %d double\n",CellCount*504);
			for(caps=1;caps<=CapsM;caps++){
				if(CellStatus[caps] != -5000){	
					for(node=1;node<=NodeM;node++){
						fprintf(fOUT, "%.7e\n", NODE_VEL[1][caps][node]);
					}
				}
			}
			
			fprintf(fOUT, "\nvz 1 %d double\n",CellCount*504);
			for(caps=1;caps<=CapsM;caps++){
				if(CellStatus[caps] != -5000){
					for(node=1;node<=NodeM;node++){
						fprintf(fOUT, "%.7e\n", NODE_VEL[2][caps][node]);
					}
				}
			}
			
			fprintf(fOUT, "\nrfx 1 %d double\n",CellCount*504);
			for(caps=1;caps<=CapsM;caps++){
				if(CellStatus[caps] != -5000){	
					for(node=1;node<=NodeM;node++){
						fprintf(fOUT, "%.7e\n", NODE_WALLREPUL[0][caps][node]);
					}
				}
			}
			
			fprintf(fOUT, "\nrfy 1 %d double\n",CellCount*504);
			for(caps=1;caps<=CapsM;caps++){
				if(CellStatus[caps] != -5000){	
					for(node=1;node<=NodeM;node++){
						fprintf(fOUT, "%.7e\n", NODE_WALLREPUL[1][caps][node]);
					}
				}
			}
			
			fprintf(fOUT, "\nrfz 1 %d double\n",CellCount*504);
			for(caps=1;caps<=CapsM;caps++){
				if(CellStatus[caps] != -5000){	
					for(node=1;node<=NodeM;node++){
						fprintf(fOUT, "%.7e\n", NODE_WALLREPUL[2][caps][node]);
					}
				}
			}
			
			fprintf(fOUT, "\naratio 1 %d double\n",CellCount*504);
			for(caps=1;caps<=CapsM;caps++){
				if(CellStatus[caps] != -5000){	
					for(node=1;node<=NodeM;node++){
						fprintf(fOUT, "%.7e\n", NODE_ARATIO[caps][node]);
					}
				}
			}
			
			fprintf(fOUT, "\nFx 1 %d double\n",CellCount*504);
			for(caps=1;caps<=CapsM;caps++){
				if(CellStatus[caps] != -5000){	
					for(node=1;node<=NodeM;node++){
						fprintf(fOUT, "%.7e\n", F_Total[0][caps][node]);
					}
				}
			}
			
			fprintf(fOUT, "\nFy 1 %d double\n",CellCount*504);
			for(caps=1;caps<=CapsM;caps++){
				if(CellStatus[caps] != -5000){	
					for(node=1;node<=NodeM;node++){
						fprintf(fOUT, "%.7e\n", F_Total[1][caps][node]);
					}
				}
			}
			
			fprintf(fOUT, "\nFz 1 %d double\n",CellCount*504);
			for(caps=1;caps<=CapsM;caps++){
				if(CellStatus[caps] != -5000){	
					for(node=1;node<=NodeM;node++){
						fprintf(fOUT, "%.7e\n", F_Total[2][caps][node]);
					}
				}
			}
			
			
			fprintf(fOUT, "\nMNoutFlag 1 %d int\n",CellCount*504);
			for(caps=1;caps<=CapsM;caps++){
				if(CellStatus[caps] != -5000){	
					for(node=1;node<=NodeM;node++){
						fprintf(fOUT, "%d\n", MembOut_NeighFlag[caps][node]);
					}
				}
			}
		fclose(fOUT);
		
		iter += outputfreq;
		
	}while(iter <= itermax);
	
	return 0;
}
