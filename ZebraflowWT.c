/*
 * ZebraflowWT.c
 * 
 * Copyright 2023 swesoe <swesoe@gmail.com>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */

#include <getopt.h>
#include <sys/time.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>


#define PI 3.14159265358979
///Thermodynamics constants
#define TEMPERATURE 300//kelvins
#define BOLTZ 1.3806488e-023//Boltzmann constant:  m^2 kg s^-2 K^-1  ==> TEMPERATURE*BOLTZ = Free energy constant (Nm)

///LBM numerical resolution [in um]
#define	DX	0.5//grid step size
#define DT 0.5//time step size

///Physical scaling
#define T_SCALE	(DT*1.e-06)
#define L_SCALE	(DX*1.e-06)
#define DENphys 1000.
#define DEN_SCALE 10.//10000.//density scaled up by 10times to keep LBM tau close to 1
#define MUDynamic 0.0012//dynamic viscosity of plasma: use the correct dynamic viscosity before scaling
#define MU (MUDynamic/DENphys/DEN_SCALE)//kinematic viscosity
#define UNITVOL (L_SCALE*L_SCALE*L_SCALE)
#define UNITMASS (UNITVOL*DEN_SCALE*DENphys)

///LBM scaling - LBM scaling is applied on top of physical scaling, these parameters control the property of the LBM transport.
///Generally, Tau should be kept around 1 to limit numerical diffusion at boundaries and regions where body forces are applied to fluid.
///Keeping Tau close to 1 can be achieved by: 1) physically scaling the density for a Stokes flow regime (Reynolds number < 1) -- take care not to change the Reynolds regime!
///                                        or 2) by reducing DT more than DX, ie. keeping MU_SCALE low -- this is computationally expensive!
#define C_LBM (L_SCALE / T_SCALE)//LBM streaming velocity
#define	CS (C_LBM/(sqrt(3.)))//LBM speed of sound
#define	EQ_A (1./(CS*CS))//LBM equilibrium ddf coefficient 1
#define	EQ_B (0.5/(CS*CS*CS*CS))//LBM equilibrium ddf coefficient 2
#define	EQ_C (-0.5/(CS*CS))//LBM equilibrium ddf coefficient 3	
#define	MU_SCALE (3.*T_SCALE/L_SCALE/L_SCALE)//LBM-scaled viscosity
#define	TauPlasma  (0.5 + MU_SCALE*MU)//BGK time-relaxation parameter for plasma phase
#define	TauCytosol (0.5 + MU_SCALE*MU*5.)//BGK time-relaxation parameter for cytoplasm phase	
#define	LBM_BFORCE_SCALE (1./3.*L_SCALE/T_SCALE/T_SCALE)//LBM-scaled body force (acceleration)
#define	LBM_PRESSURE_SCALE (CS*CS*DEN_SCALE*DENphys)//LBM-scaled pressure	
#define	VEL_SCALE (L_SCALE/T_SCALE)//LBM-scaled velocity	
// to obtain physical value after LBM calculation, multiply the LBM value by the scaling value

///LBM-BGK streaming and collision parameters
#define CM 19
#define	We1_3D	1./18.
#define	We2_3D	1./36.
#define	We3_3D	1./3.



/// The lbm grid has four domains. The ROI domain is where the zebrafish trunk network flow is modelled. The three special domains are periodic domains where RBCs are periodically re-entering the domain. 

///LBM ROI DOMAIN CARTESIAN GRID DIMENSIONS 
#define	KM	723//	  //number of lattice nodes in Z direction
#define	JM	547 //number of lattice nodes in Y direction
#define	IM	249//422  //number of lattice nodes in Z direction
///LBM ROI DOMAIN Minimum coordinates
#define Z0 -27.// Zmax = -27. + 0.5*(723-1) = 334; Z@k=KM-3 = -27 + 0.5*(720-1) = 332.5; Z@k=3 = -27 + 0.5*(3-1) = -26
#define Y0 29.// Ymax = 29. + 0.5*(547-1) = 302
#define X0 -1.// Xmax = -1 + 0.5*(249-1) = 123 


///LBMspecial zone 1: DA parent zone
#define	KM_1 121//67	  //number of x-nodes
#define	JM_1 56 //number of y-nodes
#define	IM_1 69//422  //number of z nodes
#define Z0_1 -86.// Xmax = -86. + 0.5*(121-1) = -26
#define Y0_1 210.// Ymax = 210. + 0.5*(56-1) = 237.5
#define X0_1 70.5// Zmax = 70.5 + 0.5*(69-1) = 104.5

///LBMspecial zone 2: PCV parent1 zone
#define	KM_2 121//67	  //number of x-nodes
#define	JM_2 60 //number of y-nodes
#define	IM_2 64//422  //number of z nodes
#define Z0_2 332.5// Xmax = 332.5 + 0.5*(121-1) = 392.5
#define Y0_2 232.// Ymax = 232. + 0.5*(60-1) = 261.5
#define X0_2 62.// Zmax = 62. + 0.5*(64-1) = 93.5

///LBMspecial zone 3: PCV parent2 zone
#define	KM_3 121//67	  //number of x-nodes
#define	JM_3 27 //number of y-nodes
#define	IM_3 42//422  //number of z nodes
#define Z0_3 332.5// Xmax = 332. + 0.5*(121-1) = 392.5
#define Y0_3 209.5// Ymax = 209.5 + 0.5*(27-1) = 222.5
#define X0_3 58.// Zmax = 58. + 0.5*(42-1) = 78.5

///RBC parameters
#define MasterCapsM 47//deprecated parameter
#define CapsM 3000// maximum number of RBCs declared into system memory - this must exceed actual number in circulation. As a gauge: 1000 cells populates to about 14 - 18% hematocrit in the caudal artery and vein
#define NodeM 504// Number of mesh vertices on one RBC membrane shell. 
#define EdgeM 1506//Number of mesh edges on one RBC membrane shell. The spectrin (mesh edge) here is 100 - 200 nm range. This is a coarse-grained spectrin model, actual spectrin in RBCs are 30 - 50 nm. For details into coarse-graining see the work of Fedosov (2011): https://doi.org/10.1073/pnas.110121010
#define TriM 1004//Number of triangle elements on one RBC membrane shell
#define NListM 9//maximum no. of neighbor nodes
#define NodeM_N 121//Number of mesh vertices in one RBC nucleus (nucleus is enveloped inside the RBC shell)
#define EdgeM_N 544//Number of mesh edges in one RBC nucleus

#define gamma_t 0.37306e-07//membrane viscosity - not used in the paper
#define gamma_p 1.e-07//collision dissipation parameter
#define MCV 97.5671e-018//RBC volume
#define RBC_SURF 144.916e-012//RBC membrane surface area
#define kBend 1.3856406e-19//this parameter value translates to Bending modulus 1.2e-19 J: kBend = 1.2e-019*2./sqrt(3.);
#define Es 3.e-06//N m^-1 Membrane shear modulus
#define kd 3.6e-05//local surface area incompressibility coefficient, N m^-1. Keep the value of (2Es+kd+ka) > 10*Es 
#define ka 7.7e-05//global surface area incompressibility coefficient, N m^-1 
#define kv 55.//110//22.//000//220000//3.e-004//capsule volume incompressibility coefficient, N m^-2
#define e_ratio_eq 0.3154574132//ratio of native spectrin tetramer length in CSK network to maximum extensible length of spectrin tetramer
#define MViscOn 0//Membrane vsic flag 0 is off, 1 is on

///Lumen network wall mesh
#define NodeW 1107042//Number of mesh vertices on ROI lumen wall mesh
#define CellW 2212824//Number of triangle elements on ROI wall mesh
#define NodeW_1 29282//Number of mesh vertices on DA parent lumen wall mesh
#define CellW_1 58080//Number of triangle elements on DA parent lumen wall mesh
#define NodeW_2 29040//Number of mesh vertices on CVP1 parent lumen wall mesh
#define CellW_2 57600//Number of triangle elements on CVP1 parent lumen wall mesh
#define NodeW_3 15488//Number of mesh vertices on CVP2 parent lumen wall mesh
#define CellW_3 30720//Number of triangle elements on CVP2 parent lumen wall mesh

///Fluid Structure Interaction Parameter
#define StencilWidth 2//Dirac delta width for vel interpolation and force spreading
#define WallStencilWidth 2//Dirac delta width for RBC repulsion from wall

///Solver and output control
#define SOLVERSTAGGER 50//the LBM solver updates with a higher freqeuncy than the RBC-CGSM and IBM solvers
#define iterrestart 0//31680000
#define outputfreq 10000
#define initoutfreq 1000
#define sliceoutputfreq 10000
#define itermax (iterrestart+5*720000)//simulate for 5 cardiac cycles - I recommend at least 10 cardiac cycles of flow simulation for analysis of developed flow pattern
#define iterrelax 50

///RBC membrane and nucleus vertex neighbor record list size - we keep an updated array for each particle in the simulation for quick search of its nearby neighbors during RBC-RBC interaction calculatations
#define NucMembMAX 1004
#define EulMAX 800
//#define MM 3000000//Maximum number of rows in largest file to be loaded
void LBM_COPY(void);
void LBM_STREAM(void);
void LBM_COLLISION(void);
void LBM_BC(void);
void LBM_FLUID(void);
void LBM_DOM_SETUP(void);
void LBM_BC_SETUP(void);

///Simulation system variables
int iter, iter_local, plasmaloop;
int simfailure;

///LBM variables
float X[IM+1][JM+1][KM+1], Y[IM+1][JM+1][KM+1], Z[IM+1][JM+1][KM+1];
float U[8500000], V[8500000], W[8500000], D[8500000];//macrovariables: velocities U, V and W, and density D
float F[CM+1][8500000], F_OLD[CM+1][8500000], FEQ[CM+1][8500000];//microvariables: density distribution function (ddf) F, ddf copy F_OLD, equilibirum ddf FEQ for LBM-BGK collision
int B_index[IM+1][JM+1][KM+1], Bcopy[IM+1][JM+1][KM+1], W_index[IM+1][JM+1][KM+1];//LBM grid phase property indicator B_index, copy of phase indicator
float Wdistance[IM+1][JM+1][KM+1];
float Xlattice[3][CM+1], UC[CM+1], VC[CM+1], WC[CM+1];
float Wcoeff[CM+1];//LBM state variables
float BODYFORCE[3][8500000];
int BCLIST[1200000][4], BCNLIST[1200000][3], BCMAX, DOMLIST[8500000][3], DOMMAX;
int DOMID[IM+1][JM+1][KM+1];
float Uvoid[IM+1][JM+1][KM+1], Vvoid[IM+1][JM+1][KM+1], Wvoid[IM+1][JM+1][KM+1];
int BCMAXold,DOMMAXold;

///LBMspecial zone 1: DA parent zone variables
float X_1[IM_1+1][JM_1+1][KM_1+1], Y_1[IM_1+1][JM_1+1][KM_1+1], Z_1[IM_1+1][JM_1+1][KM_1+1];
float U_1[400000], V_1[400000], W_1[400000], D_1[400000];
float F_1[CM+1][400000], F_OLD_1[CM+1][400000], FEQ_1[CM+1][400000];
int B_index_1[IM_1+1][JM_1+1][KM_1+1], W_index_1[IM_1+1][JM_1+1][KM_1+1];
float Wdistance_1[IM_1+1][JM_1+1][KM_1+1];
float BODYFORCE_1[3][400000];
int BCLIST_1[60000][4], BCNLIST_1[60000][3], BCMAX_1, DOMLIST_1[400000][3], DOMMAX_1, vBCNLIST_1[60000][3], vBCMAX_1;
int DOMID_1[IM_1+1][JM_1+1][KM_1+1];
float Uvoid_1[IM_1+1][JM_1+1][KM_1+1], Vvoid_1[IM_1+1][JM_1+1][KM_1+1], Wvoid_1[IM_1+1][JM_1+1][KM_1+1];

///LBMspecial zone 2: DPCV parent zone1 variables
float X_2[IM_2+1][JM_2+1][KM_2+1], Y_2[IM_2+1][JM_2+1][KM_2+1], Z_2[IM_2+1][JM_2+1][KM_2+1];
float U_2[400000], V_2[400000], W_2[400000], D_2[400000];
float F_2[CM+1][400000], F_OLD_2[CM+1][400000], FEQ_2[CM+1][400000];
int B_index_2[IM_2+1][JM_2+1][KM_2+1], W_index_2[IM_2+1][JM_2+1][KM_2+1];
float Wdistance_2[IM_2+1][JM_2+1][KM_2+1];
float BODYFORCE_2[3][400000];
int BCLIST_2[60000][4], BCNLIST_2[60000][3], BCMAX_2, DOMLIST_2[400000][3], DOMMAX_2, vBCNLIST_2[60000][3], vBCMAX_2;
int DOMID_2[IM_2+1][JM_2+1][KM_2+1];
float Uvoid_2[IM_2+1][JM_2+1][KM_2+1], Vvoid_2[IM_2+1][JM_2+1][KM_2+1], Wvoid_2[IM_2+1][JM_2+1][KM_2+1];

///LBMspecial zone 3: PCV parent zone2 variables
float X_3[IM_3+1][JM_3+1][KM_3+1], Y_3[IM_3+1][JM_3+1][KM_3+1], Z_3[IM_3+1][JM_3+1][KM_3+1];
float U_3[120000], V_3[120000], W_3[120000], D_3[120000];
float F_3[CM+1][120000], F_OLD_3[CM+1][120000], FEQ_3[CM+1][120000];
int B_index_3[IM_3+1][JM_3+1][KM_3+1], W_index_3[IM_3+1][JM_3+1][KM_3+1];
float Wdistance_3[IM_3+1][JM_3+1][KM_3+1];
float BODYFORCE_3[3][120000];
int BCLIST_3[25000][4], BCNLIST_3[25000][3], BCMAX_3, DOMLIST_3[120000][3], DOMMAX_3, vBCNLIST_3[25000][3], vBCMAX_3;
int DOMID_3[IM_3+1][JM_3+1][KM_3+1];
float Uvoid_3[IM_3+1][JM_3+1][KM_3+1], Vvoid_3[IM_3+1][JM_3+1][KM_3+1], Wvoid_3[IM_3+1][JM_3+1][KM_3+1];

///RBC membrane variables
float F_AREA_GLOBAL[CapsM+1], F_VOLUME[3][CapsM+1][NodeM+1], F_AGG[3][CapsM+1][NodeM+1], F_Total[3][CapsM+1][NodeM+1];
int TRINODE1[TriM+1], TRINODE2[TriM+1], TRINODE3[TriM+1], AREA_UPDATE[CapsM+1][TriM+1], EDGENODE_LIST[EdgeM+1][6+1], NEIGHNODE_LIST[NodeM+1][NListM+1], NEIGHTRI_LIST[NodeM+1][NListM+1];
double NODE_VEL[3][CapsM+1][NodeM+1], NODE_VEL_OLD[3][CapsM+1][NodeM+1], NODECOORD[3][CapsM+1][NodeM+1], NODECOORDMASTER[3][MasterCapsM+1][NodeM+1];
float NODE_AREA[CapsM+1][NodeM+1], NODE_AREA0[NodeM+1], NODE_NORMAL[3][CapsM+1][NodeM+1], TRI_AREA[CapsM+1][TriM+1], TRI_AREA_eq[TriM+1], AreaExpanRatio[CapsM+1][TriM+1], TRI_NORMAL[3][CapsM+1][TriM+1], TRI_CENTROID[3][CapsM+1][TriM+1];
float CAPS_CENTROID[3][CapsM+1], CAPS_AREA[CapsM+1];
float LMAX[EdgeM+1], Plength[EdgeM+1], kp[EdgeM+1], Theta0[EdgeM+1], SINETHETA0[EdgeM+1], COSINETHETA0[EdgeM+1];
float F_VOLUME_SCALAR[CapsM+1], CAPS_VOL[CapsM+1];
float beta;
float CAPS_MAX[3][CapsM+1],CAPS_MIN[3][CapsM+1];
int NODE_WALLNEIGH[CapsM+1][NodeM+1];
float NODE_WALLREPUL[3][CapsM+1][NodeM+1];
float AWeight1[19+1], AWeight2[19+1], AWeight3[19+1], EWeight1[4+1], EWeight2[4+1];

///RBC nucleus variables
float F_Total_N[3][CapsM+1][NodeM_N+1],NODE_WALLREPUL_N[3][CapsM+1][NodeM_N+1];
double NODE_VEL_N[3][CapsM+1][NodeM_N+1], NODE_VEL_OLD_N[3][CapsM+1][NodeM+1], NODECOORD_N[3][CapsM+1][NodeM_N+1], NODECOORDMASTER_N[3][MasterCapsM+1][NodeM_N+1];
float LMAX_N[EdgeM_N+1], Plength_N[EdgeM_N+1], kp_N[EdgeM_N+1];
int EDGENODE_LIST_N[EdgeM_N+1][2+1], NUCLNEIGH[NucMembMAX+1][CapsM+1][NodeM_N+1];

///Lumen Wall variables
float WALLNODECOORD[3][NodeW+1], WALLCELLNORMALS[3][CellW+1], WALL_CENTROID[3][CellW+1], WALL_AREA[CellW+1];
int WALLCELLNODE[3][CellW+1], WALLNODETOCELLNEIGHBORLIST[12+1][NodeW+1];

float WALLNODECOORD_1[3][NodeW_1+1], WALLCELLNORMALS_1[3][CellW_1+1], WALL_CENTROID_1[3][CellW_1+1], WALL_AREA_1[CellW_1+1];
int WALLCELLNODE_1[3][CellW_1+1], WALLNODETOCELLNEIGHBORLIST_1[6+1][NodeW_1+1];

float WALLNODECOORD_2[3][NodeW_2+1], WALLCELLNORMALS_2[3][CellW_2+1], WALL_CENTROID_2[3][CellW_2+1], WALL_AREA_2[CellW_2+1];
int WALLCELLNODE_2[3][CellW_2+1], WALLNODETOCELLNEIGHBORLIST_2[6+1][NodeW_2+1];

float WALLNODECOORD_3[3][NodeW_3+1], WALLCELLNORMALS_3[3][CellW_3+1], WALL_CENTROID_3[3][CellW_3+1], WALL_AREA_3[CellW_3+1];
int WALLCELLNODE_3[3][CellW_3+1], WALLNODETOCELLNEIGHBORLIST_3[6+1][NodeW_3+1];

///Cell population logic control
int CellStatus[CapsM+1], CAPS_REMOVE[CapsM+1];
// CellStatus:
// -5000 means not in simulation
//     0 means ROI cell that is slave to parent cell movement
//     1 means ROI cell that is independently moving
//     2 means parent cell in repeat domain interior
//     3 means parent cell at periodic crossing - start cloning new ROI cell
//     4 means parent cell at periodic crossing

float ZM_1, ZM_2, ZM_3, XM, YM, ZM, relaxfactor, relaxfactor2;
int CLONEFLAG[344+1], REGRESSFLAG[CapsM+1];//if regressflag != 0 for parent cells then it cannot be cloned; if regressflag != 0 for roi cells then they cannot be removed from domain
int NucMembUpdate[CapsM+1];
int EulNeigh[EulMAX+1][8500000], EulNeigh_1[EulMAX+1][400000], EulNeigh_2[EulMAX+1][400000], EulNeigh_3[EulMAX+1][120000];
int NUC_Neigh[CapsM+1][NodeM_N+1], NUC_NeighFlag[CapsM+1][NodeM_N+1]/*this flag will be -1 if the node has penetrated into the neighbor surface, otherwise 1 if there is a previously identified neighbor, 0 if there is a previously no identified neighbor*/;

int MembIn_Neigh[CapsM+1][NodeM+1], MembIn_NeighFlag[CapsM+1][NodeM+1]/*this flag will be -1 if the node has penetrated into the neighbor surface, otherwise 1 if there is a previously identified neighbor, 0 if there is a previously no identified neighbor*/;
int MembOut_Neigh[CapsM+1][NodeM+1], MembOut_NeighFlag[CapsM+1][NodeM+1]/*this flag will be -1 if the node has penetrated into the neighbor surface, otherwise 1 if there is a previously identified neighbor, 0 if there is a previously no identified neighbor*/;
int RBCSOLVERFLAG;
int RestartVal[IM+1][JM+1][KM+1];

int maxNodeW;
int eulmaxreport;
//There are 7451359 positive domain nodes
//There are 6424335 fluid nodes
//There are 1027024 boundary nodes

void OUTPUT(void){
	FILE *fOUT, *fIN;
	int o,i,j,k,c,cc,ii,jj,kk,l,m,iii,jjj,kkk;
	float UU,VV,WW,DD, tauLBM;
	char number[15], num[20], name[250], bulk[200], junkchar;
	float xp, yp, zp, xc, yc, zc;
	int i1, i2, j1, j2, k1, k2;
	float distance;
	float mindis;
	int bvalue;
	char zero8[] = "00000000", zero7[] = "0000000", zero6[] = "000000", zero5[]="00000", zero4[]="0000", zero3[]="000", zero2[]="00", zero1[]="0";
	float U_bc, V_bc, W_bc, D_bc, feq_bc;
	float pressure, simtime, cycletime, taucycle;
	float junkfloat;
	int junkint, ii1, jj1, kk1, ii2, jj2, kk2, ncheck;
	int caps, node, tri, n1, n2, n3, n4, tri1, tri2, edge;
	float n1coord[3], n2coord[3], n3coord[3], n4coord[3], e1v[3], e2v[3], e3v[3], Ax, Ay, Az, AreaLocal, tnormal[3];
	float e_ratio, f_n1[3], f_n2[3], f_n3[3], f_n4[3], wlc, comp, dx1_12, dy1_12, dz1_12, edge1_length;
	float betabend, costheta, sintheta, theta, sintheta0, costheta0;
	float area_diff[3], centroid_diff[3], tri1_n[3], tri2_n[3], angle_check, tri1_localcentroid[3], tri2_localcentroid[3], tri1_area, tri2_area;
	
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
			
			//for(caps=1;caps<=344;caps++){
			//caps = 144;
			
			//}
			
			printf("Iteration%d BCmax is %d\n",iter,BCMAX);
			
			strcpy(name, "SimulationOutputData/VolumeAreaHistory.dat");
			fOUT = fopen(name, "a");
				if(iter == 0) fprintf(fOUT, "Iteration#,Cell,Vol,Area,CenX,CenY,CenZ,CellState,RegressFlag\n");
				for(caps=1;caps<=CapsM;caps++){
					if(CellStatus[caps] > -5000) fprintf(fOUT, "%d,%d,%.7e,%.7e,%.7e,%.7e,%.7e,%d,%d\n",iter,caps, CAPS_VOL[caps], CAPS_AREA[caps], CAPS_CENTROID[0][caps], CAPS_CENTROID[1][caps], CAPS_CENTROID[2][caps], CellStatus[caps], REGRESSFLAG[caps]);
				}
			fclose(fOUT);
		
			
			
			strcpy(name, "SimulationOutputData/FLUIDdomain");
			strcat(name,"-");
			strcat(name,num);
			strcat(name,".csv");
			printf("Ouputting fluid at iter#%d to %s...\n", iter, name);
			fOUT = fopen(name, "w");
			fprintf(fOUT, "X,Y,Z,Bindex,U,V,W,D,Bx,By,Bz,EulNeigh,Windex\n");
				for(i=1;i<=IM;i++){
					for(j=1;j<=JM;j++){
						for(k=1;k<=KM;k++){
							if(B_index[i][j][k] >= 0){
								o = DOMID[i][j][k];
								fprintf(fOUT, "%f,%f,%f,%d,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%d,%d\n", X[i][j][k], Y[i][j][k], Z[i][j][k], B_index[i][j][k], U[o], V[o], W[o], D[o], BODYFORCE[0][o], BODYFORCE[1][o], BODYFORCE[2][o], EulNeigh[0][o], W_index[i][j][k]);
								
							}
						}
					}
				}
			fclose(fOUT);
			
			strcpy(name, "SimulationOutputData/DAparentdomain");
			strcat(name,"-");
			strcat(name,num);
			strcat(name,".csv");
			printf("Ouputting fluid at iter#%d to %s...\n", iter, name);
			fOUT = fopen(name, "w");
			fprintf(fOUT, "X,Y,Z,Bindex,U,V,W,D,Bx,By,Bz,EulNeigh,Windex\n");
				for(i=1;i<=IM_1;i++){
					for(j=1;j<=JM_1;j++){
						for(k=1;k<=KM_1;k++){
							if(B_index_1[i][j][k] >= 0){
								o = DOMID_1[i][j][k];
								fprintf(fOUT, "%f,%f,%f,%d,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%d,%d\n", X_1[i][j][k], Y_1[i][j][k], Z_1[i][j][k], B_index_1[i][j][k], U_1[o], V_1[o], W_1[o], D_1[o], BODYFORCE_1[0][o], BODYFORCE_1[1][o], BODYFORCE_1[2][o], EulNeigh_1[0][o], W_index_1[i][j][k]);
								
							}
						}
					}
				}
			fclose(fOUT);
			
			strcpy(name, "SimulationOutputData/PCV1parentdomain");
			strcat(name,"-");
			strcat(name,num);
			strcat(name,".csv");
			printf("Ouputting fluid at iter#%d to %s...\n", iter, name);
			fOUT = fopen(name, "w");
			fprintf(fOUT, "X,Y,Z,Bindex,U,V,W,D,Bx,By,Bz,EulNeigh,Windex\n");
				for(i=1;i<=IM_2;i++){
					for(j=1;j<=JM_2;j++){
						for(k=1;k<=KM_2;k++){
							if(B_index_2[i][j][k] >= 0){
								o = DOMID_2[i][j][k];
								fprintf(fOUT, "%f,%f,%f,%d,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%d,%d\n", X_2[i][j][k], Y_2[i][j][k], Z_2[i][j][k], B_index_2[i][j][k], U_2[o], V_2[o], W_2[o], D_2[o], BODYFORCE_2[0][o], BODYFORCE_2[1][o], BODYFORCE_2[2][o], EulNeigh_2[0][o], W_index_2[i][j][k]);
								
							}
						}
					}
				}
			fclose(fOUT);
			
			strcpy(name, "SimulationOutputData/PCV2aarentdomain");
			strcat(name,"-");
			strcat(name,num);
			strcat(name,".csv");
			printf("Ouputting fluid at iter#%d to %s...\n", iter, name);
			fOUT = fopen(name, "w");
			fprintf(fOUT, "X,Y,Z,Bindex,U,V,W,D,Bx,By,Bz,EulNeigh,Windex\n");
				for(i=1;i<=IM_3;i++){
					for(j=1;j<=JM_3;j++){
						for(k=1;k<=KM_3;k++){
							if(B_index_3[i][j][k] >= 0){
								o = DOMID_3[i][j][k];
								fprintf(fOUT, "%f,%f,%f,%d,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%d,%d\n", X_3[i][j][k], Y_3[i][j][k], Z_3[i][j][k], B_index_3[i][j][k], U_3[o], V_3[o], W_3[o], D_3[o], BODYFORCE_3[0][o], BODYFORCE_3[1][o], BODYFORCE_3[2][o], EulNeigh_3[0][o], W_index_3[i][j][k]);
							}
						}
					}
				}
			fclose(fOUT);
			
			strcpy(name, "SimulationOutputData/AllCellMemb");
			strcat(name,"-");
			strcat(name,num);
			strcat(name,".csv");
			printf("Ouputting fluid at iter#%d to %s...\n", iter, name);
			fOUT = fopen(name, "w");
				fprintf(fOUT, "X,Y,Z,Cell,CellStatus,Vx,Vy,Vz,WallTri,WallRepulfX,WallRepulfY,WallRepulfZ,regressflag,Aratio,Fx,Fy,Fz,FVx,FVy,FVz,OutNeighFlag\n");
				for(caps=1;caps<=CapsM;caps++){
					if(CellStatus[caps] > -5000){
						for(node=1;node<=NodeM;node++){
							fprintf(fOUT, "%f,%f,%f,%d,%d,%.7e,%.7e,%.7e,%d,%.7e,%.7e,%.7e,%d,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%d\n", NODECOORD[0][caps][node], NODECOORD[1][caps][node], NODECOORD[2][caps][node], caps, CellStatus[caps], NODE_VEL[0][caps][node], NODE_VEL[1][caps][node], NODE_VEL[2][caps][node], NODE_WALLNEIGH[caps][node],NODE_WALLREPUL[0][caps][node],NODE_WALLREPUL[1][caps][node],NODE_WALLREPUL[2][caps][node],REGRESSFLAG[caps],NODE_AREA[caps][node]/NODE_AREA0[node],F_Total[0][caps][node],F_Total[1][caps][node],F_Total[2][caps][node],F_VOLUME[0][caps][node],F_VOLUME[1][caps][node],F_VOLUME[2][caps][node],MembOut_NeighFlag[caps][node]);
							
						}
					}
				}
			fclose(fOUT);
			printf("Finished ouputting fluid at iter#%d to %s...\n", iter, name);
			
			strcpy(name, "SimulationOutputData/AllCellNucleus");
			strcat(name,"-");
			strcat(name,num);
			strcat(name,".csv");
			printf("Ouputting fluid at iter#%d to %s...\n", iter, name);
			fOUT = fopen(name, "w");
				fprintf(fOUT, "X,Y,Z,Cell,CellStatus,Vx,Vy,Vz,Fx,Fy,Fz,Neigh,Neighflag\n");
				for(caps=1;caps<=CapsM;caps++){
					if(CellStatus[caps] > -5000){
						for(node=1;node<=NodeM_N;node++){
							fprintf(fOUT, "%f,%f,%f,%d,%d,%.7e,%.7e,%.7e,%.7e,%.7e,%.7e,%d,%d\n", NODECOORD_N[0][caps][node], NODECOORD_N[1][caps][node], NODECOORD_N[2][caps][node], caps, CellStatus[caps], NODE_VEL_N[0][caps][node], NODE_VEL_N[1][caps][node], NODE_VEL_N[2][caps][node],F_Total_N[0][caps][node],F_Total_N[1][caps][node],F_Total_N[2][caps][node],NUC_Neigh[caps][node],NUC_NeighFlag[caps][node]);
						}
					}
				}
			fclose(fOUT);
			
			strcpy(name, "SimulationOutputData/AllNucEdgeMesh");
			strcat(name,"-");
			strcat(name,num);
			strcat(name,".csv");
			printf("Ouputting fluid at iter#%d to %s...\n", iter, name);
			fOUT = fopen(name, "w");
				fprintf(fOUT, "X,Y,Z,Cell,LX,LY,LZ,Stretch\n");
				for(caps=1;caps<=CapsM;caps++){
					if(CellStatus[caps] > -5000){
						for(edge=1;edge<=EdgeM_N;edge++){
							n1 = EDGENODE_LIST_N[edge][1];
							n2 = EDGENODE_LIST_N[edge][2];
															
							n1coord[0] = NODECOORD_N[0][caps][n1];
							n1coord[1] = NODECOORD_N[1][caps][n1];
							n1coord[2] = NODECOORD_N[2][caps][n1];
							
							n2coord[0] = NODECOORD_N[0][caps][n2];
							n2coord[1] = NODECOORD_N[1][caps][n2];
							n2coord[2] = NODECOORD_N[2][caps][n2];
							
							n3coord[0] = 0.5*(n1coord[0]+n2coord[0]);
							n3coord[1] = 0.5*(n1coord[1]+n2coord[1]);
							n3coord[2] = 0.5*(n1coord[2]+n2coord[2]);
							
							e1v[0] = n1coord[0] - n2coord[0];
							e1v[1] = n1coord[1] - n2coord[1];
							e1v[2] = n1coord[2] - n2coord[2];
							
							edge1_length = sqrtf(e1v[0]*e1v[0]+e1v[1]*e1v[1]+e1v[2]*e1v[2]);
							
							e_ratio = edge1_length*1.e-06/LMAX_N[edge];
							
							fprintf(fOUT, "%.7e,%.7e,%.7e,%d,%.7e,%.7e,%.7e,%.7e\n", n3coord[0],n3coord[1],n3coord[2],caps,e1v[0],e1v[1],e1v[2],e_ratio);
						}
					}
				}
			fclose(fOUT);
			
			printf("Finished ouputting fluid at iter#%d to %s...\n", iter, name);
}





void LBM_PARENTDOMAIN_SETUP(void){
	char number[15], num[20], name[250], bulk[200], junkchar;
	float junkfloat;
	int l, bval, i, j, k, m, n;
	float dval, uval, vval, wval, xval, yval, zval;
	
	FILE *fOUT, *fIN;
	
	
	for(i=1;i<=IM_1;i++){
		for(j=1;j<=JM_1;j++){
			for(k=1;k<=KM_1;k++){
				B_index_1[i][j][k] = -1;
				X_1[i][j][k] = X0_1 + (i-1)*DX;
				Y_1[i][j][k] = Y0_1 + (j-1)*DX;
				Z_1[i][j][k] = Z0_1 + (k-1)*DX;
			}
		}
	}
	
	for(i=1;i<=IM_2;i++){
		for(j=1;j<=JM_2;j++){
			for(k=1;k<=KM_2;k++){
				B_index_2[i][j][k] = -1;
				X_2[i][j][k] = X0_2 + (i-1)*DX;
				Y_2[i][j][k] = Y0_2 + (j-1)*DX;
				Z_2[i][j][k] = Z0_2 + (k-1)*DX;
			}
		}
	}
	
	for(i=1;i<=IM_3;i++){
		for(j=1;j<=JM_3;j++){
			for(k=1;k<=KM_3;k++){
				B_index_3[i][j][k] = -1;
				X_3[i][j][k] = X0_3 + (i-1)*DX;
				Y_3[i][j][k] = Y0_3 + (j-1)*DX;
				Z_3[i][j][k] = Z0_3 + (k-1)*DX;
			}
		}
	}
	///
	strcpy(name,"InputData/DAparentBoundary.csv");
	printf("loading %s file\n",name);
	fIN = fopen(name, "r");
	fgets(bulk, 180, fIN);
		for(m=1;m<=3008;m++){
			fscanf(fIN,"%d", &bval);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &junkfloat);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &junkfloat);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &junkfloat);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &junkfloat);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &xval);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &yval);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &zval);
			
			i = (int)((xval - X0_1)/DX + 1.01);
			j = (int)((yval - Y0_1)/DX + 1.01);
			k = (int)((zval - Z0_1)/DX + 1.01);
			
			B_index_1[i][j][k] = bval;
			
			X_1[i][j][k] = xval;
			Y_1[i][j][k] = yval;
			Z_1[i][j][k] = zval;
			
			if(k == KM_1 && bval > 0) B_index_1[i][j][k] = 11;
							
		}
	fclose(fIN);
	
	for(i=1;i<=IM_1;i++){
		for(j=1;j<=JM_1;j++){
			for(k=1;k<KM_1;k++){
				B_index_1[i][j][k] = B_index_1[i][j][KM_1];
				if(B_index_1[i][j][k] > 0){
					if(k > 1) B_index_1[i][j][k] = 1;
					else B_index_1[i][j][k] = 21;
				}
			}
		}
	}
	
	strcpy(name,"InputData/PCVparent1Boundary.csv");
	printf("loading %s file\n",name);
	fIN = fopen(name, "r");
	fgets(bulk, 180, fIN);
		for(m=1;m<=2918;m++){
			fscanf(fIN,"%d", &bval);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &junkfloat);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &junkfloat);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &junkfloat);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &junkfloat);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &xval);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &yval);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &zval);
			
			i = (int)((xval - X0_2)/DX + 1.01);
			j = (int)((yval - Y0_2)/DX + 1.01);
			k = (int)((zval - Z0_2)/DX + 1.01);
			
			B_index_2[i][j][k] = bval;
			
			X_2[i][j][k] = xval;
			Y_2[i][j][k] = yval;
			Z_2[i][j][k] = zval;
			
			if(k == 1 && bval > 0) B_index_2[i][j][k] = 12;
			
			
		}
	fclose(fIN);
	
	for(i=1;i<=IM_2;i++){
		for(j=1;j<=JM_2;j++){
			for(k=2;k<=KM_2;k++){
				B_index_2[i][j][k] = B_index_2[i][j][1];
				if(B_index_2[i][j][k] > 0){
					if(k < KM_2) B_index_2[i][j][k] = 1;
					else B_index_2[i][j][k] = 22;
				}
				
			}
		}
	}
	
	strcpy(name,"InputData/PCVparent2Boundary.csv");
	printf("loading %s file\n",name);
	fIN = fopen(name, "r");
	fgets(bulk, 180, fIN);
		for(m=1;m<=794;m++){
			fscanf(fIN,"%d", &bval);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &junkfloat);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &junkfloat);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &junkfloat);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &junkfloat);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &xval);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &yval);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &zval);
			
			i = (int)((xval - X0_3)/DX + 1.01);
			j = (int)((yval - Y0_3)/DX + 1.01);
			k = (int)((zval - Z0_3)/DX + 1.01);
			
			B_index_3[i][j][k] = bval;
			
			X_3[i][j][k] = xval;
			Y_3[i][j][k] = yval;
			Z_3[i][j][k] = zval;
			
			if(k == 1 && bval > 0) B_index_3[i][j][k] = 13;
			
		}
	fclose(fIN);
	
	for(i=1;i<=IM_3;i++){
		for(j=1;j<=JM_3;j++){
			for(k=2;k<=KM_3;k++){
				B_index_3[i][j][k] = B_index_3[i][j][1];
				if(B_index_3[i][j][k] > 0){
					if(k < KM_3) B_index_3[i][j][k] = 1;
					else B_index_3[i][j][k] = 23;
				}
				
			}
		}
	}
	
	///let's make the LBM parent domain, fluid and BC lists in 1D array format for OMP multithreading
	//DA parent
	DOMMAX_1 = 0;
	for(i=1;i<=IM_1;i++){
		for(j=1;j<=JM_1;j++){
			for(k=1;k<=KM_1;k++){
				DOMID_1[i][j][k] = -1;
				if(B_index_1[i][j][k] >= 0){
					DOMMAX_1 += 1;
					DOMLIST_1[DOMMAX_1][0] = i;
					DOMLIST_1[DOMMAX_1][1] = j;
					DOMLIST_1[DOMMAX_1][2] = k;
					
					DOMID_1[i][j][k] = DOMMAX_1;
				}
			}
		}
	}
	printf("There are %d positive  DA parent domain nodes\n",DOMMAX_1);
	
	//PCV parent1
	DOMMAX_2 = 0;
	for(i=1;i<=IM_2;i++){
		for(j=1;j<=JM_2;j++){
			for(k=1;k<=KM_2;k++){
				DOMID_2[i][j][k] = -1;
				if(B_index_2[i][j][k] >= 0){
					DOMMAX_2 += 1;
					DOMLIST_2[DOMMAX_2][0] = i;
					DOMLIST_2[DOMMAX_2][1] = j;
					DOMLIST_2[DOMMAX_2][2] = k;
					
					DOMID_2[i][j][k] = DOMMAX_2;
				}
			}
		}
	}
	printf("There are %d positive PCV parent1 domain nodes\n",DOMMAX_2);
	
	//PCV parent2
	DOMMAX_3 = 0;
	for(i=1;i<=IM_3;i++){
		for(j=1;j<=JM_3;j++){
			for(k=1;k<=KM_3;k++){
				DOMID_3[i][j][k] = -1;
				if(B_index_3[i][j][k] >= 0){
					DOMMAX_3 += 1;
					DOMLIST_3[DOMMAX_3][0] = i;
					DOMLIST_3[DOMMAX_3][1] = j;
					DOMLIST_3[DOMMAX_3][2] = k;
					
					DOMID_3[i][j][k] = DOMMAX_3;
				}
			}
		}
	}
	printf("There are %d positive PCV parent2 domain nodes\n",DOMMAX_3);
	
	ZM_1 = Z0_1 + (KM_1-1)*DX;
	ZM_2 = Z0_2 + (KM_2-1)*DX;
	ZM_3 = Z0_3 + (KM_3-1)*DX;
	
	printf("ZM1=%f, ZM2=%f and ZM3=%f\n",ZM_1, ZM_2, ZM_3);
	///
	
	strcpy(name,"InputData/DAparentBoundary.csv");
	printf("loading %s file\n",name);
	fIN = fopen(name, "r");
	fgets(bulk, 180, fIN);
		for(m=1;m<=3008;m++){
			fscanf(fIN,"%d", &bval);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &dval);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &uval);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &vval);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &wval);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &xval);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &yval);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &zval);
			
			if(m == 1) printf("%d, %f,  %f,  %f,  %f,  %f,  %f,  %f\n", bval, dval, uval, vval, wval, xval, yval, zval);
			
			i = (int)((xval - X0_1)/DX + 1.01);
			j = (int)((yval - Y0_1)/DX + 1.01);
			k = (int)((zval - Z0_1)/DX + 1.01);
			
			l = DOMID_1[i][j][k];
			
			B_index_1[i][j][k] = bval;
			D_1[l] = dval;
			U_1[l] = uval;
			V_1[l] = vval;
			W_1[l] = wval;
			X_1[i][j][k] = xval;
			Y_1[i][j][k] = yval;
			Z_1[i][j][k] = zval;
			
			if(k == KM_1 && bval > 0) B_index_1[i][j][k] = 11;
							
		}
	fclose(fIN);
	
	for(i=1;i<=IM_1;i++){
		for(j=1;j<=JM_1;j++){
			for(k=1;k<KM_1;k++){
				B_index_1[i][j][k] = B_index_1[i][j][KM_1];
				if(B_index_1[i][j][k] > 0){
					if(k > 1) B_index_1[i][j][k] = 1;
					else B_index_1[i][j][k] = 21;
				}
				
				l = DOMID_1[i][j][k];
				n = DOMID_1[i][j][KM_1];
				
				D_1[l] = D_1[n];
				U_1[l] = U_1[n];
				V_1[l] = V_1[n];
				W_1[l] = W_1[n];
				
			}
		}
	}
	
	strcpy(name,"InputData/PCVparent1Boundary.csv");
	printf("loading %s file\n",name);
	fIN = fopen(name, "r");
	fgets(bulk, 180, fIN);
		for(m=1;m<=2918;m++){
			fscanf(fIN,"%d", &bval);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &dval);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &uval);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &vval);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &wval);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &xval);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &yval);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &zval);
			
			i = (int)((xval - X0_2)/DX + 1.01);
			j = (int)((yval - Y0_2)/DX + 1.01);
			k = (int)((zval - Z0_2)/DX + 1.01);
			
			B_index_2[i][j][k] = bval;
			D_2[l] = dval;
			U_2[l] = uval;
			V_2[l] = vval;
			W_2[l] = wval;
			X_2[i][j][k] = xval;
			Y_2[i][j][k] = yval;
			Z_2[i][j][k] = zval;
			
			if(k == 1 && bval > 0) B_index_2[i][j][k] = 12;
			
			if(m == 1) printf("%d, %f,  %f,  %f,  %f,  %f,  %f,  %f\n", bval, dval, uval, vval, wval, xval, yval, zval);
			
		}
	fclose(fIN);
	
	for(i=1;i<=IM_2;i++){
		for(j=1;j<=JM_2;j++){
			for(k=2;k<=KM_2;k++){
				B_index_2[i][j][k] = B_index_2[i][j][1];
				if(B_index_2[i][j][k] > 0){
					if(k < KM_2) B_index_2[i][j][k] = 1;
					else B_index_2[i][j][k] = 22;
				}
				
				l = DOMID_2[i][j][k];
				n = DOMID_2[i][j][1];
				
				D_2[l] = D_2[n];
				U_2[l] = U_2[n];
				V_2[l] = V_2[n];
				W_2[l] = W_2[n];
				
			}
		}
	}
	
	strcpy(name,"InputData/PCVparent2Boundary.csv");
	printf("loading %s file\n",name);
	fIN = fopen(name, "r");
		fgets(bulk, 180, fIN);
		for(m=1;m<=794;m++){
			fscanf(fIN,"%d", &bval);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &dval);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &uval);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &vval);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &wval);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &xval);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &yval);
			fscanf(fIN,"%c", &junkchar);
			fscanf(fIN,"%f", &zval);
			
			i = (int)((xval - X0_3)/DX + 1.01);
			j = (int)((yval - Y0_3)/DX + 1.01);
			k = (int)((zval - Z0_3)/DX + 1.01);
			
			B_index_3[i][j][k] = bval;
			D_3[l] = dval;
			U_3[l] = uval;
			V_3[l] = vval;
			W_3[l] = wval;
			X_3[i][j][k] = xval;
			Y_3[i][j][k] = yval;
			Z_3[i][j][k] = zval;
			
			if(k == 1 && bval > 0) B_index_3[i][j][k] = 13;
			
			if(m == 1) printf("%d, %f,  %f,  %f,  %f,  %f,  %f,  %f\n", bval, dval, uval, vval, wval, xval, yval, zval);
			
		}
	fclose(fIN);
	
	for(i=1;i<=IM_3;i++){
		for(j=1;j<=JM_3;j++){
			for(k=2;k<=KM_3;k++){
				B_index_3[i][j][k] = B_index_3[i][j][1];
				if(B_index_3[i][j][k] > 0){
					if(k < KM_3) B_index_3[i][j][k] = 1;
					else B_index_3[i][j][k] = 23;
				}
				
				l = DOMID_3[i][j][k];
				n = DOMID_3[i][j][1];
				
				D_3[l] = D_3[n];
				U_3[l] = U_3[n];
				V_3[l] = V_3[n];
				W_3[l] = W_3[n];
				
			}
		}
	}
	
	///let's check that domains were properly generated
	strcpy(name,"LBMGridDataStructure/DAparentDomain.csv");
	printf("writing %s file\n",name);
	fOUT = fopen(name, "w");
		fprintf(fOUT,"x y z bindex\n");
			for(i=1;i<=IM_1;i++){
				for(j=1;j<=JM_1;j++){
					for(k=1;k<=KM_1;k++){
						fprintf(fOUT,"%.7e ", X_1[i][j][k]);
						fprintf(fOUT,"%.7e ", Y_1[i][j][k]);
						fprintf(fOUT,"%.7e ", Z_1[i][j][k]);
						fprintf(fOUT,"%d\n", B_index_1[i][j][k]);
						
					}
				}
			}
	fclose(fOUT);
	
	
	strcpy(name,"LBMGridDataStructure/PCVparent1Domain.csv");
	printf("writing %s file\n",name);
	fOUT = fopen(name, "w");
		fprintf(fOUT,"x y z bindex\n");
			for(i=1;i<=IM_2;i++){
				for(j=1;j<=JM_2;j++){
					for(k=1;k<=KM_2;k++){
						fprintf(fOUT,"%.7e ", X_2[i][j][k]);
						fprintf(fOUT,"%.7e ", Y_2[i][j][k]);
						fprintf(fOUT,"%.7e ", Z_2[i][j][k]);
						fprintf(fOUT,"%d\n", B_index_2[i][j][k]);
						
					}
				}
			}
	fclose(fOUT);
	
	
	strcpy(name,"LBMGridDataStructure/PCVparent2Domain.csv");
	printf("writing %s file\n",name);
	fOUT = fopen(name, "w");
		fprintf(fOUT,"x y z bindex\n");
			for(i=1;i<=IM_3;i++){
				for(j=1;j<=JM_3;j++){
					for(k=1;k<=KM_3;k++){
						fprintf(fOUT,"%.7e ", X_3[i][j][k]);
						fprintf(fOUT,"%.7e ", Y_3[i][j][k]);
						fprintf(fOUT,"%.7e ", Z_3[i][j][k]);
						fprintf(fOUT,"%d\n", B_index_3[i][j][k]);
						
					}
				}
			}
	fclose(fOUT);
	
	
	
}

void LBM_PARENTBC_SETUP(void){
	int i, j, k, ii, jj, kk, c, iii, jjj, kkk;
	int check;
	//float UU, VV, WW, DD;
	float U_bc, V_bc, W_bc, D_bc, feq_bc;
	float pressure, simtime, tauLBM;
	//int listcounter;
	
	simtime = iter*T_SCALE;
	tauLBM=TauPlasma;
	
	BCMAX_1 = 0;
	for(i=1;i<60000;i++){
		BCNLIST_1[i][0] = 0;
		BCNLIST_1[i][1] = 0;
		BCNLIST_1[i][2] = 0;
		BCLIST_1[i][0] = 0;
		BCLIST_1[i][1] = 0;
		BCLIST_1[i][2] = 0;
		BCLIST_1[i][3] = 0;
		vBCNLIST_1[i][0] = 0;
		vBCNLIST_1[i][1] = 0;
		vBCNLIST_1[i][2] = 0;
	}
	
	vBCMAX_1 = 0;
	
	for(i=1;i<=IM_1;i++){
		for(j=1;j<=JM_1;j++){
			for(k=1;k<=KM_1;k++){
				if(B_index_1[i][j][k] >= 10 || B_index_1[i][j][k] == 0){
					BCMAX_1 += 1;
					BCLIST_1[BCMAX_1][3] = B_index_1[i][j][k];
					BCLIST_1[BCMAX_1][0] = i;
					BCLIST_1[BCMAX_1][1] = j;
					BCLIST_1[BCMAX_1][2] = k;
						
					
					if(B_index_1[i][j][k] == 0 || B_index_1[i][j][k] == 21){
						
						check = 0;
						for(c=1;c<CM;c++){//try to find fluid neighbor
							//if(c == 19) continue;
							ii = i+Xlattice[0][c]; jj = j+Xlattice[1][c]; kk = k+Xlattice[2][c];
							
							if(ii > IM_1 || ii < 1 || jj > JM_1 || jj < 1 || kk > KM_1 || kk < 1) continue;
							
							if(B_index_1[ii][jj][kk] > 0 && B_index_1[ii][jj][kk] < 10){//adjacent neighbor is a fluid domain interior lattice point
								iii = ii; jjj = jj; kkk = kk;
								check = 1;
								break;
							}
						}
						
						if(check == 0){
							for(c=1;c<19;c++){//try to find inlet/outlet neighbor since interior fluid neighbor does not exist
								//if(c == 19) continue;
								ii = i+Xlattice[0][c]; jj = j+Xlattice[1][c]; kk = k+Xlattice[2][c];
								
								if(ii > IM_1 || ii < 1 || jj > JM_1 || jj < 1 || kk > KM_1 || kk < 1) continue;
								
								if(B_index_1[ii][jj][kk] > 10){//adjacent neighbor is a fluid domain B.C point
									iii = ii; jjj = jj; kkk = kk;
									check = 1;
									break;
								}
							}
						}
						
						if(check == 0){
							for(c=1;c<19;c++){//try to find wall neighbor since neither interior fluid neighbor nor inlet/outlet neighbor exists
								//if(c == 19) continue;
								ii = i+Xlattice[0][c]; jj = j+Xlattice[1][c]; kk = k+Xlattice[2][c];
								
								if(ii > IM_1 || ii < 1 || jj > JM_1 || jj < 1 || kk > KM_1 || kk < 1) continue;
								
								if(B_index_1[ii][jj][kk] == 0){//adjacent neighbor is a fluid domain B.C point
									iii = ii; jjj = jj; kkk = kk;
									check = 1;
									break;
								}
							}
						}
						
						if(check == 0){
							printf("Boundary condition error on node(%d,%d,%d) for iteration#%d... cannot find neighboring interior or wall node\n",i,j,k,iter);
						}
					}
					else if(B_index_1[i][j][k] == 11){
					
						iii = (int)((X_1[i][j][KM_1] - X0)/DX + 1.01);///this the i index of the ROI domain overlapping lattice node
						jjj = (int)((Y_1[i][j][KM_1] - Y0)/DX + 1.01);///this the j index of the ROI domain overlapping lattice node
						kkk = (int)((Z_1[i][j][KM_1] - Z0)/DX + 1.01);///this the k index of the ROI domain overlapping lattice node
					//printf("Boundary copy BC for (%d,%d,%d) on DA parent domain is referenced to (%d,%d,%d) on ROI domain\n", i,j,k,iii,jjj,kkk);	
						if(kkk < 1 || kkk > KM) printf("Boundary copy BC for (%d,%d,%d) on DA parent domain has error!\n", i,j,k);
					}
					
					if(B_index_1[i][j][k] == 21){
						vBCMAX_1 += 1;
						ii = (int)((X_1[i][j][KM_1] - X0)/DX + 1.01);
						jj = (int)((Y_1[i][j][KM_1] - Y0)/DX + 1.01);
						kk = (int)((Z_1[i][j][KM_1] - Z0)/DX + 1.01);
						
						vBCNLIST_1[BCMAX_1][0] = ii;
						vBCNLIST_1[BCMAX_1][1] = jj;
						vBCNLIST_1[BCMAX_1][2] = kk;
						
					}
						
					BCNLIST_1[BCMAX_1][0] = iii;
					BCNLIST_1[BCMAX_1][1] = jjj;
					BCNLIST_1[BCMAX_1][2] = kkk;
					
					
				}
			}
		}
	}
	
	printf("There are %d boundary nodes in DA parent\n",BCMAX_1);
	printf("There are %d special boundary nodes in DA parent\n",vBCMAX_1);
	
	BCMAX_2 = 0;
	for(i=1;i<60000;i++){
		BCNLIST_2[i][0] = 0;
		BCNLIST_2[i][1] = 0;
		BCNLIST_2[i][2] = 0;
		BCLIST_2[i][0] = 0;
		BCLIST_2[i][1] = 0;
		BCLIST_2[i][2] = 0;
		BCLIST_2[i][3] = 0;
		vBCNLIST_2[i][0] = 0;
		vBCNLIST_2[i][1] = 0;
		vBCNLIST_2[i][2] = 0;
	}
	
	vBCMAX_2 = 0;
	for(i=1;i<=IM_2;i++){
		for(j=1;j<=JM_2;j++){
			for(k=1;k<=KM_2;k++){
				if(B_index_2[i][j][k] >= 10 || B_index_2[i][j][k] == 0){
					BCMAX_2 += 1;
					BCLIST_2[BCMAX_2][3] = B_index_2[i][j][k];
					BCLIST_2[BCMAX_2][0] = i;
					BCLIST_2[BCMAX_2][1] = j;
					BCLIST_2[BCMAX_2][2] = k;
					
					if(B_index_2[i][j][k] == 0 || B_index_2[i][j][k] == 22){
						
						
						check = 0;
						for(c=1;c<CM;c++){//try to find fluid neighbor
							//if(c == 19) continue;
							ii = i+Xlattice[0][c]; jj = j+Xlattice[1][c]; kk = k+Xlattice[2][c];
							//if(ii > IM_2 || ii < 1 || jj > JM_2 || jj < 1 || kk > KM_2 || kk < 1) continue;
							if(B_index_2[ii][jj][kk] > 0 && B_index_2[ii][jj][kk] < 10){//adjacent neighbor is a fluid domain interior lattice point
								iii = ii; jjj = jj; kkk = kk;
								check = 1;
								break;
							}
						}
						
						if(check == 0){
							for(c=1;c<19;c++){//try to find inlet/outlet neighbor since interior fluid neighbor does not exist
								//if(c == 19) continue;
								ii = i+Xlattice[0][c]; jj = j+Xlattice[1][c]; kk = k+Xlattice[2][c];
								if(ii > IM_2 || ii < 1 || jj > JM_2 || jj < 1 || kk > KM_2 || kk < 1) continue;
								if(B_index_2[ii][jj][kk] > 10){//adjacent neighbor is a fluid domain B.C point
									iii = ii; jjj = jj; kkk = kk;
									check = 1;
									break;
								}
							}
						}
						
						if(check == 0){
							for(c=1;c<19;c++){//try to find wall neighbor since neither interior fluid neighbor nor inlet/outlet neighbor exists
								//if(c == 19) continue;
								ii = i+Xlattice[0][c]; jj = j+Xlattice[1][c]; kk = k+Xlattice[2][c];
								if(ii > IM_2 || ii < 1 || jj > JM_2 || jj < 1 || kk > KM_2 || kk < 1) continue;
								if(B_index_2[ii][jj][kk] == 0){//adjacent neighbor is a fluid domain B.C point
									iii = ii; jjj = jj; kkk = kk;
									check = 1;
									break;
								}
							}
						}
						
						if(check == 0){
							printf("Boundary condition error on node(%d,%d,%d) for iteration#%d... cannot find neighboring interior or wall node\n",i,j,k,iter);
						}
					}
					else if(B_index_2[i][j][k] == 12){
						
						
						iii = (int)((X_2[i][j][1] - X0)/DX + 1.01);
						jjj = (int)((Y_2[i][j][1] - Y0)/DX + 1.01);
						kkk = (int)((Z_2[i][j][1] - Z0)/DX + 1.01);
						
						if(kkk < 1 || kkk > KM) printf("Boundary copy BC for (%d,%d,%d) on PCV parent1 domain has error!\n", i,j,k);
					}
					
					if(B_index_2[i][j][k] == 22){
						vBCMAX_2 += 1;
						ii = (int)((X_2[i][j][1] - X0)/DX + 1.01);
						jj = (int)((Y_2[i][j][1] - Y0)/DX + 1.01);
						kk = (int)((Z_2[i][j][1] - Z0)/DX + 1.01);
						
						
						vBCNLIST_2[BCMAX_2][0] = ii;
						vBCNLIST_2[BCMAX_2][1] = jj;
						vBCNLIST_2[BCMAX_2][2] = kk;
					}
					
					
					BCNLIST_2[BCMAX_2][0] = iii;
					BCNLIST_2[BCMAX_2][1] = jjj;
					BCNLIST_2[BCMAX_2][2] = kkk;					
					
				}
			}
		}
	}
	
	printf("There are %d boundary nodes in PCV parent1\n",BCMAX_2);
	printf("There are %d special boundary nodes in PCV parent1\n",vBCMAX_2);
	
	BCMAX_3 = 0;
	for(i=1;i<25000;i++){
		BCNLIST_3[i][0] = 0;
		BCNLIST_3[i][1] = 0;
		BCNLIST_3[i][2] = 0;
		BCLIST_3[i][0] = 0;
		BCLIST_3[i][1] = 0;
		BCLIST_3[i][2] = 0;
		BCLIST_3[i][3] = 0;
		vBCNLIST_3[i][0] = 0;
		vBCNLIST_3[i][1] = 0;
		vBCNLIST_3[i][2] = 0;
	
	}
	
	vBCMAX_3 = 0;
	
	for(i=1;i<=IM_3;i++){
		for(j=1;j<=JM_3;j++){
			for(k=1;k<=KM_3;k++){
				if(B_index_3[i][j][k] >= 10 || B_index_3[i][j][k] == 0){
					BCMAX_3 += 1;
					BCLIST_3[BCMAX_3][3] = B_index_3[i][j][k];
					BCLIST_3[BCMAX_3][0] = i;
					BCLIST_3[BCMAX_3][1] = j;
					BCLIST_3[BCMAX_3][2] = k;
						
					if(B_index_3[i][j][k] == 23 || B_index_3[i][j][k] == 0){
						
						
						check = 0;
						for(c=1;c<CM;c++){//try to find fluid neighbor
							//if(c == 19) continue;
							ii = i+Xlattice[0][c]; jj = j+Xlattice[1][c]; kk = k+Xlattice[2][c];
							if(ii > IM_3 || ii < 1 || jj > JM_3 || jj < 1 || kk > KM_3 || kk < 1) continue;
							if(B_index_3[ii][jj][kk] > 0 && B_index_3[ii][jj][kk] < 10){//adjacent neighbor is a fluid domain interior lattice point
								iii = ii; jjj = jj; kkk = kk;
								check = 1;
								break;
							}
						}
						
						if(check == 0){
							for(c=1;c<19;c++){//try to find inlet/outlet neighbor since interior fluid neighbor does not exist
								//if(c == 19) continue;
								ii = i+Xlattice[0][c]; jj = j+Xlattice[1][c]; kk = k+Xlattice[2][c];
								if(ii > IM_3 || ii < 1 || jj > JM_3 || jj < 1 || kk > KM_3 || kk < 1) continue;
								if(B_index_3[ii][jj][kk] > 10){//adjacent neighbor is a fluid domain B.C point
									iii = ii; jjj = jj; kkk = kk;
									check = 1;
									break;
								}
							}
						}
						
						if(check == 0){
							for(c=1;c<19;c++){//try to find wall neighbor since neither interior fluid neighbor nor inlet/outlet neighbor exists
								//if(c == 19) continue;
								ii = i+Xlattice[0][c]; jj = j+Xlattice[1][c]; kk = k+Xlattice[2][c];
								if(ii > IM_3 || ii < 1 || jj > JM_3 || jj < 1 || kk > KM_3 || kk < 1) continue;
								if(B_index_3[ii][jj][kk] == 0){//adjacent neighbor is a fluid domain B.C point
									iii = ii; jjj = jj; kkk = kk;
									check = 1;
									break;
								}
							}
						}
						
						if(check == 0){
							printf("Boundary condition error on node(%d,%d,%d) for iteration#%d... cannot find neighboring interior or wall node\n",i,j,k,iter);
						}
					}
					else if(B_index_3[i][j][k] == 13){
						iii = (int)((X_3[i][j][1] - X0)/DX + 1.01);
						jjj = (int)((Y_3[i][j][1] - Y0)/DX + 1.01);
						kkk = (int)((Z_3[i][j][1] - Z0)/DX + 1.01);
						
						if(kkk < 1 || kkk > KM) printf("Boundary copy BC for (%d,%d,%d) on PCV parent2 domain has error!\n", i,j,k);
					}
					
					if(B_index_3[i][j][k] == 23){
						vBCMAX_3 += 1;
						ii = (int)((X_3[i][j][1] - X0)/DX + 1.01);
						jj = (int)((Y_3[i][j][1] - Y0)/DX + 1.01);
						kk = (int)((Z_3[i][j][1] - Z0)/DX + 1.01);
										
						vBCNLIST_3[BCMAX_3][0] = ii;
						vBCNLIST_3[BCMAX_3][1] = jj;
						vBCNLIST_3[BCMAX_3][2] = kk;
					}
					
					
					BCNLIST_3[BCMAX_3][0] = iii;
					BCNLIST_3[BCMAX_3][1] = jjj;
					BCNLIST_3[BCMAX_3][2] = kkk;	
						
				}
			}
		}
	}
	
	printf("There are %d boundary nodes in PCV parent2\n",BCMAX_3);
	printf("There are %d special boundary nodes in PCV parent2\n",vBCMAX_3);
}


void LBM_PARENTINTERPOLATESOLUTION(void){
	char number[15], num[20], name[250], bulk[200], junkchar;
	int m, l, i, j, k, imain, jmain, kmain, kc;
	float UU, VV, WW, DD, xp, yp, zp;
	FILE *fOUT, *fIN;
	
	
	for(l=1;l<=DOMMAX_1;l++){
		i = DOMLIST_1[l][0];
		j = DOMLIST_1[l][1];
		k = DOMLIST_1[l][2];
		
		imain = (int)((X_1[i][j][k] - X0)/DX + 1.01);
		jmain = (int)((Y_1[i][j][k] - Y0)/DX + 1.01);
		kmain = (int)((Z_1[i][j][k] - Z0)/DX + 1.01);
		kc = (int)((ZM_1 - Z0)/DX + 1.01);
		
		m = DOMID[imain][jmain][kc];
		
		D_1[l] = D[m];
		U_1[l] = U[m];
		V_1[l] = V[m];
		W_1[l] = W[m];
	}
	
	for(l=1;l<=DOMMAX_2;l++){
		i = DOMLIST_2[l][0];
		j = DOMLIST_2[l][1];
		k = DOMLIST_2[l][2];
		
		imain = (int)((X_2[i][j][k] - X0)/DX + 1.01);
		jmain = (int)((Y_2[i][j][k] - Y0)/DX + 1.01);
		kmain = (int)((Z_2[i][j][k] - Z0)/DX + 1.01);
		kc = (int)((Z0_2 - Z0)/DX + 1.01);
		
		m = DOMID[imain][jmain][kc];
		
		D_2[l] = D[m];
		U_2[l] = U[m];
		V_2[l] = V[m];
		W_2[l] = W[m];
	}
	
	for(l=1;l<=DOMMAX_3;l++){
		i = DOMLIST_3[l][0];
		j = DOMLIST_3[l][1];
		k = DOMLIST_3[l][2];
		
		imain = (int)((X_3[i][j][k] - X0)/DX + 1.01);
		jmain = (int)((Y_3[i][j][k] - Y0)/DX + 1.01);
		kmain = (int)((Z_3[i][j][k] - Z0)/DX + 1.01);
		kc = (int)((Z0_3 - Z0)/DX + 1.01);
		
		m = DOMID[imain][jmain][kc];
		
		D_3[l] = D[m];
		U_3[l] = U[m];
		V_3[l] = V[m];
		W_3[l] = W[m];
	}
	
	strcpy(name,"LBMGridDataStructure/DAparentDomain.csv");
	printf("writing %s file\n",name);
	fOUT = fopen(name, "w");
		fprintf(fOUT,"x y z bindex U V W D\n");
			for(i=1;i<=IM_1;i++){
				for(j=1;j<=JM_1;j++){
					for(k=1;k<=KM_1;k++){
						//if(B_index_1[i][j][k] > 0){
						l = DOMID_1[i][j][k];
							fprintf(fOUT,"%.7e ", X_1[i][j][k]);
							fprintf(fOUT,"%.7e ", Y_1[i][j][k]);
							fprintf(fOUT,"%.7e ", Z_1[i][j][k]);
							fprintf(fOUT,"%d ", B_index_1[i][j][k]);
							fprintf(fOUT,"%.7e ", U_1[l]);
							fprintf(fOUT,"%.7e ", V_1[l]);
							fprintf(fOUT,"%.7e ", W_1[l]);
							fprintf(fOUT,"%.7e\n",D_1[l]);
						//}
					}
				}
			}
	fclose(fOUT);
	
	
	strcpy(name,"LBMGridDataStructure/PCVparent1Domain.csv");
	printf("writing %s file\n",name);
	fOUT = fopen(name, "w");
		fprintf(fOUT,"x y z bindex U V W D\n");
			for(i=1;i<=IM_2;i++){
				for(j=1;j<=JM_2;j++){
					for(k=1;k<=KM_2;k++){
						//if(B_index_2[i][j][k] > 0){
						l = DOMID_2[i][j][k];
							fprintf(fOUT,"%.7e ", X_2[i][j][k]);
							fprintf(fOUT,"%.7e ", Y_2[i][j][k]);
							fprintf(fOUT,"%.7e ", Z_2[i][j][k]);
							fprintf(fOUT,"%d ", B_index_2[i][j][k]);
							fprintf(fOUT,"%.7e ", U_2[l]);
							fprintf(fOUT,"%.7e ", V_2[l]);
							fprintf(fOUT,"%.7e ", W_2[l]);
							fprintf(fOUT,"%.7e\n",D_2[l]);
						//}
					}
				}
			}
	fclose(fOUT);
	
	
	strcpy(name,"LBMGridDataStructure/PCVparent2Domain.csv");
	printf("writing %s file\n",name);
	fOUT = fopen(name, "w");
		fprintf(fOUT,"x y z bindex U V W D\n");
			for(i=1;i<=IM_3;i++){
				for(j=1;j<=JM_3;j++){
					for(k=1;k<=KM_3;k++){
						//if(B_index_3[i][j][k] > 0){
						l = DOMID_3[i][j][k];
							fprintf(fOUT,"%.7e ", X_3[i][j][k]);
							fprintf(fOUT,"%.7e ", Y_3[i][j][k]);
							fprintf(fOUT,"%.7e ", Z_3[i][j][k]);
							fprintf(fOUT,"%d ", B_index_3[i][j][k]);
							fprintf(fOUT,"%.7e ", U_3[l]);
							fprintf(fOUT,"%.7e ", V_3[l]);
							fprintf(fOUT,"%.7e ", W_3[l]);
							fprintf(fOUT,"%.7e\n",D_3[l]);
						//}
					}
				}
			}
	fclose(fOUT);
}

void LBM_DOM_SETUP(void){
	int i, j, k;
	
	DOMMAX = 0;
	for(i=1;i<=IM;i++){
		for(j=1;j<=JM;j++){
			for(k=1;k<=KM;k++){
				DOMID[i][j][k] = -1;
				if(B_index[i][j][k] >= 0){
					DOMMAX += 1;
					DOMLIST[DOMMAX][0] = i;
					DOMLIST[DOMMAX][1] = j;
					DOMLIST[DOMMAX][2] = k;
					
					DOMID[i][j][k] = DOMMAX;
				}
			}
		}
	}
	printf("There are %d positive domain nodes\n",DOMMAX);
	
}



void LBM_BC_SETUP(void){
	int i, j, k, ii, jj, kk, c, iii, jjj, kkk;
	int check;
	float U_bc, V_bc, W_bc, D_bc, feq_bc;
	float pressure, simtime, tauLBM;
	
	
	simtime = iter*T_SCALE;  
	tauLBM=TauPlasma;
	
	BCMAX = 0;
	for(i=1;i<1200000;i++){
		BCNLIST[i][0] = 0;
		BCNLIST[i][1] = 0;
		BCNLIST[i][2] = 0;
		BCLIST[i][0] = 0;
		BCLIST[i][1] = 0;
		BCLIST[i][2] = 0;
		BCLIST[i][3] = 0;
	}
	
	for(i=1;i<=IM;i++){
		for(j=1;j<=JM;j++){
			for(k=1;k<=KM;k++){
				if(B_index[i][j][k] >= 1000 || B_index[i][j][k] == 0){
					BCMAX += 1;
					BCLIST[BCMAX][3] = B_index[i][j][k];
					BCLIST[BCMAX][0] = i;
					BCLIST[BCMAX][1] = j;
					BCLIST[BCMAX][2] = k;
					
					check = 0;
					for(c=1;c<CM;c++){//try to find fluid neighbor
						
						ii = i+Xlattice[0][c]; jj = j+Xlattice[1][c]; kk = k+Xlattice[2][c];
						if(ii > IM || ii < 1 || jj > JM || jj < 1 || kk > KM || kk < 1) continue;
						if(B_index[ii][jj][kk] > 0 && B_index[ii][jj][kk] < 1000){//adjacent neighbor is a fluid domain interior lattice point
							iii = ii; jjj = jj; kkk = kk;
							check = 1;
							break;
						}
					}
					
					if(check == 0){
						for(c=1;c<19;c++){//try to find inlet/outlet neighbor since interior fluid neighbor does not exist
							
							ii = i+Xlattice[0][c]; jj = j+Xlattice[1][c]; kk = k+Xlattice[2][c];
							if(ii > IM || ii < 1 || jj > JM || jj < 1 || kk > KM || kk < 1) continue;
							if(B_index[ii][jj][kk] > 1000){//adjacent neighbor is a fluid domain B.C point
								iii = ii; jjj = jj; kkk = kk;
								check = 1;
								break;
							}
						}
					}
					
					if(check == 0){
						for(c=1;c<19;c++){//try to find wall neighbor since neither interior fluid neighbor nor inlet/outlet neighbor exists
							
							ii = i+Xlattice[0][c]; jj = j+Xlattice[1][c]; kk = k+Xlattice[2][c];
							if(ii > IM || ii < 1 || jj > JM || jj < 1 || kk > KM || kk < 1) continue;
							if(B_index[ii][jj][kk] == 0){//adjacent neighbor is a fluid domain B.C point
								iii = ii; jjj = jj; kkk = kk;
								check = 1;
								break;
							}
						}
					}
					
					if(check == 0){
						printf("Boundary condition error on node(%d,%d,%d) for iteration#%d... cannot find neighboring interior or wall node\n",i,j,k,iter);
					}
					
					BCNLIST[BCMAX][0] = iii;
					BCNLIST[BCMAX][1] = jjj;
					BCNLIST[BCMAX][2] = kkk;					
						
				}
			}
		}
	}
	
	printf("There are %d boundary nodes\n",BCMAX);
}


void CGSM_ZFISHRBC_DATA_GEN(void){
	int	caps,node,edge,tri, junk,localnode,tri1, tri2, tri3, duplicateflag, edgecounter, edge1, edge2, i;
	int localnode1, localnode2, localnode3, localnode4, neighcount, ncmax;
	float elx, ely, elz, elength, enx, eny, enz;
	char bulk[100];
	
	FILE *fOUT, *fIN;
	
	/** read in node coordinates**/
	printf("loading ZFISHRBCNODES.dat file\n");
	fIN = fopen("InputData/ZFISHRBCNODES.dat", "r");
		fgets(bulk, 30, fIN);//first line states NodeM, ignore the line, put into junk char array bulk[]
		
		
			caps = 1;
			for(node=1;node<=NodeM;node++){
				fscanf(fIN,"%lf", &NODECOORD[0][caps][node]);
				fscanf(fIN,"%lf", &NODECOORD[1][caps][node]);
				fscanf(fIN,"%lf", &NODECOORD[2][caps][node]);
			}
		
	fclose(fIN);
	
	/** read in triangle element node list **/
	printf("loading ZFISHRBCTRIANGLES.dat file\n");
	fIN = fopen("InputData/ZFISHRBCTRIANGLES.dat", "r");
		fgets(bulk, 30, fIN);//first line states TriM, ignore the line, put into junk char array bulk[]
	
		
		for(tri=1;tri<=TriM;tri++){
			fscanf(fIN,"%d", &TRINODE1[tri]);
			fscanf(fIN,"%d", &TRINODE2[tri]);
			fscanf(fIN,"%d", &TRINODE3[tri]);
		}
		
		
	fclose(fIN);
	
	
	
	/** now let's generate the edge element node list! **/
	for(edge=1;edge<=EdgeM;edge++){
		for(localnode=1;localnode<=4;localnode++){
			EDGENODE_LIST[edge][localnode] = 0;
		}
	}
	
	
	edgecounter = 1;
	caps = 1;
	for(tri=1;tri<=TriM;tri++){
		tri1=TRINODE1[tri];
		tri2=TRINODE2[tri];
		tri3=TRINODE3[tri];
				
		//1) Check edegenode list if trinode1 and trinode2 pair have already appeared
		for(edge=1;edge<=EdgeM;edge++){
			duplicateflag=0;//by default the node-pair have not appeared in the edgenode list
			if(EDGENODE_LIST[edge][1] == tri1 && EDGENODE_LIST[edge][2] == tri2){
				duplicateflag=1;
				break;
			}
			if(EDGENODE_LIST[edge][1] == tri2 && EDGENODE_LIST[edge][2] == tri1){
				duplicateflag=1;
				break;
			}
		}
		//2) If the pair are new, then assign trinode1 & trinode2 as local nodes 1 and 2 of edge: EDGENODELIST[caps][edgecounter][1] = TRINODE1[caps][tri] ...etc.
		if(duplicateflag==0){
			EDGENODE_LIST[edgecounter][1] = tri1;
			EDGENODE_LIST[edgecounter][2] = tri2;
			
		//3) Increment edge counter by 1 if instruction 2) was done
			edgecounter += 1;
		}
		
		//4) Check edgenode list if trinode2 and trinode3 pair have already appeared
		for(edge=1;edge<=EdgeM;edge++){
			duplicateflag=0;//by default the node-pair have not appeared in the edgenode list
			if(EDGENODE_LIST[edge][1] == tri2 && EDGENODE_LIST[edge][2] == tri3){
				duplicateflag=1;
				break;
			}
			if(EDGENODE_LIST[edge][1] == tri3 && EDGENODE_LIST[edge][2] == tri2){
				duplicateflag=1;
				break;
			}
		}
		//5) If the pair are new, then assign trinode2 & trinode3 as local nodes 1 and 2 of edge: EDGENODELIST[caps][edgecounter][1] = TRINODE2[caps][tri] ...etc.
		if(duplicateflag==0){
			EDGENODE_LIST[edgecounter][1] = tri2;
			EDGENODE_LIST[edgecounter][2] = tri3;
			
		//6) Increment edge counter by 1 if instruction 5) was done
			edgecounter += 1;
		}
		
		//7) Check edgenode list if trinode3 and trinode1 pair have already appeared
		for(edge=1;edge<=EdgeM;edge++){
			duplicateflag=0;//by default the node-pair have not appeared in the edgenode list
			if(EDGENODE_LIST[edge][1] == tri3 && EDGENODE_LIST[edge][2] == tri1){
				duplicateflag=1;
				break;
			}
			if(EDGENODE_LIST[edge][1] == tri1 && EDGENODE_LIST[edge][2] == tri3){
				duplicateflag=1;
				break;
			}
		}
		//8) If the pair are new, then assign trinode3 & trinode1 as local nodes 1 and 2 of edge: EDGENODELIST[caps][edgecounter][1] = TRINODE3[caps][tri] ...etc.
		if(duplicateflag==0){
			EDGENODE_LIST[edgecounter][1] = tri3;
			EDGENODE_LIST[edgecounter][2] = tri1;
			
		//9) Increment edge counter by 1 if instruction 8) was done
			edgecounter += 1;
		}
	}
	
	
	printf("number of edges = %d\n",edgecounter-1);
	
	/** now let's populate local nodes 3 and 4 of the edge member, local nodes 3 and 4 will make the vertices of the conjugate edge in the bi-triangle rhombus element **/
	//take the edge loop
	
	for(edge=1;edge<=EdgeM;edge++){
		if(EDGENODE_LIST[edge][1] == 0) break;//don't consider null members at the end of the edge list; null members exist because of the overpredicted list size
		
		edge1 = EDGENODE_LIST[edge][1];//global index of node 1 on edge
		edge2 = EDGENODE_LIST[edge][2];//global index of node 2 on edge
	for(tri=1;tri<=TriM;tri++){
			tri1=TRINODE1[tri];
			tri2=TRINODE2[tri];
			tri3=TRINODE3[tri];
		
			
			if(edge1 == tri1 && edge2 == tri2){
				if(EDGENODE_LIST[edge][3] == 0){
				EDGENODE_LIST[edge][3] = tri3;
					EDGENODE_LIST[edge][5] = tri;
			}
				else if(EDGENODE_LIST[edge][4] == 0){
				EDGENODE_LIST[edge][4] = tri3;
					EDGENODE_LIST[edge][6] = tri;
			}
				
				if(EDGENODE_LIST[edge][3] != 0 && EDGENODE_LIST[edge][4] != 0) break;
			}
			if(edge1 == tri2 && edge2 == tri1){
				if(EDGENODE_LIST[edge][3] == 0){
					EDGENODE_LIST[edge][3] = tri3;
				EDGENODE_LIST[edge][5] = tri;
					}
				else if(EDGENODE_LIST[edge][4] == 0){
					EDGENODE_LIST[edge][4] = tri3;
					EDGENODE_LIST[edge][6] = tri;
				}
				
				if(EDGENODE_LIST[edge][3] != 0 && EDGENODE_LIST[edge][4] != 0) break;
			}
			
			
			if(edge1 == tri2 && edge2 == tri3){
				if(EDGENODE_LIST[edge][3] == 0){
					EDGENODE_LIST[edge][3] = tri1;
					EDGENODE_LIST[edge][5] = tri;
				}
				else if(EDGENODE_LIST[edge][4] == 0){
					EDGENODE_LIST[edge][4] = tri1;
					EDGENODE_LIST[edge][6] = tri;
				}
				
				if(EDGENODE_LIST[edge][3] != 0 && EDGENODE_LIST[edge][4] != 0) break;
			}
			if(edge1 == tri3 && edge2 == tri2){
				if(EDGENODE_LIST[edge][3] == 0){
					EDGENODE_LIST[edge][3] = tri1;
					EDGENODE_LIST[edge][5] = tri;
				}
				else if(EDGENODE_LIST[edge][4] == 0){
					EDGENODE_LIST[edge][4] = tri1;
					EDGENODE_LIST[edge][6] = tri;
				}
				
				if(EDGENODE_LIST[edge][3] != 0 && EDGENODE_LIST[edge][4] != 0) break;
			}
			
			
			if(edge1 == tri3 && edge2 == tri1){
				if(EDGENODE_LIST[edge][3] == 0){
					EDGENODE_LIST[edge][3] = tri2;
					EDGENODE_LIST[edge][5] = tri;
				}
				else if(EDGENODE_LIST[edge][4] == 0){
					EDGENODE_LIST[edge][4] = tri2;
					EDGENODE_LIST[edge][6] = tri;
				}
				
				if(EDGENODE_LIST[edge][3] != 0 && EDGENODE_LIST[edge][4] != 0) break;
			}
			if(edge1 == tri1 && edge2 == tri3){
				if(EDGENODE_LIST[edge][3] == 0){
					EDGENODE_LIST[edge][3] = tri2;
					EDGENODE_LIST[edge][5] = tri;
				}
				else if(EDGENODE_LIST[edge][4] == 0){
					EDGENODE_LIST[edge][4] = tri2;
					EDGENODE_LIST[edge][6] = tri;
				}
				
				if(EDGENODE_LIST[edge][3] != 0 && EDGENODE_LIST[edge][4] != 0) break;	
			}
		}
	}
	
	
	
	/** output the new edge data structure **/
	fOUT = fopen("CGSMDataStructure/RBC_edgeNODES_TRIANGLES_list.dat", "w");
		fprintf(fOUT,"edge node1 node2 node3 node4 tri1 tri2 edgemidptX edgemidptY edgemidptZ edgelengthX edgelengthY edgelengthZ\n");
		
		caps = 1;
		for(edge=1;edge<=EdgeM;edge++){
			localnode1 = EDGENODE_LIST[edge][1];
			localnode2 = EDGENODE_LIST[edge][2];
			localnode3 = EDGENODE_LIST[edge][3];
			localnode4 = EDGENODE_LIST[edge][4];
			tri1 = EDGENODE_LIST[edge][5];
			tri2 = EDGENODE_LIST[edge][6];
			elx = 0.5*(NODECOORD[0][caps][localnode1]+NODECOORD[0][caps][localnode2]);
			ely = 0.5*(NODECOORD[1][caps][localnode1]+NODECOORD[1][caps][localnode2]);
			elz = 0.5*(NODECOORD[2][caps][localnode1]+NODECOORD[2][caps][localnode2]);
			elength = sqrtf((NODECOORD[0][caps][localnode1]-NODECOORD[0][caps][localnode2])*(NODECOORD[0][caps][localnode1]-NODECOORD[0][caps][localnode2]) + (NODECOORD[1][caps][localnode1]-NODECOORD[1][caps][localnode2])*(NODECOORD[1][caps][localnode1]-NODECOORD[1][caps][localnode2]) + (NODECOORD[2][caps][localnode1]-NODECOORD[2][caps][localnode2])*(NODECOORD[2][caps][localnode1]-NODECOORD[2][caps][localnode2]));
			enx = (NODECOORD[0][caps][localnode1]-NODECOORD[0][caps][localnode2]);
			eny = (NODECOORD[1][caps][localnode1]-NODECOORD[1][caps][localnode2]);
			enz = (NODECOORD[2][caps][localnode1]-NODECOORD[2][caps][localnode2]);
			fprintf(fOUT,"%d %d %d %d %d %d %d %.7e %.7e %.7e %.7e %.7e %.7e\n", edge, localnode1, localnode2, localnode3, localnode4, tri1, tri2, elx, ely, elz, enx, eny, enz);
		}
	
	fclose(fOUT);
	
	/** generate neighbourlists **/
	for(node=1;node<=NodeM;node++){
		for(i=0;i<=NListM;i++){
			NEIGHNODE_LIST[node][i] = 0;
			NEIGHTRI_LIST[node][i] = 0;
		}
	}
	
	ncmax = 0;
	for(node=1;node<=NodeM;node++){
		neighcount = 0;
		for(edge=1;edge<=EdgeM;edge++){
			localnode1 = EDGENODE_LIST[edge][1];
			localnode2 = EDGENODE_LIST[edge][2];
			if(localnode1 == node){
				neighcount += 1;
				NEIGHNODE_LIST[node][0] = neighcount;
				NEIGHNODE_LIST[node][neighcount] = localnode2;
			}
			else if(localnode2 == node){
				neighcount += 1;
				NEIGHNODE_LIST[node][0] = neighcount;
				NEIGHNODE_LIST[node][neighcount] = localnode1;
			}
			
			if(neighcount > ncmax) ncmax = neighcount;
		}
	}
	
	printf("Maximum node neighbors in the node to node neighbor list is %d\n", ncmax);
	
	fOUT = fopen("CGSMDataStructure/RBC_NODEtoNODE_list.dat", "w");
	fprintf(fOUT,"parentnode totalneighbours n1 n2 n3 n4 n5 n6 n7 n8 n9\n");
	
		for(node=1;node<=NodeM;node++){
			fprintf(fOUT,"%d %d %d %d %d %d %d %d %d %d %d\n", node, NEIGHNODE_LIST[node][0], NEIGHNODE_LIST[node][1], NEIGHNODE_LIST[node][2], NEIGHNODE_LIST[node][3], NEIGHNODE_LIST[node][4], NEIGHNODE_LIST[node][5], NEIGHNODE_LIST[node][6], NEIGHNODE_LIST[node][7], NEIGHNODE_LIST[node][8], NEIGHNODE_LIST[node][9]);
		}
	
	fclose(fOUT);
	
	for(node=1;node<=NodeM;node++){
		neighcount = 0;
		for(tri=1;tri<=TriM;tri++){
			localnode1 = TRINODE1[tri];
			localnode2 = TRINODE2[tri];
			localnode3 = TRINODE3[tri];
			if(localnode1 == node || localnode2 == node || localnode3 == node){
				neighcount += 1;
				NEIGHTRI_LIST[node][0] = neighcount;
				NEIGHTRI_LIST[node][neighcount] = tri;
			}
			
			if(neighcount > ncmax) ncmax = neighcount;
			
		}
	}
	
	printf("Maximum tri neighbors in the node to tri neighbor list is %d\n", ncmax);
	
	
	
	fOUT = fopen("CGSMDataStructure/RBC_TRIANGLEtoNODE_list.dat", "w");
	fprintf(fOUT,"parentnode totalneighbours t1 t2 t3 t4 t5 t6 t7 t8 t9\n");
	
		for(node=1;node<=NodeM;node++){
			fprintf(fOUT,"%d %d %d %d %d %d %d %d %d %d %d\n", node, NEIGHTRI_LIST[node][0], NEIGHTRI_LIST[node][1], NEIGHTRI_LIST[node][2], NEIGHTRI_LIST[node][3], NEIGHTRI_LIST[node][4], NEIGHTRI_LIST[node][5], NEIGHTRI_LIST[node][6], NEIGHTRI_LIST[node][7], NEIGHTRI_LIST[node][8], NEIGHTRI_LIST[node][9]);
		}
	
	fclose(fOUT);
	
	/** read in Nucleus edge element node list **/
	printf("loading NucleusEdges.dat file\n");
	fIN = fopen("InputData/NucleusEdges.dat", "r");
		fgets(bulk, 30, fIN);//first line states TriM, ignore the line, put into junk char array bulk[]
			for(edge=1;edge<=EdgeM_N;edge++){
				fscanf(fIN,"%d", &EDGENODE_LIST_N[edge][1]);
				fscanf(fIN,"%d", &EDGENODE_LIST_N[edge][2]);
			}
	fclose(fIN);
	
	/** read in Nucleus node coordinates **/
	printf("loading NucleusNodes.dat file\n");
	fIN = fopen("InputData/NucleusNodes.dat", "r");
		fgets(bulk, 30, fIN);//first line states TriM, ignore the line, put into junk char array bulk[]
		for(node=1;node<=NodeM_N;node++){
			fscanf(fIN,"%lf", &NODECOORD_N[0][caps][node]);
			fscanf(fIN,"%lf", &NODECOORD_N[1][caps][node]);
			fscanf(fIN,"%lf", &NODECOORD_N[2][caps][node]);
		}
	fclose(fIN);
	
	return;
}



//RBC energy relaxation for equilirbrium shape...
void CGSM_FILAMENT_PARAMETERS(void){
	int caps, edge, node1, node2, node3, node4, tri1, tri2;
	float en1[3], en2[3];//edge-node vectors
	float rhom12_length;
	float dx1_12, dy1_12, dz1_12;
	float dx2_12, dy2_12, dz2_12;
	char head1[]="WLC", And[]="-", zero5[]="00000", zero4[]="0000", zero3[]="000",zero2[]="00",zero1[]="0";
	char number[10];
	char name[250], num[13], d1[]=".csv";
	float alpha, eta, Es_local;
	float RESIDUALplength, RESIDUALplength_N, RESIDUALkp, RESIDUALkp_N, PlengthAVE, PlengthAVE_N, kpAVE, kpAVE_N;
	
	beta = 0.5*e_ratio_eq/((1.-e_ratio_eq)*(1.-e_ratio_eq)*(1.-e_ratio_eq)) - 0.25/((1.-e_ratio_eq)*(1.-e_ratio_eq)) + 0.25;
	
	caps = 1;
	RESIDUALplength = 0;
	RESIDUALkp = 0;
	PlengthAVE = 0;
	kpAVE = 0; 
	for(edge=1;edge<=EdgeM;edge++){
		
		RESIDUALplength += Plength[edge]/EdgeM;
		RESIDUALkp += kp[edge]/EdgeM;
		
		node1 = EDGENODE_LIST[edge][1];
		node2 = EDGENODE_LIST[edge][2];
		node3 = EDGENODE_LIST[edge][3];
		node4 = EDGENODE_LIST[edge][4];
		tri1 = EDGENODE_LIST[edge][5];
		tri2 = EDGENODE_LIST[edge][6];
		
		///WLC & POW
		dx1_12 = NODECOORD[0][caps][node1] - NODECOORD[0][caps][node2];
		dy1_12 = NODECOORD[1][caps][node1] - NODECOORD[1][caps][node2];
		dz1_12 = NODECOORD[2][caps][node1] - NODECOORD[2][caps][node2];
		
		dx2_12 = NODECOORD[0][caps][node2] - NODECOORD[0][caps][node1];
		dy2_12 = NODECOORD[1][caps][node2] - NODECOORD[1][caps][node1];
		dz2_12 = NODECOORD[2][caps][node2] - NODECOORD[2][caps][node1];
		
		rhom12_length = sqrtf(dx1_12*dx1_12+dy1_12*dy1_12+dz1_12*dz1_12);
		
		///this step is for setting finite csk tension at equilibrium when the zero tension length is not the native resting length of spectrin in RBC CSK network
		//rhom12_length = e_ratio_zero/e_ratio_eq*rhom12_length;//e_ratio_zero/e_ratio_eq gives the ratio of the zero tension isolated spectrin length to native CSK spectrin length
		
		///
		en1[0] = dx1_12/rhom12_length;
		en1[1] = dy1_12/rhom12_length;
		en1[2] = dz1_12/rhom12_length;
		
		en2[0] = dx2_12/rhom12_length;
		en2[1] = dy2_12/rhom12_length;
		en2[2] = dz2_12/rhom12_length;
		
		alpha = sqrtf(3.)*0.25*BOLTZ*TEMPERATURE/(rhom12_length*1.e-006);
		eta = sqrtf(3.)*0.25*(2./*<==this is m in Fedosov's thesis*/+1.)/(rhom12_length*rhom12_length*rhom12_length*1.e-018);//m=2, https://www.dam.brown.edu/dpd/lib/exe/fetch.php/wiki:phd_thesis:fedosov_thesis.pdf
		
		if(iter_local == 1) Plength[edge] = 1.e-9;
		Es_local = Es;
		kp[edge] = (Es_local - alpha*beta/Plength[edge])/eta;
		Plength[edge] = (rhom12_length*rhom12_length*1.e-012)*TEMPERATURE*BOLTZ*(e_ratio_eq - 0.25 + 1./(4.*(1.-e_ratio_eq)*(1.-e_ratio_eq)))/kp[edge];
		
		RESIDUALplength += -Plength[edge]/EdgeM;
		RESIDUALkp += -kp[edge]/EdgeM;
		
		PlengthAVE += Plength[edge]/EdgeM;
		kpAVE += kp[edge]/EdgeM; 
	}
	
	RESIDUALplength = RESIDUALplength/PlengthAVE;
	RESIDUALkp = RESIDUALkp/kpAVE;
	
	RESIDUALplength_N = 0;
	RESIDUALkp_N = 0;
	PlengthAVE_N = 0;
	kpAVE_N = 0; 
	for(edge=1;edge<=EdgeM_N;edge++){
		
		RESIDUALplength_N += Plength_N[edge]/EdgeM_N;
		RESIDUALkp_N += kp_N[edge]/EdgeM_N;
		
		node1 = EDGENODE_LIST_N[edge][1];
		node2 = EDGENODE_LIST_N[edge][2];
		
		///WLC & POW
		dx1_12 = NODECOORD_N[0][caps][node1] - NODECOORD_N[0][caps][node2];
		dy1_12 = NODECOORD_N[1][caps][node1] - NODECOORD_N[1][caps][node2];
		dz1_12 = NODECOORD_N[2][caps][node1] - NODECOORD_N[2][caps][node2];
		
		dx2_12 = NODECOORD_N[0][caps][node2] - NODECOORD_N[0][caps][node1];
		dy2_12 = NODECOORD_N[1][caps][node2] - NODECOORD_N[1][caps][node1];
		dz2_12 = NODECOORD_N[2][caps][node2] - NODECOORD_N[2][caps][node1];
		
		rhom12_length = sqrtf(dx1_12*dx1_12+dy1_12*dy1_12+dz1_12*dz1_12);
		
		
		///
		en1[0] = dx1_12/rhom12_length;
		en1[1] = dy1_12/rhom12_length;
		en1[2] = dz1_12/rhom12_length;
		
		en2[0] = dx2_12/rhom12_length;
		en2[1] = dy2_12/rhom12_length;
		en2[2] = dz2_12/rhom12_length;
		
		alpha = sqrtf(3.)*0.25*BOLTZ*TEMPERATURE/(rhom12_length*1.e-006);
		eta = sqrtf(3.)*0.25*(2.+1.)/(rhom12_length*rhom12_length*rhom12_length*1.e-018);
		
		if(iter_local == 1) Plength_N[edge] = 1.e-09;
		Es_local = Es*0.1;
		kp_N[edge] = (Es_local - alpha*beta/Plength_N[edge])/eta;
		Plength_N[edge] = (rhom12_length*rhom12_length*1.e-012)*TEMPERATURE*BOLTZ*(e_ratio_eq - 0.25 + 1./(4.*(1.-e_ratio_eq)*(1.-e_ratio_eq)))/kp_N[edge];
		
		RESIDUALplength_N += -Plength_N[edge]/EdgeM_N;
		RESIDUALkp_N += -kp_N[edge]/EdgeM_N;
		
		PlengthAVE_N += Plength_N[edge]/EdgeM_N;
		kpAVE_N += kp_N[edge]/EdgeM_N; 
	}
	RESIDUALplength_N = RESIDUALplength_N/PlengthAVE_N;
	RESIDUALkp_N = RESIDUALkp_N/kpAVE_N;
	
	return;
}


void CGSM_EQUILIBRIUM_METRICS(void){
	int caps, edge, node1, node2, node3, node4, tri1, tri2, i;
	float tri1_n[3], tri2_n[3];//triangle normal vectors
	float theta, angle_check, costheta, sintheta;//, sintheta, root, betabend, b11, b12, b22;
	float rhom12_length, rhom23_length, rhom31_length, rhom24_length,rhom41_length, semiperimeter1, semiperimeter2;//, e_ratio;
	float rhom_tri1_cen[3], rhom_tri2_cen[3], rhom_tri1_area, rhom_tri2_area;
	float sx12, sy12, sz12;
	
	float dx1_12, dy1_12, dz1_12, dx1_23, dy1_23, dz1_23, dx1_31, dy1_31, dz1_31;
	float sx1_23, sy1_23, sz1_23, sx1_31, sy1_31, sz1_31;
	float xy1_12, yz1_12, zx1_12, xy1_23, xy1_31, yz1_23, yz1_31, zx1_23, zx1_31;
	float dx2_12, dy2_12, dz2_12, dx2_23, dy2_23, dz2_23, dx2_31, dy2_31, dz2_31;
	float sx2_23, sy2_23, sz2_23, sx2_31, sy2_31, sz2_31;
	float xy2_12, yz2_12, zx2_12, xy2_23, xy2_31, yz2_23, yz2_31, zx2_23, zx2_31;
	
	float area_diff[3], centroid_diff[3];
	
	caps = 1;
	for(edge=1;edge<=EdgeM;edge++){
		
		node1 = EDGENODE_LIST[edge][1];
		node2 = EDGENODE_LIST[edge][2];
		node3 = EDGENODE_LIST[edge][3];
		node4 = EDGENODE_LIST[edge][4];
		tri1 = EDGENODE_LIST[edge][5];
		tri2 = EDGENODE_LIST[edge][6];
		
		
		///WLC & POW
		
		dx1_12 = NODECOORD[0][caps][node1] - NODECOORD[0][caps][node2];
		dy1_12 = NODECOORD[1][caps][node1] - NODECOORD[1][caps][node2];
		dz1_12 = NODECOORD[2][caps][node1] - NODECOORD[2][caps][node2];
		
		dx2_12 = NODECOORD[0][caps][node2] - NODECOORD[0][caps][node1];
		dy2_12 = NODECOORD[1][caps][node2] - NODECOORD[1][caps][node1];
		dz2_12 = NODECOORD[2][caps][node2] - NODECOORD[2][caps][node1];
		
		rhom12_length = sqrtf(dx1_12*dx1_12+dy1_12*dy1_12+dz1_12*dz1_12);
		
		
		/// 1. initializing spontaneous LMAX
		
		LMAX[edge] = 1.e-006*rhom12_length/e_ratio_eq;
			
					
		///bending
		for(i=0;i<=2;i++){
			rhom_tri1_cen[i] = (NODECOORD[i][caps][node1] + NODECOORD[i][caps][node2] + NODECOORD[i][caps][node3])/3.;//in um (microns)
			rhom_tri2_cen[i] = (NODECOORD[i][caps][node1] + NODECOORD[i][caps][node2] + NODECOORD[i][caps][node4])/3.;
			
			TRI_CENTROID[i][caps][tri1] = rhom_tri1_cen[i];//units are in um
			TRI_CENTROID[i][caps][tri2] = rhom_tri2_cen[i];
		}
		
		sx12 = NODECOORD[0][caps][node1] + NODECOORD[0][caps][node2];
		sy12 = NODECOORD[1][caps][node1] + NODECOORD[1][caps][node2];
		sz12 = NODECOORD[2][caps][node1] + NODECOORD[2][caps][node2];
		
		
		dx1_23 = NODECOORD[0][caps][node2] - NODECOORD[0][caps][node3];
		sx1_23 = NODECOORD[0][caps][node2] + NODECOORD[0][caps][node3];
		dy1_23 = NODECOORD[1][caps][node2] - NODECOORD[1][caps][node3];
		sy1_23 = NODECOORD[1][caps][node2] + NODECOORD[1][caps][node3];
		dz1_23 = NODECOORD[2][caps][node2] - NODECOORD[2][caps][node3];
		sz1_23 = NODECOORD[2][caps][node2] + NODECOORD[2][caps][node3];
		rhom23_length = sqrtf(dx1_23*dx1_23+dy1_23*dy1_23+dz1_23*dz1_23);//node 2 to centroid1
		
		dx1_31 = NODECOORD[0][caps][node3] - NODECOORD[0][caps][node1];
		sx1_31 = NODECOORD[0][caps][node3] + NODECOORD[0][caps][node1];
		dy1_31 = NODECOORD[1][caps][node3] - NODECOORD[1][caps][node1];
		sy1_31 = NODECOORD[1][caps][node3] + NODECOORD[1][caps][node1];
		dz1_31 = NODECOORD[2][caps][node3] - NODECOORD[2][caps][node1];
		sz1_31 = NODECOORD[2][caps][node3] + NODECOORD[2][caps][node1];
		rhom31_length = sqrtf(dx1_31*dx1_31+dy1_31*dy1_31+dz1_31*dz1_31);//centroid1 to node 3
		
		semiperimeter1 = 0.5*(rhom12_length + rhom23_length + rhom31_length);//units are in um
		
		
		dx2_23 = NODECOORD[0][caps][node1] - NODECOORD[0][caps][node4];
		sx2_23 = NODECOORD[0][caps][node1] + NODECOORD[0][caps][node4];
		dy2_23 = NODECOORD[1][caps][node1] - NODECOORD[1][caps][node4];
		sy2_23 = NODECOORD[1][caps][node1] + NODECOORD[1][caps][node4];
		dz2_23 = NODECOORD[2][caps][node1] - NODECOORD[2][caps][node4];
		sz2_23 = NODECOORD[2][caps][node1] + NODECOORD[2][caps][node4];
		rhom24_length = sqrtf(dx2_23*dx2_23+dy2_23*dy2_23+dz2_23*dz2_23);
		
		dx2_31 = NODECOORD[0][caps][node4] - NODECOORD[0][caps][node2];
		sx2_31 = NODECOORD[0][caps][node4] + NODECOORD[0][caps][node2];
		dy2_31 = NODECOORD[1][caps][node4] - NODECOORD[1][caps][node2];
		sy2_31 = NODECOORD[1][caps][node4] + NODECOORD[1][caps][node2];
		dz2_31 = NODECOORD[2][caps][node4] - NODECOORD[2][caps][node2];
		sz2_31 = NODECOORD[2][caps][node4] + NODECOORD[2][caps][node2];
		rhom41_length = sqrtf(dx2_31*dx2_31+dy2_31*dy2_31+dz2_31*dz2_31);
		
		semiperimeter2 = 0.5*(rhom12_length + rhom24_length + rhom41_length);//units are in um
		
		rhom_tri1_area = sqrtf(semiperimeter1*(semiperimeter1-rhom12_length)*(semiperimeter1-rhom23_length)*(semiperimeter1-rhom31_length))*1.e-12;
		rhom_tri2_area = sqrtf(semiperimeter2*(semiperimeter2-rhom12_length)*(semiperimeter2-rhom24_length)*(semiperimeter2-rhom41_length))*1.e-12;
		
		xy1_12 = dx1_12*sy12;   xy1_23 = dx1_23*sy1_23;   xy1_31 = dx1_31*sy1_31;
		yz1_12 = dy1_12*sz12;   yz1_23 = dy1_23*sz1_23;   yz1_31 = dy1_31*sz1_31;
		zx1_12 = dz1_12*sx12;   zx1_23 = dz1_23*sx1_23;   zx1_31 = dz1_31*sx1_31;
		
		xy2_12 = dx2_12*sy12;   xy2_23 = dx2_23*sy2_23;   xy2_31 = dx2_31*sy2_31;
		yz2_12 = dy2_12*sz12;   yz2_23 = dy2_23*sz2_23;   yz2_31 = dy2_31*sz2_31;
		zx2_12 = dz2_12*sx12;   zx2_23 = dz2_23*sx2_23;   zx2_31 = dz2_31*sx2_31;
		
		tri1_n[0] = 0.5*(yz1_12 + yz1_23 + yz1_31)*1.e-12/rhom_tri1_area;
		tri1_n[1] = 0.5*(zx1_12 + zx1_23 + zx1_31)*1.e-12/rhom_tri1_area;
		tri1_n[2] = 0.5*(xy1_12 + xy1_23 + xy1_31)*1.e-12/rhom_tri1_area;
		
		tri2_n[0] = 0.5*(yz2_12 + yz2_23 + yz2_31)*1.e-12/rhom_tri2_area;
		tri2_n[1] = 0.5*(zx2_12 + zx2_23 + zx2_31)*1.e-12/rhom_tri2_area;
		tri2_n[2] = 0.5*(xy2_12 + xy2_23 + xy2_31)*1.e-12/rhom_tri2_area;
		
		for(i=0;i<=2;i++){
			area_diff[i] = tri1_n[i] - tri2_n[i];
			centroid_diff[i] = (rhom_tri1_cen[i] - rhom_tri2_cen[i])*1.e-006;
		}
			
			
		costheta = tri1_n[0]*tri2_n[0] + tri1_n[1]*tri2_n[1] + tri1_n[2]*tri2_n[2];
		
		//make sure costheta does not exceed 1.0 in absolute value
		if(costheta > 1) costheta = 1.;
		else if(costheta < -1) costheta = -1.;
		
		angle_check = 0;
		for(i=0;i<=2;i++){
			angle_check += area_diff[i]*centroid_diff[i];
		}
		
		if(angle_check >= 0) sintheta = sqrtf(1.-costheta*costheta);
		else sintheta = -sqrtf(1.-costheta*costheta);
		
		theta = asinf(sintheta);
		
		
		SINETHETA0[edge] = sintheta;
		COSINETHETA0[edge] = costheta;
		Theta0[edge] = theta;
		
	}
		
	for(edge=1;edge<=EdgeM_N;edge++){
		
		node1 = EDGENODE_LIST_N[edge][1];
		node2 = EDGENODE_LIST_N[edge][2];
		
		///WLC & POW
		dx1_12 = NODECOORD_N[0][caps][node1] - NODECOORD_N[0][caps][node2];
		dy1_12 = NODECOORD_N[1][caps][node1] - NODECOORD_N[1][caps][node2];
		dz1_12 = NODECOORD_N[2][caps][node1] - NODECOORD_N[2][caps][node2];
		
		dx2_12 = NODECOORD_N[0][caps][node2] - NODECOORD_N[0][caps][node1];
		dy2_12 = NODECOORD_N[1][caps][node2] - NODECOORD_N[1][caps][node1];
		dz2_12 = NODECOORD_N[2][caps][node2] - NODECOORD_N[2][caps][node1];
		
		rhom12_length = sqrtf(dx1_12*dx1_12+dy1_12*dy1_12+dz1_12*dz1_12);
		
		
		/// 1. initializing spontaneous LMAX
		
		LMAX_N[edge] = 1.e-006*rhom12_length/e_ratio_eq;
	}
	return;
}


void OUTPUT_EQUILIBRIUM_PARAMETERS(void){
	int	tri, edge, node, caps;

	FILE *fOUT;
	char name[250], d3[]=".dat";
	
		
	strcpy(name,"CGSMDataStructure/MembEquilibriumData");
	strcat(name,d3);
	printf("writing equilibrium LMAX, Plength, kp and Theta0 into %s file\n", name);
	fOUT = fopen(name, "w");
		caps = 1;
		for(edge=1;edge<=EdgeM;edge++){
			fprintf(fOUT,"%.7e %.7e %.7e %.7e\n", LMAX[edge], Plength[edge], kp[edge], Theta0[edge]);
		}
	fclose(fOUT);
	
	
	strcpy(name,"CGSMDataStructure/MembTriAreaEq.dat");
	printf("writing TRI_AREA_eq into %s file\n", name);
	fOUT = fopen(name, "w");	
		caps = 1;
		for(tri=1;tri<=TriM;tri++){
			fprintf(fOUT,"%.7e\n", TRI_AREA_eq[tri]);
		}
	fclose(fOUT);
	printf("Outputted equilibrium membrane data successfully... \n");
	
	
	strcpy(name,"CGSMDataStructure/NucEquilibriumData");
	strcat(name,d3);
	printf("writing equilibrium LMAX, Plength, kp for nucleus filaments into %s file\n", name);
	fOUT = fopen(name, "w");
		caps = 1;
		for(edge=1;edge<=EdgeM_N;edge++){
			fprintf(fOUT,"%.7e %.7e %.7e\n", LMAX_N[edge], Plength_N[edge], kp_N[edge]);
		}
	fclose(fOUT);
	printf("Outputted nucleus equilibrium data successfully... \n");
	
	
	return;
}


int main(void){
	FILE *fOUT, *fIN;
	int v,i,j,k,c,cc,ii,jj,kk,l,m,n,o,iii,jjj,kkk, linecount, CellCount;
	float UU,VV,WW,DD, tauLBM;
	char number[15], num[20], name[250], bulk[400], junkchar;
	float xp, yp, zp, xc, yc, zc;
	int i1, i2, j1, j2, k1, k2;
	float distance;
	float mindis;
	int bvalue;
	char cin, zero8[] = "00000000", zero7[] = "0000000", zero6[] = "000000", zero5[]="00000", zero4[]="0000", zero3[]="000", zero2[]="00", zero1[]="0";
	float U_bc, V_bc, W_bc, D_bc, feq_bc, f_bc;
	float pressure, simtime, cycletime, taucycle;
	float junkfloat;
	int junkint, ii1, jj1, kk1, ii2, jj2, kk2, ncheck;
	int caps, node, tri, n1, n2, n3, n4, tri1, tri2, edge, nmax;
	float n1coord[3], n2coord[3], n3coord[3], n4coord[3], e1v[3], e2v[3], e3v[3], Ax, Ay, Az, AreaLocal, tnormal[3];
	float e_ratio, f_n1[3], f_n2[3], f_n3[3], f_n4[3], wlc, comp, dx1_12, dy1_12, dz1_12, edge1_length;// dx2_12, dy2_12, dz2_12;
	float betabend, costheta, sintheta, theta, sintheta0, costheta0;
	float area_diff[3], centroid_diff[3], tri1_n[3], tri2_n[3], angle_check, tri1_localcentroid[3], tri2_localcentroid[3], tri1_area, tri2_area;
	float dp21[3], dp24[3], dp32[3], dp41[3], dp13[3], b11, b12, b22;
	float cpdt1[3], cpdt2[3], cpdt3[3], cpdt4[3], cpdt5[3], cpdt6[3], cpdt7[3], cpdt8[3], cpdt9[3], cpdt10[3];//, cpdt11[3], cpdt12[3];
	float cpdt1v[3], cpdt4v[3], cpdt5v[3], cpdt8v[3], cpdt9v[3], cpdt12v[3], alpha1, alpha2;
	float vel_12[3], viscf1[3], viscf2[3];
	float CellVol, bfeq, weightsum, validinter[5+1][5+1][5+1], interweight[5+1][5+1][5+1], weightcount;
	int neighmax, neighnodemax, neightrimax;
	float distance1, distance2, distance3, repulf, invdis, weight1, weight2, weight3;
	int kmap, ParentID, check, nn;
	float zpmap, nzmin, nzmax;
	float zcmap, xn[3], xcen2[3], disvec[3];
	float xp1, yp1, zp1, xp2, yp2, zp2, xp3, yp3, zp3, xp4, yp4, zp4, xi, yi, zi, area1[3], area2[3], area3[3], aratio, r23[3], ri3[3], ri1[3], r21[3], r31[3];
	float xc1, yc1, zc1, xc2, yc2, zc2,dotpdt,mindis2,bx,by,bz;
	int EulNeighval, Windexval;
	float nnorm1[3], nnorm2[3];
	int nn1, nn2, nn3, region;
	int wneighmax, nodewcount, nodewall, triwall;
	float tc[16+1][3], ep[9+1][3], tp[15+1][3], rvec1[3+1][3], rvec2[3+1][3], a1mag, a2mag, a3mag;
	int nlocal[NListM+1], tlocal[NListM+1],nnn,trimax,tri_n,entri,ntri;
	float Ffluid[3],xcc1,ycc1,zcc1,xcc2,ycc2,zcc2,emindis;
	float n1coordF[3], n2coordF[3], n3coordF[3], n1coordB[3], n2coordB[3], n3coordB[3];
	int ic,jc, kc, Outneigh, Inneigh, check1, check2, check3, check4, checkn;
	float mindisIn, mindisOut,dotpdt2, VelReflect[3], ReflectProj, Vneigh[3];
	float f_a1[3], f_a2[3], f_a3[3], f_aa1[3], f_aa2[3], f_aa3[3], f_b1[3], f_b2[3], f_b3[3], f_b4[3], f_v1[3], f_v2[3], f_v3[3];
	int a1check, a2check;
	float l1,l2,l3;
	int nodenumber, cellindex, statval, nwtri, rflag, MNoutflag, MNin, MNinflag, MNout;
	float vx, vy, vz, nwrepulx, nwrepuly, nwrepulz, astrain, fx, fy, fz, fvx, fvy, fvz;
	int counter;
	float shiftx, shifty, shiftz, cellvel[3];
	
	printf("T_SCALE is %.7e, L_SCALE is %.7e, C_LBM is %.7e, CS is %.7e, EQ_A is %.7e, EQ_B is %.7e, EQ_C is %.7e\n", T_SCALE, L_SCALE, C_LBM, CS, EQ_A, EQ_B, EQ_C);
	
	printf("UnitVol is %.7e, UnitMass is %.7e\n", UNITVOL, UNITMASS);
	
	printf("MU is %.7e, MU_SCALE is %.7e, TA0 is %.7e TA02 is %.7e and the speed of sound is %.7e\n", MU, MU_SCALE, TauPlasma, TauCytosol, CS);
	
	printf("LBM scaling for the bodyforces is %.7e, pressure scaling is %.7e and the velocity scaling is %.7e\n", LBM_BFORCE_SCALE, LBM_PRESSURE_SCALE, VEL_SCALE);
	
	for(c=1;c<=CM;c++){
		if(c<=6) Wcoeff[c] = We1_3D;
		else if(c<=CM-1) Wcoeff[c] = We2_3D;
		else Wcoeff[c] = We3_3D;
	}
	
	
	CGSM_ZFISHRBC_DATA_GEN();
	
	
	iter_local = 0;
	do{
		iter_local = iter_local + 1;
				
		CGSM_FILAMENT_PARAMETERS();
		CGSM_EQUILIBRIUM_METRICS();
	}
	while(iter_local <= iterrelax);
	
	for(caps=1;caps<=CapsM;caps++){
		CellStatus[caps] = -5000;
		CAPS_REMOVE[caps] = 0;
		for(i=0;i<=2;i++){
			CAPS_MAX[i][caps] = -5000.;
			CAPS_MIN[i][caps] = 5000.;
		}
		for(node=1;node<=NodeM;node++){
			NODE_AREA[caps][node] = 0.;
			NODE_WALLNEIGH[caps][node] = 0;
			MembIn_NeighFlag[caps][node] = 0;
			MembOut_NeighFlag[caps][node] = 0;
			for(i=0;i<=2;i++){
				NODE_WALLREPUL[i][caps][node] = 0;
				NODE_NORMAL[i][caps][node] = 0.;
				F_VOLUME[i][caps][node] = 0.;
				F_AGG[i][caps][node] = 0.;
				F_Total[i][caps][node] = 0.;
				NODE_VEL[i][caps][node] = 0.;
				NODE_VEL_OLD[i][caps][node] = 0.;
				
				if(NODECOORD[i][caps][node] > CAPS_MAX[i][caps]) CAPS_MAX[i][caps] = NODECOORD[i][caps][node];
				if(NODECOORD[i][caps][node] < CAPS_MIN[i][caps]) CAPS_MIN[i][caps] = NODECOORD[i][caps][node];
			}
		}
		for(node=1;node<=NodeM_N;node++){
			NUC_NeighFlag[caps][node] = 0;
			for(i=0;i<=2;i++){
				NODE_WALLREPUL_N[i][caps][node] = 0;
				F_Total_N[i][caps][node] = 0.;
				NODE_VEL_N[i][caps][node] = 0.;
				NODE_VEL_OLD_N[i][caps][node] = 0.;
				
				if(NODECOORD_N[i][caps][node] > CAPS_MAX[i][caps]) CAPS_MAX[i][caps] = NODECOORD_N[i][caps][node];
				if(NODECOORD_N[i][caps][node] < CAPS_MIN[i][caps]) CAPS_MIN[i][caps] = NODECOORD_N[i][caps][node];
			}
		}
	}
	
	caps = 1;
	for(tri=1;tri<=TriM;tri++){
					
		n1 = TRINODE1[tri];
		n2 = TRINODE2[tri];
		n3 = TRINODE3[tri];
		
		
		for(i=0;i<=2;i++){
			n1coord[i] = NODECOORD[i][caps][n1];
			n2coord[i] = NODECOORD[i][caps][n2];
			n3coord[i] = NODECOORD[i][caps][n3];
			TRI_CENTROID[i][caps][tri] = (n1coord[i] + n2coord[i] + n3coord[i])/3.;
		}
		e1v[0] = n2coord[0] - n1coord[0];
		e1v[1] = n2coord[1] - n1coord[1];
		e1v[2] = n2coord[2] - n1coord[2];
		
		e2v[0] = n3coord[0] - n2coord[0];
		e2v[1] = n3coord[1] - n2coord[1];
		e2v[2] = n3coord[2] - n2coord[2];
		
		e3v[0] = n1coord[0] - n3coord[0];
		e3v[1] = n1coord[1] - n3coord[1];
		e3v[2] = n1coord[2] - n3coord[2];
		
		
		Ax = 0.5*(e1v[1]*e3v[2] - e1v[2]*e3v[1])*1.e-012;
		Ay = 0.5*(e1v[2]*e3v[0] - e1v[0]*e3v[2])*1.e-012;
		Az = 0.5*(e1v[0]*e3v[1] - e1v[1]*e3v[0])*1.e-012;
		AreaLocal = sqrtf(Ax*Ax+Ay*Ay+Az*Az);
		
		TRI_AREA[caps][tri] = AreaLocal;
		
		
		TRI_AREA_eq[tri] = AreaLocal;
		
					
		AreaExpanRatio[caps][tri] = AreaLocal/TRI_AREA_eq[tri];
		
		
		tnormal[0] = -Ax/AreaLocal;
		tnormal[1] = -Ay/AreaLocal;
		tnormal[2] = -Az/AreaLocal;
		
		TRI_NORMAL[0][caps][tri] = tnormal[0];
		TRI_NORMAL[1][caps][tri] = tnormal[1];
		TRI_NORMAL[2][caps][tri] = tnormal[2];
			
		NODE_AREA[caps][n1] += AreaLocal/3.;
		NODE_AREA[caps][n2] += AreaLocal/3.;
		NODE_AREA[caps][n3] += AreaLocal/3.; 
		
		NODE_NORMAL[0][caps][n1] += tnormal[0]*AreaLocal/3.;
		NODE_NORMAL[1][caps][n1] += tnormal[1]*AreaLocal/3.;
		NODE_NORMAL[2][caps][n1] += tnormal[2]*AreaLocal/3.;
		
		NODE_NORMAL[0][caps][n2] += tnormal[0]*AreaLocal/3.;
		NODE_NORMAL[1][caps][n2] += tnormal[1]*AreaLocal/3.;
		NODE_NORMAL[2][caps][n2] += tnormal[2]*AreaLocal/3.;
		
		NODE_NORMAL[0][caps][n3] += tnormal[0]*AreaLocal/3.;
		NODE_NORMAL[1][caps][n3] += tnormal[1]*AreaLocal/3.;
		NODE_NORMAL[2][caps][n3] += tnormal[2]*AreaLocal/3.;
		
	}
			
	for(node=1; node<=NodeM; node++){
		NODE_NORMAL[0][caps][node] = NODE_NORMAL[0][caps][node]/NODE_AREA[caps][node];
		NODE_NORMAL[1][caps][node] = NODE_NORMAL[1][caps][node]/NODE_AREA[caps][node];
		NODE_NORMAL[2][caps][node] = NODE_NORMAL[2][caps][node]/NODE_AREA[caps][node];
					
		NODE_AREA0[node]= NODE_AREA[caps][node];
	}
		
	OUTPUT_EQUILIBRIUM_PARAMETERS();
	
	printf("1\n");
	
	for(caps=1;caps<=CapsM;caps++){
		CellStatus[caps] = -5000;
	}
	
	//main i,j,k directions
	Xlattice[0][1]= 1; Xlattice[1][1]= 0; Xlattice[2][1]= 0;
	Xlattice[0][2]=-1; Xlattice[1][2]= 0; Xlattice[2][2]= 0;
	Xlattice[0][3]= 0; Xlattice[1][3]= 1; Xlattice[2][3]= 0;
	Xlattice[0][4]= 0; Xlattice[1][4]=-1; Xlattice[2][4]= 0;
	Xlattice[0][5]= 0; Xlattice[1][5]= 0; Xlattice[2][5]= 1;
	Xlattice[0][6]= 0; Xlattice[1][6]= 0; Xlattice[2][6]=-1;
	
	//xy diagonals
	Xlattice[0][7]= 1; Xlattice[1][7]= 1; Xlattice[2][7]= 0;
	Xlattice[0][8]=-1; Xlattice[1][8]=-1; Xlattice[2][8]= 0;
	Xlattice[0][9]=-1; Xlattice[1][9]= 1; Xlattice[2][9]= 0;
	Xlattice[0][10]= 1; Xlattice[1][10]=-1; Xlattice[2][10]= 0;
	
	//zx diagonals
	Xlattice[0][11]= 1; Xlattice[1][11]= 0; Xlattice[2][11]= 1;
	Xlattice[0][12]=-1; Xlattice[1][12]= 0; Xlattice[2][12]=-1;
	Xlattice[0][13]=-1; Xlattice[1][13]= 0; Xlattice[2][13]= 1;
	Xlattice[0][14]= 1; Xlattice[1][14]= 0; Xlattice[2][14]=-1;
	
	//yz diagonals
	Xlattice[0][15]= 0; Xlattice[1][15]= 1; Xlattice[2][15]= 1;
	Xlattice[0][16]= 0; Xlattice[1][16]=-1; Xlattice[2][16]=-1;
	Xlattice[0][17]= 0; Xlattice[1][17]=-1; Xlattice[2][17]= 1;
	Xlattice[0][18]= 0; Xlattice[1][18]= 1; Xlattice[2][18]=-1;
	
	//zerodirection
	Xlattice[0][19]= 0; Xlattice[1][19]= 0; Xlattice[2][19]= 0;
	
	printf("2\n");
	for(c=1;c<=CM;c++){
		UC[c]=(float)Xlattice[0][c]*C_LBM;
		VC[c]=(float)Xlattice[1][c]*C_LBM;
		WC[c]=(float)Xlattice[2][c]*C_LBM;
	}
	
	printf("3\n");
	//intialize coordinates accordingly and for all field values to be zero
	for(i=1;i<=IM;i++){
		for(j=1;j<=JM;j++){
			for(k=1;k<=KM;k++){
				X[i][j][k] = X0 + (i-1)*DX;// x-coordinates from logical i position
				Y[i][j][k] = Y0 + (j-1)*DX;//y-coord from logical j position
				Z[i][j][k] = Z0 + (k-1)*DX;
				
				RestartVal[i][j][k] = 0;
				B_index[i][j][k] = -1;
				
				Uvoid[i][j][k] = 0; Vvoid[i][j][k] = 0; Wvoid[i][j][k] = 0;
				
			}
		}
	}
	XM = X0 + (IM-1)*DX;
	YM = Y0 + (JM-1)*DX;
	ZM = Z0 + (KM-1)*DX;
	
	printf("4\n");
	for(i=1;i<=IM_1;i++){
		for(j=1;j<=JM_1;j++){
			for(k=1;k<=KM_1;k++){
				
				B_index_1[i][j][k] = -1;
				
				Uvoid_1[i][j][k] = 0; Vvoid_1[i][j][k] = 0; Wvoid_1[i][j][k] = 0;
			}
		}
	}
	
	for(i=1;i<=IM_2;i++){
		for(j=1;j<=JM_2;j++){
			for(k=1;k<=KM_2;k++){
				
				B_index_2[i][j][k] = -1;
				
				Uvoid_2[i][j][k] = 0; Vvoid_2[i][j][k] = 0; Wvoid_2[i][j][k] = 0;
			}
		}
	}
	
	for(i=1;i<=IM_3;i++){
		for(j=1;j<=JM_3;j++){
			for(k=1;k<=KM_3;k++){
				
				B_index_3[i][j][k] = -1;
				
				Uvoid_3[i][j][k] = 0; Vvoid_3[i][j][k] = 0; Wvoid_3[i][j][k] = 0;
			}
		}
	}
	
	printf("5\n");
	for(l=0;l<8500000;l++){
		U[l] = 0; V[l] = 0; W[l] = 0; D[l] = 1.;
		BODYFORCE[0][l] = 0;
		BODYFORCE[1][l] = 0;
		BODYFORCE[2][l] = 0;
		for(c=1;c<=CM;c++){
			FEQ[c][l] = D[l]*Wcoeff[c]*(1. + EQ_A*(U[l]*UC[c] + V[l]*VC[c] + W[l]*WC[c]) + EQ_B*(U[l]*UC[c] + V[l]*VC[c] + W[l]*WC[c])*(U[l]*UC[c] + V[l]*VC[c] + W[l]*WC[c]) + EQ_C*(U[l]*U[l] + V[l]*V[l] + W[l]*W[l]));
			F[c][l] = FEQ[c][l];
			F_OLD[c][l] = F[c][l];
		}
	}
	
	linecount = 0;
	strcpy(name, "InputData/FLUIDdomainBindex.csv");
		// Open the file 
		printf("Opening file %s", name); 
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
	
	
	strcpy(name, "InputData/FLUIDdomainBindex.csv");
	printf("Opening %s \n", name);
	fIN = fopen(name, "r");
		for(m=1;m<=linecount;m++){
			if(m==1) fgets(bulk, 180, fIN);
			else{
				fscanf(fIN,"%d",&bvalue);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&xp);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&yp);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&zp);
				
				i = (int)((xp-X0)/DX + 1.01);
				j = (int)((yp-Y0)/DX + 1.01);
				k = (int)((zp-Z0)/DX + 1.01);
				
				B_index[i][j][k] = bvalue;	
				
			}			
		}
	fclose(fIN);
	
	printf("Setup CFD domain\n");
	LBM_DOM_SETUP();
	printf("Setup BC\n");
	LBM_BC_SETUP();
	
	/// the start files for fluid domain is too big for GitHub so it was split into ten files
	linecount = 0;
	strcpy(name, "InputData/FLUIDdomainStart.aa");
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
	
	
	strcpy(name, "InputData/FLUIDdomainStart.aa");
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
				fscanf(fIN,"%f",&bx);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&by);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&bz);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%d",&EulNeighval);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%d",&Windexval);
				if(bvalue >= 0){
					i = (int)((xp-X0)/DX + 1.01);
					j = (int)((yp-Y0)/DX + 1.01);
					k = (int)((zp-Z0)/DX + 1.01);
					
					l = DOMID[i][j][k];
					
					RestartVal[i][j][k] = 1;
					
					D[l] = DD;
					U[l] = UU;
					V[l] = VV;
					W[l] = WW;
					BODYFORCE[0][l] = bx;
					BODYFORCE[1][l] = by;
					BODYFORCE[2][l] = bz;
					EulNeigh[0][l] = EulNeighval;
					
					W_index[i][j][k] = Windexval;
					
				}
			}			
		}
	fclose(fIN);
	
	linecount = 0;
	strcpy(name, "InputData/FLUIDdomainStart.ab");
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
	
	
	strcpy(name, "InputData/FLUIDdomainStart.ab");
	printf("Opening %s \n", name);
	fIN = fopen(name, "r");
		for(m=1;m<=linecount;m++){
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
				fscanf(fIN,"%f",&bx);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&by);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&bz);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%d",&EulNeighval);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%d",&Windexval);
				if(bvalue >= 0){
					i = (int)((xp-X0)/DX + 1.01);
					j = (int)((yp-Y0)/DX + 1.01);
					k = (int)((zp-Z0)/DX + 1.01);
					
					l = DOMID[i][j][k];
					
					RestartVal[i][j][k] = 1;
					
					D[l] = DD;
					U[l] = UU;
					V[l] = VV;
					W[l] = WW;
					BODYFORCE[0][l] = bx;
					BODYFORCE[1][l] = by;
					BODYFORCE[2][l] = bz;
					EulNeigh[0][l] = EulNeighval;
					
					W_index[i][j][k] = Windexval;
					
				}
						
		}
	fclose(fIN);
	
	linecount = 0;
	strcpy(name, "InputData/FLUIDdomainStart.ac");
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
	
	
	strcpy(name, "InputData/FLUIDdomainStart.ac");
	printf("Opening %s \n", name);
	fIN = fopen(name, "r");
		for(m=1;m<=linecount;m++){
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
				fscanf(fIN,"%f",&bx);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&by);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&bz);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%d",&EulNeighval);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%d",&Windexval);
				if(bvalue >= 0){
					i = (int)((xp-X0)/DX + 1.01);
					j = (int)((yp-Y0)/DX + 1.01);
					k = (int)((zp-Z0)/DX + 1.01);
					
					l = DOMID[i][j][k];
					
					RestartVal[i][j][k] = 1;
					
					D[l] = DD;
					U[l] = UU;
					V[l] = VV;
					W[l] = WW;
					BODYFORCE[0][l] = bx;
					BODYFORCE[1][l] = by;
					BODYFORCE[2][l] = bz;
					EulNeigh[0][l] = EulNeighval;
					
					W_index[i][j][k] = Windexval;
					
				}
						
		}
	fclose(fIN);
	
	linecount = 0;
	strcpy(name, "InputData/FLUIDdomainStart.ad");
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
	
	
	strcpy(name, "InputData/FLUIDdomainStart.ad");
	printf("Opening %s \n", name);
	fIN = fopen(name, "r");
		for(m=1;m<=linecount;m++){
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
				fscanf(fIN,"%f",&bx);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&by);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&bz);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%d",&EulNeighval);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%d",&Windexval);
				if(bvalue >= 0){
					i = (int)((xp-X0)/DX + 1.01);
					j = (int)((yp-Y0)/DX + 1.01);
					k = (int)((zp-Z0)/DX + 1.01);
					
					l = DOMID[i][j][k];
					
					RestartVal[i][j][k] = 1;
					
					D[l] = DD;
					U[l] = UU;
					V[l] = VV;
					W[l] = WW;
					BODYFORCE[0][l] = bx;
					BODYFORCE[1][l] = by;
					BODYFORCE[2][l] = bz;
					EulNeigh[0][l] = EulNeighval;
					
					W_index[i][j][k] = Windexval;
					
				}
						
		}
	fclose(fIN);
	
	linecount = 0;
	strcpy(name, "InputData/FLUIDdomainStart.ae");
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
	
	
	strcpy(name, "InputData/FLUIDdomainStart.ae");
	printf("Opening %s \n", name);
	fIN = fopen(name, "r");
		for(m=1;m<=linecount;m++){
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
				fscanf(fIN,"%f",&bx);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&by);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&bz);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%d",&EulNeighval);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%d",&Windexval);
				if(bvalue >= 0){
					i = (int)((xp-X0)/DX + 1.01);
					j = (int)((yp-Y0)/DX + 1.01);
					k = (int)((zp-Z0)/DX + 1.01);
					
					l = DOMID[i][j][k];
					
					RestartVal[i][j][k] = 1;
					
					D[l] = DD;
					U[l] = UU;
					V[l] = VV;
					W[l] = WW;
					BODYFORCE[0][l] = bx;
					BODYFORCE[1][l] = by;
					BODYFORCE[2][l] = bz;
					EulNeigh[0][l] = EulNeighval;
					
					W_index[i][j][k] = Windexval;
					
				}
						
		}
	fclose(fIN);
	
	linecount = 0;
	strcpy(name, "InputData/FLUIDdomainStart.af");
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
	
	
	strcpy(name, "InputData/FLUIDdomainStart.af");
	printf("Opening %s \n", name);
	fIN = fopen(name, "r");
		for(m=1;m<=linecount;m++){
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
				fscanf(fIN,"%f",&bx);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&by);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&bz);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%d",&EulNeighval);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%d",&Windexval);
				if(bvalue >= 0){
					i = (int)((xp-X0)/DX + 1.01);
					j = (int)((yp-Y0)/DX + 1.01);
					k = (int)((zp-Z0)/DX + 1.01);
					
					l = DOMID[i][j][k];
					
					RestartVal[i][j][k] = 1;
					
					D[l] = DD;
					U[l] = UU;
					V[l] = VV;
					W[l] = WW;
					BODYFORCE[0][l] = bx;
					BODYFORCE[1][l] = by;
					BODYFORCE[2][l] = bz;
					EulNeigh[0][l] = EulNeighval;
					
					W_index[i][j][k] = Windexval;
					
				}
						
		}
	fclose(fIN);
	
	linecount = 0;
	strcpy(name, "InputData/FLUIDdomainStart.ag");
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
	
	
	strcpy(name, "InputData/FLUIDdomainStart.ag");
	printf("Opening %s \n", name);
	fIN = fopen(name, "r");
		for(m=1;m<=linecount;m++){
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
				fscanf(fIN,"%f",&bx);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&by);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&bz);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%d",&EulNeighval);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%d",&Windexval);
				if(bvalue >= 0){
					i = (int)((xp-X0)/DX + 1.01);
					j = (int)((yp-Y0)/DX + 1.01);
					k = (int)((zp-Z0)/DX + 1.01);
					
					l = DOMID[i][j][k];
					
					RestartVal[i][j][k] = 1;
					
					D[l] = DD;
					U[l] = UU;
					V[l] = VV;
					W[l] = WW;
					BODYFORCE[0][l] = bx;
					BODYFORCE[1][l] = by;
					BODYFORCE[2][l] = bz;
					EulNeigh[0][l] = EulNeighval;
					
					W_index[i][j][k] = Windexval;
					
				}
						
		}
	fclose(fIN);
	
	linecount = 0;
	strcpy(name, "InputData/FLUIDdomainStart.ah");
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
	
	
	strcpy(name, "InputData/FLUIDdomainStart.ah");
	printf("Opening %s \n", name);
	fIN = fopen(name, "r");
		for(m=1;m<=linecount;m++){
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
				fscanf(fIN,"%f",&bx);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&by);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&bz);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%d",&EulNeighval);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%d",&Windexval);
				if(bvalue >= 0){
					i = (int)((xp-X0)/DX + 1.01);
					j = (int)((yp-Y0)/DX + 1.01);
					k = (int)((zp-Z0)/DX + 1.01);
					
					l = DOMID[i][j][k];
					
					RestartVal[i][j][k] = 1;
					
					D[l] = DD;
					U[l] = UU;
					V[l] = VV;
					W[l] = WW;
					BODYFORCE[0][l] = bx;
					BODYFORCE[1][l] = by;
					BODYFORCE[2][l] = bz;
					EulNeigh[0][l] = EulNeighval;
					
					W_index[i][j][k] = Windexval;
					
				}
						
		}
	fclose(fIN);
	
	linecount = 0;
	strcpy(name, "InputData/FLUIDdomainStart.ai");
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
	
	
	strcpy(name, "InputData/FLUIDdomainStart.ai");
	printf("Opening %s \n", name);
	fIN = fopen(name, "r");
		for(m=1;m<=linecount;m++){
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
				fscanf(fIN,"%f",&bx);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&by);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&bz);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%d",&EulNeighval);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%d",&Windexval);
				if(bvalue >= 0){
					i = (int)((xp-X0)/DX + 1.01);
					j = (int)((yp-Y0)/DX + 1.01);
					k = (int)((zp-Z0)/DX + 1.01);
					
					l = DOMID[i][j][k];
					
					RestartVal[i][j][k] = 1;
					
					D[l] = DD;
					U[l] = UU;
					V[l] = VV;
					W[l] = WW;
					BODYFORCE[0][l] = bx;
					BODYFORCE[1][l] = by;
					BODYFORCE[2][l] = bz;
					EulNeigh[0][l] = EulNeighval;
					
					W_index[i][j][k] = Windexval;
					
				}
						
		}
	fclose(fIN);
	
	linecount = 0;
	strcpy(name, "InputData/FLUIDdomainStart.aj");
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
	
	
	strcpy(name, "InputData/FLUIDdomainStart.aj");
	printf("Opening %s \n", name);
	fIN = fopen(name, "r");
		for(m=1;m<=linecount;m++){
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
				fscanf(fIN,"%f",&bx);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&by);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&bz);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%d",&EulNeighval);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%d",&Windexval);
				if(bvalue >= 0){
					i = (int)((xp-X0)/DX + 1.01);
					j = (int)((yp-Y0)/DX + 1.01);
					k = (int)((zp-Z0)/DX + 1.01);
					
					l = DOMID[i][j][k];
					
					RestartVal[i][j][k] = 1;
					
					D[l] = DD;
					U[l] = UU;
					V[l] = VV;
					W[l] = WW;
					BODYFORCE[0][l] = bx;
					BODYFORCE[1][l] = by;
					BODYFORCE[2][l] = bz;
					EulNeigh[0][l] = EulNeighval;
					
					W_index[i][j][k] = Windexval;
					
				}
						
		}
	fclose(fIN);
	/// //////////////////////////
	
	LBM_PARENTDOMAIN_SETUP();
	LBM_PARENTBC_SETUP();
	LBM_PARENTINTERPOLATESOLUTION();
	
	
	linecount = 0;
	strcpy(name, "InputData/DAparentdomainStart.csv");
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
	
	strcpy(name, "InputData/DAparentdomainStart.csv");
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
				fscanf(fIN,"%f",&bx);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&by);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&bz);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%d",&EulNeighval);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%d",&Windexval);
				
				if(bvalue >= 0){
					i = (int)((xp-X0_1)/DX + 1.01);
					j = (int)((yp-Y0_1)/DX + 1.01);
					k = (int)((zp-Z0_1)/DX + 1.01);
					
					B_index_1[i][j][k] = bvalue;
					
					l = DOMID_1[i][j][k];
					
					
					D_1[l] = DD;
					U_1[l] = UU;
					V_1[l] = VV;
					W_1[l] = WW;
					BODYFORCE_1[0][l] = bx;
					BODYFORCE_1[1][l] = by;
					BODYFORCE_1[2][l] = bz;
					
					EulNeigh_1[0][l] = EulNeighval;
					
					W_index_1[i][j][k] = Windexval;
				}
			}			
		}
	fclose(fIN);
	
	linecount = 0;
	strcpy(name, "InputData/PCV1parentdomainStart.csv");
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
		
	strcpy(name, "InputData/PCV1parentdomainStart.csv");
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
				fscanf(fIN,"%f",&bx);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&by);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&bz);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%d",&EulNeighval);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%d",&Windexval);
				
				if(bvalue >= 0){
					i = (int)((xp-X0_2)/DX + 1.01);
					j = (int)((yp-Y0_2)/DX + 1.01);
					k = (int)((zp-Z0_2)/DX + 1.01);
					
					B_index_2[i][j][k] = bvalue;
					
					l = DOMID_2[i][j][k];
					
					D_2[l] = DD;
					U_2[l] = UU;
					V_2[l] = VV;
					W_2[l] = WW;
					BODYFORCE_2[0][l] = bx;
					BODYFORCE_2[1][l] = by;
					BODYFORCE_2[2][l] = bz;
					
					EulNeigh_2[0][l] = EulNeighval;
					
					W_index_2[i][j][k] = Windexval;
				}
			}			
		}
	fclose(fIN);
	
	linecount = 0;
	strcpy(name, "InputData/PCV2parentdomainStart.csv");
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
	
	strcpy(name, "InputData/PCV2parentdomainStart.csv");
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
				fscanf(fIN,"%f",&bx);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&by);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&bz);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%d",&EulNeighval);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%d",&Windexval);
				
				if(bvalue >= 0){
					i = (int)((xp-X0_3)/DX + 1.01);
					j = (int)((yp-Y0_3)/DX + 1.01);
					k = (int)((zp-Z0_3)/DX + 1.01);
					
					B_index_3[i][j][k] = bvalue;
					
					l = DOMID_3[i][j][k];
					
					D_3[l] = DD;
					U_3[l] = UU;
					V_3[l] = VV;
					W_3[l] = WW;
					BODYFORCE_3[0][l] = bx;
					BODYFORCE_3[1][l] = by;
					BODYFORCE_3[2][l] = bz;
					
					EulNeigh_3[0][l] = EulNeighval;
					
					W_index_3[i][j][k] = Windexval;
				}
			}			
		}
	fclose(fIN);
	
	
	
	/// main loop for LBM
	iter = 0;
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
	
	
	
	printf("Initializing near wall flags for ROI lattice nodes to zero\n");
	///Set-up NearWallDomain
	for(k=1;k<=KM;k++){
		for(j=1;j<=JM;j++){
			for(i=1;i<=IM;i++){
				W_index[i][j][k] = 0;
				Wdistance[i][j][k] = 5000.;
			}
		}
	}
	
	printf("Initializing near wall flags for DAP lattice nodes to zero\n");
	for(k=1;k<=KM_1;k++){
		for(j=1;j<=JM_1;j++){
			for(i=1;i<=IM_1;i++){
				W_index_1[i][j][k] = 0;
				Wdistance_1[i][j][k] = 5000.;
			}
		}
	}
	
	printf("Initializing near wall flags for PCV1 lattice nodes to zero\n");
	for(k=1;k<=KM_2;k++){
		for(j=1;j<=JM_2;j++){
			for(i=1;i<=IM_2;i++){
				W_index_2[i][j][k] = 0;
				Wdistance_2[i][j][k] = 5000.;
			}
		}
	}
	
	printf("Initializing near wall flags for PCV2 lattice nodes to zero\n");
	for(k=1;k<=KM_3;k++){
		for(j=1;j<=JM_3;j++){
			for(i=1;i<=IM_3;i++){
				W_index_3[i][j][k] = 0;
				Wdistance_3[i][j][k] = 5000.;
			}
		}
	}
	
	
	/// ///////////////////////
	if(NodeW%3 == 0){
		counter = NodeW/3;
		check = 0;
	}
	else{
		counter = NodeW/3 + 1;
		check = NodeW%3;
	}
	
	strcpy(name, "InputData/WallMeshNodeCoord.dat");	
	printf("Opening %s \n", name);
	fIN = fopen(name, "r");
		fgets(bulk, 180, fIN);
		for(m=1;m<=counter;m++){
			l = (m-1)*3+1;
			fscanf(fIN,"%f",&WALLNODECOORD[0][l]);
			fscanf(fIN,"%f",&WALLNODECOORD[1][l]);
			fscanf(fIN,"%f",&WALLNODECOORD[2][l]);
			
			
			///populate the eulerian neighborlist W_index that stores the node of the wall cell closest to lattice node
			xc = WALLNODECOORD[0][l];
			yc = WALLNODECOORD[1][l];
			zc = WALLNODECOORD[2][l];
							
			i = (int)((xc-X0)/DX + 1.01);    j = (int)((yc-Y0)/DX + 1.01);    k = (int)((zc-Z0)/DX + 1.01);
			i1 = i - 4; i2 = i + 4;		j1 = j - 4; j2 = j + 4;		k1 = k - 4; k2 = k + 4;
			
			if(i1 < 1) i1 = 1; if(i2 > IM) i2 = IM;
			if(j1 < 1) j1 = 1; if(j2 > JM) j2 = JM;
			if(k1 < 1) k1 = 1; if(k2 > KM) k2 = KM;
			
			for(i=i1;i<=i2;i++){
				for(j=j1;j<=j2;j++){
					for(k=k1;k<=k2;k++){
						distance = sqrtf((xc-X[i][j][k])*(xc-X[i][j][k])+(yc-Y[i][j][k])*(yc-Y[i][j][k])+(zc-Z[i][j][k])*(zc-Z[i][j][k]));
						
						if(distance < Wdistance[i][j][k]){
							Wdistance[i][j][k] = distance;
							W_index[i][j][k] = l;
						}
					}
				}
			}
			
			if(m == counter && check == 1) continue;
			
			l = (m-1)*3+2;
			fscanf(fIN,"%f",&WALLNODECOORD[0][l]);
			fscanf(fIN,"%f",&WALLNODECOORD[1][l]);
			fscanf(fIN,"%f",&WALLNODECOORD[2][l]);
			
			
			///populate the eulerian neighborlist W_index that stores the node of the wall cell closest to lattice node
			xc = WALLNODECOORD[0][l];
			yc = WALLNODECOORD[1][l];
			zc = WALLNODECOORD[2][l];
							
			i = (int)((xc-X0)/DX + 1.01);    j = (int)((yc-Y0)/DX + 1.01);    k = (int)((zc-Z0)/DX + 1.01);
			i1 = i - 4; i2 = i + 4;		j1 = j - 4; j2 = j + 4;		k1 = k - 4; k2 = k + 4;
			
			if(i1 < 1) i1 = 1; if(i2 > IM) i2 = IM;
			if(j1 < 1) j1 = 1; if(j2 > JM) j2 = JM;
			if(k1 < 1) k1 = 1; if(k2 > KM) k2 = KM;
			
			for(i=i1;i<=i2;i++){
				for(j=j1;j<=j2;j++){
					for(k=k1;k<=k2;k++){
						distance = sqrtf((xc-X[i][j][k])*(xc-X[i][j][k])+(yc-Y[i][j][k])*(yc-Y[i][j][k])+(zc-Z[i][j][k])*(zc-Z[i][j][k]));
						
						if(distance < Wdistance[i][j][k]){
							Wdistance[i][j][k] = distance;
							W_index[i][j][k] = l;
						}
					}
				}
			}
			
			if(m == counter && check == 2) continue;
			
			l = (m-1)*3+3;
			fscanf(fIN,"%f",&WALLNODECOORD[0][l]);
			fscanf(fIN,"%f",&WALLNODECOORD[1][l]);
			fscanf(fIN,"%f",&WALLNODECOORD[2][l]);
			
			
			///populate the eulerian neighborlist W_index that stores the node of the wall cell closest to lattice node
			xc = WALLNODECOORD[0][l];
			yc = WALLNODECOORD[1][l];
			zc = WALLNODECOORD[2][l];
							
			i = (int)((xc-X0)/DX + 1.01);    j = (int)((yc-Y0)/DX + 1.01);    k = (int)((zc-Z0)/DX + 1.01);
			i1 = i - 4; i2 = i + 4;		j1 = j - 4; j2 = j + 4;		k1 = k - 4; k2 = k + 4;
			
			if(i1 < 1) i1 = 1; if(i2 > IM) i2 = IM;
			if(j1 < 1) j1 = 1; if(j2 > JM) j2 = JM;
			if(k1 < 1) k1 = 1; if(k2 > KM) k2 = KM;
			
			for(i=i1;i<=i2;i++){
				for(j=j1;j<=j2;j++){
					for(k=k1;k<=k2;k++){
						distance = sqrtf((xc-X[i][j][k])*(xc-X[i][j][k])+(yc-Y[i][j][k])*(yc-Y[i][j][k])+(zc-Z[i][j][k])*(zc-Z[i][j][k]));
						
						if(distance < Wdistance[i][j][k]){
							Wdistance[i][j][k] = distance;
							W_index[i][j][k] = l;
						}
					}
				}
			}
		}		
	fclose(fIN);
	
	
	//if(CellW%3 == 0){
		//counter = CellW/3;
		//check = 0;
	//}
	//else{
		//counter = CellW/3 + 1;
		//check = CellW%3;
	//}
		
	strcpy(name, "InputData/WallMeshCells.dat");	
	printf("Opening %s \n", name);
	fIN = fopen(name, "r");
		fgets(bulk, 180, fIN);
		for(m=1;m<=CellW;m++){
			fscanf(fIN,"%d",&junkint);
			fscanf(fIN,"%d",&i);
			fscanf(fIN,"%d",&j);
			fscanf(fIN,"%d",&k);
			WALLCELLNODE[0][m] = i+1;
			WALLCELLNODE[1][m] = j+1;
			WALLCELLNODE[2][m] = k+1;
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
			fscanf(fIN,"%f",&WALLCELLNORMALS[0][m]);
			fscanf(fIN,"%f",&WALLCELLNORMALS[1][m]);
			fscanf(fIN,"%f",&WALLCELLNORMALS[2][m]);
		}	
	fclose(fIN);
	
	strcpy(name, "InputData/WallMeshCellNormals.ab");	
	printf("Opening %s \n", name);
	fIN = fopen(name, "r");
		for(m=counter+1;m<=CellW;m++){
			fscanf(fIN,"%f",&WALLCELLNORMALS[0][m]);
			fscanf(fIN,"%f",&WALLCELLNORMALS[1][m]);
			fscanf(fIN,"%f",&WALLCELLNORMALS[2][m]);
		}	
	fclose(fIN);
	/// //////////////
	
	///initialize Wall mesh connectivity list
	for(m=1;m<=NodeW;m++){
		for(i=0;i<=12;i++){
			WALLNODETOCELLNEIGHBORLIST[i][m] = 0;
		}
	}
	
	///generate Wall mesh connectivity list
	maxNodeW = 0;
	for(m=1;m<=CellW;m++){
		
		n1 = WALLCELLNODE[0][m];
		WALLNODETOCELLNEIGHBORLIST[0][n1] += 1;
		check = WALLNODETOCELLNEIGHBORLIST[0][n1];
		if(check > 24) printf("List size exceeded.. ");
		WALLNODETOCELLNEIGHBORLIST[check][n1] = m;
		
		if(check > maxNodeW) maxNodeW = check;
		
		n2 = WALLCELLNODE[1][m];
		WALLNODETOCELLNEIGHBORLIST[0][n2] += 1;
		check = WALLNODETOCELLNEIGHBORLIST[0][n2];
		if(check > 24) printf("List size exceeded.. ");
		WALLNODETOCELLNEIGHBORLIST[check][n2] = m;
		
		if(check > maxNodeW) maxNodeW = check;
		
		n3 = WALLCELLNODE[2][m];
		WALLNODETOCELLNEIGHBORLIST[0][n3] += 1;
		check = WALLNODETOCELLNEIGHBORLIST[0][n3];
		if(check > 24) printf("List size exceeded.. ");
		WALLNODETOCELLNEIGHBORLIST[check][n3] = m;
		
		if(check > maxNodeW) maxNodeW = check;
	}
	
	
	
	printf("\n ROI wall mesh max connectivity is %d\n", maxNodeW);
	
	
	strcpy(name, "WallMeshDataStructure/WallNodeToTriConnectivity.csv");
	printf("Finished ouputting node to tri connectivity to %s...\n", name);
	fOUT = fopen(name, "w");
		fprintf(fOUT, "X Y Z TotalTri T1 T2 T3 T4 T5 T6 T7 T8 T9 T10 T11 T12\n");
		for(m=1;m<=NodeW;m++){
			fprintf(fOUT, "%.7e %.7e %.7e %d %d %d %d %d %d %d %d %d %d %d %d %d\n", WALLNODECOORD[0][m],WALLNODECOORD[1][m],WALLNODECOORD[2][m],WALLNODETOCELLNEIGHBORLIST[0][m],WALLNODETOCELLNEIGHBORLIST[1][m],WALLNODETOCELLNEIGHBORLIST[2][m],WALLNODETOCELLNEIGHBORLIST[3][m],WALLNODETOCELLNEIGHBORLIST[4][m],WALLNODETOCELLNEIGHBORLIST[5][m],WALLNODETOCELLNEIGHBORLIST[6][m],WALLNODETOCELLNEIGHBORLIST[7][m],WALLNODETOCELLNEIGHBORLIST[8][m],WALLNODETOCELLNEIGHBORLIST[9][m],WALLNODETOCELLNEIGHBORLIST[10][m],WALLNODETOCELLNEIGHBORLIST[11][m],WALLNODETOCELLNEIGHBORLIST[12][m]);
		}
	fclose(fOUT);
	printf("Finished ouputting node to tri connectivity to %s...\n", name);	
	
	strcpy(name, "WallMeshDataStructure/WallTriNormals.csv");
	printf("Finished ouputting node to tri connectivity to %s...\n", name);
	fOUT = fopen(name, "w");
		fprintf(fOUT, "X Y Z Nx Ny Nz\n");
		for(tri=1;tri<=CellW;tri++){
			n1 = WALLCELLNODE[0][tri];
			n2 = WALLCELLNODE[1][tri];
			n3 = WALLCELLNODE[2][tri];
			
			xc = (WALLNODECOORD[0][n1]+WALLNODECOORD[0][n2]+WALLNODECOORD[0][n3])/3.;
			yc = (WALLNODECOORD[1][n1]+WALLNODECOORD[1][n2]+WALLNODECOORD[1][n3])/3.;
			zc = (WALLNODECOORD[2][n1]+WALLNODECOORD[2][n2]+WALLNODECOORD[2][n3])/3.;
			
			fprintf(fOUT, "%.7e %.7e %.7e %.7e %.7e %.7e\n", xc,yc,zc,WALLCELLNORMALS[0][tri],WALLCELLNORMALS[1][tri],WALLCELLNORMALS[2][tri]);
		}
	fclose(fOUT);
	printf("Finished ouputting node to tri connectivity to %s...\n", name);	
	
	///DA wall setup
	if(NodeW_1%3 == 0){
		counter = NodeW_1/3;
		check = 0;
	}
	else{
		counter = NodeW_1/3 + 1;
		check = NodeW_1%3;
	}
	
	strcpy(name, "InputData/DAWallMeshNodeCoord.dat");	
	printf("Opening %s \n", name);
	fIN = fopen(name, "r");
		fgets(bulk, 180, fIN);
		for(m=1;m<=counter;m++){
			l = (m-1)*3+1;
			fscanf(fIN,"%f",&WALLNODECOORD_1[0][l]);
			fscanf(fIN,"%f",&WALLNODECOORD_1[1][l]);
			fscanf(fIN,"%f",&WALLNODECOORD_1[2][l]);
			
			xc = WALLNODECOORD_1[0][l];
			yc = WALLNODECOORD_1[1][l];
			zc = WALLNODECOORD_1[2][l];
							
			i = (int)((xc-X0_1)/DX + 1.01);    j = (int)((yc-Y0_1)/DX + 1.01);    k = (int)((zc-Z0_1)/DX + 1.01);
			i1 = i - 4; i2 = i + 4;		j1 = j - 4; j2 = j + 4;		k1 = k - 4; k2 = k + 4;
			
			if(i1 < 1) i1 = 1; if(i2 > IM_1) i2 = IM_1;
			if(j1 < 1) j1 = 1; if(j2 > JM_1) j2 = JM_1;
			if(k1 < 1) k1 = 1; if(k2 > KM_1) k2 = KM_1;
			
			for(i=i1;i<=i2;i++){
				for(j=j1;j<=j2;j++){
					for(k=k1;k<=k2;k++){
						distance = sqrtf((xc-X_1[i][j][k])*(xc-X_1[i][j][k])+(yc-Y_1[i][j][k])*(yc-Y_1[i][j][k])+(zc-Z_1[i][j][k])*(zc-Z_1[i][j][k]));
						
						if(distance < Wdistance_1[i][j][k]){
							Wdistance_1[i][j][k] = distance;
							W_index_1[i][j][k] = l;
						}
					}
				}
			}
			
			if(m == counter && check == 1) continue;
			l = (m-1)*3+2;
			fscanf(fIN,"%f",&WALLNODECOORD_1[0][l]);
			fscanf(fIN,"%f",&WALLNODECOORD_1[1][l]);
			fscanf(fIN,"%f",&WALLNODECOORD_1[2][l]);
			
			xc = WALLNODECOORD_1[0][l];
			yc = WALLNODECOORD_1[1][l];
			zc = WALLNODECOORD_1[2][l];
							
			i = (int)((xc-X0_1)/DX + 1.01);    j = (int)((yc-Y0_1)/DX + 1.01);    k = (int)((zc-Z0_1)/DX + 1.01);
			i1 = i - 4; i2 = i + 4;		j1 = j - 4; j2 = j + 4;		k1 = k - 4; k2 = k + 4;
			
			if(i1 < 1) i1 = 1; if(i2 > IM_1) i2 = IM_1;
			if(j1 < 1) j1 = 1; if(j2 > JM_1) j2 = JM_1;
			if(k1 < 1) k1 = 1; if(k2 > KM_1) k2 = KM_1;
			
			for(i=i1;i<=i2;i++){
				for(j=j1;j<=j2;j++){
					for(k=k1;k<=k2;k++){
						distance = sqrtf((xc-X_1[i][j][k])*(xc-X_1[i][j][k])+(yc-Y_1[i][j][k])*(yc-Y_1[i][j][k])+(zc-Z_1[i][j][k])*(zc-Z_1[i][j][k]));
						
						if(distance < Wdistance_1[i][j][k]){
							Wdistance_1[i][j][k] = distance;
							W_index_1[i][j][k] = l;
						}
					}
				}
			}
			
			if(m == counter && check == 2) continue;
			l = (m-1)*3+3;
			fscanf(fIN,"%f",&WALLNODECOORD_1[0][l]);
			fscanf(fIN,"%f",&WALLNODECOORD_1[1][l]);
			fscanf(fIN,"%f",&WALLNODECOORD_1[2][l]);
			
			xc = WALLNODECOORD_1[0][l];
			yc = WALLNODECOORD_1[1][l];
			zc = WALLNODECOORD_1[2][l];
							
			i = (int)((xc-X0_1)/DX + 1.01);    j = (int)((yc-Y0_1)/DX + 1.01);    k = (int)((zc-Z0_1)/DX + 1.01);
			i1 = i - 4; i2 = i + 4;		j1 = j - 4; j2 = j + 4;		k1 = k - 4; k2 = k + 4;
			
			if(i1 < 1) i1 = 1; if(i2 > IM_1) i2 = IM_1;
			if(j1 < 1) j1 = 1; if(j2 > JM_1) j2 = JM_1;
			if(k1 < 1) k1 = 1; if(k2 > KM_1) k2 = KM_1;
			
			for(i=i1;i<=i2;i++){
				for(j=j1;j<=j2;j++){
					for(k=k1;k<=k2;k++){
						distance = sqrtf((xc-X_1[i][j][k])*(xc-X_1[i][j][k])+(yc-Y_1[i][j][k])*(yc-Y_1[i][j][k])+(zc-Z_1[i][j][k])*(zc-Z_1[i][j][k]));
						
						if(distance < Wdistance_1[i][j][k]){
							Wdistance_1[i][j][k] = distance;
							W_index_1[i][j][k] = l;
						}
					}
				}
			}
		}		
	fclose(fIN);
	
	if(CellW_1%3 == 0){
		counter = CellW_1/3;
		check = 0;
	}
	else{
		counter = CellW_1/3 + 1;
		check = CellW_1%3;
	}
	
	strcpy(name, "InputData/DAWallMeshCells.dat");	
	printf("Opening %s \n", name);
	fIN = fopen(name, "r");
		fgets(bulk, 180, fIN);
		for(m=1;m<=CellW_1;m++){
			fscanf(fIN,"%d",&junkint);
			fscanf(fIN,"%d",&i);
			fscanf(fIN,"%d",&j);
			fscanf(fIN,"%d",&k);
			WALLCELLNODE_1[0][m] = i+1;
			WALLCELLNODE_1[1][m] = j+1;
			WALLCELLNODE_1[2][m] = k+1;
		}		
	fclose(fIN);
	
	strcpy(name, "InputData/DAWallMeshCellNormals.dat");	
	printf("Opening %s \n", name);
	fIN = fopen(name, "r");
		fgets(bulk, 180, fIN);
		fgets(bulk, 180, fIN);
		for(m=1;m<=counter;m++){
			l = (m-1)*3+1;
			fscanf(fIN,"%f",&WALLCELLNORMALS_1[0][l]);
			fscanf(fIN,"%f",&WALLCELLNORMALS_1[1][l]);
			fscanf(fIN,"%f",&WALLCELLNORMALS_1[2][l]);
			
			if(m == counter && check == 1) continue;
			
			l = (m-1)*3+2;
			fscanf(fIN,"%f",&WALLCELLNORMALS_1[0][l]);
			fscanf(fIN,"%f",&WALLCELLNORMALS_1[1][l]);
			fscanf(fIN,"%f",&WALLCELLNORMALS_1[2][l]);
			
			if(m == counter && check == 2) continue;
			
			l = (m-1)*3+3;
			fscanf(fIN,"%f",&WALLCELLNORMALS_1[0][l]);
			fscanf(fIN,"%f",&WALLCELLNORMALS_1[1][l]);
			fscanf(fIN,"%f",&WALLCELLNORMALS_1[2][l]);
		}		
	fclose(fIN);
	
	///initialize Wall mesh connectivity list
	for(m=1;m<=NodeW_1;m++){
		for(i=0;i<=6;i++){
			WALLNODETOCELLNEIGHBORLIST_1[i][m] = 0;
		}
	}
	
	///generate Wall mesh connectivity list
	maxNodeW = 0;
	for(m=1;m<=CellW_1;m++){
		
		n1 = WALLCELLNODE_1[0][m];
		WALLNODETOCELLNEIGHBORLIST_1[0][n1] += 1;
		check = WALLNODETOCELLNEIGHBORLIST_1[0][n1];
		if(check > 10) printf("List size exceeded.. ");
		WALLNODETOCELLNEIGHBORLIST_1[check][n1] = m;
		
		if(check > maxNodeW) maxNodeW = check;
		
		n2 = WALLCELLNODE_1[1][m];
		WALLNODETOCELLNEIGHBORLIST_1[0][n2] += 1;
		check = WALLNODETOCELLNEIGHBORLIST_1[0][n2];
		if(check > 10) printf("List size exceeded.. ");
		WALLNODETOCELLNEIGHBORLIST_1[check][n2] = m;
		
		if(check > maxNodeW) maxNodeW = check;
		
		n3 = WALLCELLNODE_1[2][m];
		WALLNODETOCELLNEIGHBORLIST_1[0][n3] += 1;
		check = WALLNODETOCELLNEIGHBORLIST_1[0][n3];
		if(check > 10) printf("List size exceeded.. ");
		WALLNODETOCELLNEIGHBORLIST_1[check][n3] = m;
		
		if(check > maxNodeW) maxNodeW = check;
	}
	printf("\n DA wall mesh max connectivity is %d\n", maxNodeW);
	
	
	
	strcpy(name, "WallMeshDataStructure/DAWallNodeToTriConnectivity.csv");
	printf("Finished ouputting node to tri connectivity to %s...\n", name);
	fOUT = fopen(name, "w");
		fprintf(fOUT, "X Y Z TotalTri T1 T2 T3 T4 T5 T6\n");
		for(m=1;m<=NodeW_1;m++){
			fprintf(fOUT, "%.7e %.7e %.7e %d %d %d %d %d %d %d\n",WALLNODECOORD_1[0][m],WALLNODECOORD_1[1][m],WALLNODECOORD_1[2][m],WALLNODETOCELLNEIGHBORLIST_1[0][m],WALLNODETOCELLNEIGHBORLIST_1[1][m],WALLNODETOCELLNEIGHBORLIST_1[2][m],WALLNODETOCELLNEIGHBORLIST_1[3][m],WALLNODETOCELLNEIGHBORLIST_1[4][m],WALLNODETOCELLNEIGHBORLIST_1[5][m],WALLNODETOCELLNEIGHBORLIST_1[6][m]);
		}
	fclose(fOUT);
	printf("Finished ouputting node to tri connectivity to %s...\n", name);
	
	
	strcpy(name, "WallMeshDataStructure/DAWallTriNormals.csv");
	printf("Finished ouputting node to tri connectivity to %s...\n", name);
	fOUT = fopen(name, "w");
		fprintf(fOUT, "X Y Z Nx Ny Nz\n");
		for(tri=1;tri<=CellW_1;tri++){
			n1 = WALLCELLNODE_1[0][tri];
			n2 = WALLCELLNODE_1[1][tri];
			n3 = WALLCELLNODE_1[2][tri];
			
			xc = (WALLNODECOORD_1[0][n1]+WALLNODECOORD_1[0][n2]+WALLNODECOORD_1[0][n3])/3.;
			yc = (WALLNODECOORD_1[1][n1]+WALLNODECOORD_1[1][n2]+WALLNODECOORD_1[1][n3])/3.;
			zc = (WALLNODECOORD_1[2][n1]+WALLNODECOORD_1[2][n2]+WALLNODECOORD_1[2][n3])/3.;
			
			fprintf(fOUT, "%.7e %.7e %.7e %.7e %.7e %.7e\n", xc,yc,zc,WALLCELLNORMALS_1[0][tri],WALLCELLNORMALS_1[1][tri],WALLCELLNORMALS_1[2][tri]);
		}
	fclose(fOUT);
	printf("Finished ouputting node to tri connectivity to %s...\n", name);	
	
	
	///PCV1 wall setup
	if(NodeW_2%3 == 0){
		counter = NodeW_2/3;
		check = 0;
	}
	else{
		counter = NodeW_2/3 + 1;
		check = NodeW_2%3;
	}
	
	strcpy(name, "InputData/PCV1WallMeshNodeCoord.dat");	
	printf("Opening %s \n", name);
	fIN = fopen(name, "r");
		fgets(bulk, 180, fIN);
		for(m=1;m<=counter;m++){
			l = (m-1)*3+1;
			fscanf(fIN,"%f",&WALLNODECOORD_2[0][l]);
			fscanf(fIN,"%f",&WALLNODECOORD_2[1][l]);
			fscanf(fIN,"%f",&WALLNODECOORD_2[2][l]);
			
			xc = WALLNODECOORD_2[0][l];
			yc = WALLNODECOORD_2[1][l];
			zc = WALLNODECOORD_2[2][l];
							
			i = (int)((xc-X0_2)/DX + 1.01);    j = (int)((yc-Y0_2)/DX + 1.01);    k = (int)((zc-Z0_2)/DX + 1.01);
			i1 = i - 4; i2 = i + 4;		j1 = j - 4; j2 = j + 4;		k1 = k - 4; k2 = k + 4;
			
			if(i1 < 1) i1 = 1; if(i2 > IM_2) i2 = IM_2;
			if(j1 < 1) j1 = 1; if(j2 > JM_2) j2 = JM_2;
			if(k1 < 1) k1 = 1; if(k2 > KM_2) k2 = KM_2;
			
			for(i=i1;i<=i2;i++){
				for(j=j1;j<=j2;j++){
					for(k=k1;k<=k2;k++){
						distance = sqrtf((xc-X_2[i][j][k])*(xc-X_2[i][j][k])+(yc-Y_2[i][j][k])*(yc-Y_2[i][j][k])+(zc-Z_2[i][j][k])*(zc-Z_2[i][j][k]));
						
						if(distance < Wdistance_2[i][j][k]){
							Wdistance_2[i][j][k] = distance;
							W_index_2[i][j][k] = l;
						}
					}
				}
			}
			
			if(m == counter && check == 1) continue;
			l = (m-1)*3+2;
			fscanf(fIN,"%f",&WALLNODECOORD_2[0][l]);
			fscanf(fIN,"%f",&WALLNODECOORD_2[1][l]);
			fscanf(fIN,"%f",&WALLNODECOORD_2[2][l]);
			
			xc = WALLNODECOORD_2[0][l];
			yc = WALLNODECOORD_2[1][l];
			zc = WALLNODECOORD_2[2][l];
							
			i = (int)((xc-X0_2)/DX + 1.01);    j = (int)((yc-Y0_2)/DX + 1.01);    k = (int)((zc-Z0_2)/DX + 1.01);
			i1 = i - 4; i2 = i + 4;		j1 = j - 4; j2 = j + 4;		k1 = k - 4; k2 = k + 4;
			
			if(i1 < 1) i1 = 1; if(i2 > IM_2) i2 = IM_2;
			if(j1 < 1) j1 = 1; if(j2 > JM_2) j2 = JM_2;
			if(k1 < 1) k1 = 1; if(k2 > KM_2) k2 = KM_2;
			
			for(i=i1;i<=i2;i++){
				for(j=j1;j<=j2;j++){
					for(k=k1;k<=k2;k++){
						distance = sqrtf((xc-X_2[i][j][k])*(xc-X_2[i][j][k])+(yc-Y_2[i][j][k])*(yc-Y_2[i][j][k])+(zc-Z_2[i][j][k])*(zc-Z_2[i][j][k]));
						
						if(distance < Wdistance_2[i][j][k]){
							Wdistance_2[i][j][k] = distance;
							W_index_2[i][j][k] = l;
						}
					}
				}
			}
			
			if(m == counter && check == 2) continue;
			l = (m-1)*3+3;
			fscanf(fIN,"%f",&WALLNODECOORD_2[0][l]);
			fscanf(fIN,"%f",&WALLNODECOORD_2[1][l]);
			fscanf(fIN,"%f",&WALLNODECOORD_2[2][l]);
			
			xc = WALLNODECOORD_2[0][l];
			yc = WALLNODECOORD_2[1][l];
			zc = WALLNODECOORD_2[2][l];
							
			i = (int)((xc-X0_2)/DX + 1.01);    j = (int)((yc-Y0_2)/DX + 1.01);    k = (int)((zc-Z0_2)/DX + 1.01);
			i1 = i - 4; i2 = i + 4;		j1 = j - 4; j2 = j + 4;		k1 = k - 4; k2 = k + 4;
			
			if(i1 < 1) i1 = 1; if(i2 > IM_2) i2 = IM_2;
			if(j1 < 1) j1 = 1; if(j2 > JM_2) j2 = JM_2;
			if(k1 < 1) k1 = 1; if(k2 > KM_2) k2 = KM_2;
			
			for(i=i1;i<=i2;i++){
				for(j=j1;j<=j2;j++){
					for(k=k1;k<=k2;k++){
						distance = sqrtf((xc-X_2[i][j][k])*(xc-X_2[i][j][k])+(yc-Y_2[i][j][k])*(yc-Y_2[i][j][k])+(zc-Z_2[i][j][k])*(zc-Z_2[i][j][k]));
						
						if(distance < Wdistance_2[i][j][k]){
							Wdistance_2[i][j][k] = distance;
							W_index_2[i][j][k] = l;
						}
					}
				}
			}
		}	
	fclose(fIN);
	
	if(CellW_2%3 == 0){
		counter = CellW_2/3;
		check = 0;
	}
	else{
		counter = CellW_2/3 + 1;
		check = CellW_2%3;
	}
	
	strcpy(name, "InputData/PCV1WallMeshCells.dat");	
	printf("Opening %s \n", name);
	fIN = fopen(name, "r");
		fgets(bulk, 180, fIN);
		for(m=1;m<=CellW_2;m++){
			fscanf(fIN,"%d",&junkint);
			fscanf(fIN,"%d",&i);
			fscanf(fIN,"%d",&j);
			fscanf(fIN,"%d",&k);
			WALLCELLNODE_2[0][m] = i+1;
			WALLCELLNODE_2[1][m] = j+1;
			WALLCELLNODE_2[2][m] = k+1;
		}		
	fclose(fIN);
	
	strcpy(name, "InputData/PCV1WallMeshCellNormals.dat");	
	printf("Opening %s \n", name);
	fIN = fopen(name, "r");
		fgets(bulk, 180, fIN);
		fgets(bulk, 180, fIN);
		for(m=1;m<=counter;m++){
			l = (m-1)*3+1;
			fscanf(fIN,"%f",&WALLCELLNORMALS_2[0][l]);
			fscanf(fIN,"%f",&WALLCELLNORMALS_2[1][l]);
			fscanf(fIN,"%f",&WALLCELLNORMALS_2[2][l]);
			
			if(m == counter && check == 1) continue;
			
			l = (m-1)*3+2;
			fscanf(fIN,"%f",&WALLCELLNORMALS_2[0][l]);
			fscanf(fIN,"%f",&WALLCELLNORMALS_2[1][l]);
			fscanf(fIN,"%f",&WALLCELLNORMALS_2[2][l]);
			
			if(m == counter && check == 2) continue;
			
			l = (m-1)*3+3;
			fscanf(fIN,"%f",&WALLCELLNORMALS_2[0][l]);
			fscanf(fIN,"%f",&WALLCELLNORMALS_2[1][l]);
			fscanf(fIN,"%f",&WALLCELLNORMALS_2[2][l]);
		}		
	fclose(fIN);
	
	///initialize Wall mesh connectivity list
	for(m=1;m<=NodeW_2;m++){
		for(i=0;i<=6;i++){
			WALLNODETOCELLNEIGHBORLIST_2[i][m] = 0;
		}
	}
	
	///generate Wall mesh connectivity list
	maxNodeW = 0;
	for(m=1;m<=CellW_2;m++){
		
		n1 = WALLCELLNODE_2[0][m];
		WALLNODETOCELLNEIGHBORLIST_2[0][n1] += 1;
		check = WALLNODETOCELLNEIGHBORLIST_2[0][n1];
		if(check > 10) printf("List size exceeded.. ");
		WALLNODETOCELLNEIGHBORLIST_2[check][n1] = m;
		
		if(check > maxNodeW) maxNodeW = check;
		
		n2 = WALLCELLNODE_2[1][m];
		WALLNODETOCELLNEIGHBORLIST_2[0][n2] += 1;
		check = WALLNODETOCELLNEIGHBORLIST_2[0][n2];
		if(check > 10) printf("List size exceeded.. ");
		WALLNODETOCELLNEIGHBORLIST_2[check][n2] = m;
		
		if(check > maxNodeW) maxNodeW = check;
		
		n3 = WALLCELLNODE_2[2][m];
		WALLNODETOCELLNEIGHBORLIST_2[0][n3] += 1;
		check = WALLNODETOCELLNEIGHBORLIST_2[0][n3];
		if(check > 10) printf("List size exceeded.. ");
		WALLNODETOCELLNEIGHBORLIST_2[check][n3] = m;
		
		if(check > maxNodeW) maxNodeW = check;
	}
	printf("\n PCV1 wall mesh max connectivity is %d\n", maxNodeW);
	
		
	
	strcpy(name, "WallMeshDataStructure/PCV1WallNodeToTriConnectivity.csv");
	printf("Finished ouputting node to tri connectivity to %s...\n", name);
	fOUT = fopen(name, "w");
		fprintf(fOUT, "X Y Z TotalTri T1 T2 T3 T4 T5 T6\n");
		for(m=1;m<=NodeW_2;m++){
			fprintf(fOUT, "%.7e %.7e %.7e %d %d %d %d %d %d %d\n",WALLNODECOORD_2[0][m],WALLNODECOORD_2[1][m],WALLNODECOORD_2[2][m],WALLNODETOCELLNEIGHBORLIST_2[0][m],WALLNODETOCELLNEIGHBORLIST_2[1][m],WALLNODETOCELLNEIGHBORLIST_2[2][m],WALLNODETOCELLNEIGHBORLIST_2[3][m],WALLNODETOCELLNEIGHBORLIST_2[4][m],WALLNODETOCELLNEIGHBORLIST_2[5][m],WALLNODETOCELLNEIGHBORLIST_2[6][m]);
		}
	fclose(fOUT);
	printf("Finished ouputting node to tri connectivity to %s...\n", name);
	
	
	strcpy(name, "WallMeshDataStructure/PCV1WallTriNormals.csv");
	printf("Finished ouputting node to tri connectivity to %s...\n", name);
	fOUT = fopen(name, "w");
		fprintf(fOUT, "X Y Z Nx Ny Nz\n");
		for(tri=1;tri<=CellW_2;tri++){
			n1 = WALLCELLNODE_2[0][tri];
			n2 = WALLCELLNODE_2[1][tri];
			n3 = WALLCELLNODE_2[2][tri];
			
			xc = (WALLNODECOORD_2[0][n1]+WALLNODECOORD_2[0][n2]+WALLNODECOORD_2[0][n3])/3.;
			yc = (WALLNODECOORD_2[1][n1]+WALLNODECOORD_2[1][n2]+WALLNODECOORD_2[1][n3])/3.;
			zc = (WALLNODECOORD_2[2][n1]+WALLNODECOORD_2[2][n2]+WALLNODECOORD_2[2][n3])/3.;
			
			fprintf(fOUT, "%.7e %.7e %.7e %.7e %.7e %.7e\n", xc,yc,zc,WALLCELLNORMALS_2[0][tri],WALLCELLNORMALS_2[1][tri],WALLCELLNORMALS_2[2][tri]);
		}
	fclose(fOUT);
	printf("Finished ouputting node to tri connectivity to %s...\n", name);	
	
	///PCV2 wall setup
	if(NodeW_3%3 == 0){
		counter = NodeW_3/3;
		check = 0;
	}
	else{
		counter = NodeW_3/3 + 1;
		check = NodeW_3%3;
	}
	
	strcpy(name, "InputData/PCV2WallMeshNodeCoord.dat");	
	printf("Opening %s \n", name);
	fIN = fopen(name, "r");
		fgets(bulk, 180, fIN);
		for(m=1;m<=counter;m++){
			l = (m-1)*3+1;
			fscanf(fIN,"%f",&WALLNODECOORD_3[0][l]);
			fscanf(fIN,"%f",&WALLNODECOORD_3[1][l]);
			fscanf(fIN,"%f",&WALLNODECOORD_3[2][l]);
			
			xc = WALLNODECOORD_3[0][l];
			yc = WALLNODECOORD_3[1][l];
			zc = WALLNODECOORD_3[2][l];
							
			i = (int)((xc-X0_3)/DX + 1.01);    j = (int)((yc-Y0_3)/DX + 1.01);    k = (int)((zc-Z0_3)/DX + 1.01);
			i1 = i - 4; i2 = i + 4;		j1 = j - 4; j2 = j + 4;		k1 = k - 4; k2 = k + 4;
			
			if(i1 < 1) i1 = 1; if(i2 > IM_3) i2 = IM_3;
			if(j1 < 1) j1 = 1; if(j2 > JM_3) j2 = JM_3;
			if(k1 < 1) k1 = 1; if(k2 > KM_3) k2 = KM_3;
			
			for(i=i1;i<=i2;i++){
				for(j=j1;j<=j2;j++){
					for(k=k1;k<=k2;k++){
						distance = sqrtf((xc-X_3[i][j][k])*(xc-X_3[i][j][k])+(yc-Y_3[i][j][k])*(yc-Y_3[i][j][k])+(zc-Z_3[i][j][k])*(zc-Z_3[i][j][k]));
						
						if(distance < Wdistance_3[i][j][k]){
							Wdistance_3[i][j][k] = distance;
							W_index_3[i][j][k] = l;
						}
					}
				}
			}
			
			if(m == counter && check == 1) continue;
			l = (m-1)*3+2;
			fscanf(fIN,"%f",&WALLNODECOORD_3[0][l]);
			fscanf(fIN,"%f",&WALLNODECOORD_3[1][l]);
			fscanf(fIN,"%f",&WALLNODECOORD_3[2][l]);
			
			xc = WALLNODECOORD_3[0][l];
			yc = WALLNODECOORD_3[1][l];
			zc = WALLNODECOORD_3[2][l];
							
			i = (int)((xc-X0_3)/DX + 1.01);    j = (int)((yc-Y0_3)/DX + 1.01);    k = (int)((zc-Z0_3)/DX + 1.01);
			i1 = i - 4; i2 = i + 4;		j1 = j - 4; j2 = j + 4;		k1 = k - 4; k2 = k + 4;
			
			if(i1 < 1) i1 = 1; if(i2 > IM_3) i2 = IM_3;
			if(j1 < 1) j1 = 1; if(j2 > JM_3) j2 = JM_3;
			if(k1 < 1) k1 = 1; if(k2 > KM_3) k2 = KM_3;
			
			for(i=i1;i<=i2;i++){
				for(j=j1;j<=j2;j++){
					for(k=k1;k<=k2;k++){
						distance = sqrtf((xc-X_3[i][j][k])*(xc-X_3[i][j][k])+(yc-Y_3[i][j][k])*(yc-Y_3[i][j][k])+(zc-Z_3[i][j][k])*(zc-Z_3[i][j][k]));
						
						if(distance < Wdistance_3[i][j][k]){
							Wdistance_3[i][j][k] = distance;
							W_index_3[i][j][k] = l;
						}
					}
				}
			}
			
			if(m == counter && check == 2) continue;
			l = (m-1)*3+3;
			fscanf(fIN,"%f",&WALLNODECOORD_3[0][l]);
			fscanf(fIN,"%f",&WALLNODECOORD_3[1][l]);
			fscanf(fIN,"%f",&WALLNODECOORD_3[2][l]);
			
			xc = WALLNODECOORD_3[0][l];
			yc = WALLNODECOORD_3[1][l];
			zc = WALLNODECOORD_3[2][l];
							
			i = (int)((xc-X0_3)/DX + 1.01);    j = (int)((yc-Y0_3)/DX + 1.01);    k = (int)((zc-Z0_3)/DX + 1.01);
			i1 = i - 4; i2 = i + 4;		j1 = j - 4; j2 = j + 4;		k1 = k - 4; k2 = k + 4;
			
			if(i1 < 1) i1 = 1; if(i2 > IM_3) i2 = IM_3;
			if(j1 < 1) j1 = 1; if(j2 > JM_3) j2 = JM_3;
			if(k1 < 1) k1 = 1; if(k2 > KM_3) k2 = KM_3;
			
			for(i=i1;i<=i2;i++){
				for(j=j1;j<=j2;j++){
					for(k=k1;k<=k2;k++){
						distance = sqrtf((xc-X_3[i][j][k])*(xc-X_3[i][j][k])+(yc-Y_3[i][j][k])*(yc-Y_3[i][j][k])+(zc-Z_3[i][j][k])*(zc-Z_3[i][j][k]));
						
						if(distance < Wdistance_3[i][j][k]){
							Wdistance_3[i][j][k] = distance;
							W_index_3[i][j][k] = l;
						}
					}
				}
			}
		}	
	fclose(fIN);
	
	if(CellW_3%3 == 0){
		counter = CellW_3/3;
		check = 0;
	}
	else{
		counter = CellW_3/3 + 1;
		check = CellW_3%3;
	}
	
	
	strcpy(name, "InputData/PCV2WallMeshCells.dat");	
	printf("Opening %s \n", name);
	fIN = fopen(name, "r");
		fgets(bulk, 180, fIN);
		for(m=1;m<=CellW_3;m++){
			fscanf(fIN,"%d",&junkint);
			fscanf(fIN,"%d",&i);
			fscanf(fIN,"%d",&j);
			fscanf(fIN,"%d",&k);
			WALLCELLNODE_3[0][m] = i+1;
			WALLCELLNODE_3[1][m] = j+1;
			WALLCELLNODE_3[2][m] = k+1;
		}		
	fclose(fIN);
	
	
	strcpy(name, "InputData/PCV2WallMeshCellNormals.dat");	
	printf("Opening %s \n", name);
	fIN = fopen(name, "r");
		fgets(bulk, 180, fIN);
		fgets(bulk, 180, fIN);
		for(m=1;m<=counter;m++){
			l = (m-1)*3+1;
			fscanf(fIN,"%f",&WALLCELLNORMALS_3[0][l]);
			fscanf(fIN,"%f",&WALLCELLNORMALS_3[1][l]);
			fscanf(fIN,"%f",&WALLCELLNORMALS_3[2][l]);
			
			if(m == counter && check == 1) continue;
			
			l = (m-1)*3+2;
			fscanf(fIN,"%f",&WALLCELLNORMALS_3[0][l]);
			fscanf(fIN,"%f",&WALLCELLNORMALS_3[1][l]);
			fscanf(fIN,"%f",&WALLCELLNORMALS_3[2][l]);
			
			if(m == counter && check == 2) continue;
			
			l = (m-1)*3+3;
			fscanf(fIN,"%f",&WALLCELLNORMALS_3[0][l]);
			fscanf(fIN,"%f",&WALLCELLNORMALS_3[1][l]);
			fscanf(fIN,"%f",&WALLCELLNORMALS_3[2][l]);
		}		
	fclose(fIN);
	
	///initialize Wall mesh connectivity list
	for(m=1;m<=NodeW_3;m++){
		for(i=0;i<=6;i++){
			WALLNODETOCELLNEIGHBORLIST_3[i][m] = 0;
		}
	}
	
	///generate Wall mesh connectivity list
	maxNodeW = 0;
	for(m=1;m<=CellW_3;m++){
		
		n1 = WALLCELLNODE_3[0][m];
		WALLNODETOCELLNEIGHBORLIST_3[0][n1] += 1;
		check = WALLNODETOCELLNEIGHBORLIST_3[0][n1];
		if(check > 6) printf("List size exceeded.. ");
		WALLNODETOCELLNEIGHBORLIST_3[check][n1] = m;
		
		if(check > maxNodeW) maxNodeW = check;
		
		n2 = WALLCELLNODE_3[1][m];
		WALLNODETOCELLNEIGHBORLIST_3[0][n2] += 1;
		check = WALLNODETOCELLNEIGHBORLIST_3[0][n2];
		if(check > 6) printf("List size exceeded.. ");
		WALLNODETOCELLNEIGHBORLIST_3[check][n2] = m;
		
		if(check > maxNodeW) maxNodeW = check;
		
		n3 = WALLCELLNODE_3[2][m];
		WALLNODETOCELLNEIGHBORLIST_3[0][n3] += 1;
		check = WALLNODETOCELLNEIGHBORLIST_3[0][n3];
		if(check > 6) printf("List size exceeded.. ");
		WALLNODETOCELLNEIGHBORLIST_3[check][n3] = m;
		
		if(check > maxNodeW) maxNodeW = check;
	}
	printf("\n PCV2 wall mesh max connectivity is %d\n", maxNodeW);
	
	strcpy(name, "WallMeshDataStructure/PCV2WallNodeToTriConnectivity.csv");
	printf("Finished ouputting node to tri connectivity to %s...\n", name);
	fOUT = fopen(name, "w");
		fprintf(fOUT, "X Y Z TotalTri T1 T2 T3 T4 T5 T6\n");
		for(m=1;m<=NodeW_3;m++){
			fprintf(fOUT, "%.7e %.7e %.7e %d %d %d %d %d %d %d\n",WALLNODECOORD_3[0][m],WALLNODECOORD_3[1][m],WALLNODECOORD_3[2][m],WALLNODETOCELLNEIGHBORLIST_3[0][m],WALLNODETOCELLNEIGHBORLIST_3[1][m],WALLNODETOCELLNEIGHBORLIST_3[2][m],WALLNODETOCELLNEIGHBORLIST_3[3][m],WALLNODETOCELLNEIGHBORLIST_3[4][m],WALLNODETOCELLNEIGHBORLIST_3[5][m],WALLNODETOCELLNEIGHBORLIST_3[6][m]);
		}
	fclose(fOUT);
	printf("Finished ouputting node to tri connectivity to %s...\n", name);
	
	
	strcpy(name, "WallMeshDataStructure/PCV2WallTriNormals.csv");
	printf("Finished ouputting node to tri connectivity to %s...\n", name);
	fOUT = fopen(name, "w");
		fprintf(fOUT, "X Y Z Nx Ny Nz\n");
		for(tri=1;tri<=CellW_3;tri++){
			n1 = WALLCELLNODE_3[0][tri];
			n2 = WALLCELLNODE_3[1][tri];
			n3 = WALLCELLNODE_3[2][tri];
			
			xc = (WALLNODECOORD_3[0][n1]+WALLNODECOORD_3[0][n2]+WALLNODECOORD_3[0][n3])/3.;
			yc = (WALLNODECOORD_3[1][n1]+WALLNODECOORD_3[1][n2]+WALLNODECOORD_3[1][n3])/3.;
			zc = (WALLNODECOORD_3[2][n1]+WALLNODECOORD_3[2][n2]+WALLNODECOORD_3[2][n3])/3.;
			
			fprintf(fOUT, "%.7e %.7e %.7e %.7e %.7e %.7e\n", xc,yc,zc,WALLCELLNORMALS_3[0][tri],WALLCELLNORMALS_3[1][tri],WALLCELLNORMALS_3[2][tri]);
		}
	fclose(fOUT);
	printf("Finished ouputting node to tri connectivity to %s...\n", name);	

	/// ///////////////////////////////
	

	printf("Calculating ROI domain lumen wall areas\n");
	for(tri=1;tri<=CellW;tri++){
		n1 = WALLCELLNODE[0][tri];
		n2 = WALLCELLNODE[1][tri];
		n3 = WALLCELLNODE[2][tri];
		
		for(i=0;i<=2;i++){
			n1coord[i] = WALLNODECOORD[i][n1];
			n2coord[i] = WALLNODECOORD[i][n2];
			n3coord[i] = WALLNODECOORD[i][n3];
		}
							
		for(i=0;i<=2;i++) WALL_CENTROID[i][tri] = (n1coord[i] + n2coord[i] + n3coord[i])/3.;
							
		e1v[0] = n2coord[0] - n1coord[0];
		e1v[1] = n2coord[1] - n1coord[1];
		e1v[2] = n2coord[2] - n1coord[2];
							
		e2v[0] = n3coord[0] - n2coord[0];
		e2v[1] = n3coord[1] - n2coord[1];
		e2v[2] = n3coord[2] - n2coord[2];
							
		e3v[0] = n1coord[0] - n3coord[0];
		e3v[1] = n1coord[1] - n3coord[1];
		e3v[2] = n1coord[2] - n3coord[2];
					
		Ax = 0.5*(e1v[1]*e3v[2] - e1v[2]*e3v[1])*1.e-012;
		Ay = 0.5*(e1v[2]*e3v[0] - e1v[0]*e3v[2])*1.e-012;
		Az = 0.5*(e1v[0]*e3v[1] - e1v[1]*e3v[0])*1.e-012;
		AreaLocal = sqrtf(Ax*Ax+Ay*Ay+Az*Az);
		
		WALL_AREA[tri] = AreaLocal;													
	}
	
	printf("Calculating DA parent domain lumen wall areas\n");
	for(tri=1;tri<=CellW_1;tri++){
		n1 = WALLCELLNODE_1[0][tri];
		n2 = WALLCELLNODE_1[1][tri];
		n3 = WALLCELLNODE_1[2][tri];
		
		for(i=0;i<=2;i++){
			n1coord[i] = WALLNODECOORD_1[i][n1];
			n2coord[i] = WALLNODECOORD_1[i][n2];
			n3coord[i] = WALLNODECOORD_1[i][n3];
		}
							
		for(i=0;i<=2;i++) WALL_CENTROID_1[i][tri] = (n1coord[i] + n2coord[i] + n3coord[i])/3.;
							
		e1v[0] = n2coord[0] - n1coord[0];
		e1v[1] = n2coord[1] - n1coord[1];
		e1v[2] = n2coord[2] - n1coord[2];
							
		e2v[0] = n3coord[0] - n2coord[0];
		e2v[1] = n3coord[1] - n2coord[1];
		e2v[2] = n3coord[2] - n2coord[2];
							
		e3v[0] = n1coord[0] - n3coord[0];
		e3v[1] = n1coord[1] - n3coord[1];
		e3v[2] = n1coord[2] - n3coord[2];
					
		Ax = 0.5*(e1v[1]*e3v[2] - e1v[2]*e3v[1])*1.e-012;
		Ay = 0.5*(e1v[2]*e3v[0] - e1v[0]*e3v[2])*1.e-012;
		Az = 0.5*(e1v[0]*e3v[1] - e1v[1]*e3v[0])*1.e-012;
		AreaLocal = sqrtf(Ax*Ax+Ay*Ay+Az*Az);
		
		WALL_AREA_1[tri] = AreaLocal;													
	}
	
	printf("Calculating PCV parent1 domain lumen wall areas\n");
	for(tri=1;tri<=CellW_2;tri++){
		n1 = WALLCELLNODE_2[0][tri];
		n2 = WALLCELLNODE_2[1][tri];
		n3 = WALLCELLNODE_2[2][tri];
		
		for(i=0;i<=2;i++){
			n1coord[i] = WALLNODECOORD_2[i][n1];
			n2coord[i] = WALLNODECOORD_2[i][n2];
			n3coord[i] = WALLNODECOORD_2[i][n3];
		}
							
		for(i=0;i<=2;i++) WALL_CENTROID_2[i][tri] = (n1coord[i] + n2coord[i] + n3coord[i])/3.;
							
		e1v[0] = n2coord[0] - n1coord[0];
		e1v[1] = n2coord[1] - n1coord[1];
		e1v[2] = n2coord[2] - n1coord[2];
							
		e2v[0] = n3coord[0] - n2coord[0];
		e2v[1] = n3coord[1] - n2coord[1];
		e2v[2] = n3coord[2] - n2coord[2];
							
		e3v[0] = n1coord[0] - n3coord[0];
		e3v[1] = n1coord[1] - n3coord[1];
		e3v[2] = n1coord[2] - n3coord[2];
					
		Ax = 0.5*(e1v[1]*e3v[2] - e1v[2]*e3v[1])*1.e-012;
		Ay = 0.5*(e1v[2]*e3v[0] - e1v[0]*e3v[2])*1.e-012;
		Az = 0.5*(e1v[0]*e3v[1] - e1v[1]*e3v[0])*1.e-012;
		AreaLocal = sqrtf(Ax*Ax+Ay*Ay+Az*Az);
		
		WALL_AREA_2[tri] = AreaLocal;													
	}
	
	printf("Calculating PCV parent2 domain lumen wall areas\n");
	for(tri=1;tri<=CellW_3;tri++){
		n1 = WALLCELLNODE_3[0][tri];
		n2 = WALLCELLNODE_3[1][tri];
		n3 = WALLCELLNODE_3[2][tri];
		
		for(i=0;i<=2;i++){
			n1coord[i] = WALLNODECOORD_3[i][n1];
			n2coord[i] = WALLNODECOORD_3[i][n2];
			n3coord[i] = WALLNODECOORD_3[i][n3];
		}
							
		for(i=0;i<=2;i++) WALL_CENTROID_3[i][tri] = (n1coord[i] + n2coord[i] + n3coord[i])/3.;
							
		e1v[0] = n2coord[0] - n1coord[0];
		e1v[1] = n2coord[1] - n1coord[1];
		e1v[2] = n2coord[2] - n1coord[2];
							
		e2v[0] = n3coord[0] - n2coord[0];
		e2v[1] = n3coord[1] - n2coord[1];
		e2v[2] = n3coord[2] - n2coord[2];
							
		e3v[0] = n1coord[0] - n3coord[0];
		e3v[1] = n1coord[1] - n3coord[1];
		e3v[2] = n1coord[2] - n3coord[2];
					
		Ax = 0.5*(e1v[1]*e3v[2] - e1v[2]*e3v[1])*1.e-012;
		Ay = 0.5*(e1v[2]*e3v[0] - e1v[0]*e3v[2])*1.e-012;
		Az = 0.5*(e1v[0]*e3v[1] - e1v[1]*e3v[0])*1.e-012;
		AreaLocal = sqrtf(Ax*Ax+Ay*Ay+Az*Az);
		
		WALL_AREA_3[tri] = AreaLocal;													
	}
	
	
	
	
	
	for(i=1;i<=IM;i++){
		for(j=1;j<=JM;j++){
			for(k=3;k<=KM-2;k++){
				if(W_index[i][j][k] > 0){
					if(B_index[i][j][k] <= 0){
						mindis = 5000.;
						nodewall = W_index[i][j][k];
												
						neighmax = WALLNODETOCELLNEIGHBORLIST[0][nodewall];
						for(l=1;l<=neighmax;l++){
							tri = WALLNODETOCELLNEIGHBORLIST[l][nodewall];
							n1 = WALLCELLNODE[0][tri]; n2 = WALLCELLNODE[1][tri]; n3 = WALLCELLNODE[2][tri];
							xp1 = WALLNODECOORD[0][n1];
							yp1 = WALLNODECOORD[1][n1];
							zp1 = WALLNODECOORD[2][n1];
							distance = sqrtf((xp-xp1)*(xp-xp1) + (yp-yp1)*(yp-yp1) + (zp-zp1)*(zp-zp1));
							if(distance < mindis){
								mindis = distance;
								triwall = tri;
							}
													
							xp2 = WALLNODECOORD[0][n2];
							yp2 = WALLNODECOORD[1][n2];
							zp2 = WALLNODECOORD[2][n2];
							distance = sqrtf((xp-xp2)*(xp-xp2) + (yp-yp2)*(yp-yp2) + (zp-zp2)*(zp-zp2));
							if(distance < mindis){
								mindis = distance;
								triwall = tri;
							}
													
							xp3 = WALLNODECOORD[0][n3];
							yp3 = WALLNODECOORD[1][n3];
							zp3 = WALLNODECOORD[2][n3];
							distance = sqrtf((xp-xp3)*(xp-xp3) + (yp-yp3)*(yp-yp3) + (zp-zp3)*(zp-zp3));
							if(distance < mindis){
								mindis = distance;
								triwall = tri;
							}
												
							xp4 = WALL_CENTROID[0][tri];
							yp4 = WALL_CENTROID[1][tri];
							zp4 = WALL_CENTROID[2][tri];
							distance = sqrtf((xp-xp4)*(xp-xp4) + (yp-yp4)*(yp-yp4) + (zp-zp4)*(zp-zp4));
							if(distance < mindis){
								mindis = distance;
								triwall = tri;
							}	
						}
						
						e1v[0] = WALLCELLNORMALS[0][triwall];
						e1v[1] = WALLCELLNORMALS[1][triwall];
						e1v[2] = WALLCELLNORMALS[2][triwall];
						
						if(Wdistance[i][j][k] <= 0.7) weight1 = 0;
						else if(Wdistance[i][j][k] <= 3.) weight1 = (Wdistance[i][j][k]-0.7)/(3.-0.7);
						else weight1 = 1.;
						
						Uvoid[i][j][k] = (weight1*1.e-03)/VEL_SCALE*WALLCELLNORMALS[0][triwall];
						Vvoid[i][j][k] = (weight1*1.e-03)/VEL_SCALE*WALLCELLNORMALS[1][triwall];
						Wvoid[i][j][k] = (weight1*1.e-03)/VEL_SCALE*WALLCELLNORMALS[2][triwall];
					}
				}
			}
		}		
	}
	
	
	///open special pleural zone
	xp1 = 74.23347;
	yp1 = 241.36881;
	zp1 = 319.79235;
	for(i=1;i<=IM;i++){
		for(j=1;j<=JM;j++){
			for(k=3;k<=KM-2;k++){
				xp = X[i][j][k]; yp = Y[i][j][k]; zp = Z[i][j][k];
				distance = sqrtf((xp-xp1)*(xp-xp1) + (yp-yp1)*(yp-yp1) + (zp-zp1)*(zp-zp1));
				
				if(distance < 3.){
					if(distance > 0.001){
						e1v[0] = xp-xp1/distance;
						e1v[1] = yp-yp1/distance;
						e1v[2] = zp-zp1/distance;
					}
					else{
						e1v[0] = 0;
						e1v[1] = 0;
						e1v[2] = 1;
					}
					
					
					weight1 = (distance)/(3.);
					
						
						Uvoid[i][j][k] = (weight1*1.e-03)/e1v[0];
						Vvoid[i][j][k] = (weight1*1.e-03)/e1v[1];
						Wvoid[i][j][k] = (weight1*1.e-03)/e1v[2];
					
					
				}
				
				
				
			}
		}
	}
	
	
	///
	for(i=1;i<=IM_1;i++){
		for(j=1;j<=JM_1;j++){
			for(k=1;k<=KM_1;k++){
				if(W_index_1[i][j][k] > 0){
					if(B_index_1[i][j][k] <= 0){
						mindis = 5000.;
						nodewall = W_index_1[i][j][k];
												
						neighmax = WALLNODETOCELLNEIGHBORLIST_1[0][nodewall];
						for(l=1;l<=neighmax;l++){
							tri = WALLNODETOCELLNEIGHBORLIST_1[l][nodewall];
							n1 = WALLCELLNODE_1[0][tri]; n2 = WALLCELLNODE_1[1][tri]; n3 = WALLCELLNODE_1[2][tri];
							xp1 = WALLNODECOORD_1[0][n1];
							yp1 = WALLNODECOORD_1[1][n1];
							zp1 = WALLNODECOORD_1[2][n1];
							distance = sqrtf((xp-xp1)*(xp-xp1) + (yp-yp1)*(yp-yp1) + (zp-zp1)*(zp-zp1));
							if(distance < mindis){
								mindis = distance;
								triwall = tri;
							}
													
							xp2 = WALLNODECOORD_1[0][n2];
							yp2 = WALLNODECOORD_1[1][n2];
							zp2 = WALLNODECOORD_1[2][n2];
							distance = sqrtf((xp-xp2)*(xp-xp2) + (yp-yp2)*(yp-yp2) + (zp-zp2)*(zp-zp2));
							if(distance < mindis){
								mindis = distance;
								triwall = tri;
							}
													
							xp3 = WALLNODECOORD_1[0][n3];
							yp3 = WALLNODECOORD_1[1][n3];
							zp3 = WALLNODECOORD_1[2][n3];
							distance = sqrtf((xp-xp3)*(xp-xp3) + (yp-yp3)*(yp-yp3) + (zp-zp3)*(zp-zp3));
							if(distance < mindis){
								mindis = distance;
								triwall = tri;
							}
												
							xp4 = WALL_CENTROID_1[0][tri];
							yp4 = WALL_CENTROID_1[1][tri];
							zp4 = WALL_CENTROID_1[2][tri];
							distance = sqrtf((xp-xp4)*(xp-xp4) + (yp-yp4)*(yp-yp4) + (zp-zp4)*(zp-zp4));
							if(distance < mindis){
								mindis = distance;
								triwall = tri;
							}	
						}
						
						e1v[0] = WALLCELLNORMALS_1[0][triwall];
						e1v[1] = WALLCELLNORMALS_1[1][triwall];
						e1v[2] = WALLCELLNORMALS_1[2][triwall];
						
						if(Wdistance_1[i][j][k] <= 2.) weight1 = (Wdistance_1[i][j][k])/(2.);
						else weight1 = 1.;
						
						Uvoid_1[i][j][k] = (weight1*1.e-03)/VEL_SCALE*WALLCELLNORMALS_1[0][triwall];
						Vvoid_1[i][j][k] = (weight1*1.e-03)/VEL_SCALE*WALLCELLNORMALS_1[1][triwall];
						Wvoid_1[i][j][k] = (weight1*1.e-03)/VEL_SCALE*WALLCELLNORMALS_1[2][triwall];
					}
				}
			}
		}		
	}
	
	for(i=1;i<=IM_2;i++){
		for(j=1;j<=JM_2;j++){
			for(k=1;k<=KM_2;k++){
				if(W_index_2[i][j][k] > 0){
					if(B_index_2[i][j][k] <= 0){
						mindis = 5000.;
						nodewall = W_index_2[i][j][k];
												
						neighmax = WALLNODETOCELLNEIGHBORLIST_2[0][nodewall];
						for(l=1;l<=neighmax;l++){
							tri = WALLNODETOCELLNEIGHBORLIST_2[l][nodewall];
							n1 = WALLCELLNODE_2[0][tri]; n2 = WALLCELLNODE_2[1][tri]; n3 = WALLCELLNODE_2[2][tri];
							xp1 = WALLNODECOORD_2[0][n1];
							yp1 = WALLNODECOORD_2[1][n1];
							zp1 = WALLNODECOORD_2[2][n1];
							distance = sqrtf((xp-xp1)*(xp-xp1) + (yp-yp1)*(yp-yp1) + (zp-zp1)*(zp-zp1));
							if(distance < mindis){
								mindis = distance;
								triwall = tri;
							}
													
							xp2 = WALLNODECOORD_2[0][n2];
							yp2 = WALLNODECOORD_2[1][n2];
							zp2 = WALLNODECOORD_2[2][n2];
							distance = sqrtf((xp-xp2)*(xp-xp2) + (yp-yp2)*(yp-yp2) + (zp-zp2)*(zp-zp2));
							if(distance < mindis){
								mindis = distance;
								triwall = tri;
							}
													
							xp3 = WALLNODECOORD_2[0][n3];
							yp3 = WALLNODECOORD_2[1][n3];
							zp3 = WALLNODECOORD_2[2][n3];
							distance = sqrtf((xp-xp3)*(xp-xp3) + (yp-yp3)*(yp-yp3) + (zp-zp3)*(zp-zp3));
							if(distance < mindis){
								mindis = distance;
								triwall = tri;
							}
												
							xp4 = WALL_CENTROID_2[0][tri];
							yp4 = WALL_CENTROID_2[1][tri];
							zp4 = WALL_CENTROID_2[2][tri];
							distance = sqrtf((xp-xp4)*(xp-xp4) + (yp-yp4)*(yp-yp4) + (zp-zp4)*(zp-zp4));
							if(distance < mindis){
								mindis = distance;
								triwall = tri;
							}	
						}
						
						e1v[0] = WALLCELLNORMALS_2[0][triwall];
						e1v[1] = WALLCELLNORMALS_2[1][triwall];
						e1v[2] = WALLCELLNORMALS_2[2][triwall];
						
						if(Wdistance_2[i][j][k] <= 2.) weight1 = (Wdistance_2[i][j][k])/(2.);
						else weight1 = 1.;
						
						Uvoid_2[i][j][k] = (weight1*1.e-03)/VEL_SCALE*WALLCELLNORMALS_2[0][triwall];
						Vvoid_2[i][j][k] = (weight1*1.e-03)/VEL_SCALE*WALLCELLNORMALS_2[1][triwall];
						Wvoid_2[i][j][k] = (weight1*1.e-03)/VEL_SCALE*WALLCELLNORMALS_2[2][triwall];
					}
				}
			}
		}		
	}
	
	for(i=1;i<=IM_3;i++){
		for(j=1;j<=JM_3;j++){
			for(k=1;k<=KM_3;k++){
				if(W_index_3[i][j][k] > 0){
					if(B_index_3[i][j][k] <= 0){
						mindis = 5000.;
						nodewall = W_index_3[i][j][k];
												
						neighmax = WALLNODETOCELLNEIGHBORLIST_3[0][nodewall];
						for(l=1;l<=neighmax;l++){
							tri = WALLNODETOCELLNEIGHBORLIST_3[l][nodewall];
							n1 = WALLCELLNODE_3[0][tri]; n2 = WALLCELLNODE_3[1][tri]; n3 = WALLCELLNODE_3[2][tri];
							xp1 = WALLNODECOORD_3[0][n1];
							yp1 = WALLNODECOORD_3[1][n1];
							zp1 = WALLNODECOORD_3[2][n1];
							distance = sqrtf((xp-xp1)*(xp-xp1) + (yp-yp1)*(yp-yp1) + (zp-zp1)*(zp-zp1));
							if(distance < mindis){
								mindis = distance;
								triwall = tri;
							}
													
							xp2 = WALLNODECOORD_3[0][n2];
							yp2 = WALLNODECOORD_3[1][n2];
							zp2 = WALLNODECOORD_3[2][n2];
							distance = sqrtf((xp-xp2)*(xp-xp2) + (yp-yp2)*(yp-yp2) + (zp-zp2)*(zp-zp2));
							if(distance < mindis){
								mindis = distance;
								triwall = tri;
							}
													
							xp3 = WALLNODECOORD_3[0][n3];
							yp3 = WALLNODECOORD_3[1][n3];
							zp3 = WALLNODECOORD_3[2][n3];
							distance = sqrtf((xp-xp3)*(xp-xp3) + (yp-yp3)*(yp-yp3) + (zp-zp3)*(zp-zp3));
							if(distance < mindis){
								mindis = distance;
								triwall = tri;
							}
												
							xp4 = WALL_CENTROID_3[0][tri];
							yp4 = WALL_CENTROID_3[1][tri];
							zp4 = WALL_CENTROID_3[2][tri];
							distance = sqrtf((xp-xp4)*(xp-xp4) + (yp-yp4)*(yp-yp4) + (zp-zp4)*(zp-zp4));
							if(distance < mindis){
								mindis = distance;
								triwall = tri;
							}	
						}
						
						e1v[0] = WALLCELLNORMALS_3[0][triwall];
						e1v[1] = WALLCELLNORMALS_3[1][triwall];
						e1v[2] = WALLCELLNORMALS_3[2][triwall];
						
						if(Wdistance_3[i][j][k] <= 2.) weight1 = (Wdistance_3[i][j][k])/(2.);
						else weight1 = 1.;
						
						Uvoid_3[i][j][k] = (weight1*1.e-03)/VEL_SCALE*WALLCELLNORMALS_3[0][triwall];
						Vvoid_3[i][j][k] = (weight1*1.e-03)/VEL_SCALE*WALLCELLNORMALS_3[1][triwall];
						Wvoid_3[i][j][k] = (weight1*1.e-03)/VEL_SCALE*WALLCELLNORMALS_3[2][triwall];
					}
				}
			}
		}		
	}
	
	strcpy(name, "LBMGridDataStructure/DOMID.csv");
	printf("Ouputting fluid at iter#%d to %s...\n", iter, name);
	fOUT = fopen(name, "w");
		fprintf(fOUT, "x,y,z,I,J,K,1DarrayID\n");
		for(k=1;k<=KM;k++){
			for(j=1;j<=JM;j++){
				for(i=1;i<=IM;i++){
					xp = X[i][j][k];
					yp = Y[i][j][k];
					zp = Z[i][j][k];
					
					if(B_index[i][j][k] >= 0)fprintf(fOUT, "%.7e, %.7e, %.7e, %d, %d, %d, %d\n", xp, yp, zp, i, j, k, DOMID[i][j][k]);
				}
			}
		}		
	fclose(fOUT);
	
	strcpy(name, "LBMGridDataStructure/DOMID_1.csv");
	printf("Ouputting fluid at iter#%d to %s...\n", iter, name);
	fOUT = fopen(name, "w");
		fprintf(fOUT, "x,y,z,I,J,K,1DarrayID\n");
		for(k=1;k<=KM_1;k++){
			for(j=1;j<=JM_1;j++){
				for(i=1;i<=IM_1;i++){
					xp = X_1[i][j][k];
					yp = Y_1[i][j][k];
					zp = Z_1[i][j][k];
					
					if(B_index[i][j][k] >= 0)fprintf(fOUT, "%.7e, %.7e, %.7e, %d, %d, %d, %d\n", xp, yp, zp, i, j, k, DOMID_1[i][j][k]);
				}
			}
		}
	fclose(fOUT);
	
	strcpy(name, "LBMGridDataStructure/DOMID_2.csv");
	printf("Ouputting fluid at iter#%d to %s...\n", iter, name);
	fOUT = fopen(name, "w");
		fprintf(fOUT, "x,y,z,I,J,K,1DarrayID\n");
		for(k=1;k<=KM_2;k++){
			for(j=1;j<=JM_2;j++){
				for(i=1;i<=IM_2;i++){
					xp = X_2[i][j][k];
					yp = Y_2[i][j][k];
					zp = Z_2[i][j][k];
					
					if(B_index[i][j][k] >= 0)fprintf(fOUT, "%.7e, %.7e, %.7e, %d, %d, %d, %d\n", xp, yp, zp, i, j, k, DOMID_2[i][j][k]);
				}
			}
		}
	fclose(fOUT);
	
	strcpy(name, "LBMGridDataStructure/DOMID_3.csv");
	printf("Ouputting fluid at iter#%d to %s...\n", iter, name);
	fOUT = fopen(name, "w");
		fprintf(fOUT, "x,y,z,I,J,K,1DarrayID\n");
		for(k=1;k<=KM_3;k++){
			for(j=1;j<=JM_3;j++){
				for(i=1;i<=IM_3;i++){
					xp = X_3[i][j][k];
					yp = Y_3[i][j][k];
					zp = Z_3[i][j][k];
					
					if(B_index[i][j][k] >= 0)fprintf(fOUT, "%.7e, %.7e, %.7e, %d, %d, %d, %d\n", xp, yp, zp, i, j, k, DOMID_3[i][j][k]);
				}
			}
		}
	fclose(fOUT);	
	
	
	
	
	
	
	strcpy(name, "LBMGridDataStructure/LumenBlock");
	strcat(name,".vtk");
	printf("Ouputting fluid at iter#%d to %s...\n", iter, name);
	fOUT = fopen(name, "w");
		fprintf(fOUT, "# vtk DataFile Version 4.0\n");
		fprintf(fOUT, "vtk output\n");
		fprintf(fOUT, "ASCII\n");
		fprintf(fOUT, "DATASET STRUCTURED_GRID\n");
		fprintf(fOUT, "DIMENSIONS %d %d %d\n",IM,JM,KM);
		fprintf(fOUT, "POINTS %d double\n",IM*JM*KM);
		
		for(k=1;k<=KM;k++){
			for(j=1;j<=JM;j++){
				for(i=1;i<=IM;i++){
					fprintf(fOUT, "%f %f %f\n", X[i][j][k], Y[i][j][k], Z[i][j][k]);
				}
			}
		}
		
		fprintf(fOUT, "\nPOINT_DATA %d\n",IM*JM*KM);
		fprintf(fOUT, "FIELD FieldData 4\n");
		fprintf(fOUT, "Bindex 1 %d int\n",IM*JM*KM);
		for(k=1;k<=KM;k++){
			for(j=1;j<=JM;j++){
				for(i=1;i<=IM;i++){
					fprintf(fOUT, "%d\n", B_index[i][j][k]);
				}
			}
		}
		
		fprintf(fOUT, "\nUvoid 1 %d float\n",IM*JM*KM);
		for(k=1;k<=KM;k++){
			for(j=1;j<=JM;j++){
				for(i=1;i<=IM;i++){
					fprintf(fOUT, "%.7e\n", Uvoid[i][j][k]);
				}
			}
		}
		
		fprintf(fOUT, "\nVvoid 1 %d float\n",IM*JM*KM);
		for(k=1;k<=KM;k++){
			for(j=1;j<=JM;j++){
				for(i=1;i<=IM;i++){
					fprintf(fOUT, "%.7e\n", Vvoid[i][j][k]);
				}
			}
		}
		
		fprintf(fOUT, "\nWvoid 1 %d float\n",IM*JM*KM);
		for(k=1;k<=KM;k++){
			for(j=1;j<=JM;j++){
				for(i=1;i<=IM;i++){
					fprintf(fOUT, "%.7e\n", Wvoid[i][j][k]);
				}
			}
		}
		
	fclose(fOUT);
	printf("Finished ouputting fluid at iter#%d to %s...\n", iter, name);
	
	printf("Initializing RBCnode list members in ROI lattice nodes to zero\n");
	for(l=0;l<8500000;l++){
		for(c=0;c<=EulMAX;c++){
			EulNeigh[c][l] = 0;
		}
	}
	
	printf("Initializing RBCnode list members in DAP lattice nodes to zero\n");
	for(l=0;l<400000;l++){
		for(c=0;c<=EulMAX;c++){
			EulNeigh_1[c][l] = 0;
		}
	}
	
	printf("Initializing RBCnode list members in PCV1 lattice nodes to zero\n");
	for(l=0;l<400000;l++){
		for(c=0;c<=EulMAX;c++){
			EulNeigh_2[c][l] = 0;
		}
	}
	
	printf("Initializing RBCnode list members in PCV2 lattice nodes to zero\n");
	for(l=0;l<120000;l++){
		for(c=0;c<=EulMAX;c++){
			EulNeigh_3[c][l] = 0;
		}
	}
	
	
	printf("Initializing microfluxes in ROI lattice nodes to maxwellian equilibrium\n");
	for(l=1;l<=DOMMAX;l++){
		
		for(c=1;c<=CM;c++){
			FEQ[c][l] = D[l]*Wcoeff[c]*(1. + EQ_A*(U[l]*UC[c] + V[l]*VC[c] + W[l]*WC[c]) + EQ_B*(U[l]*UC[c] + V[l]*VC[c] + W[l]*WC[c])*(U[l]*UC[c] + V[l]*VC[c] + W[l]*WC[c]) + EQ_C*(U[l]*U[l] + V[l]*V[l] + W[l]*W[l]));
			F[c][l] = FEQ[c][l];
			F_OLD[c][l] = F[c][l];
		}
		
	}
	
	printf("Initializing microfluxes in DAP lattice nodes to maxwellian equilibrium\n");
	for(l=1;l<=DOMMAX_1;l++){
		
		for(c=1;c<=CM;c++){
			FEQ_1[c][l] = D_1[l]*Wcoeff[c]*(1. + EQ_A*(U_1[l]*UC[c] + V_1[l]*VC[c] + W_1[l]*WC[c]) + EQ_B*(U_1[l]*UC[c] + V_1[l]*VC[c] + W_1[l]*WC[c])*(U_1[l]*UC[c] + V_1[l]*VC[c] + W_1[l]*WC[c]) + EQ_C*(U_1[l]*U_1[l] + V_1[l]*V_1[l] + W_1[l]*W_1[l]));
			F_1[c][l] = FEQ_1[c][l];
			F_OLD_1[c][l] = F_1[c][l];
		}
		
	}
	
	printf("Initializing microfluxes in PCV1 lattice nodes to maxwellian equilibrium\n");
	for(l=1;l<=DOMMAX_2;l++){
		
		for(c=1;c<=CM;c++){
			FEQ_2[c][l] = D_2[l]*Wcoeff[c]*(1. + EQ_A*(U_2[l]*UC[c] + V_2[l]*VC[c] + W_2[l]*WC[c]) + EQ_B*(U_2[l]*UC[c] + V_2[l]*VC[c] + W_2[l]*WC[c])*(U_2[l]*UC[c] + V_2[l]*VC[c] + W_2[l]*WC[c]) + EQ_C*(U_2[l]*U_2[l] + V_2[l]*V_2[l] + W_2[l]*W_2[l]));
			F_2[c][l] = FEQ_2[c][l];
			F_OLD_2[c][l] = F_2[c][l];
		}
		
	}
	
	printf("Initializing microfluxes in PCV2 lattice nodes to maxwellian equilibrium\n");
	for(l=1;l<=DOMMAX_3;l++){
		
		for(c=1;c<=CM;c++){
			FEQ_3[c][l] = D_3[l]*Wcoeff[c]*(1. + EQ_A*(U_3[l]*UC[c] + V_3[l]*VC[c] + W_3[l]*WC[c]) + EQ_B*(U_3[l]*UC[c] + V_3[l]*VC[c] + W_3[l]*WC[c])*(U_3[l]*UC[c] + V_3[l]*VC[c] + W_3[l]*WC[c]) + EQ_C*(U_3[l]*U_3[l] + V_3[l]*V_3[l] + W_3[l]*W_3[l]));
			F_3[c][l] = FEQ_3[c][l];
			F_OLD_3[c][l] = F_3[c][l];
		}
		
	}
	
	printf("Initializing RBC parameters to default\n");
	for(caps=1;caps<=CapsM;caps++){
		CellStatus[caps] = -5000;
		REGRESSFLAG[caps] = 0;
		NucMembUpdate[caps] = 0;
		for(node=1;node<=NodeM;node++){
			NODE_WALLNEIGH[caps][node] = 0;
			MembIn_NeighFlag[caps][node] = 0;
			MembOut_NeighFlag[caps][node] = 0;
								
			for(i=0;i<=2;i++){
				NODE_WALLREPUL[i][caps][node] = 0;
				NODECOORD[i][caps][node] = -5000.;
				F_Total[i][caps][node] = 0.;
				NODE_VEL[i][caps][node] = 0.;
			}
		}
								
		for(node=1;node<=NodeM_N;node++){
			NUC_NeighFlag[caps][node] = 0;
			for(i=0;i<=2;i++){
				NODECOORD_N[i][caps][node] = -5000.;
				F_Total_N[i][caps][node] = 0.;
				NODE_VEL_N[i][caps][node] = 0.;
			}
		}
	}
	
	/// The cell membrane data input file was too large for Github, so it was split into two files ///////
	CellCount = 0;
	strcpy(name, "InputData/AllCellMembStart.aa");
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
	
	strcpy(name, "InputData/AllCellMembStart.aa");
		printf("Opening %s \n", name);
		fIN = fopen(name, "r");
		fgets(bulk, 180, fIN);
			for(m=1;m<=NodeM*CellCount;m++){
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
				
				
				NODE_WALLNEIGH[cellindex][nodenumber] = nwtri;
				NODE_WALLREPUL[0][cellindex][nodenumber] = nwrepulx;
				NODE_WALLREPUL[1][cellindex][nodenumber] = nwrepuly;
				NODE_WALLREPUL[2][cellindex][nodenumber] = nwrepulz;
				
				MembIn_Neigh[cellindex][nodenumber] = MNin;
				MembIn_NeighFlag[cellindex][nodenumber] = MNinflag;
				
				MembOut_Neigh[cellindex][nodenumber] = MNout;
				MembOut_NeighFlag[cellindex][nodenumber] = MNoutflag;
				
				
				if(m != 1 && (m%NodeM == 0)) nodenumber = 0;
			}
		fclose(fIN);
		
	CellCount = 0;
	strcpy(name, "InputData/AllCellMembStart.ab");
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
		
	
	CellCount = CellCount/504;
		
	printf("The file %s has %d number of cells\n ", name, CellCount); 
	
	strcpy(name, "InputData/AllCellMembStart.ab");
	printf("Opening %s \n", name);
	fIN = fopen(name, "r");
		for(m=1;m<=NodeM*CellCount;m++){
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
			
			
			NODE_WALLNEIGH[cellindex][nodenumber] = nwtri;
			NODE_WALLREPUL[0][cellindex][nodenumber] = nwrepulx;
			NODE_WALLREPUL[1][cellindex][nodenumber] = nwrepuly;
			NODE_WALLREPUL[2][cellindex][nodenumber] = nwrepulz;
			
			MembIn_Neigh[cellindex][nodenumber] = MNin;
			MembIn_NeighFlag[cellindex][nodenumber] = MNinflag;
			
			MembOut_Neigh[cellindex][nodenumber] = MNout;
			MembOut_NeighFlag[cellindex][nodenumber] = MNoutflag;
			
				
			if(m != 1 && (m%NodeM == 0)) nodenumber = 0;
		}
	fclose(fIN);
	
	/// //////
	
	strcpy(name, "InputData/AllCellNucStart.csv");
	printf("Opening %s \n", name);
		fIN = fopen(name, "r");
		fgets(bulk, 180, fIN);
			for(m=1;m<=NodeM_N*CellCount;m++){	
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
				
				fscanf(fIN,"%f",&fx);//vx
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&fy);//vy
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%f",&fz);//vz
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%d",&MNin);
				fscanf(fIN,"%c",&junkchar);
				fscanf(fIN,"%d",&MNinflag);
				
				
				NODECOORD_N[0][cellindex][nodenumber] = xp;
				NODECOORD_N[1][cellindex][nodenumber] = yp;
				NODECOORD_N[2][cellindex][nodenumber] = zp;
				
				NODE_VEL_N[0][cellindex][nodenumber] = vx;
				NODE_VEL_N[1][cellindex][nodenumber] = vy;
				NODE_VEL_N[2][cellindex][nodenumber] = vz;
				
				F_Total_N[0][cellindex][nodenumber] = fx;
				F_Total_N[1][cellindex][nodenumber] = fy;
				F_Total_N[2][cellindex][nodenumber] = fz;
				
				
				NUC_Neigh[cellindex][nodenumber] = MNin;
				NUC_NeighFlag[cellindex][nodenumber] = MNinflag;
				
				if(m != 1 && (m%NodeM_N == 0)) nodenumber = 0;
			}
		fclose(fIN);
	
	//Remove bad red blood cells
	caps = 772;
	CellStatus[caps] = -5000;
	
	caps = 1065;
	CellStatus[caps] = -5000;	
	
	caps = 362;
	CellStatus[caps] = -5000;
	
	caps = 1162;
	CellStatus[caps] = -5000;
	
	for(caps=1;caps<=CapsM;caps++){
		if(CellStatus[caps] == -5000){//Let's reposition all nodes to null domain
			REGRESSFLAG[caps] = 0;
			NucMembUpdate[caps] = 0;
			for(node=1;node<=NodeM;node++){
						
				NODE_WALLNEIGH[caps][node] = 0;
				MembIn_NeighFlag[caps][node] = 0;
				MembOut_NeighFlag[caps][node] = 0;
				
				for(i=0;i<=2;i++){
					NODE_WALLREPUL[i][caps][node] = 0;
					NODECOORD[i][caps][node] = -5000.;
					F_Total[i][caps][node] = 0.;
					NODE_VEL[i][caps][node] = 0.;
				}
			}
			
			for(node=1;node<=NodeM_N;node++){
				NUC_NeighFlag[caps][node] = 0;
				for(i=0;i<=2;i++){
					NODECOORD_N[i][caps][node] = -5000.;
					F_Total_N[i][caps][node] = 0.;
					NODE_VEL_N[i][caps][node] = 0.;
				}
			}
		}
	}
	
	
	printf("Starting do loop..\n");
	
	
	plasmaloop = 0;
	iter = iterrestart-1;
	
	do{
		simfailure = 0;
		iter += 1;
		
		RBCSOLVERFLAG = 0;
					
		if(iter%SOLVERSTAGGER == 0) RBCSOLVERFLAG = 1;//means update forces
		
		relaxfactor = 1.;
		relaxfactor2 = 1.;
		
		///reset eulerian list storing neighboring RBC node ids to zero 
		if(RBCSOLVERFLAG == 1){
			#pragma omp parallel private(i, j, k, l, c)
			{
				#pragma omp for schedule(static,5000) nowait
				for(l=1;l<=DOMMAX;l++){
					
					for(c=0;c<=EulMAX;c++){
						EulNeigh[c][l] = 0;
					}
				}
			}
			
			#pragma omp parallel private(i, j, k, l, c)
			{
				#pragma omp for schedule(static,5000) nowait
				for(l=1;l<=DOMMAX_1;l++){
					
					for(c=0;c<=EulMAX;c++){
						EulNeigh_1[c][l] = 0;
					}
				}
			}
			
			#pragma omp parallel private(i, j, k, l, c)
			{
				#pragma omp for schedule(static,5000) nowait
				for(l=1;l<=DOMMAX_2;l++){
					
					for(c=0;c<=EulMAX;c++){
						EulNeigh_2[c][l] = 0;
					}
				}
			}
			
			#pragma omp parallel private(i, j, k, l, c)
			{
				#pragma omp for schedule(static,2000) nowait
				for(l=1;l<=DOMMAX_3;l++){
					
					for(c=0;c<=EulMAX;c++){
						EulNeigh_3[c][l] = 0;
					}
				}
			}
		}
		
		if(RBCSOLVERFLAG == 1){
			///CGSM_NODAL_FORCE_ZERO
			#pragma omp parallel private(caps, node)
			{
				#pragma omp for schedule(static,50) nowait
				for(caps=1;caps<=CapsM;caps++){
					if(CellStatus[caps] > -5000){	
						CAPS_REMOVE[caps] = 0;
						for(i=0;i<=2;i++){
							CAPS_MAX[i][caps] = -5000.;
							CAPS_MIN[i][caps] = 5000.;
						}
						for(node=1;node<=NodeM;node++){
							NODE_AREA[caps][node] = 0.;
							NODE_WALLNEIGH[caps][node] = 0;
							for(i=0;i<=2;i++){
								NODE_WALLREPUL[i][caps][node] = 0;
								NODE_NORMAL[i][caps][node] = 0.;
								F_VOLUME[i][caps][node] = 0.;
								F_AGG[i][caps][node] = 0.;
								F_Total[i][caps][node] = 0.;
								NODE_VEL[i][caps][node] = 0;
								
								if(NODECOORD[i][caps][node] > CAPS_MAX[i][caps]) CAPS_MAX[i][caps] = NODECOORD[i][caps][node];
								if(NODECOORD[i][caps][node] < CAPS_MIN[i][caps]) CAPS_MIN[i][caps] = NODECOORD[i][caps][node];
							}
							
						}
						for(node=1;node<=NodeM_N;node++){
							for(i=0;i<=2;i++){
								NODE_WALLREPUL_N[i][caps][node] = 0;
								F_Total_N[i][caps][node] = 0.;
								NODE_VEL_N[i][caps][node] = 0;
								
								if(NODECOORD_N[i][caps][node] > CAPS_MAX[i][caps]) CAPS_MAX[i][caps] = NODECOORD_N[i][caps][node];
								if(NODECOORD_N[i][caps][node] < CAPS_MIN[i][caps]) CAPS_MIN[i][caps] = NODECOORD_N[i][caps][node];
							}
						}
					}
				}
			}
		}
		
		if(RBCSOLVERFLAG == 1){
			///CGSM_AREA
			#pragma omp parallel private(v,l1,l2,l3,caps,tri,n1,n2,n3,i,n1coord,n2coord,n3coord,e1v,e2v,e3v,Ax,Ay,Az,AreaLocal,tnormal,mindis,xp,yp,zp,j,k,i1,i2,j1,j2,k1,k2,distance,ii,jj,kk,l,neightrimax,check,kkk,node)
			{
				#pragma omp for schedule(static,50) nowait
				for(caps=1;caps<=CapsM;caps++){
					if(CellStatus[caps] > -5000){
						for(tri=1;tri<=TriM;tri++){
							
							
							n1 = TRINODE1[tri];
							n2 = TRINODE2[tri];
							n3 = TRINODE3[tri];
							
							
								for(i=0;i<=2;i++){
									n1coord[i] = NODECOORD[i][caps][n1];
									n2coord[i] = NODECOORD[i][caps][n2];
									n3coord[i] = NODECOORD[i][caps][n3];
									
								}
							
							if(CellStatus[caps] > 3){
								if(CAPS_MAX[2][caps] - CAPS_MIN[2][caps] > 30.){//cell is crossing the periodic boundary
									if(caps <= 188){
										
											if(NODECOORD[2][caps][n1] < 0.5*(Z0_1+ZM_1)) n1coord[2] += (ZM_1-Z0_1);
											
											if(NODECOORD[2][caps][n2] < 0.5*(Z0_1+ZM_1)) n2coord[2] += (ZM_1-Z0_1);
											
											if(NODECOORD[2][caps][n3] < 0.5*(Z0_1+ZM_1)) n3coord[2] += (ZM_1-Z0_1);
										
										
									}
									else if(caps >= 189 && caps <= 329){
										
											if(NODECOORD[2][caps][n1] > 0.5*(Z0_2+ZM_2)) n1coord[2] += -(ZM_2-Z0_2);
											
											if(NODECOORD[2][caps][n2] > 0.5*(Z0_2+ZM_2)) n2coord[2] += -(ZM_2-Z0_2);
											
											if(NODECOORD[2][caps][n3] > 0.5*(Z0_2+ZM_2)) n3coord[2] += -(ZM_2-Z0_2);
										
									}
									else if(caps >= 330 && caps <= 344){
										
											if(NODECOORD[2][caps][n1] > 0.5*(Z0_3+ZM_3)) n1coord[2] += -(ZM_3-Z0_3);
											
											if(NODECOORD[2][caps][n2] > 0.5*(Z0_3+ZM_3)) n2coord[2] += -(ZM_3-Z0_3);
											
											if(NODECOORD[2][caps][n3] > 0.5*(Z0_3+ZM_3)) n3coord[2] += -(ZM_3-Z0_3);
										
									}
								}
							}
							
							for(i=0;i<=2;i++) TRI_CENTROID[i][caps][tri] = (n1coord[i] + n2coord[i] + n3coord[i])/3.;
							
							e1v[0] = n2coord[0] - n1coord[0];
							e1v[1] = n2coord[1] - n1coord[1];
							e1v[2] = n2coord[2] - n1coord[2];
							
							e2v[0] = n3coord[0] - n2coord[0];
							e2v[1] = n3coord[1] - n2coord[1];
							e2v[2] = n3coord[2] - n2coord[2];
							
							e3v[0] = n1coord[0] - n3coord[0];
							e3v[1] = n1coord[1] - n3coord[1];
							e3v[2] = n1coord[2] - n3coord[2];
							
							
							Ax = 0.5*(e1v[1]*e3v[2] - e1v[2]*e3v[1])*1.e-012;
							Ay = 0.5*(e1v[2]*e3v[0] - e1v[0]*e3v[2])*1.e-012;
							Az = 0.5*(e1v[0]*e3v[1] - e1v[1]*e3v[0])*1.e-012;
							AreaLocal = sqrtf(Ax*Ax+Ay*Ay+Az*Az);
							
							TRI_AREA[caps][tri] = AreaLocal;
							
							
							AreaExpanRatio[caps][tri] = AreaLocal/TRI_AREA_eq[tri];
							
							
							tnormal[0] = -Ax/AreaLocal;
							tnormal[1] = -Ay/AreaLocal;
							tnormal[2] = -Az/AreaLocal;
							
							TRI_NORMAL[0][caps][tri] = tnormal[0];
							TRI_NORMAL[1][caps][tri] = tnormal[1];
							TRI_NORMAL[2][caps][tri] = tnormal[2];
							
							NODE_AREA[caps][n1] += AreaLocal/3.;
							NODE_AREA[caps][n2] += AreaLocal/3.;
							NODE_AREA[caps][n3] += AreaLocal/3.; 
							
							NODE_NORMAL[0][caps][n1] += tnormal[0]*AreaLocal/3.;
							NODE_NORMAL[1][caps][n1] += tnormal[1]*AreaLocal/3.;
							NODE_NORMAL[2][caps][n1] += tnormal[2]*AreaLocal/3.;
							
							NODE_NORMAL[0][caps][n2] += tnormal[0]*AreaLocal/3.;
							NODE_NORMAL[1][caps][n2] += tnormal[1]*AreaLocal/3.;
							NODE_NORMAL[2][caps][n2] += tnormal[2]*AreaLocal/3.;
							
							NODE_NORMAL[0][caps][n3] += tnormal[0]*AreaLocal/3.;
							NODE_NORMAL[1][caps][n3] += tnormal[1]*AreaLocal/3.;
							NODE_NORMAL[2][caps][n3] += tnormal[2]*AreaLocal/3.;
							
							///update eulerian grid
							mindis = 5000.;
							xp = TRI_CENTROID[0][caps][tri]; yp = TRI_CENTROID[1][caps][tri]; zp = TRI_CENTROID[2][caps][tri];
							
							///ROI domain
							if(caps > 344){
								i = (int)((xp-X0)/DX + 1.01);    j = (int)((yp-Y0)/DX + 1.01);    k = (int)((zp-Z0)/DX + 1.01);
								i1 = i-4; i2 = i+4;		j1 = j-4; j2 = j+4;		k1 = k-4; k2 = k+4;
								
								if(i1 < 1) i1 = 1; if(i2 > IM) i2 = IM;
								if(j1 < 1) j1 = 1; if(j2 > JM) j2 = JM;
								if(k1 < 3) k1 = 3; if(k2 > KM-3) k2 = KM-3;
								
								check = 0;
								
								if(zp < Z0 || zp > ZM) check = 1;
								
								
								if(check == 0){
									for(i=i1;i<=i2;i++){
										for(j=j1;j<=j2;j++){
											for(k=k1;k<=k2;k++){
												if(B_index[i][j][k] >= 0){
													distance = sqrtf((X[i][j][k]-xp)*(X[i][j][k]-xp) + (Y[i][j][k]-yp)*(Y[i][j][k]-yp) + (Z[i][j][k]-zp)*(Z[i][j][k]-zp));
													
													if(distance < mindis){
														mindis = distance;
														ii = i; 	jj = j; 	kk =k;
													}
												}
											}
										}
									}
									
									if(mindis < 1.){
										l = DOMID[ii][jj][kk];
										#pragma omp atomic capture
										{ 
											neightrimax = EulNeigh[0][l];
											EulNeigh[0][l] += 1; 
											
										}
										
										neightrimax += 1;
										if(neightrimax >= EulMAX) printf("Eulerian cell list has exceeded size!\n");
										
										///This should be unique to each thread in the OMP process since neightrimax values will never overlap across the multithreads
										EulNeigh[neightrimax][l] = tri + caps*TriM;
									}
								}
							}
							else if(caps <= 188){
								if(CellStatus[caps] != 2 && zp >= ZM_1) zp = zp - (ZM_1-Z0_1);//keep point within eulerian domain
								 
								i = (int)((xp-X0_1)/DX + 1.01);    j = (int)((yp-Y0_1)/DX + 1.01);    k = (int)((zp-Z0_1)/DX + 1.01);
								i1 = i-4; i2 = i+4;		j1 = j-4; j2 = j+4;		k1 = k-4; k2 = k+4;
								
								if(i1 < 1) i1 = 1; if(i2 > IM_1) i2 = IM_1;
								if(j1 < 1) j1 = 1; if(j2 > JM_1) j2 = JM_1;
								
								for(i=i1;i<=i2;i++){
									for(j=j1;j<=j2;j++){
										for(k=k1;k<=k2;k++){
											check = 0;
											kkk = k;
											if(k > KM_1){
												kkk = k-KM_1+1;
												check = 1;//you need to forward translate the neighbor location
											}
											else if(k < 1){
												kkk = k+KM_1-1;
												check = 2;//you need to backward translate the neighbor location
											}
											
											if(B_index_1[i][j][kkk] >= 0){
												if(check == 0) distance = sqrtf((X_1[i][j][kkk]-xp)*(X_1[i][j][kkk]-xp) + (Y_1[i][j][kkk]-yp)*(Y_1[i][j][kkk]-yp) + (Z_1[i][j][kkk]-zp)*(Z_1[i][j][kkk]-zp));
												else if(check == 1) distance = sqrtf((X_1[i][j][kkk]-xp)*(X_1[i][j][kkk]-xp) + (Y_1[i][j][kkk]-yp)*(Y_1[i][j][kkk]-yp) + (Z_1[i][j][kkk]+(ZM_1-Z0_1)-zp)*(Z_1[i][j][kkk]+(ZM_1-Z0_1)-zp));
												else if(check == 2) distance = sqrtf((X_1[i][j][kkk]-xp)*(X_1[i][j][kkk]-xp) + (Y_1[i][j][kkk]-yp)*(Y_1[i][j][kkk]-yp) + (Z_1[i][j][kkk]-(ZM_1-Z0_1)-zp)*(Z_1[i][j][kkk]-(ZM_1-Z0_1)-zp));
												
												if(distance < mindis){
													mindis = distance;
													ii = i; 	jj = j; 	kk =kkk;
												}
											}
										}
									}
								}
								
								l = DOMID_1[ii][jj][kk];
								#pragma omp atomic capture
								{ neightrimax = EulNeigh_1[0][l]; EulNeigh_1[0][l] += 1; }
								neightrimax += 1;
								
								if(neightrimax >= EulMAX) printf("Eulerian cell list in DAparent has exceeded size!\n");		
														
								EulNeigh_1[neightrimax][l] = tri + caps*TriM;
							}
							else if(caps >= 189 && caps <= 329){
								if(CellStatus[caps] != 2 && zp <= Z0_2) zp = zp + (ZM_2-Z0_2);
								
								i = (int)((xp-X0_2)/DX + 1.01);    j = (int)((yp-Y0_2)/DX + 1.01);    k = (int)((zp-Z0_2)/DX + 1.01);
								i1 = i-4; i2 = i+4;		j1 = j-4; j2 = j+4;		k1 = k-4; k2 = k+4;
								
								if(i1 < 1) i1 = 1; if(i2 > IM_2) i2 = IM_2;
								if(j1 < 1) j1 = 1; if(j2 > JM_2) j2 = JM_2;
								for(i=i1;i<=i2;i++){
									for(j=j1;j<=j2;j++){
										for(k=k1;k<=k2;k++){
											check = 0;
											kkk = k;
											if(k > KM_2){
												kkk = k-KM_2+1;
												check = 1;//you need to forward translate the neighbor location
											}
											else if(k < 1){
												kkk = k+KM_2-1;
												check = 2;//you need to backward translate the neighbor location
											}
											
											if(B_index_2[i][j][kkk] >= 0){
												if(check == 0) distance = sqrtf((X_2[i][j][kkk]-xp)*(X_2[i][j][kkk]-xp) + (Y_2[i][j][kkk]-yp)*(Y_2[i][j][kkk]-yp) + (Z_2[i][j][kkk]-zp)*(Z_2[i][j][kkk]-zp));
												else if(check == 1) distance = sqrtf((X_2[i][j][kkk]-xp)*(X_2[i][j][kkk]-xp) + (Y_2[i][j][kkk]-yp)*(Y_2[i][j][kkk]-yp) + (Z_2[i][j][kkk]+(ZM_2-Z0_2)-zp)*(Z_2[i][j][kkk]+(ZM_2-Z0_2)-zp));
												else if(check == 2) distance = sqrtf((X_2[i][j][kkk]-xp)*(X_2[i][j][kkk]-xp) + (Y_2[i][j][kkk]-yp)*(Y_2[i][j][kkk]-yp) + (Z_2[i][j][kkk]-(ZM_2-Z0_2)-zp)*(Z_2[i][j][kkk]-(ZM_2-Z0_2)-zp));
												
												if(distance < mindis){
													mindis = distance;
													ii = i; 	jj = j; 	kk =kkk;
												}
											}
										}
									}
								}
								
								l = DOMID_2[ii][jj][kk];
								#pragma omp atomic capture
								{ neightrimax = EulNeigh_2[0][l]; EulNeigh_2[0][l] += 1; }
								neightrimax += 1;
								
								if(neightrimax >= EulMAX) printf("Eulerian cell list in PCV1parent has exceeded size!\n");
								
								EulNeigh_2[neightrimax][l] = tri + caps*TriM;
							}
							else if(caps >= 330 && caps <= 344){
								if(CellStatus[caps] != 2 && zp <= Z0_3) zp = zp + (ZM_3-Z0_3);
								
								i = (int)((xp-X0_3)/DX + 1.01);    j = (int)((yp-Y0_3)/DX + 1.01);    k = (int)((zp-Z0_3)/DX + 1.01);
								i1 = i-4; i2 = i+4;		j1 = j-4; j2 = j+4;		k1 = k-4; k2 = k+4;
								
								if(i1 < 1) i1 = 1; if(i2 > IM_3) i2 = IM_3;
								if(j1 < 1) j1 = 1; if(j2 > JM_3) j2 = JM_3;
								
								for(i=i1;i<=i2;i++){
									for(j=j1;j<=j2;j++){
										for(k=k1;k<=k2;k++){
											check = 0;
											kkk = k;
											if(k > KM_3){
												kkk = k-KM_3+1;
												check = 1;//you need to forward translate the neighbor location
											}
											else if(k < 1){
												kkk = k+KM_3-1;
												check = 2;//you need to backward translate the neighbor location
											}
											
											if(B_index_3[i][j][kkk] >= 0){
												if(check == 0) distance = sqrtf((X_3[i][j][kkk]-xp)*(X_3[i][j][kkk]-xp) + (Y_3[i][j][kkk]-yp)*(Y_3[i][j][kkk]-yp) + (Z_3[i][j][kkk]-zp)*(Z_3[i][j][kkk]-zp));
												else if(check == 1) distance = sqrtf((X_3[i][j][kkk]-xp)*(X_3[i][j][kkk]-xp) + (Y_3[i][j][kkk]-yp)*(Y_3[i][j][kkk]-yp) + (Z_3[i][j][kkk]+(ZM_3-Z0_3)-zp)*(Z_3[i][j][kkk]+(ZM_3-Z0_3)-zp));
												else if(check == 2) distance = sqrtf((X_3[i][j][kkk]-xp)*(X_3[i][j][kkk]-xp) + (Y_3[i][j][kkk]-yp)*(Y_3[i][j][kkk]-yp) + (Z_3[i][j][kkk]-(ZM_3-Z0_3)-zp)*(Z_3[i][j][kkk]-(ZM_3-Z0_3)-zp));
												
												if(distance < mindis){
													mindis = distance;
													ii = i; 	jj = j; 	kk =kkk;
												}
											}
										}
									}
								}
								
								l = DOMID_3[ii][jj][kk];
								#pragma omp atomic capture
								{ neightrimax = EulNeigh_3[0][l]; EulNeigh_3[0][l] += 1; }
								neightrimax += 1;
								
								if(neightrimax >= EulMAX) printf("Eulerian cell list in PCV2parent has exceeded size!\n");
																
								EulNeigh_3[neightrimax][l] = tri + caps*TriM;
								
							}
							
							//apply steric repulsion forces that prevent triangle elements from deforming into slivers (high aspect ratio triangles that can corrupt area and bend angle calculation)
							for(v=0;v<=2;v++){
								e1v[v] = n1coord[v] - TRI_CENTROID[v][caps][tri];
								e2v[v] = n2coord[v] - TRI_CENTROID[v][caps][tri];
								e3v[v] = n3coord[v] - TRI_CENTROID[v][caps][tri];
							}
							
							l1 = sqrtf(e1v[0]*e1v[0] + e1v[1]*e1v[1] + e1v[2]*e1v[2]);
							l2 = sqrtf(e2v[0]*e2v[0] + e2v[1]*e2v[1] + e2v[2]*e2v[2]);
							l3 = sqrtf(e3v[0]*e3v[0] + e3v[1]*e3v[1] + e3v[2]*e3v[2]);
							
							for(v=0;v<=2;v++){
								e1v[v] = e1v[v]/l1;
								e2v[v] = e2v[v]/l2;
								e3v[v] = e3v[v]/l3;
							}
							
							if(l1 < 0.05){
								if(l1 > 0.001){
									for(v=0;v<=2;v++) F_Total[v][caps][n1] += relaxfactor2*1.e-14/l1*e1v[v];
								}
								else{
									for(v=0;v<=2;v++) F_Total[v][caps][n1] += relaxfactor2*1.e-11*e1v[v];
								}
							}
							
							if(l2 < 0.05){
								if(l2 > 0.001){
									for(v=0;v<=2;v++) F_Total[v][caps][n2] += relaxfactor2*1.e-14/l2*e2v[v];
								}
								else{
									for(v=0;v<=2;v++) F_Total[v][caps][n2] += relaxfactor2*1.e-11*e2v[v];
								}
							}
							
							if(l3 < 0.05){
								if(l3 > 0.001){
									for(v=0;v<=2;v++) F_Total[v][caps][n3] += relaxfactor2*1.e-14/l3*e3v[v];
								}
								else{
									for(v=0;v<=2;v++) F_Total[v][caps][n3] += relaxfactor2*1.e-11*e3v[v];
								}
							}
						}
					
					
						for(node=1; node<=NodeM; node++){
							NODE_NORMAL[0][caps][node] = NODE_NORMAL[0][caps][node]/NODE_AREA[caps][node];
							NODE_NORMAL[1][caps][node] = NODE_NORMAL[1][caps][node]/NODE_AREA[caps][node];
							NODE_NORMAL[2][caps][node] = NODE_NORMAL[2][caps][node]/NODE_AREA[caps][node];
						}
					}
				}
			}
		}
		
		///LBM Code /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
		///LBM_STREAM() & LBM_COLLISION() combined;  /////////////////////////////////////////////////////////////////////////////////////////////////////////
		#pragma omp parallel private(bfeq, i, j, k, ii, jj, kk, l, m, tauLBM, cc, c, UU, VV, WW, DD)
		{
			#pragma omp for schedule(static,5000) nowait
			for(l=1;l<=DOMMAX;l++){
				i = DOMLIST[l][0];
				j = DOMLIST[l][1];
				k = DOMLIST[l][2];
				
				if(B_index[i][j][k] > 0 && B_index[i][j][k] < 10){
					tauLBM=0.5+MU_SCALE*MU*B_index[i][j][k];
					
					///stream
					for(c=1;c<CM;c++){
						if(c%2 == 0) cc = c-1;//cc is the opposing lattice direction to c
						else cc = c+1;
						ii = i+Xlattice[0][cc]; jj = j+Xlattice[1][cc]; kk = k+Xlattice[2][cc];
						
						m = DOMID[ii][jj][kk];
						F[c][l]=F_OLD[c][m];
					}
					c=19; cc=19;
						ii = i+Xlattice[0][cc]; jj = j+Xlattice[1][cc]; kk = k+Xlattice[2][cc];
						
						m = DOMID[ii][jj][kk];
						F[c][l]=F_OLD[c][m];
					///macro update
					UU=0.0;	VV=0.0; WW=0.0;	DD=0.0;
					for(c=1;c<=CM;c++){
						DD=DD+F[c][l]; UU=UU+F[c][l]*UC[c]; VV=VV+F[c][l]*VC[c]; WW=WW+F[c][l]*WC[c];
					}
					
					D[l]=DD;
					U[l]=(UU+DT*BODYFORCE[0][l]/2.)/DD;
					V[l]=(VV+DT*BODYFORCE[1][l]/2.)/DD;
					W[l]=(WW+DT*BODYFORCE[2][l]/2.)/DD;
					
					///collision
					for(c=1;c<=CM;c++){
						FEQ[c][l] = D[l]*Wcoeff[c]*(1. + EQ_A*(U[l]*UC[c] + V[l]*VC[c] + W[l]*WC[c]) + EQ_B*(U[l]*UC[c] + V[l]*VC[c] + W[l]*WC[c])*(U[l]*UC[c] + V[l]*VC[c] + W[l]*WC[c]) + EQ_C*(U[l]*U[l] + V[l]*V[l] + W[l]*W[l]));
						bfeq = (1.-0.5/tauLBM)*Wcoeff[c]*(((UC[c]-U[l])*BODYFORCE[0][l]+(VC[c]-V[l])*BODYFORCE[1][l]+(WC[c]-W[l])*BODYFORCE[2][l])/CS/CS+(U[l]*UC[c] + V[l]*VC[c] + W[l]*WC[c])*(UC[c]*BODYFORCE[0][l] + VC[c]*BODYFORCE[1][l] + WC[c]*BODYFORCE[2][l])/CS/CS/CS/CS);
						F[c][l]=F[c][l]*(1.0-1.0/tauLBM)/*this is Fneq*/+FEQ[c][l]/tauLBM+DT*bfeq;
					}
				}
			}
		}
		
		
		#pragma omp parallel private(bfeq, i, j, k, ii, jj, kk, l, m, tauLBM, cc, c, UU, VV, WW, DD)
		{
			#pragma omp for schedule(static,5000) nowait
			for(l=1;l<=DOMMAX_1;l++){
				i = DOMLIST_1[l][0];
				j = DOMLIST_1[l][1];
				k = DOMLIST_1[l][2];
				
				if(B_index_1[i][j][k] > 0 && B_index_1[i][j][k] < 10){
					tauLBM=0.5+MU_SCALE*MU*B_index_1[i][j][k];
					
					///stream
					for(c=1;c<CM;c++){
						if(c%2 == 0) cc = c-1;//cc is the opposing lattice direction to c
						else cc = c+1;
						ii = i+Xlattice[0][cc]; jj = j+Xlattice[1][cc]; kk = k+Xlattice[2][cc];
						
						m = DOMID_1[ii][jj][kk];
						F_1[c][l]=F_OLD_1[c][m];
					}
					c=19; cc=19;
						ii = i+Xlattice[0][cc]; jj = j+Xlattice[1][cc]; kk = k+Xlattice[2][cc];
						
						m = DOMID_1[ii][jj][kk];
						F_1[c][l]=F_OLD_1[c][m];
					///macro update
					UU=0.0;	VV=0.0; WW=0.0;	DD=0.0;
					for(c=1;c<=CM;c++){
						DD=DD+F_1[c][l]; UU=UU+F_1[c][l]*UC[c]; VV=VV+F_1[c][l]*VC[c]; WW=WW+F_1[c][l]*WC[c];
					}
					
					D_1[l]=DD;
					U_1[l]=(UU+DT*BODYFORCE_1[0][l]/2.)/DD;
					V_1[l]=(VV+DT*BODYFORCE_1[1][l]/2.)/DD;
					W_1[l]=(WW+DT*BODYFORCE_1[2][l]/2.)/DD;
					
					///collision
					for(c=1;c<=CM;c++){
						FEQ_1[c][l] = D_1[l]*Wcoeff[c]*(1. + EQ_A*(U_1[l]*UC[c] + V_1[l]*VC[c] + W_1[l]*WC[c]) + EQ_B*(U_1[l]*UC[c] + V_1[l]*VC[c] + W_1[l]*WC[c])*(U_1[l]*UC[c] + V_1[l]*VC[c] + W_1[l]*WC[c]) + EQ_C*(U_1[l]*U_1[l] + V_1[l]*V_1[l] + W_1[l]*W_1[l]));
						bfeq = (1.-0.5/tauLBM)*Wcoeff[c]*(((UC[c]-U_1[l])*BODYFORCE_1[0][l]+(VC[c]-V_1[l])*BODYFORCE_1[1][l]+(WC[c]-W_1[l])*BODYFORCE_1[2][l])/CS/CS+(U_1[l]*UC[c] + V_1[l]*VC[c] + W_1[l]*WC[c])*(UC[c]*BODYFORCE_1[0][l] + VC[c]*BODYFORCE_1[1][l] + WC[c]*BODYFORCE_1[2][l])/CS/CS/CS/CS);
						F_1[c][l]=F_1[c][l]*(1.0-1.0/tauLBM)/*this is Fneq*/+FEQ_1[c][l]/tauLBM+DT*bfeq;
					}
				}
			}
		}
		
		
		#pragma omp parallel private(bfeq, i, j, k, ii, jj, kk, l, m, tauLBM, cc, c, UU, VV, WW, DD)
		{
			#pragma omp for schedule(static,5000) nowait
			for(l=1;l<=DOMMAX_2;l++){
				i = DOMLIST_2[l][0];
				j = DOMLIST_2[l][1];
				k = DOMLIST_2[l][2];
				
				if(B_index_2[i][j][k] > 0 && B_index_2[i][j][k] < 10){
					tauLBM=0.5+MU_SCALE*MU*B_index_2[i][j][k];
					
					///stream
					for(c=1;c<CM;c++){
						if(c%2 == 0) cc = c-1;//cc is the opposing lattice direction to c
						else cc = c+1;
						ii = i+Xlattice[0][cc]; jj = j+Xlattice[1][cc]; kk = k+Xlattice[2][cc];
						
						m = DOMID_2[ii][jj][kk];
						F_2[c][l]=F_OLD_2[c][m];
					}
					c=19; cc=19;
						ii = i+Xlattice[0][cc]; jj = j+Xlattice[1][cc]; kk = k+Xlattice[2][cc];
						
						m = DOMID_2[ii][jj][kk];
						F_2[c][l]=F_OLD_2[c][m];
					///macro update
					UU=0.0;	VV=0.0; WW=0.0;	DD=0.0;
					for(c=1;c<=CM;c++){
						DD=DD+F_2[c][l]; UU=UU+F_2[c][l]*UC[c]; VV=VV+F_2[c][l]*VC[c]; WW=WW+F_2[c][l]*WC[c];
					}
					
					D_2[l]=DD;
					U_2[l]=(UU+DT*BODYFORCE_2[0][l]/2.)/DD;
					V_2[l]=(VV+DT*BODYFORCE_2[1][l]/2.)/DD;
					W_2[l]=(WW+DT*BODYFORCE_2[2][l]/2.)/DD;
					
					///collision
					for(c=1;c<=CM;c++){
						FEQ_2[c][l] = D_2[l]*Wcoeff[c]*(1. + EQ_A*(U_2[l]*UC[c] + V_2[l]*VC[c] + W_2[l]*WC[c]) + EQ_B*(U_2[l]*UC[c] + V_2[l]*VC[c] + W_2[l]*WC[c])*(U_2[l]*UC[c] + V_2[l]*VC[c] + W_2[l]*WC[c]) + EQ_C*(U_2[l]*U_2[l] + V_2[l]*V_2[l] + W_2[l]*W_2[l]));
						bfeq = (1.-0.5/tauLBM)*Wcoeff[c]*(((UC[c]-U_2[l])*BODYFORCE_2[0][l]+(VC[c]-V_2[l])*BODYFORCE_2[1][l]+(WC[c]-W_2[l])*BODYFORCE_2[2][l])/CS/CS+(U_2[l]*UC[c] + V_2[l]*VC[c] + W_2[l]*WC[c])*(UC[c]*BODYFORCE_2[0][l] + VC[c]*BODYFORCE_2[1][l] + WC[c]*BODYFORCE_2[2][l])/CS/CS/CS/CS);
						F_2[c][l]=F_2[c][l]*(1.0-1.0/tauLBM)/*this is Fneq*/+FEQ_2[c][l]/tauLBM+DT*bfeq;
					}
				}
			}
		}
		
		
		#pragma omp parallel private(bfeq, i, j, k, ii, jj, kk, l, m, tauLBM, cc, c, UU, VV, WW, DD)
		{
			#pragma omp for schedule(static,2000) nowait
			for(l=1;l<=DOMMAX_3;l++){
				i = DOMLIST_3[l][0];
				j = DOMLIST_3[l][1];
				k = DOMLIST_3[l][2];
				
				if(B_index_3[i][j][k] > 0 && B_index_3[i][j][k] < 10){
					tauLBM=0.5+MU_SCALE*MU*B_index_3[i][j][k];
					
					///stream
					for(c=1;c<CM;c++){
						if(c%2 == 0) cc = c-1;//cc is the opposing lattice direction to c
						else cc = c+1;
						ii = i+Xlattice[0][cc]; jj = j+Xlattice[1][cc]; kk = k+Xlattice[2][cc];
						
						m = DOMID_3[ii][jj][kk];
						F_3[c][l]=F_OLD_3[c][m];
					}
					c=19; cc=19;
						ii = i+Xlattice[0][cc]; jj = j+Xlattice[1][cc]; kk = k+Xlattice[2][cc];
						
						m = DOMID_3[ii][jj][kk];
						F_3[c][l]=F_OLD_3[c][m];
					///macro update
					UU=0.0;	VV=0.0; WW=0.0;	DD=0.0;
					for(c=1;c<=CM;c++){
						DD=DD+F_3[c][l]; UU=UU+F_3[c][l]*UC[c]; VV=VV+F_3[c][l]*VC[c]; WW=WW+F_3[c][l]*WC[c];
					}
					
					D_3[l]=DD;
					U_3[l]=(UU+DT*BODYFORCE_3[0][l]/2.)/DD;
					V_3[l]=(VV+DT*BODYFORCE_3[1][l]/2.)/DD;
					W_3[l]=(WW+DT*BODYFORCE_3[2][l]/2.)/DD;
					
					///collision
					for(c=1;c<=CM;c++){
						FEQ_3[c][l] = D_3[l]*Wcoeff[c]*(1. + EQ_A*(U_3[l]*UC[c] + V_3[l]*VC[c] + W_3[l]*WC[c]) + EQ_B*(U_3[l]*UC[c] + V_3[l]*VC[c] + W_3[l]*WC[c])*(U_3[l]*UC[c] + V_3[l]*VC[c] + W_3[l]*WC[c]) + EQ_C*(U_3[l]*U_3[l] + V_3[l]*V_3[l] + W_3[l]*W_3[l]));
						bfeq = (1.-0.5/tauLBM)*Wcoeff[c]*(((UC[c]-U_3[l])*BODYFORCE_3[0][l]+(VC[c]-V_3[l])*BODYFORCE_3[1][l]+(WC[c]-W_3[l])*BODYFORCE_3[2][l])/CS/CS+(U_3[l]*UC[c] + V_3[l]*VC[c] + W_3[l]*WC[c])*(UC[c]*BODYFORCE_3[0][l] + VC[c]*BODYFORCE_3[1][l] + WC[c]*BODYFORCE_3[2][l])/CS/CS/CS/CS);
						F_3[c][l]=F_3[c][l]*(1.0-1.0/tauLBM)/*this is Fneq*/+FEQ_3[c][l]/tauLBM+DT*bfeq;
					}
				}
			}
		}
		/// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		
		/// LBM_BC(); //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		#pragma omp parallel private(c, i, j, k, iii, jjj, kkk, l, m, n, simtime, cycletime, taucycle, tauLBM, pressure, D_bc, U_bc, V_bc, W_bc, feq_bc, bfeq)
		{
			#pragma omp for schedule(static,2000)  nowait
			for(l=1;l<=BCMAX;l++){
				iii = BCNLIST[l][0];
				jjj = BCNLIST[l][1];
				kkk = BCNLIST[l][2];
				i = BCLIST[l][0];
				j = BCLIST[l][1];
				k = BCLIST[l][2];
				
				simtime = iter*T_SCALE;
				tauLBM=0.5+MU_SCALE*MU;	
				
				m = DOMID[iii][jjj][kkk];
				n = DOMID[i][j][k];
				
				if(BCLIST[l][3] == 0 || BCLIST[l][3] == 4998){	
					
					D_bc = D[m];
					U_bc = 0;
					V_bc = 0;
					W_bc = 0;
				}
				else{
					
					if(BCLIST[l][3] == 1001){//pressure boundary on DA anterior
						cycletime = fmod(simtime,0.36);
						taucycle = cycletime/0.36;
						if(taucycle <= 0.23){
							pressure = 120.+54.+82.*sinf(3.14159265359*taucycle/0.46)*sinf(3.14159265359*taucycle/0.46);
						}
						else{
							pressure = 120.+54.+82.*(1.+sinf(2*3.14159265359*(taucycle+0.54)/1.54 - 3.14159265359/2.))/2.;
						}
						
					}
					else if(BCLIST[l][3] == 1002){//pressure boundary on DA posterior
						cycletime = fmod(simtime-0.06*0.36,0.36);
						if(cycletime < 0) cycletime += 0.36;
						taucycle = cycletime/0.36;
						if(taucycle <= 0.23){
							pressure = 120.+48.+68.*sinf(3.14159265359*taucycle/0.46)*sinf(3.14159265359*taucycle/0.46);
						}
						else{
							pressure = 120.+48.+68.*(1.+sinf(2*3.14159265359*(taucycle+0.54)/1.54 - 3.14159265359/2.))/2.;
						}
						
					}
					else if(BCLIST[l][3] == 1003){//pressure boundary on PCV anterior
						cycletime = fmod(simtime-0.2*0.36,0.36);
						if(cycletime < 0) cycletime += 0.36;
						taucycle = cycletime/0.36;
						if(taucycle <= 0.3){
							pressure = 120.+11.+23.*sinf(3.14159265359*taucycle/0.6);
						}
						else{
							pressure = 120.+11.+23.*(1.+sinf(2*3.14159265359*(taucycle+0.4)/1.4 - 3.14159265359/2.))/2.;
						}
					}
					else if(BCLIST[l][3] == 1004){//pressure boundary on CVP posterior main 
						cycletime = fmod(simtime-0.2*0.36,0.36);
						if(cycletime < 0) cycletime += 0.36;
						taucycle = cycletime/0.36;
						if(taucycle <= 0.3){
							pressure = 120.+18.+27.*sinf(3.14159265359*taucycle/0.6);
						}
						else{
							pressure = 120.+18.+27.*(1.+sinf(2*3.14159265359*(taucycle+0.4)/1.4 - 3.14159265359/2.))/2.;
						}
					}
					else if(BCLIST[l][3] == 1005){//pressure boundary on CVP posterior sub branch
						cycletime = fmod(simtime-0.2*0.36,0.36);
						if(cycletime < 0) cycletime += 0.36;
						taucycle = cycletime/0.36;
						if(taucycle <= 0.3){
							pressure = 120.+23.+27.*sinf(3.14159265359*taucycle/0.6);
						}
						else{
							pressure = 120.+23.+27.*(1.+sinf(2*3.14159265359*(taucycle+0.4)/1.4 - 3.14159265359/2.))/2.;
						}
					}
					
					else if(BCLIST[l][3] == 1006){//pressure boundary on DLAV anterior
						cycletime = fmod(simtime-0.2*0.36,0.36);
						if(cycletime < 0) cycletime += 0.36;
						taucycle = cycletime/0.36;
						if(taucycle <= 0.23){
							pressure = 120.+29.+66.*sinf(3.14159265359*taucycle/0.46)*sinf(3.14159265359*taucycle/0.46);
						}
						else{
							pressure = 120.+29.+66.*(1.+sinf(2*3.14159265359*(taucycle+0.54)/1.54 - 3.14159265359/2.))/2.;
						}
					}
					else if(BCLIST[l][3] == 1008){//pressure boundary on DLAV anterior
						cycletime = fmod(simtime-0.2*0.36,0.36);
						if(cycletime < 0) cycletime += 0.36;
						taucycle = cycletime/0.36;
						if(taucycle <= 0.23){
							pressure = 120.+29.+66.*sinf(3.14159265359*taucycle/0.46)*sinf(3.14159265359*taucycle/0.46);
						}
						else{
							pressure = 120.+29.+66.*(1.+sinf(2*3.14159265359*(taucycle+0.54)/1.54 - 3.14159265359/2.))/2.;
						}
					}
					else if(BCLIST[l][3] == 1007){//pressure boundary on DLAV posterior
						cycletime = fmod(simtime-0.2*0.36,0.36);
						if(cycletime < 0) cycletime += 0.36;
						taucycle = cycletime/0.36;
						if(taucycle <= 0.23){
							pressure = 120.+28.5+60.*sinf(3.14159265359*taucycle/0.46)*sinf(3.14159265359*taucycle/0.46);
						}
						else{
							pressure = 120.+28.5+60.*(1.+sinf(2*3.14159265359*(taucycle+0.54)/1.54 - 3.14159265359/2.))/2.;
						}
					}
					else if(BCLIST[l][3] == 1009){//pressure boundary on DLAV posterior
						cycletime = fmod(simtime-0.2*0.36,0.36);
						if(cycletime < 0) cycletime += 0.36;
						taucycle = cycletime/0.36;
						if(taucycle <= 0.23){
							pressure = 120.+28.5+60.*sinf(3.14159265359*taucycle/0.46)*sinf(3.14159265359*taucycle/0.46);
						}
						else{
							
							pressure = 120.+28.5+60.*(1.+sinf(2*3.14159265359*(taucycle+0.54)/1.54 - 3.14159265359/2.))/2.;
						}
					}
					
					D_bc = pressure/LBM_PRESSURE_SCALE + 1.;
					U_bc = U[m];
					V_bc = V[m];
					W_bc = W[m];
					
					
				}
					
				for(c=1;c<=CM;c++){
					feq_bc = D_bc*Wcoeff[c]*(1. + EQ_A*(U_bc*UC[c] + V_bc*VC[c] + W_bc*WC[c]) + EQ_B*(U_bc*UC[c] + V_bc*VC[c] + W_bc*WC[c])*(U_bc*UC[c] + V_bc*VC[c] + W_bc*WC[c]) + EQ_C*(U_bc*U_bc + V_bc*V_bc + W_bc*W_bc));
					F[c][n] = feq_bc;
				}	
			}
		}
		
		#pragma omp parallel private(v, c, i, j, k, ii, jj, kk, iii, jjj, kkk, l, m, n, o, simtime, cycletime, taucycle, tauLBM, pressure, D_bc, U_bc, V_bc, W_bc, feq_bc, bfeq, f_bc)
		{
			#pragma omp for schedule(static,1000)  nowait
			for(l=1;l<=BCMAX_1;l++){
				iii = BCNLIST_1[l][0];
				jjj = BCNLIST_1[l][1];
				kkk = BCNLIST_1[l][2];
				
				ii = vBCNLIST_1[l][0];
				jj = vBCNLIST_1[l][1];
				kk = vBCNLIST_1[l][2];
				
				i = BCLIST_1[l][0];
				j = BCLIST_1[l][1];
				k = BCLIST_1[l][2];
				
				v = DOMID[iii][jjj][kkk];
				m = DOMID_1[iii][jjj][kkk];
				n = DOMID_1[i][j][k];
				o = DOMID[ii][jj][kk];
				
				tauLBM=0.5+MU_SCALE*MU;	

				if(BCLIST_1[l][3] == 0){//wall bc
					D_bc = D_1[m];
					U_bc = 0;
					V_bc = 0;
					W_bc = 0;
				}
				else if(BCLIST_1[l][3] == 21){//upstream bc
					D_bc = D_1[m];
					U_bc = U[o];
					V_bc = V[o];
					W_bc = W[o];
					
				}
				else{//downstream bc: take from ROI domain
					D_bc = D[v];
					U_bc = U[v];
					V_bc = V[v];
					W_bc = W[v];
					
				}
			
				for(c=1;c<=CM;c++){
					if(BCLIST_1[l][3] == 0 || BCLIST_1[l][3] == 21) f_bc = F_OLD_1[c][m];
					else f_bc = F_OLD[c][v];
					
					feq_bc = D_bc*Wcoeff[c]*(1. + EQ_A*(U_bc*UC[c] + V_bc*VC[c] + W_bc*WC[c]) + EQ_B*(U_bc*UC[c] + V_bc*VC[c] + W_bc*WC[c])*(U_bc*UC[c] + V_bc*VC[c] + W_bc*WC[c]) + EQ_C*(U_bc*U_bc + V_bc*V_bc + W_bc*W_bc));
					F_1[c][n] = feq_bc;
				}	
			}
		}
		
		
		#pragma omp parallel private(v, c, i, j, k, ii, jj, kk, iii, jjj, kkk, l, m, n, o, simtime, cycletime, taucycle, tauLBM, pressure, D_bc, U_bc, V_bc, W_bc, feq_bc, f_bc)
		{
			#pragma omp for schedule(static,1000)  nowait
			for(l=1;l<=BCMAX_2;l++){
				iii = BCNLIST_2[l][0];
				jjj = BCNLIST_2[l][1];
				kkk = BCNLIST_2[l][2];
				
				ii = vBCNLIST_2[l][0];
				jj = vBCNLIST_2[l][1];
				kk = vBCNLIST_2[l][2];
				
				i = BCLIST_2[l][0];
				j = BCLIST_2[l][1];
				k = BCLIST_2[l][2];
				
				v = DOMID[iii][jjj][kkk];
				m = DOMID_2[iii][jjj][kkk];
				n = DOMID_2[i][j][k];
				o = DOMID[ii][jj][kk];
				
				tauLBM=0.5+MU_SCALE*MU;

				if(BCLIST_2[l][3] == 0){
					D_bc = D_2[m];
					U_bc = 0;
					V_bc = 0;
					W_bc = 0;
				}
				else if(BCLIST_2[l][3] == 22){
					D_bc = D_2[m];
					U_bc = U[o];
					V_bc = V[o];
					W_bc = W[o];
				}
				else{
					D_bc = D[v];
					U_bc = U[v];
					V_bc = V[v];
					W_bc = W[v];	
				}
			
				for(c=1;c<=CM;c++){
					if(BCLIST_2[l][3] == 0 || BCLIST_2[l][3] == 22) f_bc = F_OLD_2[c][m];
					else f_bc = F_OLD[c][v];
					
					feq_bc = D_bc*Wcoeff[c]*(1. + EQ_A*(U_bc*UC[c] + V_bc*VC[c] + W_bc*WC[c]) + EQ_B*(U_bc*UC[c] + V_bc*VC[c] + W_bc*WC[c])*(U_bc*UC[c] + V_bc*VC[c] + W_bc*WC[c]) + EQ_C*(U_bc*U_bc + V_bc*V_bc + W_bc*W_bc));
					F_2[c][n] = feq_bc;
				}	
			}
		}
		
		
		#pragma omp parallel private(v, c, i, j, k, ii, jj, kk, iii, jjj, kkk, l, m, n, o, simtime, cycletime, taucycle, tauLBM, pressure, D_bc, U_bc, V_bc, W_bc, feq_bc, f_bc)
		{
			#pragma omp for schedule(static,400)  nowait
			for(l=1;l<=BCMAX_3;l++){
				iii = BCNLIST_3[l][0];
				jjj = BCNLIST_3[l][1];
				kkk = BCNLIST_3[l][2];
				
				ii = vBCNLIST_3[l][0];
				jj = vBCNLIST_3[l][1];
				kk = vBCNLIST_3[l][2];
				
				i = BCLIST_3[l][0];
				j = BCLIST_3[l][1];
				k = BCLIST_3[l][2];
				
				v = DOMID[iii][jjj][kkk];
				m = DOMID_3[iii][jjj][kkk];
				n = DOMID_3[i][j][k];
				o = DOMID[ii][jj][kk];
				
				tauLBM=0.5+MU_SCALE*MU;

				if(BCLIST_3[l][3] == 0){
					D_bc = D_3[m];
					U_bc = 0;
					V_bc = 0;
					W_bc = 0;
				}
				else if(BCLIST_3[l][3] == 23){
					D_bc = D_3[m];
					U_bc = U[o];
					V_bc = V[o];
					W_bc = W[o];
				}
				else{
					D_bc = D[v];
					U_bc = U[v];
					V_bc = V[v];
					W_bc = W[v];	
				}
			
				for(c=1;c<=CM;c++){
					if(BCLIST_3[l][3] == 0 || BCLIST_3[l][3] == 23) f_bc = F_OLD_3[c][m];
					else f_bc = F_OLD[c][v];
					
					feq_bc = D_bc*Wcoeff[c]*(1. + EQ_A*(U_bc*UC[c] + V_bc*VC[c] + W_bc*WC[c]) + EQ_B*(U_bc*UC[c] + V_bc*VC[c] + W_bc*WC[c])*(U_bc*UC[c] + V_bc*VC[c] + W_bc*WC[c]) + EQ_C*(U_bc*U_bc + V_bc*V_bc + W_bc*W_bc));
					F_3[c][n] = feq_bc;
				}	
			}
		}
		
		/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		
		///LBM_FLUID & LBM_COPY(); ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		#pragma omp parallel private(l, c, UU, VV, WW, DD)
		{
			#pragma omp for schedule(static,5000) nowait
			for(l=1;l<=DOMMAX;l++){			
				///macro update
				UU=0.0;	VV=0.0; WW=0.0;	DD=0.0;
				for(c=1;c<=CM;c++){
					DD=DD+F[c][l]; UU=UU+F[c][l]*UC[c]; VV=VV+F[c][l]*VC[c]; WW=WW+F[c][l]*WC[c];
				}
				
				D[l]=DD;
				U[l]=UU/DD;
				V[l]=VV/DD;
				W[l]=WW/DD;
													
				//Perform FEQ update
				for(c=1;c<=CM;c++){
					FEQ[c][l] = D[l]*Wcoeff[c]*(1. + EQ_A*(U[l]*UC[c] + V[l]*VC[c] + W[l]*WC[c]) + EQ_B*(U[l]*UC[c] + V[l]*VC[c] + W[l]*WC[c])*(U[l]*UC[c] + V[l]*VC[c] + W[l]*WC[c]) + EQ_C*(U[l]*U[l] + V[l]*V[l] + W[l]*W[l]));
					F_OLD[c][l] = F[c][l];
				}
								
				if(RBCSOLVERFLAG == 1){
					for(c=0;c<=2;c++) BODYFORCE[c][l] = 0;
				}
			}
		}
		
		
		#pragma omp parallel private(l, c, UU, VV, WW, DD)
		{
			#pragma omp for schedule(static,5000) nowait
			for(l=1;l<=DOMMAX_1;l++){			
				///macro update
				UU=0.0;	VV=0.0; WW=0.0;	DD=0.0;
				for(c=1;c<=CM;c++){
					DD=DD+F_1[c][l]; UU=UU+F_1[c][l]*UC[c]; VV=VV+F_1[c][l]*VC[c]; WW=WW+F_1[c][l]*WC[c];
				}
				
				D_1[l]=DD;
				U_1[l]=UU/DD;
				V_1[l]=VV/DD;
				W_1[l]=WW/DD;
													
				//Perform FEQ update
				for(c=1;c<=CM;c++){
					FEQ_1[c][l] = D_1[l]*Wcoeff[c]*(1. + EQ_A*(U_1[l]*UC[c] + V_1[l]*VC[c] + W_1[l]*WC[c]) + EQ_B*(U_1[l]*UC[c] + V_1[l]*VC[c] + W_1[l]*WC[c])*(U_1[l]*UC[c] + V_1[l]*VC[c] + W_1[l]*WC[c]) + EQ_C*(U_1[l]*U_1[l] + V_1[l]*V_1[l] + W_1[l]*W_1[l]));
					F_OLD_1[c][l] = F_1[c][l];
				}
								
				if(RBCSOLVERFLAG == 1){
					for(c=0;c<=2;c++) BODYFORCE_1[c][l] = 0;
				}
			}
		}
		
		
		#pragma omp parallel private(l, c, UU, VV, WW, DD)
		{
			#pragma omp for schedule(static,5000) nowait
			for(l=1;l<=DOMMAX_2;l++){			
				///macro update
				UU=0.0;	VV=0.0; WW=0.0;	DD=0.0;
				for(c=1;c<=CM;c++){
					DD=DD+F_2[c][l]; UU=UU+F_2[c][l]*UC[c]; VV=VV+F_2[c][l]*VC[c]; WW=WW+F_2[c][l]*WC[c];
				}
				
				D_2[l]=DD;
				U_2[l]=UU/DD;
				V_2[l]=VV/DD;
				W_2[l]=WW/DD;
											
				//Perform FEQ update
				for(c=1;c<=CM;c++){
					FEQ_2[c][l] = D_2[l]*Wcoeff[c]*(1. + EQ_A*(U_2[l]*UC[c] + V_2[l]*VC[c] + W_2[l]*WC[c]) + EQ_B*(U_2[l]*UC[c] + V_2[l]*VC[c] + W_2[l]*WC[c])*(U_2[l]*UC[c] + V_2[l]*VC[c] + W_2[l]*WC[c]) + EQ_C*(U_2[l]*U_2[l] + V_2[l]*V_2[l] + W_2[l]*W_2[l]));
					F_OLD_2[c][l] = F_2[c][l];
				}
				
				if(RBCSOLVERFLAG == 1){				
					for(c=0;c<=2;c++) BODYFORCE_2[c][l] = 0;
				}
			}
		}
		
		
		#pragma omp parallel private(l, c, UU, VV, WW, DD)
		{
			#pragma omp for schedule(static,2000) nowait
			for(l=1;l<=DOMMAX_3;l++){			
				///macro update
				UU=0.0;	VV=0.0; WW=0.0;	DD=0.0;
				for(c=1;c<=CM;c++){
					DD=DD+F_3[c][l]; UU=UU+F_3[c][l]*UC[c]; VV=VV+F_3[c][l]*VC[c]; WW=WW+F_3[c][l]*WC[c];
				}
				
				D_3[l]=DD;
				U_3[l]=UU/DD;
				V_3[l]=VV/DD;
				W_3[l]=WW/DD;
										
				//Perform FEQ update
				for(c=1;c<=CM;c++){
					FEQ_3[c][l] = D_3[l]*Wcoeff[c]*(1. + EQ_A*(U_3[l]*UC[c] + V_3[l]*VC[c] + W_3[l]*WC[c]) + EQ_B*(U_3[l]*UC[c] + V_3[l]*VC[c] + W_3[l]*WC[c])*(U_3[l]*UC[c] + V_3[l]*VC[c] + W_3[l]*WC[c]) + EQ_C*(U_3[l]*U_3[l] + V_3[l]*V_3[l] + W_3[l]*W_3[l]));
					F_OLD_3[c][l] = F_3[c][l];
				}
								
				if(RBCSOLVERFLAG == 1){
					for(c=0;c<=2;c++) BODYFORCE_3[c][l] = 0;
				}
			}
		}
		
		///LBM end ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		///RBC model start  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		///Cell CENTROID and Area calculation
		
		if(RBCSOLVERFLAG == 1){
			#pragma omp parallel private(caps, i, node, tri, xp, yp, zp,nzmax, nzmin)
			{
				#pragma omp for schedule(static,50)  nowait
				for(caps=1;caps<=CapsM;caps++){
					if(CellStatus[caps] > -5000){
						CAPS_AREA[caps] = 0.;
						for(i=0;i<=2;i++){
							CAPS_CENTROID[i][caps] = 0.;
						}
						
						nzmax = CAPS_MAX[2][caps];
						nzmin = CAPS_MIN[2][caps];
						
						for(node=1;node<=NodeM;node++){
							xp = NODECOORD[0][caps][node];
							yp = NODECOORD[1][caps][node];
							zp = NODECOORD[2][caps][node];
							
							if(CellStatus[caps] > 3){
								if((nzmax - nzmin) > 30.){//this cell is crossing periodic BC
									if(caps <= 188){
										
										if(zp < 0.5*(Z0_1+ZM_1)) zp += (ZM_1-Z0_1);
										
									}
									else if(caps >= 189 && caps <= 329){
										
										if(zp > 0.5*(Z0_2+ZM_2)) zp += -(ZM_2-Z0_2);
										
									}
									else if(caps >= 330 && caps <= 344){
										
										if(zp > 0.5*(Z0_3+ZM_3)) zp += -(ZM_3-Z0_3);
										
									}
								}
							}
							
							CAPS_CENTROID[0][caps] += xp/(float)NodeM;
							CAPS_CENTROID[1][caps] += yp/(float)NodeM;
							CAPS_CENTROID[2][caps] += zp/(float)NodeM;
							
						} 
						
						if(caps <= 188)
						{
							if((nzmax - nzmin) < 30.){//this cell is not crossing periodic BC
								if(CAPS_CENTROID[2][caps] > (Z0_1+0.25*(ZM_1-Z0_1)) && CAPS_CENTROID[2][caps] < (Z0_1+0.75*(ZM_1-Z0_1))) REGRESSFLAG[caps] = 0;//okay to clone when crossing cloning line
							}
						}
						else if(caps >= 189 && caps <= 329)
						{
							if((nzmax - nzmin) < 30.){//this cell is not crossing periodic BC
								if(CAPS_CENTROID[2][caps] > (Z0_2+0.25*(ZM_2-Z0_2)) && CAPS_CENTROID[2][caps] < (Z0_2+0.75*(ZM_2-Z0_2))) REGRESSFLAG[caps] = 0;//okay to clone when crossing cloning line
							}
						}
						else if(caps >= 330 && caps <= 344)
						{
							if((nzmax - nzmin) < 30.){//this cell is not crossing periodic BC
								if(CAPS_CENTROID[2][caps] > (Z0_3+0.25*(ZM_3-Z0_3)) && CAPS_CENTROID[2][caps] < (Z0_3+0.75*(ZM_3-Z0_3))) REGRESSFLAG[caps] = 0;//okay to clone when crossing cloning line
							}
						}
						else{
							if(CAPS_CENTROID[2][caps] > (Z0+0.25*(ZM-Z0)) && CAPS_CENTROID[2][caps] < (Z0+0.75*(ZM-Z0))) REGRESSFLAG[caps] = 0;//okay for removal when exiting roi
							if(CAPS_CENTROID[1][caps] < 100.) REGRESSFLAG[caps] = 0;//okay for removal when exiting roi
						}
						
						for(tri=1;tri<=TriM;tri++){
							CAPS_AREA[caps] += TRI_AREA[caps][tri];
							AREA_UPDATE[caps][tri] = 0;
						}
						
						if(CellStatus[caps] > 3000) F_AREA_GLOBAL[caps] = 0;
						else F_AREA_GLOBAL[caps] = -ka*relaxfactor*(CAPS_AREA[caps] - RBC_SURF)/RBC_SURF;
					}
				}
			}
			
		}
		
		///CGSM_LOCAL_FORCES(); - refer to Fedosov"s thesis https://www.dam.brown.edu/dpd/lib/exe/fetch.php/wiki:phd_thesis:fedosov_thesis.pdf and Maung Ye's thesis http://scholarbank.nus.edu.sg/handle/10635/125238 for the formulations in the following section 
		
		if(RBCSOLVERFLAG == 1){
			#pragma omp parallel private(nzmin, nzmax, a1check, a2check, check1, check2, check3, check4, f_a1, f_a2, f_a3, f_aa1, f_aa2, f_aa3, f_b1, f_b2, f_b3, f_b4, f_v1, f_v2, f_v3, check, caps, edge, n1, n2, n3, n4, tri1, tri2, f_n1, f_n2, f_n3, f_n4, i, n1coord, n2coord, n3coord, n4coord, dx1_12, dy1_12, dz1_12, edge1_length, e1v, e2v, e_ratio, wlc, comp, tri1_area, tri2_area, tri1_n, tri2_n, area_diff, centroid_diff, costheta, angle_check, sintheta, theta, sintheta0, costheta0, betabend, dp21, dp24, dp32, dp41, dp13, cpdt1, cpdt2, cpdt3, cpdt4, cpdt5, cpdt6, cpdt7, cpdt8, cpdt9, cpdt10, tri1_localcentroid, tri2_localcentroid, cpdt1v, cpdt5v, cpdt9v, cpdt4v, cpdt8v, cpdt12v, b11, b12, b22, alpha1, alpha2, vel_12, viscf1, viscf2)
			{
				#pragma omp for schedule(static,50)  nowait
				for(caps=1;caps<=CapsM;caps++){
					if(CellStatus[caps] > -4998){
						nzmax = CAPS_MAX[2][caps];
						nzmin = CAPS_MIN[2][caps];
						
						for(edge=1;edge<=EdgeM;edge++){
							check = 0;
							
							n1 = EDGENODE_LIST[edge][1];
							n2 = EDGENODE_LIST[edge][2];
							n3 = EDGENODE_LIST[edge][3];
							n4 = EDGENODE_LIST[edge][4];
							tri1 = EDGENODE_LIST[edge][5];
							tri2 = EDGENODE_LIST[edge][6];
							
							for(i=0;i<=2;i++){
								f_n1[i] = 0.; f_n2[i] = 0.; f_n3[i] = 0.; f_n4[i] = 0.;
								f_v1[i] = 0.; f_v2[i] = 0.; f_v3[i] = 0.;
								f_a1[i] = 0.; f_a2[i] = 0.; f_a3[i] = 0.;
								f_aa1[i] = 0.; f_aa2[i] = 0.; f_aa3[i] = 0.;
								f_b1[i] = 0.; f_b2[i] = 0.; f_b3[i] = 0.; f_b4[i] = 0.;
							}
							
															
							n1coord[0] = NODECOORD[0][caps][n1];
							n1coord[1] = NODECOORD[1][caps][n1];
							n1coord[2] = NODECOORD[2][caps][n1];
							
							n2coord[0] = NODECOORD[0][caps][n2];
							n2coord[1] = NODECOORD[1][caps][n2];
							n2coord[2] = NODECOORD[2][caps][n2];
							
							n3coord[0] = NODECOORD[0][caps][n3];
							n3coord[1] = NODECOORD[1][caps][n3];
							n3coord[2] = NODECOORD[2][caps][n3];
							
							n4coord[0] = NODECOORD[0][caps][n4];
							n4coord[1] = NODECOORD[1][caps][n4];
							n4coord[2] = NODECOORD[2][caps][n4];
							
							if(CellStatus[caps] > 3000){//CA and CV exit candidates
								if(n1coord[2] > Z0_2 || n2coord[2] > Z0_2 || n3coord[2] > Z0_2 || n4coord[2] > Z0_2 || n1coord[2] < ZM_1 || n2coord[2] < ZM_1 || n3coord[2] < ZM_1 || n4coord[2] < ZM_1) check = 1;
							}
							
							if(check == 0){
								if(CellStatus[caps] > 3){
									if((nzmax - nzmin) > 30.){//this cell is crossing periodic BC
										if(caps <= 188){
											
											if(n1coord[2] < 0.5*(Z0_1+ZM_1)) n1coord[2] += (ZM_1-Z0_1);
											if(n2coord[2] < 0.5*(Z0_1+ZM_1)) n2coord[2] += (ZM_1-Z0_1);
											if(n3coord[2] < 0.5*(Z0_1+ZM_1)) n3coord[2] += (ZM_1-Z0_1);
											if(n4coord[2] < 0.5*(Z0_1+ZM_1)) n4coord[2] += (ZM_1-Z0_1);
											
										}
										else if(caps >= 189 && caps <= 329){
											
											if(n1coord[2] > 0.5*(Z0_2+ZM_2)) n1coord[2] += -(ZM_2-Z0_2);
											if(n2coord[2] > 0.5*(Z0_2+ZM_2)) n2coord[2] += -(ZM_2-Z0_2);
											if(n3coord[2] > 0.5*(Z0_2+ZM_2)) n3coord[2] += -(ZM_2-Z0_2);
											if(n4coord[2] > 0.5*(Z0_2+ZM_2)) n4coord[2] += -(ZM_2-Z0_2);
											
										}
										else if(caps >= 330 && caps <= 344){
											
											if(n1coord[2] > 0.5*(Z0_3+ZM_3)) n1coord[2] += -(ZM_3-Z0_3);
											if(n2coord[2] > 0.5*(Z0_3+ZM_3)) n2coord[2] += -(ZM_3-Z0_3);
											if(n3coord[2] > 0.5*(Z0_3+ZM_3)) n3coord[2] += -(ZM_3-Z0_3);
											if(n4coord[2] > 0.5*(Z0_3+ZM_3)) n4coord[2] += -(ZM_3-Z0_3);
											
										}
									}
								}
								
								
								///
								///WLC & POW
								
								dx1_12 = n1coord[0] - n2coord[0];
								dy1_12 = n1coord[1] - n2coord[1];
								dz1_12 = n1coord[2] - n2coord[2];
								
								
								edge1_length = sqrtf(dx1_12*dx1_12+dy1_12*dy1_12+dz1_12*dz1_12);
								
								e1v[0] = dx1_12/edge1_length;//pointing outwards (tension) is positive direction
								e1v[1] = dy1_12/edge1_length;
								e1v[2] = dz1_12/edge1_length;
								
								e2v[0] = -dx1_12/edge1_length;
								e2v[1] = -dy1_12/edge1_length;
								e2v[2] = -dz1_12/edge1_length;
								
								e_ratio = edge1_length*1.e-06/LMAX[edge];
								
								if(e_ratio > 0.95){
									e_ratio = 0.95;
								}
								if(e_ratio < 0.02435){
									
									edge1_length = 0.02435*LMAX[edge]/1.e-06;
								}
								wlc = -TEMPERATURE*BOLTZ/Plength[edge]*(e_ratio - 0.25 + 1./(4.*(1.-e_ratio)*(1.-e_ratio)));
								comp = kp[edge]*(1./(edge1_length*edge1_length*1.e-012));
								
								
								
								f_n1[0] = relaxfactor2*e1v[0]*wlc;
								f_n1[1] = relaxfactor2*e1v[1]*wlc;
								f_n1[2] = relaxfactor2*e1v[2]*wlc;
								f_n1[0] = f_n1[0] + relaxfactor2*e1v[0]*comp;
								f_n1[1] = f_n1[1] + relaxfactor2*e1v[1]*comp;
								f_n1[2] = f_n1[2] + relaxfactor2*e1v[2]*comp;
								
								
								f_n2[0] = relaxfactor2*e2v[0]*wlc;
								f_n2[1] = relaxfactor2*e2v[1]*wlc;
								f_n2[2] = relaxfactor2*e2v[2]*wlc;
								f_n2[0] = f_n2[0] + relaxfactor2*e2v[0]*comp;
								f_n2[1] = f_n2[1] + relaxfactor2*e2v[1]*comp;
								f_n2[2] = f_n2[2] + relaxfactor2*e2v[2]*comp;
								
								
								///Geometrics 	
								tri1_area = TRI_AREA[caps][tri1];
								tri2_area = TRI_AREA[caps][tri2];
								
								
								a1check = 0;
								a2check = 0;
									
								if(tri1_area < 0.1*TRI_AREA_eq[tri1]){
									a1check = 1;//problematic element, flag it for min area threshold exceeded
									tri1_area = 0.1*TRI_AREA_eq[tri1];
								}
								else if(tri1_area > 10.*TRI_AREA_eq[tri1]){
									a1check = 2;//problematic element, flag it for max area threshold exceeded
									tri1_area = 10.*TRI_AREA_eq[tri1];
								}
								
								if(tri2_area < 0.1*TRI_AREA_eq[tri2]){
									a2check = 1;//problematic element, flag it for min area threshold exceeded
									tri2_area = 0.1*TRI_AREA_eq[tri2];
								}
								else if(tri2_area > 10.*TRI_AREA_eq[tri2]){
									a2check = 2;//problematic element, flag it for max area threshold exceeded
									tri2_area = 10.*TRI_AREA_eq[tri2];
								}
								
								
								
								for(i=0;i<=2;i++){
									tri1_n[i] = TRI_NORMAL[i][caps][tri1];
									tri2_n[i] = TRI_NORMAL[i][caps][tri2];
								}
									
								
								for(i=0;i<=2;i++){
									area_diff[i] = TRI_NORMAL[i][caps][tri1] - TRI_NORMAL[i][caps][tri2];
									centroid_diff[i] = (TRI_CENTROID[i][caps][tri1] - TRI_CENTROID[i][caps][tri2])*1.e-006;
								}
								
								
								costheta = TRI_NORMAL[0][caps][tri1]*TRI_NORMAL[0][caps][tri2] + TRI_NORMAL[1][caps][tri1]*TRI_NORMAL[1][caps][tri2] + TRI_NORMAL[2][caps][tri1]*TRI_NORMAL[2][caps][tri2];
								
								//make sure costheta does not exceed 1.0 in absolute value
								if(costheta > 1) costheta = 1.;
								else if(costheta < -1) costheta = -1.;
								
								angle_check = 0;
								for(i=0;i<=2;i++){
									angle_check += area_diff[i]*centroid_diff[i];
								}
										
								if(angle_check >= 0) sintheta = sqrtf(1.-costheta*costheta);
								else sintheta = -sqrtf(1.-costheta*costheta);
								
								theta = asinf(sintheta);
								
								sintheta0 = SINETHETA0[edge];
								costheta0 = COSINETHETA0[edge];
								
								
								
								///Bending
								
								if(fabsf(costheta) > 0.99)betabend = 0;
								else betabend = kBend*(sintheta*costheta0 - costheta*sintheta0)/sintheta;
								
								for(i=0;i<=2;i++){
									dp21[i] = 1.e-6*(n2coord[i] - n1coord[i]);
									dp24[i] = 1.e-6*(n2coord[i] - n4coord[i]);
									dp32[i] = 1.e-6*(n3coord[i] - n2coord[i]);
									dp41[i] = 1.e-6*(n4coord[i] - n1coord[i]);
									dp13[i] = 1.e-6*(n1coord[i] - n3coord[i]);
								}
										
								//define vector cross products
								cpdt1[0] = tri1_area*(tri1_n[1]*dp21[2] - tri1_n[2]*dp21[1]);
								cpdt1[1] = tri1_area*(tri1_n[2]*dp21[0] - tri1_n[0]*dp21[2]);
								cpdt1[2] = tri1_area*(tri1_n[0]*dp21[1] - tri1_n[1]*dp21[0]);
								
								cpdt2[0] = tri2_area*(tri2_n[1]*dp21[2] - tri2_n[2]*dp21[1]);
								cpdt2[1] = tri2_area*(tri2_n[2]*dp21[0] - tri2_n[0]*dp21[2]);
								cpdt2[2] = tri2_area*(tri2_n[0]*dp21[1] - tri2_n[1]*dp21[0]);
								
								cpdt3[0] = tri1_area*(tri1_n[1]*dp24[2] - tri1_n[2]*dp24[1]);
								cpdt3[1] = tri1_area*(tri1_n[2]*dp24[0] - tri1_n[0]*dp24[2]);
								cpdt3[2] = tri1_area*(tri1_n[0]*dp24[1] - tri1_n[1]*dp24[0]);
								
								cpdt4[0] = tri2_area*(tri2_n[1]*dp24[2] - tri2_n[2]*dp24[1]);
								cpdt4[1] = tri2_area*(tri2_n[2]*dp24[0] - tri2_n[0]*dp24[2]);
								cpdt4[2] = tri2_area*(tri2_n[0]*dp24[1] - tri2_n[1]*dp24[0]);
								
								cpdt5[0] = tri1_area*(tri1_n[1]*dp32[2] - tri1_n[2]*dp32[1]);
								cpdt5[1] = tri1_area*(tri1_n[2]*dp32[0] - tri1_n[0]*dp32[2]);
								cpdt5[2] = tri1_area*(tri1_n[0]*dp32[1] - tri1_n[1]*dp32[0]);
								
								cpdt6[0] = tri2_area*(tri2_n[1]*dp32[2] - tri2_n[2]*dp32[1]);
								cpdt6[1] = tri2_area*(tri2_n[2]*dp32[0] - tri2_n[0]*dp32[2]);
								cpdt6[2] = tri2_area*(tri2_n[0]*dp32[1] - tri2_n[1]*dp32[0]);
								
								cpdt7[0] = tri1_area*(tri1_n[1]*dp41[2] - tri1_n[2]*dp41[1]);
								cpdt7[1] = tri1_area*(tri1_n[2]*dp41[0] - tri1_n[0]*dp41[2]);
								cpdt7[2] = tri1_area*(tri1_n[0]*dp41[1] - tri1_n[1]*dp41[0]);
								
								cpdt8[0] = tri2_area*(tri2_n[1]*dp41[2] - tri2_n[2]*dp41[1]);
								cpdt8[1] = tri2_area*(tri2_n[2]*dp41[0] - tri2_n[0]*dp41[2]);
								cpdt8[2] = tri2_area*(tri2_n[0]*dp41[1] - tri2_n[1]*dp41[0]);
								
								cpdt9[0] = tri1_area*(tri1_n[1]*dp13[2] - tri1_n[2]*dp13[1]);
								cpdt9[1] = tri1_area*(tri1_n[2]*dp13[0] - tri1_n[0]*dp13[2]);
								cpdt9[2] = tri1_area*(tri1_n[0]*dp13[1] - tri1_n[1]*dp13[0]);
								
								cpdt10[0] = tri2_area*(tri2_n[1]*dp13[2] - tri2_n[2]*dp13[1]);
								cpdt10[1] = tri2_area*(tri2_n[2]*dp13[0] - tri2_n[0]*dp13[2]);
								cpdt10[2] = tri2_area*(tri2_n[0]*dp13[1] - tri2_n[1]*dp13[0]);
								
								for(i=0;i<=2;i++){
									tri1_localcentroid[i] = (TRI_CENTROID[i][caps][tri1] - CAPS_CENTROID[i][caps])*1.e-006;
									tri2_localcentroid[i] = (TRI_CENTROID[i][caps][tri2] - CAPS_CENTROID[i][caps])*1.e-006;
								}
										
								cpdt1v[0] = (tri1_localcentroid[1]*dp21[2] - tri1_localcentroid[2]*dp21[1]);
								cpdt1v[1] = (tri1_localcentroid[2]*dp21[0] - tri1_localcentroid[0]*dp21[2]);
								cpdt1v[2] = (tri1_localcentroid[0]*dp21[1] - tri1_localcentroid[1]*dp21[0]);
								
								cpdt5v[0] = (tri1_localcentroid[1]*dp32[2] - tri1_localcentroid[2]*dp32[1]);
								cpdt5v[1] = (tri1_localcentroid[2]*dp32[0] - tri1_localcentroid[0]*dp32[2]);
								cpdt5v[2] = (tri1_localcentroid[0]*dp32[1] - tri1_localcentroid[1]*dp32[0]);
								
								cpdt9v[0] = (tri1_localcentroid[1]*dp13[2] - tri1_localcentroid[2]*dp13[1]);
								cpdt9v[1] = (tri1_localcentroid[2]*dp13[0] - tri1_localcentroid[0]*dp13[2]);
								cpdt9v[2] = (tri1_localcentroid[0]*dp13[1] - tri1_localcentroid[1]*dp13[0]);
								
								
								cpdt4v[0] = (tri2_localcentroid[1]*dp24[2] - tri2_localcentroid[2]*dp24[1]);
								cpdt4v[1] = (tri2_localcentroid[2]*dp24[0] - tri2_localcentroid[0]*dp24[2]);
								cpdt4v[2] = (tri2_localcentroid[0]*dp24[1] - tri2_localcentroid[1]*dp24[0]);
								
								cpdt8v[0] = (tri2_localcentroid[1]*dp41[2] - tri2_localcentroid[2]*dp41[1]);
								cpdt8v[1] = (tri2_localcentroid[2]*dp41[0] - tri2_localcentroid[0]*dp41[2]);
								cpdt8v[2] = (tri2_localcentroid[0]*dp41[1] - tri2_localcentroid[1]*dp41[0]);
								
								cpdt12v[0] = -(tri2_localcentroid[1]*dp21[2] - tri2_localcentroid[2]*dp21[1]);
								cpdt12v[1] = -(tri2_localcentroid[2]*dp21[0] - tri2_localcentroid[0]*dp21[2]);
								cpdt12v[2] = -(tri2_localcentroid[0]*dp21[1] - tri2_localcentroid[1]*dp21[0]);
								
								b11 = -betabend*costheta/(tri1_area*tri1_area);
								b12 = betabend/(tri1_area*tri2_area);
								b22 = -betabend*costheta/(tri2_area*tri2_area);
								
								for(i=0;i<=2;i++){
									f_b1[i] = (b11*cpdt5[i] + b12*(cpdt3[i]+cpdt6[i]) + b22*cpdt4[i]);
									f_b2[i] = (b11*cpdt9[i] + b12*(cpdt7[i]+cpdt10[i]) + b22*cpdt8[i]);
									f_b3[i] = (b11*cpdt1[i] + b12*cpdt2[i]);
									f_b4[i] = (b12*-cpdt1[i] + b22*-cpdt2[i]);
								}
								
								check1 = 0; check2 = 0; check3 = 0; check4 = 0;
								
								for(i=0;i<=2;i++){
									if(isnan(f_b1[i]) == 1 || isinf(f_b1[i]) == 1 || fabsf(f_b1[i]) > 1.e-10){
										check1 = 1;
										printf("Iter%d, Bending problem on edge%d between tri%d and tri%d on RBC%d, node1 forces over the limit, cell status is %d regression flag is %d\n", iter, edge, tri1, tri2, caps, CellStatus[caps], REGRESSFLAG[caps]);
										CAPS_REMOVE[caps] += 1;
										simfailure += 1;
									}
									if(isnan(f_b2[i]) == 1 || isinf(f_b2[i]) == 1 || fabsf(f_b2[i]) > 1.e-10){
										check2 = 1;
										printf("Iter%d, Bending problem on edge%d between tri%d and tri%d on RBC%d, node2 forces over the limit, cell status is %d regression flag is %d\n", iter, edge, tri1, tri2, caps, CellStatus[caps], REGRESSFLAG[caps]);
										CAPS_REMOVE[caps] += 1;
										simfailure += 1;
									}
									if(isnan(f_b3[i]) == 1 || isinf(f_b3[i]) == 1 || fabsf(f_b3[i]) > 1.e-10){
										check3 = 1;
										printf("Iter%d, Bending problem on edge%d between tri%d and tri%d on RBC%d, node3 forces over the limit, cell status is %d regression flag is %d\n", iter, edge, tri1, tri2, caps, CellStatus[caps], REGRESSFLAG[caps]);
										CAPS_REMOVE[caps] += 1;
										simfailure += 1;
									}
									if(isnan(f_b4[i]) == 1 || isinf(f_b4[i]) == 1 || fabsf(f_b4[i]) > 1.e-10){
										check4 = 1;
										printf("Iter%d, Bending problem on edge%d between tri%d and tri%d on RBC%d, node4 forces over the limit, cell status is %d regression flag is %d\n", iter, edge, tri1, tri2, caps, CellStatus[caps], REGRESSFLAG[caps]);
										CAPS_REMOVE[caps] += 1;
										simfailure += 1;
									}
								}
								
								if(check1 == 0){
									for(i=0;i<=2;i++) f_n1[i] += f_b1[i];
								}
								if(check2 == 0){
									for(i=0;i<=2;i++) f_n2[i] += f_b2[i];
								}
								if(check3 == 0){
									for(i=0;i<=2;i++) f_n3[i] += f_b3[i];
								}
								if(check4 == 0){
									for(i=0;i<=2;i++) f_n4[i] += f_b4[i];
								}
								
									
								///area compression & volume compression
								if(AREA_UPDATE[caps][tri1] == 0){
									AREA_UPDATE[caps][tri1] = 1;
									alpha1 = -kd*relaxfactor*(tri1_area - TRI_AREA_eq[tri1])/(4.*TRI_AREA_eq[tri1]*tri1_area);
									
									for(i=0;i<=2;i++){
										f_a1[i] = alpha1*cpdt5[i];
										f_a2[i] = alpha1*cpdt9[i];
										f_a3[i] = alpha1*cpdt1[i];
										
										f_aa1[i] = 0.25*F_AREA_GLOBAL[caps]/tri1_area*cpdt5[i];
										f_aa2[i] = 0.25*F_AREA_GLOBAL[caps]/tri1_area*cpdt9[i];
										f_aa3[i] = 0.25*F_AREA_GLOBAL[caps]/tri1_area*cpdt1[i];
										
										f_v1[i] = 1./6.*(tri1_area*tri1_n[i]/3. + cpdt5v[i]);
										f_v2[i] = 1./6.*(tri1_area*tri1_n[i]/3. + cpdt9v[i]);
										f_v3[i] = 1./6.*(tri1_area*tri1_n[i]/3. + cpdt1v[i]);
										
									}
									
									check1 = 0; check2 = 0; check3 = 0;
								
									for(i=0;i<=2;i++){
										if(isnan(f_a1[i]) == 1 || isinf(f_a1[i]) == 1 || fabsf(f_a1[i]) > 1.e-10){
											check1 = 1;
											printf("Iter%d, local area problem on tri%d on RBC%d, node1 forces over the limit, cell status is %d regression flag is %d\n", iter, tri1, caps, CellStatus[caps], REGRESSFLAG[caps]);
											CAPS_REMOVE[caps] += 1;
											simfailure += 1;
										}
										if(isnan(f_a2[i]) == 1 || isinf(f_a2[i]) == 1 || fabsf(f_a2[i]) > 1.e-10){
											check2 = 1;
											printf("Iter%d, local area problem on tri%d on RBC%d, node2 forces over the limit, cell status is %d regression flag is %d\n", iter, tri1, caps, CellStatus[caps], REGRESSFLAG[caps]);
											CAPS_REMOVE[caps] += 1;
											simfailure += 1;
										}
										if(isnan(f_a3[i]) == 1 || isinf(f_a3[i]) == 1 || fabsf(f_a3[i]) > 1.e-10){
											check3 = 1;
											printf("Iter%d, local area problem on tri%d on RBC%d, node3 forces over the limit, cell status is %d regression flag is %d\n", iter, tri1, caps, CellStatus[caps], REGRESSFLAG[caps]);
											CAPS_REMOVE[caps] += 1;
											simfailure += 1;
										}
									}
									
									if(check1 == 0){
										for(i=0;i<=2;i++) f_n1[i] += f_a1[i];
									}
									if(check2 == 0){
										for(i=0;i<=2;i++) f_n2[i] += f_a2[i];
									}
									if(check3 == 0){
										for(i=0;i<=2;i++) f_n3[i] += f_a3[i];
									}
									
									
									check1 = 0; check2 = 0; check3 = 0;
								
									for(i=0;i<=2;i++){
										if(isnan(f_aa1[i]) == 1 || isinf(f_aa1[i]) == 1 || fabsf(f_aa1[i]) > 1.e-10){
											check1 = 1;
											printf("Iter%d, global area problem on tri%d on RBC%d, node1 forces over the limit, cell status is %d regression flag is %d\n", iter, tri1, caps, CellStatus[caps], REGRESSFLAG[caps]);
											CAPS_REMOVE[caps] += 1;
											simfailure += 1;
										}
										if(isnan(f_aa2[i]) == 1 || isinf(f_aa2[i]) == 1 || fabsf(f_aa2[i]) > 1.e-10){
											check2 = 1;
											printf("Iter%d, global area problem on tri%d on RBC%d, node2 forces over the limit, cell status is %d regression flag is %d\n", iter, tri1, caps, CellStatus[caps], REGRESSFLAG[caps]);
											CAPS_REMOVE[caps] += 1;
											simfailure += 1;
										}
										if(isnan(f_aa3[i]) == 1 || isinf(f_aa3[i]) == 1 || fabsf(f_aa3[i]) > 1.e-10){
											check3 = 1;
											printf("Iter%d, global area problem on tri%d on RBC%d, node3 forces over the limit, cell status is %d regression flag is %d\n", iter, tri1, caps, CellStatus[caps], REGRESSFLAG[caps]);
											CAPS_REMOVE[caps] += 1;
											simfailure += 1;
										}
									}
									
									if(check1 == 0){
										for(i=0;i<=2;i++) f_n1[i] += f_aa1[i];
									}
									if(check2 == 0){
										for(i=0;i<=2;i++) f_n2[i] += f_aa2[i];
									}
									if(check3 == 0){
										for(i=0;i<=2;i++) f_n3[i] += f_aa3[i];
									}
									
									
									check1 = 0; check2 = 0; check3 = 0;
								
									for(i=0;i<=2;i++){
										if(isnan(f_v1[i]) == 1 || isinf(f_v1[i]) == 1 || fabsf(f_v1[i]) > 1.e-10){
											check1 = 1;
											printf("Iter%d, vol problem on tri%d on RBC%d, node1 forces over the limit, cell status is %d regression flag is %d\n", iter, tri1, caps, CellStatus[caps], REGRESSFLAG[caps]);
											CAPS_REMOVE[caps] += 1;
											simfailure += 1;
										}
										if(isnan(f_v2[i]) == 1 || isinf(f_v2[i]) == 1 || fabsf(f_v2[i]) > 1.e-10){
											check2 = 1;
											printf("Iter%d, vol problem on tri%d on RBC%d, node2 forces over the limit, cell status is %d regression flag is %d\n", iter, tri1, caps, CellStatus[caps], REGRESSFLAG[caps]);
											CAPS_REMOVE[caps] += 1;
											simfailure += 1;
										}
										if(isnan(f_v3[i]) == 1 || isinf(f_v3[i]) == 1 || fabsf(f_v3[i]) > 1.e-10){
											check3 = 1;
											printf("Iter%d, vol problem on tri%d on RBC%d, node3 forces over the limit, cell status is %d regression flag is %d\n", iter, tri1, caps, CellStatus[caps], REGRESSFLAG[caps]);
											CAPS_REMOVE[caps] += 1;
											simfailure += 1;
										}
									}
									
									if(check1 == 0){
										for(i=0;i<=2;i++) F_VOLUME[i][caps][n1] += f_v1[i];
									}
									if(check2 == 0){
										for(i=0;i<=2;i++) F_VOLUME[i][caps][n2] += f_v2[i];
									}
									if(check3 == 0){
										for(i=0;i<=2;i++) F_VOLUME[i][caps][n3] += f_v3[i];
									}
								}
								
								
									
								if(AREA_UPDATE[caps][tri2] == 0){
									AREA_UPDATE[caps][tri2] = 1;
									alpha2 = -kd*relaxfactor*(tri2_area - TRI_AREA_eq[tri2])/(4.*TRI_AREA_eq[tri2]*tri2_area);
									
									for(i=0;i<=2;i++){
										f_a1[i] = alpha2*cpdt4[i];
										f_a2[i] = alpha2*cpdt8[i];
										f_a3[i] = alpha2*-cpdt2[i];
									
										f_aa1[i] = 0.25*F_AREA_GLOBAL[caps]/tri2_area*cpdt4[i];
										f_aa2[i] = 0.25*F_AREA_GLOBAL[caps]/tri2_area*cpdt8[i];
										f_aa3[i] = 0.25*F_AREA_GLOBAL[caps]/tri2_area*-cpdt2[i];
													
										f_v1[i] = 1./6.*(tri2_area*tri2_n[i]/3. + cpdt4v[i]);
										f_v2[i] = 1./6.*(tri2_area*tri2_n[i]/3. + cpdt8v[i]);
										f_v3[i] = 1./6.*(tri2_area*tri2_n[i]/3. + cpdt12v[i]);
									}
									
									check1 = 0; check2 = 0; check3 = 0;
								
									for(i=0;i<=2;i++){
										if(isnan(f_a1[i]) == 1 || isinf(f_a1[i]) == 1 || fabsf(f_a1[i]) > 1.e-10){
											check1 = 1;
											simfailure += 1;
											printf("Iter%d, local area problem on tri%d on RBC%d, node1 forces over the limit, cell status is %d regression flag is %d\n", iter, tri2, caps, CellStatus[caps], REGRESSFLAG[caps]);
										}
										if(isnan(f_a2[i]) == 1 || isinf(f_a2[i]) == 1 || fabsf(f_a2[i]) > 1.e-10){
											check2 = 1;
											simfailure += 1;
											printf("Iter%d, local area problem on tri%d on RBC%d, node2 forces over the limit, cell status is %d regression flag is %d\n", iter, tri2, caps, CellStatus[caps], REGRESSFLAG[caps]);
										}
										if(isnan(f_a3[i]) == 1 || isinf(f_a3[i]) == 1 || fabsf(f_a3[i]) > 1.e-10){
											check3 = 1;
											simfailure += 1;
											printf("Iter%d, local area problem on tri%d on RBC%d, node3 forces over the limit, cell status is %d regression flag is %d\n", iter, tri2, caps, CellStatus[caps], REGRESSFLAG[caps]);
										}
									}
									
									if(check1 == 0){
										for(i=0;i<=2;i++) f_n1[i] += f_a1[i];
									}
									if(check2 == 0){
										for(i=0;i<=2;i++) f_n2[i] += f_a2[i];
									}
									if(check3 == 0){
										for(i=0;i<=2;i++) f_n4[i] += f_a3[i];
									}
									
									
									check1 = 0; check2 = 0; check3 = 0;
								
									for(i=0;i<=2;i++){
										if(isnan(f_aa1[i]) == 1 || isinf(f_aa1[i]) == 1 || fabsf(f_aa1[i]) > 1.e-10){
											check1 = 1;
											simfailure += 1;
											printf("Iter%d, global area problem on tri%d on RBC%d, node1 forces over the limit, cell status is %d regression flag is %d\n", iter, tri2, caps, CellStatus[caps], REGRESSFLAG[caps]);
										}
										if(isnan(f_aa2[i]) == 1 || isinf(f_aa2[i]) == 1 || fabsf(f_aa2[i]) > 1.e-10){
											check2 = 1;
											simfailure += 1;
											printf("Iter%d, global area problem on tri%d on RBC%d, node2 forces over the limit, cell status is %d regression flag is %d\n", iter, tri2, caps, CellStatus[caps], REGRESSFLAG[caps]);
										}
										if(isnan(f_aa3[i]) == 1 || isinf(f_aa3[i]) == 1 || fabsf(f_aa3[i]) > 1.e-10){
											check3 = 1;
											simfailure += 1;
											printf("Iter%d, global area problem on tri%d on RBC%d, node3 forces over the limit, cell status is %d regression flag is %d\n", iter, tri2, caps, CellStatus[caps], REGRESSFLAG[caps]);
										}
									}
									
									if(check1 == 0){
										for(i=0;i<=2;i++) f_n1[i] += f_aa1[i];
									}
									if(check2 == 0){
										for(i=0;i<=2;i++) f_n2[i] += f_aa2[i];
									}
									if(check3 == 0){
										for(i=0;i<=2;i++) f_n4[i] += f_aa3[i];
									}
									
									check1 = 0; check2 = 0; check3 = 0;
								
									for(i=0;i<=2;i++){
										if(isnan(f_v1[i]) == 1 || isinf(f_v1[i]) == 1 || fabsf(f_v1[i]) > 1.e-10){
											check1 = 1;
											simfailure += 1;
											printf("Iter%d, vol problem on tri%d on RBC%d, node1 forces over the limit, cell status is %d regression flag is %d\n", iter, tri2, caps, CellStatus[caps], REGRESSFLAG[caps]);
										}
										if(isnan(f_v2[i]) == 1 || isinf(f_v2[i]) == 1 || fabsf(f_v2[i]) > 1.e-10){
											check2 = 1;
											simfailure += 1;
											printf("Iter%d, vol problem on tri%d on RBC%d, node1 forces over the limit, cell status is %d regression flag is %d\n", iter, tri2, caps, CellStatus[caps], REGRESSFLAG[caps]);
										}
										if(isnan(f_v3[i]) == 1 || isinf(f_v3[i]) == 1 || fabsf(f_v3[i]) > 1.e-10){
											check3 = 1;
											simfailure += 1;
											printf("Iter%d, vol problem on tri%d on RBC%d, node1 forces over the limit, cell status is %d regression flag is %d\n", iter, tri2, caps, CellStatus[caps], REGRESSFLAG[caps]);
										}
									}
									
									if(check1 == 0){
										for(i=0;i<=2;i++) F_VOLUME[i][caps][n1] += f_v1[i];
									}
									if(check2 == 0){
										for(i=0;i<=2;i++) F_VOLUME[i][caps][n2] += f_v2[i];
									}
									if(check3 == 0){
										for(i=0;i<=2;i++) F_VOLUME[i][caps][n4] += f_v3[i];
									}
									
								}
								
								
								///membrane viscous force
								if(MViscOn == 1){//this value is 0, the membrane viscosity feature is not included in the paper
									for(i=0;i<=2;i++) vel_12[i] = NODE_VEL_OLD[i][caps][n1]-NODE_VEL_OLD[i][caps][n2];
									
									
									for(i=0;i<=2;i++){
										viscf1[i] = gamma_t*(-vel_12[i] - 1./3.*((vel_12[0]*e1v[0]) + (vel_12[1]*e1v[1]) + (vel_12[2]*e1v[2]))*e1v[i]);
										viscf2[i] = gamma_t*(vel_12[i] - 1./3.*((-vel_12[0]*e2v[0]) + (-vel_12[1]*e2v[1]) + (-vel_12[2]*e2v[2]))*e2v[i]);
									}
								
									for(i=0;i<=2;i++){
										F_Total[i][caps][n1] += relaxfactor*viscf1[i];
										F_Total[i][caps][n2] += relaxfactor*viscf2[i];
									}
								}
								
								for(i=0;i<=2;i++){
									F_Total[i][caps][n1] += f_n1[i];
									F_Total[i][caps][n2] += f_n2[i];
									F_Total[i][caps][n3] += f_n3[i];
									F_Total[i][caps][n4] += f_n4[i];
									
								}
							}
						}
					
						for(edge=1;edge<=EdgeM_N;edge++){
							n1 = EDGENODE_LIST_N[edge][1];
							n2 = EDGENODE_LIST_N[edge][2];
							
							check = 0;
							for(i=0;i<=2;i++){
								f_n1[i] = 0.; f_n2[i] = 0.;
							}
															
							n1coord[0] = NODECOORD_N[0][caps][n1];
							n1coord[1] = NODECOORD_N[1][caps][n1];
							n1coord[2] = NODECOORD_N[2][caps][n1];
							
							n2coord[0] = NODECOORD_N[0][caps][n2];
							n2coord[1] = NODECOORD_N[1][caps][n2];
							n2coord[2] = NODECOORD_N[2][caps][n2];
							
							if(CellStatus[caps] > 3000){//CA and CV exit candidates
								if(n1coord[2] > Z0_2 || n2coord[2] > Z0_2 || n1coord[2] < ZM_1 || n2coord[2] < ZM_1) check = 1;
							}
							
							
							if(check == 0){
								if(CellStatus[caps] > 3){
									if((nzmax - nzmin) > 30.){//this cell is crossing periodic BC
										if(caps <= 188){
											if(n1coord[2] < 0.5*(Z0_1+ZM_1)) n1coord[2] += (ZM_1-Z0_1);
											if(n2coord[2] < 0.5*(Z0_1+ZM_1)) n2coord[2] += (ZM_1-Z0_1);
										}
										else if(caps >= 189 && caps <= 329){
											if(n1coord[2] > 0.5*(Z0_2+ZM_2)) n1coord[2] += -(ZM_2-Z0_2);
											if(n2coord[2] > 0.5*(Z0_2+ZM_2)) n2coord[2] += -(ZM_2-Z0_2);
										}
										else if(caps >= 330 && caps <= 344){
											if(n1coord[2] > 0.5*(Z0_3+ZM_3)) n1coord[2] += -(ZM_3-Z0_3);
											if(n2coord[2] > 0.5*(Z0_3+ZM_3)) n2coord[2] += -(ZM_3-Z0_3);
										}
									}
								}
								
								///
								///WLC & POW
								dx1_12 = n1coord[0] - n2coord[0];
								dy1_12 = n1coord[1] - n2coord[1];
								dz1_12 = n1coord[2] - n2coord[2];
								
								edge1_length = sqrtf(dx1_12*dx1_12+dy1_12*dy1_12+dz1_12*dz1_12);
								
								e1v[0] = dx1_12/edge1_length;//pointing outwards (tension) is positive direction
								e1v[1] = dy1_12/edge1_length;
								e1v[2] = dz1_12/edge1_length;
								
								e2v[0] = -dx1_12/edge1_length;
								e2v[1] = -dy1_12/edge1_length;
								e2v[2] = -dz1_12/edge1_length;
								
								e_ratio = edge1_length*1.e-06/LMAX_N[edge];
								
								if(e_ratio > 0.95) e_ratio = 0.95;
								if(e_ratio < 0.02435) edge1_length = 0.02435*LMAX_N[edge]/1.e-06;
								
								wlc = -TEMPERATURE*BOLTZ/Plength_N[edge]*(e_ratio - 0.25 + 1./(4.*(1.-e_ratio)*(1.-e_ratio)));
								comp = kp_N[edge]*(1./(edge1_length*edge1_length*1.e-012));
								
								f_n1[0] = relaxfactor2*e1v[0]*wlc;
								f_n1[1] = relaxfactor2*e1v[1]*wlc;
								f_n1[2] = relaxfactor2*e1v[2]*wlc;
								f_n1[0] = f_n1[0] + relaxfactor2*e1v[0]*comp;
								f_n1[1] = f_n1[1] + relaxfactor2*e1v[1]*comp;
								f_n1[2] = f_n1[2] + relaxfactor2*e1v[2]*comp;
								
								f_n2[0] = relaxfactor2*e2v[0]*wlc;
								f_n2[1] = relaxfactor2*e2v[1]*wlc;
								f_n2[2] = relaxfactor2*e2v[2]*wlc;
								f_n2[0] = f_n2[0] + relaxfactor2*e2v[0]*comp;
								f_n2[1] = f_n2[1] + relaxfactor2*e2v[1]*comp;
								f_n2[2] = f_n2[2] + relaxfactor2*e2v[2]*comp;
								
								for(i=0;i<=2;i++){
									F_Total_N[i][caps][n1] += f_n1[i];
									F_Total_N[i][caps][n2] += f_n2[i];
								}
							}
						}
						
						if(caps > 344 && CAPS_REMOVE[caps] > 500){
							printf("Iter%d, RBC%d is being removed from the domain for bad geometries\n", iter, caps);
							CellStatus[caps] = -5000;
							
							REGRESSFLAG[caps] = 0;
							NucMembUpdate[caps] = 0;
							for(node=1;node<=NodeM;node++){
								
								NODE_WALLNEIGH[caps][node] = 0;
								MembIn_NeighFlag[caps][node] = 0;
								MembOut_NeighFlag[caps][node] = 0;
								
								for(i=0;i<=2;i++){
									NODE_WALLREPUL[i][caps][node] = 0;
									NODECOORD[i][caps][node] = -5000.;
									F_Total[i][caps][node] = 0.;
									NODE_VEL[i][caps][node] = 0.;
								}
							}
							
							for(node=1;node<=NodeM_N;node++){
								NUC_NeighFlag[caps][node] = 0;
								for(i=0;i<=2;i++){
									NODECOORD_N[i][caps][node] = -5000.;
									F_Total_N[i][caps][node] = 0.;
									NODE_VEL_N[i][caps][node] = 0.;
								}
							}
						}
					
					}
				}
			}
			
		}
		/// end of membrane force calculation 
		
		///CGSM_VOLUME_CALC();
		
		if(RBCSOLVERFLAG == 1){
			#pragma omp parallel private(caps, tri, node, i, CellVol)
			{
				#pragma omp for schedule(static,50)  nowait
				for(caps=1;caps<=CapsM;caps++){
					if(CellStatus[caps] > -5000 && CellStatus[caps] < 3000){
						CellVol = MCV;
						CAPS_VOL[caps] = 0.;
						for(tri=1;tri<=TriM;tri++){
							CAPS_VOL[caps] = CAPS_VOL[caps] + TRI_AREA[caps][tri]*(TRI_CENTROID[0][caps][tri]*TRI_NORMAL[0][caps][tri] + TRI_CENTROID[1][caps][tri]*TRI_NORMAL[1][caps][tri] + TRI_CENTROID[2][caps][tri]*TRI_NORMAL[2][caps][tri])*1.e-006;
						}
						CAPS_VOL[caps] = CAPS_VOL[caps]/3.;
						F_VOLUME_SCALAR[caps] = -kv*(CAPS_VOL[caps] - CellVol)/CellVol;
						
						for(node=1;node<=NodeM;node++){
							for(i=0;i<=2;i++){
								F_VOLUME[i][caps][node] = relaxfactor*F_VOLUME_SCALAR[caps]*F_VOLUME[i][caps][node];
								F_Total[i][caps][node] += F_VOLUME[i][caps][node];
							}
						}
					}
				}
			}
		}
		
		
		///RBC velocity interpolation and Wall repulsion velocity
		/// Parent cells
		
		if(RBCSOLVERFLAG == 1){
			#pragma omp parallel private(ReflectProj,VelReflect,o,Ffluid,caps,node,i,xp,yp,zp,j,k,i1,j1,k1,i2,j2,k2,weightsum,weightcount,bvalue,mindis,triwall,xc,yc,zc,kmap,validinter,distance,interweight,wneighmax,m,nodewall,neighmax,l,tri,n1,n2,n3,xp1,yp1,zp1,repulf,v,xp2,yp2,zp2,xp3,yp3,zp3,xp4,yp4,zp4,zpmap)
			{
				#pragma omp for schedule(static,5)  nowait
				for(caps=1;caps<=344;caps++){
					if(CellStatus[caps] == 2){//Parent cells in parent interior domains
						//DAparent
						if(caps <= 188){
							for(node=1;node<=NodeM;node++){
								
								///interpolate velocities on membrane from surrounding fluid
								xp = NODECOORD[0][caps][node]; yp = NODECOORD[1][caps][node]; zp = NODECOORD[2][caps][node];
								
								i = (int)((xp-X0_1)/DX + 1.01);
								j = (int)((yp-Y0_1)/DX + 1.01);
								k = (int)((zp-Z0_1)/DX + 1.01);
								
								i1 = i-4; i2 = i+4;
								j1 = j-4; j2 = j+4;
								k1 = k-4; k2 = k+4;
								
								if(i1 < 1) i1 = 1; if(i2> IM_1) i2 = IM_1;
								if(j1 < 1) j1 = 1; if(j2> JM_1) j2 = JM_1;
								
								weightsum = 0; weightcount = 0; bvalue = 0; mindis = 1.5; triwall = 0;
								for(i=i1;i<=i2;i++){
									for(j=j1;j<=j2;j++){
										for(k=k1;k<=k2;k++){
											
											xc = X0_1 + (i-1)*DX;
											yc = Y0_1 + (j-1)*DX;
											zc = Z0_1 + (k-1)*DX;
											
											if(k < 1) kmap = k + (KM_1-1);
											else if(k > KM_1) kmap = k - (KM_1-1);
											else kmap = k;
											
											o = DOMID_1[i][j][kmap];
											
											if(B_index_1[i][j][kmap] > 0){
												
												distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
												if(distance <= StencilWidth*DX){
													weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
													NODE_VEL[0][caps][node] += weightsum*(VEL_SCALE*U_1[o]);
													NODE_VEL[1][caps][node] += weightsum*(VEL_SCALE*V_1[o]);
													NODE_VEL[2][caps][node] += weightsum*(VEL_SCALE*W_1[o]);
													
												}
											}
											else if(W_index_1[i][j][kmap] > 0){
												distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
												if(distance <= WallStencilWidth*DX){
													weightsum = (1.+cosf(PI*fabsf(xc-xp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(WallStencilWidth*DX)))/((WallStencilWidth*2)*(WallStencilWidth*2)*(WallStencilWidth*2));
													NODE_VEL[0][caps][node] += weightsum*(VEL_SCALE*Uvoid_1[i][j][kmap]);
													NODE_VEL[1][caps][node] += weightsum*(VEL_SCALE*Vvoid_1[i][j][kmap]);
													NODE_VEL[2][caps][node] += weightsum*(VEL_SCALE*Wvoid_1[i][j][kmap]);
													
												}
											}
											
											if(W_index_1[i][j][kmap] > 0){
												nodewall = W_index_1[i][j][kmap];
												
												neighmax = WALLNODETOCELLNEIGHBORLIST_1[0][nodewall];
														
												for(l=1;l<=neighmax;l++){
													tri = WALLNODETOCELLNEIGHBORLIST_1[l][nodewall];
													n1 = WALLCELLNODE_1[0][tri]; n2 = WALLCELLNODE_1[1][tri]; n3 = WALLCELLNODE_1[2][tri];
													xp1 = WALLNODECOORD_1[0][n1];
													yp1 = WALLNODECOORD_1[1][n1];
													zp1 = WALLNODECOORD_1[2][n1];
													distance = sqrtf((xp-xp1)*(xp-xp1) + (yp-yp1)*(yp-yp1) + (zp-zp1)*(zp-zp1));
													if(distance < mindis){
														mindis = distance;
														triwall = tri;
													}
													
													xp2 = WALLNODECOORD_1[0][n2];
													yp2 = WALLNODECOORD_1[1][n2];
													zp2 = WALLNODECOORD_1[2][n2];
													distance = sqrtf((xp-xp2)*(xp-xp2) + (yp-yp2)*(yp-yp2) + (zp-zp2)*(zp-zp2));
													if(distance < mindis){
														mindis = distance;
														triwall = tri;
													}
													
													xp3 = WALLNODECOORD_1[0][n3];
													yp3 = WALLNODECOORD_1[1][n3];
													zp3 = WALLNODECOORD_1[2][n3];
													distance = sqrtf((xp-xp3)*(xp-xp3) + (yp-yp3)*(yp-yp3) + (zp-zp3)*(zp-zp3));
													if(distance < mindis){
														mindis = distance;
														triwall = tri;
													}
													
													xp4 = WALL_CENTROID_1[0][tri];
													yp4 = WALL_CENTROID_1[1][tri];
													zp4 = WALL_CENTROID_1[2][tri];
													distance = sqrtf((xp-xp4)*(xp-xp4) + (yp-yp4)*(yp-yp4) + (zp-zp4)*(zp-zp4));
													if(distance < mindis){
														mindis = distance;
														triwall = tri;
													}	
												}
												
											}
										}
									}
								}
								
								if(triwall > 0){
									NODE_WALLNEIGH[caps][node] = triwall+2400000;
									n1 = WALLCELLNODE_1[0][triwall]; n2 = WALLCELLNODE_1[1][triwall]; n3 = WALLCELLNODE_1[2][triwall];
									xp1 = WALLNODECOORD_1[0][n1];
									yp1 = WALLNODECOORD_1[1][n1];
									zp1 = WALLNODECOORD_1[2][n1];
									
									distance = (xp-xp1)*WALLCELLNORMALS_1[0][triwall] + (yp-yp1)*WALLCELLNORMALS_1[1][triwall] + (zp-zp1)*WALLCELLNORMALS_1[2][triwall];
									
									if(distance < 0.25){
										if(fabsf(distance) < 1.0){
											ReflectProj = 0;
											repulf = (0.25 - distance)/1.0;
											if(distance < -0.75) repulf = 1.;
											
											for(v=0;v<=2;v++){
												NODE_WALLREPUL[v][caps][node] += repulf*(1.e-04)*WALLCELLNORMALS_1[v][triwall];
											}
										}
									}
								}	
							}
						
							for(node=1;node<=NodeM_N;node++){
								
								///interpolate velocities on membrane from surrounding fluid
								xp = NODECOORD_N[0][caps][node]; yp = NODECOORD_N[1][caps][node]; zp = NODECOORD_N[2][caps][node];
								
								i = (int)((xp-X0_1)/DX + 1.01);
								j = (int)((yp-Y0_1)/DX + 1.01);
								k = (int)((zp-Z0_1)/DX + 1.01);
								
								
								i1 = i-4; i2 = i+4;
								j1 = j-4; j2 = j+4;
								k1 = k-4; k2 = k+4;
								
								if(i1 < 1) i1 = 1; if(i2> IM_1) i2 = IM_1;
								if(j1 < 1) j1 = 1; if(j2> JM_1) j2 = JM_1;
								weightsum = 0; weightcount = 0; bvalue = 0;
								for(i=i1;i<=i2;i++){
									for(j=j1;j<=j2;j++){
										for(k=k1;k<=k2;k++){
											
											xc = X0_1 + (i-1)*DX;
											yc = Y0_1 + (j-1)*DX;
											zc = Z0_1 + (k-1)*DX;
											
											if(k < 1) kmap = k + (KM_1-1);
											else if(k > KM_1) kmap = k - (KM_1-1);
											else kmap = k;
											
											o = DOMID_1[i][j][kmap];
											
											if(B_index_1[i][j][kmap] > 0){
												
												distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
												if(distance <= StencilWidth*DX){
													weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
													NODE_VEL_N[0][caps][node] += weightsum*(VEL_SCALE*U_1[o]);
													NODE_VEL_N[1][caps][node] += weightsum*(VEL_SCALE*V_1[o]);
													NODE_VEL_N[2][caps][node] += weightsum*(VEL_SCALE*W_1[o]);
												}
											}
											else if(W_index_1[i][j][kmap] > 0){
												distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
												if(distance <= WallStencilWidth*DX){
													weightsum = (1.+cosf(PI*fabsf(xc-xp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(WallStencilWidth*DX)))/((WallStencilWidth*2)*(WallStencilWidth*2)*(WallStencilWidth*2));
													NODE_VEL_N[0][caps][node] += weightsum*(VEL_SCALE*Uvoid_1[i][j][kmap]);
													NODE_VEL_N[1][caps][node] += weightsum*(VEL_SCALE*Vvoid_1[i][j][kmap]);
													NODE_VEL_N[2][caps][node] += weightsum*(VEL_SCALE*Wvoid_1[i][j][kmap]);
													
												}
											}
										}
									}
								}
							}
						
						}
						//PCV1 parent
						if(caps >= 189 && caps <= 329){
							for(node=1;node<=NodeM;node++){
								
								///interpolate velocities on membrane from surrounding fluid
								xp = NODECOORD[0][caps][node]; yp = NODECOORD[1][caps][node]; zp = NODECOORD[2][caps][node];
								
								i = (int)((xp-X0_2)/DX + 1.01);
								j = (int)((yp-Y0_2)/DX + 1.01);
								k = (int)((zp-Z0_2)/DX + 1.01);
								
								
								i1 = i-4; i2 = i+4;
								j1 = j-4; j2 = j+4;
								k1 = k-4; k2 = k+4;
								
								if(i1 < 1) i1 = 1; if(i2> IM_2) i2 = IM_2;
								if(j1 < 1) j1 = 1; if(j2> JM_2) j2 = JM_2;
								weightsum = 0; weightcount = 0; bvalue = 0; mindis = 1.5; triwall = 0;
								for(i=i1;i<=i2;i++){
									for(j=j1;j<=j2;j++){
										for(k=k1;k<=k2;k++){
											xc = X0_2 + (i-1)*DX;
											yc = Y0_2 + (j-1)*DX;
											zc = Z0_2 + (k-1)*DX;
											
											if(k < 1) kmap = k + (KM_2-1);
											else if(k > KM_2) kmap = k - (KM_2-1);
											else kmap = k;
											
											o = DOMID_2[i][j][kmap];
											
											if(B_index_2[i][j][kmap] > 0){
												
												distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
												if(distance <= StencilWidth*DX){
													weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
													NODE_VEL[0][caps][node] += weightsum*(VEL_SCALE*U_2[o]);
													NODE_VEL[1][caps][node] += weightsum*(VEL_SCALE*V_2[o]);
													NODE_VEL[2][caps][node] += weightsum*(VEL_SCALE*W_2[o]);
												}
											}
											else if(W_index_2[i][j][kmap] > 0){
												distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
												if(distance <= WallStencilWidth*DX){
													weightsum = (1.+cosf(PI*fabsf(xc-xp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(WallStencilWidth*DX)))/((WallStencilWidth*2)*(WallStencilWidth*2)*(WallStencilWidth*2));
													NODE_VEL[0][caps][node] += weightsum*(VEL_SCALE*Uvoid_2[i][j][kmap]);
													NODE_VEL[1][caps][node] += weightsum*(VEL_SCALE*Vvoid_2[i][j][kmap]);
													NODE_VEL[2][caps][node] += weightsum*(VEL_SCALE*Wvoid_2[i][j][kmap]);
													
												}
											}
											
											
											if(W_index_2[i][j][kmap] > 0){
												nodewall = W_index_2[i][j][kmap];
												
												neighmax = WALLNODETOCELLNEIGHBORLIST_2[0][nodewall];
												for(l=1;l<=neighmax;l++){
													tri = WALLNODETOCELLNEIGHBORLIST_2[l][nodewall];
													n1 = WALLCELLNODE_2[0][tri]; n2 = WALLCELLNODE_2[1][tri]; n3 = WALLCELLNODE_2[2][tri];
													xp1 = WALLNODECOORD_2[0][n1];
													yp1 = WALLNODECOORD_2[1][n1];
													zp1 = WALLNODECOORD_2[2][n1];
													distance = sqrtf((xp-xp1)*(xp-xp1) + (yp-yp1)*(yp-yp1) + (zp-zp1)*(zp-zp1));
													if(distance < mindis){
														mindis = distance;
														triwall = tri;
													}
													
													xp2 = WALLNODECOORD_2[0][n2];
													yp2 = WALLNODECOORD_2[1][n2];
													zp2 = WALLNODECOORD_2[2][n2];
													distance = sqrtf((xp-xp2)*(xp-xp2) + (yp-yp2)*(yp-yp2) + (zp-zp2)*(zp-zp2));
													if(distance < mindis){
														mindis = distance;
														triwall = tri;
													}
													
													xp3 = WALLNODECOORD_2[0][n3];
													yp3 = WALLNODECOORD_2[1][n3];
													zp3 = WALLNODECOORD_2[2][n3];
													distance = sqrtf((xp-xp3)*(xp-xp3) + (yp-yp3)*(yp-yp3) + (zp-zp3)*(zp-zp3));
													if(distance < mindis){
														mindis = distance;
														triwall = tri;
													}
													
													xp4 = WALL_CENTROID_2[0][tri];
													yp4 = WALL_CENTROID_2[1][tri];
													zp4 = WALL_CENTROID_2[2][tri];
													distance = sqrtf((xp-xp4)*(xp-xp4) + (yp-yp4)*(yp-yp4) + (zp-zp4)*(zp-zp4));
													if(distance < mindis){
														mindis = distance;
														triwall = tri;
													}	
												}
												
											}
										}
									}
								}
								
								if(triwall > 0){
									NODE_WALLNEIGH[caps][node] = triwall+3000000;
									n1 = WALLCELLNODE_2[0][triwall]; n2 = WALLCELLNODE_2[1][triwall]; n3 = WALLCELLNODE_2[2][triwall];
									xp1 = WALLNODECOORD_2[0][n1];
									yp1 = WALLNODECOORD_2[1][n1];
									zp1 = WALLNODECOORD_2[2][n1];
									
									distance = (xp-xp1)*WALLCELLNORMALS_2[0][triwall] + (yp-yp1)*WALLCELLNORMALS_2[1][triwall] + (zp-zp1)*WALLCELLNORMALS_2[2][triwall];
									
									if(distance < 0.25){
										if(fabsf(distance) < 1.0){
											ReflectProj = 0;
											repulf = (0.25 - distance)/1.0;
											if(distance < -0.75) repulf = 1.;
											
											for(v=0;v<=2;v++){
												NODE_WALLREPUL[v][caps][node] += repulf*(1.e-04)*WALLCELLNORMALS_2[v][triwall];
											}
										}
									}
								}
								
							}
						
							for(node=1;node<=NodeM_N;node++){
								
								///interpolate velocities on membrane from surrounding fluid
								xp = NODECOORD_N[0][caps][node]; yp = NODECOORD_N[1][caps][node]; zp = NODECOORD_N[2][caps][node];
								
								i = (int)((xp-X0_2)/DX + 1.01);
								j = (int)((yp-Y0_2)/DX + 1.01);
								k = (int)((zp-Z0_2)/DX + 1.01);
								
								
								i1 = i-4; i2 = i+4;
								j1 = j-4; j2 = j+4;
								k1 = k-4; k2 = k+4;
								
								if(i1 < 1) i1 = 1; if(i2> IM_2) i2 = IM_2;
								if(j1 < 1) j1 = 1; if(j2> JM_2) j2 = JM_2;
								weightsum = 0; weightcount = 0; bvalue = 0;
								for(i=i1;i<=i2;i++){
									for(j=j1;j<=j2;j++){
										for(k=k1;k<=k2;k++){
											xc = X0_2 + (i-1)*DX;
											yc = Y0_2 + (j-1)*DX;
											zc = Z0_2 + (k-1)*DX;
											
											if(k < 1) kmap = k + (KM_2-1);
											else if(k > KM_2) kmap = k - (KM_2-1);
											else kmap = k;
											
											o = DOMID_2[i][j][kmap];
											
											if(B_index_2[i][j][kmap] > 0){
												
												distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
												if(distance <= StencilWidth*DX){
													weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
													NODE_VEL_N[0][caps][node] += weightsum*(VEL_SCALE*U_2[o]);
													NODE_VEL_N[1][caps][node] += weightsum*(VEL_SCALE*V_2[o]);
													NODE_VEL_N[2][caps][node] += weightsum*(VEL_SCALE*W_2[o]);
												}
											}
											else if(W_index_2[i][j][kmap] > 0){
												distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
												if(distance <= WallStencilWidth*DX){
													weightsum = (1.+cosf(PI*fabsf(xc-xp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(WallStencilWidth*DX)))/((WallStencilWidth*2)*(WallStencilWidth*2)*(WallStencilWidth*2));
													NODE_VEL_N[0][caps][node] += weightsum*(VEL_SCALE*Uvoid_2[i][j][kmap]);
													NODE_VEL_N[1][caps][node] += weightsum*(VEL_SCALE*Vvoid_2[i][j][kmap]);
													NODE_VEL_N[2][caps][node] += weightsum*(VEL_SCALE*Wvoid_2[i][j][kmap]);
													
												}
											}
										}
									}
								}
								
								
								
							}
						
						}
						//PCV2 parent
						if(caps >= 330 && caps <= 344){
							for(node=1;node<=NodeM;node++){
								
								///interpolate velocities on membrane from surrounding fluid
								xp = NODECOORD[0][caps][node]; yp = NODECOORD[1][caps][node]; zp = NODECOORD[2][caps][node];
								
								i = (int)((xp-X0_3)/DX + 1.01);
								j = (int)((yp-Y0_3)/DX + 1.01);
								k = (int)((zp-Z0_3)/DX + 1.01);
								
								
								i1 = i-4; i2 = i+4;
								j1 = j-4; j2 = j+4;
								k1 = k-4; k2 = k+4;
								
								if(i1 < 1) i1 = 1; if(i2> IM_3) i2 = IM_3;
								if(j1 < 1) j1 = 1; if(j2> JM_3) j2 = JM_3;
								weightsum = 0; weightcount = 0; bvalue = 0; mindis = 1.5; triwall = 0;
								for(i=i1;i<=i2;i++){
									for(j=j1;j<=j2;j++){
										for(k=k1;k<=k2;k++){
											xc = X0_3 + (i-1)*DX;
											yc = Y0_3 + (j-1)*DX;
											zc = Z0_3 + (k-1)*DX;
											
											if(k < 1) kmap = k + (KM_3-1);
											else if(k > KM_3) kmap = k - (KM_3-1);
											else kmap = k;
										
											o = DOMID_3[i][j][kmap];
											
											if(B_index_3[i][j][kmap] > 0){
												
												distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
												if(distance <= StencilWidth*DX){
													weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
													NODE_VEL[0][caps][node] += weightsum*(VEL_SCALE*U_3[o]);
													NODE_VEL[1][caps][node] += weightsum*(VEL_SCALE*V_3[o]);
													NODE_VEL[2][caps][node] += weightsum*(VEL_SCALE*W_3[o]);
												}
											}
											else if(W_index_3[i][j][kmap] > 0){
												distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
												if(distance <= WallStencilWidth*DX){
													weightsum = (1.+cosf(PI*fabsf(xc-xp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(WallStencilWidth*DX)))/((WallStencilWidth*2)*(WallStencilWidth*2)*(WallStencilWidth*2));
													NODE_VEL[0][caps][node] += weightsum*(VEL_SCALE*Uvoid_3[i][j][kmap]);
													NODE_VEL[1][caps][node] += weightsum*(VEL_SCALE*Vvoid_3[i][j][kmap]);
													NODE_VEL[2][caps][node] += weightsum*(VEL_SCALE*Wvoid_3[i][j][kmap]);
													
												}
											}
											
											if(W_index_3[i][j][kmap] > 0){
												
												nodewall = W_index_3[i][j][kmap];//bvalue is the BC list index number
												
												neighmax = WALLNODETOCELLNEIGHBORLIST_3[0][nodewall];
												
												for(l=1;l<=neighmax;l++){
													tri = WALLNODETOCELLNEIGHBORLIST_3[l][nodewall];
													n1 = WALLCELLNODE_3[0][tri]; n2 = WALLCELLNODE_3[1][tri]; n3 = WALLCELLNODE_3[2][tri];
													xp1 = WALLNODECOORD_3[0][n1];
													yp1 = WALLNODECOORD_3[1][n1];
													zp1 = WALLNODECOORD_3[2][n1];
													distance = sqrtf((xp-xp1)*(xp-xp1) + (yp-yp1)*(yp-yp1) + (zp-zp1)*(zp-zp1));
													if(distance < mindis){
														mindis = distance;
														triwall = tri;
													}
													
													xp2 = WALLNODECOORD_3[0][n2];
													yp2 = WALLNODECOORD_3[1][n2];
													zp2 = WALLNODECOORD_3[2][n2];
													distance = sqrtf((xp-xp2)*(xp-xp2) + (yp-yp2)*(yp-yp2) + (zp-zp2)*(zp-zp2));
													if(distance < mindis){
														mindis = distance;
														triwall = tri;
													}
													
													xp3 = WALLNODECOORD_3[0][n3];
													yp3 = WALLNODECOORD_3[1][n3];
													zp3 = WALLNODECOORD_3[2][n3];
													distance = sqrtf((xp-xp3)*(xp-xp3) + (yp-yp3)*(yp-yp3) + (zp-zp3)*(zp-zp3));
													if(distance < mindis){
														mindis = distance;
														triwall = tri;
													}
													
													xp4 = WALL_CENTROID_3[0][tri];
													yp4 = WALL_CENTROID_3[1][tri];
													zp4 = WALL_CENTROID_3[2][tri];
													distance = sqrtf((xp-xp4)*(xp-xp4) + (yp-yp4)*(yp-yp4) + (zp-zp4)*(zp-zp4));
													if(distance < mindis){
														mindis = distance;
														triwall = tri;
													}
												}
											}
										}
									}
								}
								
								if(triwall > 0){
									NODE_WALLNEIGH[caps][node] = triwall+3600000;
									n1 = WALLCELLNODE_3[0][triwall]; n2 = WALLCELLNODE_3[1][triwall]; n3 = WALLCELLNODE_3[2][triwall];
									xp1 = WALLNODECOORD_3[0][n1];
									yp1 = WALLNODECOORD_3[1][n1];
									zp1 = WALLNODECOORD_3[2][n1];
									
									distance = (xp-xp1)*WALLCELLNORMALS_3[0][triwall] + (yp-yp1)*WALLCELLNORMALS_3[1][triwall] + (zp-zp1)*WALLCELLNORMALS_3[2][triwall];
									
									if(distance < 0.25){
										if(fabsf(distance) < 1.0){
											ReflectProj = 0;
											repulf = (0.25 - distance)/1.0;
											if(distance < -0.75) repulf = 1.;
											
											for(v=0;v<=2;v++){
												NODE_WALLREPUL[v][caps][node] += repulf*(1.e-04)*WALLCELLNORMALS_3[v][triwall];
											}
										}
									}
								}
							}
						
							for(node=1;node<=NodeM_N;node++){
								
								///interpolate velocities on membrane from surrounding fluid
								xp = NODECOORD_N[0][caps][node]; yp = NODECOORD_N[1][caps][node]; zp = NODECOORD_N[2][caps][node];
								
								i = (int)((xp-X0_3)/DX + 1.01);
								j = (int)((yp-Y0_3)/DX + 1.01);
								k = (int)((zp-Z0_3)/DX + 1.01);
								
								
								i1 = i-4; i2 = i+4;
								j1 = j-4; j2 = j+4;
								k1 = k-4; k2 = k+4;
								
								if(i1 < 1) i1 = 1; if(i2> IM_3) i2 = IM_3;
								if(j1 < 1) j1 = 1; if(j2> JM_3) j2 = JM_3;
								weightsum = 0; weightcount = 0; bvalue = 0;
								for(i=i1;i<=i2;i++){
									for(j=j1;j<=j2;j++){
										for(k=k1;k<=k2;k++){
											xc = X0_3 + (i-1)*DX;
											yc = Y0_3 + (j-1)*DX;
											zc = Z0_3 + (k-1)*DX;
											
											if(k < 1) kmap = k + (KM_3-1);
											else if(k > KM_3) kmap = k - (KM_3-1);
											else kmap = k;
											
											o = DOMID_3[i][j][kmap];
											
											if(B_index_3[i][j][kmap] > 0){
												distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
												if(distance <= StencilWidth*DX){
													weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
													NODE_VEL_N[0][caps][node] += weightsum*(VEL_SCALE*U_3[o]);
													NODE_VEL_N[1][caps][node] += weightsum*(VEL_SCALE*V_3[o]);
													NODE_VEL_N[2][caps][node] += weightsum*(VEL_SCALE*W_3[o]);
												}
											}
											else if(W_index_3[i][j][kmap] > 0){
												distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
												if(distance <= WallStencilWidth*DX){
													weightsum = (1.+cosf(PI*fabsf(xc-xp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(WallStencilWidth*DX)))/((WallStencilWidth*2)*(WallStencilWidth*2)*(WallStencilWidth*2));
													NODE_VEL_N[0][caps][node] += weightsum*(VEL_SCALE*Uvoid_3[i][j][kmap]);
													NODE_VEL_N[1][caps][node] += weightsum*(VEL_SCALE*Vvoid_3[i][j][kmap]);
													NODE_VEL_N[2][caps][node] += weightsum*(VEL_SCALE*Wvoid_3[i][j][kmap]);
													
												}
											}
											
										}
									}
								}
								
							}
						}
					}
					else if(CellStatus[caps] > 2){//Parent cells crossing periodic boundaries
						if(caps <= 188){
							for(node=1;node<=NodeM;node++){
								
								///interpolate velocities on membrane from surrounding fluid
								xp = NODECOORD[0][caps][node]; yp = NODECOORD[1][caps][node]; zp = NODECOORD[2][caps][node];
								
								i = (int)((xp-X0_1)/DX + 1.01);
								j = (int)((yp-Y0_1)/DX + 1.01);
								k = (int)((zp-Z0_1)/DX + 1.01);
											
								if(zp < 0.5*(Z0_1 + ZM_1)){
									zpmap = zp+(ZM_1-Z0_1);
								}
								else zpmap = zp;
																
								i = (int)((xp-X0_1)/DX + 1.01);
								j = (int)((yp-Y0_1)/DX + 1.01);
								k = (int)((zpmap-Z0_1)/DX + 1.01);
								
								i1 = i-4; i2 = i+4;
								j1 = j-4; j2 = j+4;
								k1 = k-4; k2 = k+4;
								if(i1 < 1) i1 = 1; if(i2 > IM_1) i2 = IM_1;
								if(j1 < 1) j1 = 1; if(j2 > JM_1) j2 = JM_1;
								
								weightsum = 0; weightcount = 0; bvalue = 0; mindis = 1.5; triwall = 0;
								for(i=i1;i<=i2;i++){
									for(j=j1;j<=j2;j++){
										for(k=k1;k<=k2;k++){
											xc = X0_1 + (i-1)*DX;
											yc = Y0_1 + (j-1)*DX;
											zc = Z0_1 + (k-1)*DX;
											
											if(k > KM_1){
												kmap = k - (KM_1-1);
											}
											else if(k < 1) kmap = k + (KM_1-1);
											else kmap = k;
											
											o = DOMID_1[i][j][kmap];
											
											if(B_index_1[i][j][kmap] > 0){
												
												distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zpmap)*(zc-zpmap));
												if(distance <= StencilWidth*DX){
													weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
													NODE_VEL[0][caps][node] += weightsum*(VEL_SCALE*U_1[o]);
													NODE_VEL[1][caps][node] += weightsum*(VEL_SCALE*V_1[o]);
													NODE_VEL[2][caps][node] += weightsum*(VEL_SCALE*W_1[o]);
												}
											}
											else if(W_index_1[i][j][kmap] > 0){
												distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
												if(distance <= WallStencilWidth*DX){
													weightsum = (1.+cosf(PI*fabsf(xc-xp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(WallStencilWidth*DX)))/((WallStencilWidth*2)*(WallStencilWidth*2)*(WallStencilWidth*2));
													NODE_VEL[0][caps][node] += weightsum*(VEL_SCALE*Uvoid_1[i][j][kmap]);
													NODE_VEL[1][caps][node] += weightsum*(VEL_SCALE*Vvoid_1[i][j][kmap]);
													NODE_VEL[2][caps][node] += weightsum*(VEL_SCALE*Wvoid_1[i][j][kmap]);
												}
											}
											
											if(W_index_1[i][j][kmap] > 0){
												
												nodewall = W_index_1[i][j][kmap];
												neighmax = WALLNODETOCELLNEIGHBORLIST_1[0][nodewall];
												
												for(l=1;l<=neighmax;l++){
													tri = WALLNODETOCELLNEIGHBORLIST_1[l][nodewall];
													n1 = WALLCELLNODE_1[0][tri]; n2 = WALLCELLNODE_1[1][tri]; n3 = WALLCELLNODE_1[2][tri];
													xp1 = WALLNODECOORD_1[0][n1];
													yp1 = WALLNODECOORD_1[1][n1];
													zp1 = WALLNODECOORD_1[2][n1];
													distance = sqrtf((xp-xp1)*(xp-xp1) + (yp-yp1)*(yp-yp1) + (zp-zp1)*(zp-zp1));
													if(distance < mindis){
														mindis = distance;
														triwall = tri;
													}
													
													xp2 = WALLNODECOORD_1[0][n2];
													yp2 = WALLNODECOORD_1[1][n2];
													zp2 = WALLNODECOORD_1[2][n2];
													distance = sqrtf((xp-xp2)*(xp-xp2) + (yp-yp2)*(yp-yp2) + (zp-zp2)*(zp-zp2));
													if(distance < mindis){
														mindis = distance;
														triwall = tri;
													}
													
													xp3 = WALLNODECOORD_1[0][n3];
													yp3 = WALLNODECOORD_1[1][n3];
													zp3 = WALLNODECOORD_1[2][n3];
													distance = sqrtf((xp-xp3)*(xp-xp3) + (yp-yp3)*(yp-yp3) + (zp-zp3)*(zp-zp3));
													if(distance < mindis){
														mindis = distance;
														triwall = tri;
													}
													
													xp4 = WALL_CENTROID_1[0][tri];
													yp4 = WALL_CENTROID_1[1][tri];
													zp4 = WALL_CENTROID_1[2][tri];
													distance = sqrtf((xp-xp4)*(xp-xp4) + (yp-yp4)*(yp-yp4) + (zp-zp4)*(zp-zp4));
													if(distance < mindis){
														mindis = distance;
														triwall = tri;
													}
															
												}
											}
										}
									}
								}
								
								if(triwall > 0){
									NODE_WALLNEIGH[caps][node] = triwall+2400000;
									n1 = WALLCELLNODE_1[0][triwall]; n2 = WALLCELLNODE_1[1][triwall]; n3 = WALLCELLNODE_1[2][triwall];
									xp1 = WALLNODECOORD_1[0][n1];
									yp1 = WALLNODECOORD_1[1][n1];
									zp1 = WALLNODECOORD_1[2][n1];
									
									distance = (xp-xp1)*WALLCELLNORMALS_1[0][triwall] + (yp-yp1)*WALLCELLNORMALS_1[1][triwall] + (zp-zp1)*WALLCELLNORMALS_1[2][triwall];
									
									if(distance < 0.25){
										if(fabsf(distance) < 1.0){
											ReflectProj = 0;
											repulf = (0.25 - distance)/1.0;
											if(distance < -0.75) repulf = 1.;
											
											for(v=0;v<=2;v++){
												NODE_WALLREPUL[v][caps][node] += repulf*(1.e-04)*WALLCELLNORMALS_1[v][triwall];
											}
										}
									}
								}
								
							}
							
							for(node=1;node<=NodeM_N;node++){
								
								///interpolate velocities on membrane from surrounding fluid
								xp = NODECOORD_N[0][caps][node]; yp = NODECOORD_N[1][caps][node]; zp = NODECOORD_N[2][caps][node];
								
								if(zp < 0.5*(Z0_1 + ZM_1)){
									zpmap = zp+(ZM_1-Z0_1);
								}
								else zpmap = zp;
																
								i = (int)((xp-X0_1)/DX + 1.01);
								j = (int)((yp-Y0_1)/DX + 1.01);
								k = (int)((zpmap-Z0_1)/DX + 1.01);
								
								i1 = i-4; i2 = i+4;
								j1 = j-4; j2 = j+4;
								k1 = k-4; k2 = k+4;
								if(i1 < 1) i1 = 1; if(i2 > IM_1) i2 = IM_1;
								if(j1 < 1) j1 = 1; if(j2 > JM_1) j2 = JM_1;
								
								weightsum = 0; weightcount = 0; bvalue = 0;
								for(i=i1;i<=i2;i++){
									for(j=j1;j<=j2;j++){
										for(k=k1;k<=k2;k++){
											xc = X0_1 + (i-1)*DX;
											yc = Y0_1 + (j-1)*DX;
											zc = Z0_1 + (k-1)*DX;
											
											if(k > KM_1) kmap = k - (KM_1-1);
											else if(k < 1) kmap = k + (KM_1-1);
											else kmap = k;
											
											o = DOMID_1[i][j][kmap];
											
											if(B_index_1[i][j][kmap] > 0){
												distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zpmap)*(zc-zpmap));
												
												if(distance <= StencilWidth*DX){
													weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
													NODE_VEL_N[0][caps][node] += weightsum*(VEL_SCALE*U_1[o]);
													NODE_VEL_N[1][caps][node] += weightsum*(VEL_SCALE*V_1[o]);
													NODE_VEL_N[2][caps][node] += weightsum*(VEL_SCALE*W_1[o]);
												}
											}
											else if(W_index_1[i][j][kmap] > 0){
												distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
												if(distance <= WallStencilWidth*DX){
													weightsum = (1.+cosf(PI*fabsf(xc-xp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(WallStencilWidth*DX)))/((WallStencilWidth*2)*(WallStencilWidth*2)*(WallStencilWidth*2));
													NODE_VEL_N[0][caps][node] += weightsum*(VEL_SCALE*Uvoid_1[i][j][kmap]);
													NODE_VEL_N[1][caps][node] += weightsum*(VEL_SCALE*Vvoid_1[i][j][kmap]);
													NODE_VEL_N[2][caps][node] += weightsum*(VEL_SCALE*Wvoid_1[i][j][kmap]);
												}
											}
										}
									}
								}
								
								
							}
							 
						}
						if(caps >= 189 && caps <= 329){
							for(node=1;node<=NodeM;node++){
								
								///interpolate velocities on membrane from surrounding fluid
								xp = NODECOORD[0][caps][node]; yp = NODECOORD[1][caps][node]; zp = NODECOORD[2][caps][node];
								
								i = (int)((xp-X0_2)/DX + 1.01);
								j = (int)((yp-Y0_2)/DX + 1.01);
								k = (int)((zp-Z0_2)/DX + 1.01);
								
								if(zp > 0.5*(Z0_2 + ZM_2)){
									zpmap = zp-(ZM_2-Z0_2);
								}
								else zpmap = zp;
								
								i = (int)((xp-X0_2)/DX + 1.01);
								j = (int)((yp-Y0_2)/DX + 1.01);
								k = (int)((zpmap-Z0_2)/DX + 1.01);
								
								i1 = i-4; i2 = i+4;
								j1 = j-4; j2 = j+4;
								k1 = k-4; k2 = k+4;
								if(i1 < 1) i1 = 1; if(i2 > IM_2) i2 = IM_2;
								if(j1 < 1) j1 = 1; if(j2 > JM_2) j2 = JM_2;
								
								weightsum = 0; weightcount = 0; bvalue = 0; mindis = 1.5; triwall = 0;
								for(i=i1;i<=i2;i++){
									for(j=j1;j<=j2;j++){
										for(k=k1;k<=k2;k++){
											xc = X0_2 + (i-1)*DX;
											yc = Y0_2 + (j-1)*DX;
											zc = Z0_2 + (k-1)*DX;
											
											if(k < 1) kmap = k + (KM_2-1);
											else if(k > KM_2) kmap = k - (KM_2-1);
											else kmap = k;
											
											o = DOMID_2[i][j][kmap];
											
											if(B_index_2[i][j][kmap] > 0){
												distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zpmap)*(zc-zpmap));
												
												if(distance <= StencilWidth*DX){
													weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
													NODE_VEL[0][caps][node] += weightsum*(VEL_SCALE*U_2[o]);
													NODE_VEL[1][caps][node] += weightsum*(VEL_SCALE*V_2[o]);
													NODE_VEL[2][caps][node] += weightsum*(VEL_SCALE*W_2[o]);
												}
											}
											else if(W_index_2[i][j][kmap] > 0){
												distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
												if(distance <= WallStencilWidth*DX){
													weightsum = (1.+cosf(PI*fabsf(xc-xp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(WallStencilWidth*DX)))/((WallStencilWidth*2)*(WallStencilWidth*2)*(WallStencilWidth*2));
													NODE_VEL[0][caps][node] += weightsum*(VEL_SCALE*Uvoid_2[i][j][kmap]);
													NODE_VEL[1][caps][node] += weightsum*(VEL_SCALE*Vvoid_2[i][j][kmap]);
													NODE_VEL[2][caps][node] += weightsum*(VEL_SCALE*Wvoid_2[i][j][kmap]);
												}
											}
											
											if(W_index_2[i][j][kmap] > 0){
												nodewall = W_index_2[i][j][kmap];
												
												neighmax = WALLNODETOCELLNEIGHBORLIST_2[0][nodewall];
														
												for(l=1;l<=neighmax;l++){
													tri = WALLNODETOCELLNEIGHBORLIST_2[l][nodewall];
													n1 = WALLCELLNODE_2[0][tri]; n2 = WALLCELLNODE_2[1][tri]; n3 = WALLCELLNODE_2[2][tri];
													xp1 = WALLNODECOORD_2[0][n1];
													yp1 = WALLNODECOORD_2[1][n1];
													zp1 = WALLNODECOORD_2[2][n1];
													distance = sqrtf((xp-xp1)*(xp-xp1) + (yp-yp1)*(yp-yp1) + (zp-zp1)*(zp-zp1));
													if(distance < mindis){
														mindis = distance;
														triwall = tri;
													}
													
													xp2 = WALLNODECOORD_2[0][n2];
													yp2 = WALLNODECOORD_2[1][n2];
													zp2 = WALLNODECOORD_2[2][n2];
													distance = sqrtf((xp-xp2)*(xp-xp2) + (yp-yp2)*(yp-yp2) + (zp-zp2)*(zp-zp2));
													if(distance < mindis){
														mindis = distance;
														triwall = tri;
													}
													
													xp3 = WALLNODECOORD_2[0][n3];
													yp3 = WALLNODECOORD_2[1][n3];
													zp3 = WALLNODECOORD_2[2][n3];
													distance = sqrtf((xp-xp3)*(xp-xp3) + (yp-yp3)*(yp-yp3) + (zp-zp3)*(zp-zp3));
													if(distance < mindis){
														mindis = distance;
														triwall = tri;
													}
													
													xp4 = WALL_CENTROID_2[0][tri];
													yp4 = WALL_CENTROID_2[1][tri];
													zp4 = WALL_CENTROID_2[2][tri];
													distance = sqrtf((xp-xp4)*(xp-xp4) + (yp-yp4)*(yp-yp4) + (zp-zp4)*(zp-zp4));
													if(distance < mindis){
														mindis = distance;
														triwall = tri;
													}
															
												}
											}
										}
									}
								}
								
								if(triwall > 0){
									NODE_WALLNEIGH[caps][node] = triwall+3000000;
									n1 = WALLCELLNODE_2[0][triwall]; n2 = WALLCELLNODE_2[1][triwall]; n3 = WALLCELLNODE_2[2][triwall];
									xp1 = WALLNODECOORD_2[0][n1];
									yp1 = WALLNODECOORD_2[1][n1];
									zp1 = WALLNODECOORD_2[2][n1];
									
									distance = (xp-xp1)*WALLCELLNORMALS_2[0][triwall] + (yp-yp1)*WALLCELLNORMALS_2[1][triwall] + (zp-zp1)*WALLCELLNORMALS_2[2][triwall];
									
									if(distance < 0.25){
										if(fabsf(distance) < 1.0){
											ReflectProj = 0;
											repulf = (0.25 - distance)/1.0;
											if(distance < -0.75) repulf = 1.;
											
											for(v=0;v<=2;v++){
												NODE_WALLREPUL[v][caps][node] += repulf*(1.e-04)*WALLCELLNORMALS_2[v][triwall];
											}
										}
									}
								}
							}
							
							for(node=1;node<=NodeM_N;node++){
								
								///interpolate velocities on membrane from surrounding fluid
								xp = NODECOORD_N[0][caps][node]; yp = NODECOORD_N[1][caps][node]; zp = NODECOORD_N[2][caps][node];
								
								
								if(zp > 0.5*(Z0_2 + ZM_2)){
									zpmap = zp-(ZM_2-Z0_2);
								}
								else zpmap = zp;
								
								i = (int)((xp-X0_2)/DX + 1.01);
								j = (int)((yp-Y0_2)/DX + 1.01);
								k = (int)((zpmap-Z0_2)/DX + 1.01);
								
								i1 = i-4; i2 = i+4;
								j1 = j-4; j2 = j+4;
								k1 = k-4; k2 = k+4;
								if(i1 < 1) i1 = 1; if(i2 > IM_2) i2 = IM_2;
								if(j1 < 1) j1 = 1; if(j2 > JM_2) j2 = JM_2;
								
								weightsum = 0; weightcount = 0; bvalue = 0;
								for(i=i1;i<=i2;i++){
									for(j=j1;j<=j2;j++){
										for(k=k1;k<=k2;k++){
											xc = X0_2 + (i-1)*DX;
											yc = Y0_2 + (j-1)*DX;
											zc = Z0_2 + (k-1)*DX;
											
											if(k < 1) kmap = k + (KM_2-1);
											else if(k > KM_2) kmap = k - (KM_2-1);
											else kmap = k;
											
											o = DOMID_2[i][j][kmap];
											
											if(B_index_2[i][j][kmap] > 0){
												distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zpmap)*(zc-zpmap));
												
												if(distance <= StencilWidth*DX){
													weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
													NODE_VEL_N[0][caps][node] += weightsum*(VEL_SCALE*U_2[o]);
													NODE_VEL_N[1][caps][node] += weightsum*(VEL_SCALE*V_2[o]);
													NODE_VEL_N[2][caps][node] += weightsum*(VEL_SCALE*W_2[o]);
												}
											}
											else if(W_index_2[i][j][kmap] > 0){
												distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
												if(distance <= WallStencilWidth*DX){
													weightsum = (1.+cosf(PI*fabsf(xc-xp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(WallStencilWidth*DX)))/((WallStencilWidth*2)*(WallStencilWidth*2)*(WallStencilWidth*2));
													NODE_VEL_N[0][caps][node] += weightsum*(VEL_SCALE*Uvoid_2[i][j][kmap]);
													NODE_VEL_N[1][caps][node] += weightsum*(VEL_SCALE*Vvoid_2[i][j][kmap]);
													NODE_VEL_N[2][caps][node] += weightsum*(VEL_SCALE*Wvoid_2[i][j][kmap]);
												}
											}
										}
									}
								}
								
							}
						}
						if(caps >= 330 && caps <= 344){
							for(node=1;node<=NodeM;node++){
								
								///interpolate velocities on membrane from surrounding fluid
								triwall = 0;
								xp = NODECOORD[0][caps][node]; yp = NODECOORD[1][caps][node]; zp = NODECOORD[2][caps][node];
								
								i = (int)((xp-X0_3)/DX + 1.01);
								j = (int)((yp-Y0_3)/DX + 1.01);
								k = (int)((zp-Z0_3)/DX + 1.01);
								
								if(zp > 0.5*(Z0_3 + ZM_3)){
									zpmap = zp-(ZM_3-Z0_3);
								}
								else zpmap = zp;
								
								i = (int)((xp-X0_3)/DX + 1.01);
								j = (int)((yp-Y0_3)/DX + 1.01);
								k = (int)((zpmap-Z0_3)/DX + 1.01);
								
								i1 = i-4; i2 = i+4;
								j1 = j-4; j2 = j+4;
								k1 = k-4; k2 = k+4;
								if(i1 < 1) i1 = 1; if(i2 > IM_3) i2 = IM_3;
								if(j1 < 1) j1 = 1; if(j2 > JM_3) j2 = JM_3;
								
								weightsum = 0; weightcount = 0; bvalue = 0; mindis = 5000.;
								for(i=i1;i<=i2;i++){
									for(j=j1;j<=j2;j++){
										for(k=k1;k<=k2;k++){
											xc = X0_3 + (i-1)*DX;
											yc = Y0_3 + (j-1)*DX;
											zc = Z0_3 + (k-1)*DX;
											
											if(k < 1) kmap = k + (KM_3-1);
											else if(k > KM_3) kmap = k - (KM_3-1);
											else kmap = k;
											
											o = DOMID_3[i][j][kmap];
											
											if(B_index_3[i][j][kmap] > 0){
												distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zpmap)*(zc-zpmap));
												
												if(distance <= StencilWidth*DX){
													weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
													NODE_VEL[0][caps][node] += weightsum*(VEL_SCALE*U_3[o]);
													NODE_VEL[1][caps][node] += weightsum*(VEL_SCALE*V_3[o]);
													NODE_VEL[2][caps][node] += weightsum*(VEL_SCALE*W_3[o]);
												}
											}
											else if(W_index_3[i][j][kmap] > 0){
												distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
												if(distance <= WallStencilWidth*DX){
													weightsum = (1.+cosf(PI*fabsf(xc-xp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(WallStencilWidth*DX)))/((WallStencilWidth*2)*(WallStencilWidth*2)*(WallStencilWidth*2));
													NODE_VEL[0][caps][node] += weightsum*(VEL_SCALE*Uvoid_3[i][j][kmap]);
													NODE_VEL[1][caps][node] += weightsum*(VEL_SCALE*Vvoid_3[i][j][kmap]);
													NODE_VEL[2][caps][node] += weightsum*(VEL_SCALE*Wvoid_3[i][j][kmap]);
												}
											}
											
											if(W_index_3[i][j][kmap] > 0){
												
												nodewall = W_index_3[i][j][kmap];
												
												neighmax = WALLNODETOCELLNEIGHBORLIST_3[0][nodewall];
														
												for(l=1;l<=neighmax;l++){
													tri = WALLNODETOCELLNEIGHBORLIST_3[l][nodewall];
													n1 = WALLCELLNODE_3[0][tri]; n2 = WALLCELLNODE_3[1][tri]; n3 = WALLCELLNODE_3[2][tri];
													xp1 = WALLNODECOORD_3[0][n1];
													yp1 = WALLNODECOORD_3[1][n1];
													zp1 = WALLNODECOORD_3[2][n1];
													distance = sqrtf((xp-xp1)*(xp-xp1) + (yp-yp1)*(yp-yp1) + (zp-zp1)*(zp-zp1));
													if(distance < mindis){
														mindis = distance;
														triwall = tri;
													}
													
													xp2 = WALLNODECOORD_3[0][n2];
													yp2 = WALLNODECOORD_3[1][n2];
													zp2 = WALLNODECOORD_3[2][n2];
													distance = sqrtf((xp-xp2)*(xp-xp2) + (yp-yp2)*(yp-yp2) + (zp-zp2)*(zp-zp2));
													if(distance < mindis){
														mindis = distance;
														triwall = tri;
													}
													
													xp3 = WALLNODECOORD_3[0][n3];
													yp3 = WALLNODECOORD_3[1][n3];
													zp3 = WALLNODECOORD_3[2][n3];
													distance = sqrtf((xp-xp3)*(xp-xp3) + (yp-yp3)*(yp-yp3) + (zp-zp3)*(zp-zp3));
													if(distance < mindis){
														mindis = distance;
														triwall = tri;
													}
													
													xp4 = WALL_CENTROID_3[0][tri];
													yp4 = WALL_CENTROID_3[1][tri];
													zp4 = WALL_CENTROID_3[2][tri];
													distance = sqrtf((xp-xp4)*(xp-xp4) + (yp-yp4)*(yp-yp4) + (zp-zp4)*(zp-zp4));
													if(distance < mindis){
														mindis = distance;
														triwall = tri;
													}
													
												}
											}
										}
									}
								}
								
								if(triwall > 0){
									NODE_WALLNEIGH[caps][node] = triwall+3600000;
									n1 = WALLCELLNODE_3[0][triwall]; n2 = WALLCELLNODE_3[1][triwall]; n3 = WALLCELLNODE_3[2][triwall];
									xp1 = WALLNODECOORD_3[0][n1];
									yp1 = WALLNODECOORD_3[1][n1];
									zp1 = WALLNODECOORD_3[2][n1];
									
									distance = (xp-xp1)*WALLCELLNORMALS_3[0][triwall] + (yp-yp1)*WALLCELLNORMALS_3[1][triwall] + (zp-zp1)*WALLCELLNORMALS_3[2][triwall];
									
									if(distance < 0.25){
										if(fabsf(distance) < 1.0){
											ReflectProj = 0;
											repulf = (0.25 - distance)/1.0;
											if(distance < -0.75) repulf = 1.;
											
											for(v=0;v<=2;v++){
												NODE_WALLREPUL[v][caps][node] += repulf*(1.e-04)*WALLCELLNORMALS_3[v][triwall];
											}
										}
									}
								}	
							}
														
							for(node=1;node<=NodeM_N;node++){
								
								///interpolate velocities on membrane from surrounding fluid
								xp = NODECOORD_N[0][caps][node]; yp = NODECOORD_N[1][caps][node]; zp = NODECOORD_N[2][caps][node];
								
								if(zp > 0.5*(Z0_3 + ZM_3)){
									zpmap = zp-(ZM_3-Z0_3);
								}
								else zpmap = zp;
								
								i = (int)((xp-X0_3)/DX + 1.01);
								j = (int)((yp-Y0_3)/DX + 1.01);
								k = (int)((zpmap-Z0_3)/DX + 1.01);
								
								i1 = i-4; i2 = i+4;
								j1 = j-4; j2 = j+4;
								k1 = k-4; k2 = k+4;
								if(i1 < 1) i1 = 1; if(i2 > IM_3) i2 = IM_3;
								if(j1 < 1) j1 = 1; if(j2 > JM_3) j2 = JM_3;
								
								weightsum = 0; weightcount = 0; bvalue = 0;
								for(i=i1;i<=i2;i++){
									for(j=j1;j<=j2;j++){
										for(k=k1;k<=k2;k++){
											xc = X0_3 + (i-1)*DX;
											yc = Y0_3 + (j-1)*DX;
											zc = Z0_3 + (k-1)*DX;
											
											validinter[i-i1][j-j1][k-k1] = 0;
											
											if(k < 1) kmap = k + (KM_3-1);
											else if(k > KM_3) kmap = k - (KM_3-1);
											else kmap = k;
											
											o = DOMID_3[i][j][kmap];
											
											if(B_index_3[i][j][kmap] > 0){
												distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zpmap)*(zc-zpmap));
												
												if(distance <= StencilWidth*DX){
													weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
													NODE_VEL_N[0][caps][node] += weightsum*(VEL_SCALE*U_3[o]);
													NODE_VEL_N[1][caps][node] += weightsum*(VEL_SCALE*V_3[o]);
													NODE_VEL_N[2][caps][node] += weightsum*(VEL_SCALE*W_3[o]);
												}
											}
											else if(W_index_3[i][j][kmap] > 0){
												distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
												if(distance <= WallStencilWidth*DX){
													weightsum = (1.+cosf(PI*fabsf(xc-xp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(WallStencilWidth*DX)))/((WallStencilWidth*2)*(WallStencilWidth*2)*(WallStencilWidth*2));
													NODE_VEL_N[0][caps][node] += weightsum*(VEL_SCALE*Uvoid_3[i][j][kmap]);
													NODE_VEL_N[1][caps][node] += weightsum*(VEL_SCALE*Vvoid_3[i][j][kmap]);
													NODE_VEL_N[2][caps][node] += weightsum*(VEL_SCALE*Wvoid_3[i][j][kmap]);
												}
											}
										}
									}
								}
								
							}
						}
					}
				}
			}
		}
		
		///RBC velocity interpolation and Wall repulsion velocity
		/// ROI cells
		
		if(RBCSOLVERFLAG == 1){
			#pragma omp parallel private(ReflectProj,o,wneighmax,nodewall,triwall,region,n1,n2,n3,xp1,yp1,zp1,xp2,yp2,zp2,xp3,yp3,zp3,xp4,yp4,zp4,xi,yi,zi,r23,ri3,r21,ri1,r31,area1,area2,area3,aratio,v,l,neighmax,tri,mindis2,xc2,yc2,zc2,xc1,yc1,zc1,iii,jjj,kkk,mindis,bvalue,dotpdt,m,caps,node,xp,yp,zp,i,j,k,i1,j1,k1,i2,j2,k2,weightsum,validinter,xc,yc,zc,ii,jj,kk,distance,interweight,repulf,ParentID,check) 
			{
				#pragma omp for schedule(static,50)  nowait
				for(caps=345;caps<=CapsM;caps++){
					if(CellStatus[caps] > 0){//Cloned ROI cells released for independent movement
						for(node=1;node<=NodeM;node++){
							
							///interpolate velocities on membrane from surrounding fluid
							xp = NODECOORD[0][caps][node]; yp = NODECOORD[1][caps][node]; zp = NODECOORD[2][caps][node];
							
							if(CellStatus[caps] == 1 || (zp > ZM_1 && zp < Z0_2)){
								i = (int)((xp-X0)/DX + 1.01);
								j = (int)((yp-Y0)/DX + 1.01);
								k = (int)((zp-Z0)/DX + 1.01);
								
								i1 = i-4; i2 = i+4;
								j1 = j-4; j2 = j+4;
								k1 = k-4; k2 = k+4;
								
								if(i1 < 1) i1 = 1; if(i2> IM) i2 = IM;
								if(j1 < 1) j1 = 1; if(j2> JM) j2 = JM;
								
								weightsum = 0; bvalue = 0; mindis = 1.5; region = 5000; triwall = 0;
								for(i=i1;i<=i2;i++){
									for(j=j1;j<=j2;j++){
										for(k=k1;k<=k2;k++){
											xc = X0 + (i-1)*DX;
											yc = Y0 + (j-1)*DX;
											zc = Z0 + (k-1)*DX;
											
											if(k < 3 && REGRESSFLAG[caps] >= -188 && REGRESSFLAG[caps] < 0){
												ii = (int)((xc-X0_1)/DX + 1.01);
												jj = (int)((yc-Y0_1)/DX + 1.01);
												kk = (int)((zc-Z0_1)/DX + 1.01);
												
												o = DOMID_1[ii][jj][kk];
												
												if(B_index_1[ii][jj][kk] > 0){
													distance = sqrtf((xc-xp)*(xc-xp) + (yc-yp)*(yc-yp) + (zc-zp)*(zc-zp));
													if(distance <= StencilWidth*DX){
														weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
														NODE_VEL[0][caps][node] += weightsum*(VEL_SCALE*U_1[o]);
														NODE_VEL[1][caps][node] += weightsum*(VEL_SCALE*V_1[o]);
														NODE_VEL[2][caps][node] += weightsum*(VEL_SCALE*W_1[o]);
													}
												}
												else if(W_index_1[ii][jj][kk] > 0){
													distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
													if(distance <= WallStencilWidth*DX){
														weightsum = (1.+cosf(PI*fabsf(xc-xp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(WallStencilWidth*DX)))/((WallStencilWidth*2)*(WallStencilWidth*2)*(WallStencilWidth*2));
														NODE_VEL[0][caps][node] += weightsum*(VEL_SCALE*Uvoid_1[ii][jj][kk]);
														NODE_VEL[1][caps][node] += weightsum*(VEL_SCALE*Vvoid_1[ii][jj][kk]);
														NODE_VEL[2][caps][node] += weightsum*(VEL_SCALE*Wvoid_1[ii][jj][kk]);
													}
												}
												
												
												if(W_index_1[ii][jj][kk] > 0){
													nodewall = W_index_1[ii][jj][kk];
													
													neighmax = WALLNODETOCELLNEIGHBORLIST_1[0][nodewall];
													
													if(nodewall <= NodeW_1){	
														for(l=1;l<=neighmax;l++){
															tri = WALLNODETOCELLNEIGHBORLIST_1[l][nodewall];
															n1 = WALLCELLNODE_1[0][tri]; n2 = WALLCELLNODE_1[1][tri]; n3 = WALLCELLNODE_1[2][tri];
															xp1 = WALLNODECOORD_1[0][n1];
															yp1 = WALLNODECOORD_1[1][n1];
															zp1 = WALLNODECOORD_1[2][n1];
															distance = sqrtf((xp-xp1)*(xp-xp1) + (yp-yp1)*(yp-yp1) + (zp-zp1)*(zp-zp1));
															if(distance < mindis){
																region = 1;
																mindis = distance;
																triwall = tri;
															}
															
															xp2 = WALLNODECOORD_1[0][n2];
															yp2 = WALLNODECOORD_1[1][n2];
															zp2 = WALLNODECOORD_1[2][n2];
															distance = sqrtf((xp-xp2)*(xp-xp2) + (yp-yp2)*(yp-yp2) + (zp-zp2)*(zp-zp2));
															if(distance < mindis){
																region = 1;
																mindis = distance;
																triwall = tri;
															}
															
															xp3 = WALLNODECOORD_1[0][n3];
															yp3 = WALLNODECOORD_1[1][n3];
															zp3 = WALLNODECOORD_1[2][n3];
															distance = sqrtf((xp-xp3)*(xp-xp3) + (yp-yp3)*(yp-yp3) + (zp-zp3)*(zp-zp3));
															if(distance < mindis){
																region = 1;
																mindis = distance;
																triwall = tri;
															}
															
															xp4 = WALL_CENTROID_1[0][tri];
															yp4 = WALL_CENTROID_1[1][tri];
															zp4 = WALL_CENTROID_1[2][tri];
															distance = sqrtf((xp-xp4)*(xp-xp4) + (yp-yp4)*(yp-yp4) + (zp-zp4)*(zp-zp4));
															if(distance < mindis){
																region = 1;
																mindis = distance;
																triwall = tri;
															}
															
														}
													}
												}
											}
											else if(k > KM-3 && REGRESSFLAG[caps] >= -329 && REGRESSFLAG[caps] <= -189){
												ii = (int)((xc-X0_2)/DX + 1.01);
												jj = (int)((yc-Y0_2)/DX + 1.01);
												kk = (int)((zc-Z0_2)/DX + 1.01);
												
												o = DOMID_2[ii][jj][kk];
												
												if(B_index_2[ii][jj][kk] > 0){
													distance = sqrtf((xc-xp)*(xc-xp) + (yc-yp)*(yc-yp) + (zc-zp)*(zc-zp));
													if(distance <= StencilWidth*DX){
														weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
														NODE_VEL[0][caps][node] += weightsum*(VEL_SCALE*U_2[o]);
														NODE_VEL[1][caps][node] += weightsum*(VEL_SCALE*V_2[o]);
														NODE_VEL[2][caps][node] += weightsum*(VEL_SCALE*W_2[o]);
													}
												}
												else if(W_index_2[ii][jj][kk] > 0){
													distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
													if(distance <= WallStencilWidth*DX){
														weightsum = (1.+cosf(PI*fabsf(xc-xp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(WallStencilWidth*DX)))/((WallStencilWidth*2)*(WallStencilWidth*2)*(WallStencilWidth*2));
														NODE_VEL[0][caps][node] += weightsum*(VEL_SCALE*Uvoid_2[ii][jj][kk]);
														NODE_VEL[1][caps][node] += weightsum*(VEL_SCALE*Vvoid_2[ii][jj][kk]);
														NODE_VEL[2][caps][node] += weightsum*(VEL_SCALE*Wvoid_2[ii][jj][kk]);
													}
												}
												
												if(W_index_2[ii][jj][kk] > 0){
													nodewall = W_index_2[ii][jj][kk];
													
													neighmax = WALLNODETOCELLNEIGHBORLIST_2[0][nodewall];
													
													if(nodewall <= NodeW_2){		
														for(l=1;l<=neighmax;l++){
															tri = WALLNODETOCELLNEIGHBORLIST_2[l][nodewall];
															n1 = WALLCELLNODE_2[0][tri]; n2 = WALLCELLNODE_2[1][tri]; n3 = WALLCELLNODE_2[2][tri];
															xp1 = WALLNODECOORD_2[0][n1];
															yp1 = WALLNODECOORD_2[1][n1];
															zp1 = WALLNODECOORD_2[2][n1];
															distance = sqrtf((xp-xp1)*(xp-xp1) + (yp-yp1)*(yp-yp1) + (zp-zp1)*(zp-zp1));
															if(distance < mindis){
																region = 2;
																mindis = distance;
																triwall = tri;
															}
															
															xp2 = WALLNODECOORD_2[0][n2];
															yp2 = WALLNODECOORD_2[1][n2];
															zp2 = WALLNODECOORD_2[2][n2];
															distance = sqrtf((xp-xp2)*(xp-xp2) + (yp-yp2)*(yp-yp2) + (zp-zp2)*(zp-zp2));
															if(distance < mindis){
																region = 2;
																mindis = distance;
																triwall = tri;
															}
															
															xp3 = WALLNODECOORD_2[0][n3];
															yp3 = WALLNODECOORD_2[1][n3];
															zp3 = WALLNODECOORD_2[2][n3];
															distance = sqrtf((xp-xp3)*(xp-xp3) + (yp-yp3)*(yp-yp3) + (zp-zp3)*(zp-zp3));
															if(distance < mindis){
																region = 2;
																mindis = distance;
																triwall = tri;
															}
															
															xp4 = WALL_CENTROID_2[0][tri];
															yp4 = WALL_CENTROID_2[1][tri];
															zp4 = WALL_CENTROID_2[2][tri];
															distance = sqrtf((xp-xp4)*(xp-xp4) + (yp-yp4)*(yp-yp4) + (zp-zp4)*(zp-zp4));
															if(distance < mindis){
																region = 2;
																mindis = distance;
																triwall = tri;
															}
															
														}
													}
												}
											}
											else if(k > KM-3 && REGRESSFLAG[caps] >= -344 && REGRESSFLAG[caps] <= -330){
												ii = (int)((xc-X0_3)/DX + 1.01);
												jj = (int)((yc-Y0_3)/DX + 1.01);
												kk = (int)((zc-Z0_3)/DX + 1.01);
												
												o = DOMID_3[ii][jj][kk];
												
												if(B_index_3[ii][jj][kk] > 0){
													distance = sqrtf((xc-xp)*(xc-xp) + (yc-yp)*(yc-yp) + (zc-zp)*(zc-zp));
													if(distance <= StencilWidth*DX){
														weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
														NODE_VEL[0][caps][node] += weightsum*(VEL_SCALE*U_3[o]);
														NODE_VEL[1][caps][node] += weightsum*(VEL_SCALE*V_3[o]);
														NODE_VEL[2][caps][node] += weightsum*(VEL_SCALE*W_3[o]);
													}
												}
												else if(W_index_3[ii][jj][kk] > 0){
													distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
													if(distance <= WallStencilWidth*DX){
														weightsum = (1.+cosf(PI*fabsf(xc-xp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(WallStencilWidth*DX)))/((WallStencilWidth*2)*(WallStencilWidth*2)*(WallStencilWidth*2));
														NODE_VEL[0][caps][node] += weightsum*(VEL_SCALE*Uvoid_3[ii][jj][kk]);
														NODE_VEL[1][caps][node] += weightsum*(VEL_SCALE*Vvoid_3[ii][jj][kk]);
														NODE_VEL[2][caps][node] += weightsum*(VEL_SCALE*Wvoid_3[ii][jj][kk]);
													}
												}
												
												if(W_index_3[ii][jj][kk] > 0){
													nodewall = W_index_3[ii][jj][kk];
													
													neighmax = WALLNODETOCELLNEIGHBORLIST_3[0][nodewall];
													
													if(nodewall <= NodeW_3){		
														for(l=1;l<=neighmax;l++){
															tri = WALLNODETOCELLNEIGHBORLIST_3[l][nodewall];
															n1 = WALLCELLNODE_3[0][tri]; n2 = WALLCELLNODE_3[1][tri]; n3 = WALLCELLNODE_3[2][tri];
															xp1 = WALLNODECOORD_3[0][n1];
															yp1 = WALLNODECOORD_3[1][n1];
															zp1 = WALLNODECOORD_3[2][n1];
															distance = sqrtf((xp-xp1)*(xp-xp1) + (yp-yp1)*(yp-yp1) + (zp-zp1)*(zp-zp1));
															if(distance < mindis){
																region = 3;
																mindis = distance;
																triwall = tri;
															}
															
															xp2 = WALLNODECOORD_3[0][n2];
															yp2 = WALLNODECOORD_3[1][n2];
															zp2 = WALLNODECOORD_3[2][n2];
															distance = sqrtf((xp-xp2)*(xp-xp2) + (yp-yp2)*(yp-yp2) + (zp-zp2)*(zp-zp2));
															if(distance < mindis){
																region = 3;
																mindis = distance;
																triwall = tri;
															}
															
															xp3 = WALLNODECOORD_3[0][n3];
															yp3 = WALLNODECOORD_3[1][n3];
															zp3 = WALLNODECOORD_3[2][n3];
															distance = sqrtf((xp-xp3)*(xp-xp3) + (yp-yp3)*(yp-yp3) + (zp-zp3)*(zp-zp3));
															if(distance < mindis){
																region = 3;
																mindis = distance;
																triwall = tri;
															}
															
															xp4 = WALL_CENTROID_3[0][tri];
															yp4 = WALL_CENTROID_3[1][tri];
															zp4 = WALL_CENTROID_3[2][tri];
															distance = sqrtf((xp-xp4)*(xp-xp4) + (yp-yp4)*(yp-yp4) + (zp-zp4)*(zp-zp4));
															if(distance < mindis){
																region = 3;
																mindis = distance;
																triwall = tri;
															}
															
														}
													}
												}
											}
											
											if(k >= 3 && k <= KM-3){
												ii = i;
												jj = j;
												kk = k;
												
												o = DOMID[ii][jj][kk];
												
												if(B_index[ii][jj][kk] > 0 && B_index[ii][jj][kk] < 4998){
													distance = sqrtf((xc-xp)*(xc-xp) + (yc-yp)*(yc-yp) + (zc-zp)*(zc-zp));
													if(distance <= StencilWidth*DX){
														weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
														NODE_VEL[0][caps][node] += weightsum*(VEL_SCALE*U[o]);
														NODE_VEL[1][caps][node] += weightsum*(VEL_SCALE*V[o]);
														NODE_VEL[2][caps][node] += weightsum*(VEL_SCALE*W[o]);
													}
												}
												else if(B_index[ii][jj][kk] == 4998){//plugged cell location - this is deprecated - it is no longer used 
													distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
													if(distance <= StencilWidth*DX){
														weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
														NODE_VEL[0][caps][node] += weightsum*(VEL_SCALE*Uvoid[ii][jj][kk]);
														NODE_VEL[1][caps][node] += weightsum*(VEL_SCALE*Vvoid[ii][jj][kk]);
														NODE_VEL[2][caps][node] += weightsum*(VEL_SCALE*Wvoid[ii][jj][kk]);
													}
												}
												else if(W_index[ii][jj][kk] > 0){
													distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
													if(distance <= WallStencilWidth*DX){
														weightsum = (1.+cosf(PI*fabsf(xc-xp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(WallStencilWidth*DX)))/((WallStencilWidth*2)*(WallStencilWidth*2)*(WallStencilWidth*2));
														NODE_VEL[0][caps][node] += weightsum*(VEL_SCALE*Uvoid[ii][jj][kk]);
														NODE_VEL[1][caps][node] += weightsum*(VEL_SCALE*Vvoid[ii][jj][kk]);
														NODE_VEL[2][caps][node] += weightsum*(VEL_SCALE*Wvoid[ii][jj][kk]);
													}
												}
												
												if(W_index[ii][jj][kk] > 0){
													nodewall = W_index[ii][jj][kk];
													
													neighmax = WALLNODETOCELLNEIGHBORLIST[0][nodewall];
													
													if(nodewall <= NodeW){		
														for(l=1;l<=neighmax;l++){
															tri = WALLNODETOCELLNEIGHBORLIST[l][nodewall];
															n1 = WALLCELLNODE[0][tri]; n2 = WALLCELLNODE[1][tri]; n3 = WALLCELLNODE[2][tri];
															xp1 = WALLNODECOORD[0][n1];
															yp1 = WALLNODECOORD[1][n1];
															zp1 = WALLNODECOORD[2][n1];
															distance = sqrtf((xp-xp1)*(xp-xp1) + (yp-yp1)*(yp-yp1) + (zp-zp1)*(zp-zp1));
															if(distance < mindis){
																region = 0;
																mindis = distance;
																triwall = tri;
															}
															
															xp2 = WALLNODECOORD[0][n2];
															yp2 = WALLNODECOORD[1][n2];
															zp2 = WALLNODECOORD[2][n2];
															distance = sqrtf((xp-xp2)*(xp-xp2) + (yp-yp2)*(yp-yp2) + (zp-zp2)*(zp-zp2));
															if(distance < mindis){
																region = 0;
																mindis = distance;
																triwall = tri;
															}
															
															xp3 = WALLNODECOORD[0][n3];
															yp3 = WALLNODECOORD[1][n3];
															zp3 = WALLNODECOORD[2][n3];
															distance = sqrtf((xp-xp3)*(xp-xp3) + (yp-yp3)*(yp-yp3) + (zp-zp3)*(zp-zp3));
															if(distance < mindis){
																region = 0;
																mindis = distance;
																triwall = tri;
															}
															
															xp4 = WALL_CENTROID[0][tri];
															yp4 = WALL_CENTROID[1][tri];
															zp4 = WALL_CENTROID[2][tri];
															distance = sqrtf((xp-xp4)*(xp-xp4) + (yp-yp4)*(yp-yp4) + (zp-zp4)*(zp-zp4));
															if(distance < mindis){
																region = 0;
																mindis = distance;
																triwall = tri;
															}
															
														}
													}
												}
											}
										}
									}
								}
								
								///Calculate wall repulsion from nearby wall triangle elements
								if(triwall > 0){
									
									if(region == 0){
										NODE_WALLNEIGH[caps][node] = triwall;
										n1 = WALLCELLNODE[0][triwall]; n2 = WALLCELLNODE[1][triwall]; n3 = WALLCELLNODE[2][triwall];
										xp1 = WALLNODECOORD[0][n1];
										yp1 = WALLNODECOORD[1][n1];
										zp1 = WALLNODECOORD[2][n1];
										
										distance = (xp-xp1)*WALLCELLNORMALS[0][triwall] + (yp-yp1)*WALLCELLNORMALS[1][triwall] + (zp-zp1)*WALLCELLNORMALS[2][triwall];
										
										if(distance < 0.25){
											if(fabsf(distance) < 1.0){
												ReflectProj = 0;
												repulf = (0.25 - distance)/1.0;
												if(distance < -0.75) repulf = 1.;
												
												for(v=0;v<=2;v++){
													NODE_WALLREPUL[v][caps][node] += repulf*(1.e-04)*WALLCELLNORMALS[v][triwall];		
												}
											}
										}
									}
									if(region == 1){
										NODE_WALLNEIGH[caps][node] = triwall+2400000;
										n1 = WALLCELLNODE_1[0][triwall]; n2 = WALLCELLNODE_1[1][triwall]; n3 = WALLCELLNODE_1[2][triwall];
										xp1 = WALLNODECOORD_1[0][n1];
										yp1 = WALLNODECOORD_1[1][n1];
										zp1 = WALLNODECOORD_1[2][n1];
										
										distance = (xp-xp1)*WALLCELLNORMALS_1[0][triwall] + (yp-yp1)*WALLCELLNORMALS_1[1][triwall] + (zp-zp1)*WALLCELLNORMALS_1[2][triwall];
										
										if(distance < 0.25){
											if(fabsf(distance) < 1.0){
												ReflectProj = 0;
												repulf = (0.25 - distance)/1.0;
												if(distance < -0.75) repulf = 1.;
												
												for(v=0;v<=2;v++){	
													NODE_WALLREPUL[v][caps][node] += repulf*(1.e-04)*WALLCELLNORMALS_1[v][triwall];		
												}
											}
										}
									}
									if(region == 2){
										NODE_WALLNEIGH[caps][node] = triwall+3000000;
										n1 = WALLCELLNODE_2[0][triwall]; n2 = WALLCELLNODE_2[1][triwall]; n3 = WALLCELLNODE_2[2][triwall];
										xp1 = WALLNODECOORD_2[0][n1];
										yp1 = WALLNODECOORD_2[1][n1];
										zp1 = WALLNODECOORD_2[2][n1];
										
										distance = (xp-xp1)*WALLCELLNORMALS_2[0][triwall] + (yp-yp1)*WALLCELLNORMALS_2[1][triwall] + (zp-zp1)*WALLCELLNORMALS_2[2][triwall];
										
										if(distance < 0.25){
											if(fabsf(distance) < 1.0){
												ReflectProj = 0;
												repulf = (0.25 - distance)/1.0;
												if(distance < -0.75) repulf = 1.;
												
												for(v=0;v<=2;v++){
													NODE_WALLREPUL[v][caps][node] += repulf*(1.e-04)*WALLCELLNORMALS_2[v][triwall];
												}
											}
										}
									}
									if(region == 3){
										NODE_WALLNEIGH[caps][node] = triwall+3600000;
										n1 = WALLCELLNODE_3[0][triwall]; n2 = WALLCELLNODE_3[1][triwall]; n3 = WALLCELLNODE_3[2][triwall];
										xp1 = WALLNODECOORD_3[0][n1];
										yp1 = WALLNODECOORD_3[1][n1];
										zp1 = WALLNODECOORD_3[2][n1];
										
										distance = (xp-xp1)*WALLCELLNORMALS_3[0][triwall] + (yp-yp1)*WALLCELLNORMALS_3[1][triwall] + (zp-zp1)*WALLCELLNORMALS_3[2][triwall];
										
										if(distance < 0.25){
											if(fabsf(distance) < 1.0){
												ReflectProj = 0;
												repulf = (0.25 - distance)/1.0;
												if(distance < -0.75) repulf = 1.;
												
												for(v=0;v<=2;v++){
													NODE_WALLREPUL[v][caps][node] += repulf*(1.e-04)*WALLCELLNORMALS_3[v][triwall];	
												}
											}
										}
									}
								}
							}
						}
						
						
						
						
						for(node=1;node<=NodeM_N;node++){
								
							///interpolate velocities on membrane from surrounding fluid
							xp = NODECOORD_N[0][caps][node]; yp = NODECOORD_N[1][caps][node]; zp = NODECOORD_N[2][caps][node];
							
							if(CellStatus[caps] == 1 || (zp > ZM_1 && zp < Z0_2)){
								i = (int)((xp-X0)/DX + 1.01);
								j = (int)((yp-Y0)/DX + 1.01);
								k = (int)((zp-Z0)/DX + 1.01);
								
								
								i1 = i-4; i2 = i+4;
								j1 = j-4; j2 = j+4;
								k1 = k-4; k2 = k+4;
								
								if(i1 < 1) i1 = 1; if(i2> IM) i2 = IM;
								if(j1 < 1) j1 = 1; if(j2> JM) j2 = JM;
								
								weightsum = 0; bvalue = 0; 
								for(i=i1;i<=i2;i++){
									for(j=j1;j<=j2;j++){
										for(k=k1;k<=k2;k++){
											xc = X0 + (i-1)*DX;
											yc = Y0 + (j-1)*DX;
											zc = Z0 + (k-1)*DX;
											
											if(k < 3 && REGRESSFLAG[caps] >= -188 && REGRESSFLAG[caps] < 0){
												ii = (int)((xc-X0_1)/DX + 1.01);
												jj = (int)((yc-Y0_1)/DX + 1.01);
												kk = (int)((zc-Z0_1)/DX + 1.01);
												
												o = DOMID_1[ii][jj][kk];
												
												if(B_index_1[ii][jj][kk] > 0){
													distance = sqrtf((xc-xp)*(xc-xp) + (yc-yp)*(yc-yp) + (zc-zp)*(zc-zp));
													if(distance <= StencilWidth*DX){
														weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
														NODE_VEL_N[0][caps][node] += weightsum*(VEL_SCALE*U_1[o]);
														NODE_VEL_N[1][caps][node] += weightsum*(VEL_SCALE*V_1[o]);
														NODE_VEL_N[2][caps][node] += weightsum*(VEL_SCALE*W_1[o]);
													}
												}
												else if(W_index_1[ii][jj][kk] > 0){
													distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
													if(distance <= WallStencilWidth*DX){
														weightsum = (1.+cosf(PI*fabsf(xc-xp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(WallStencilWidth*DX)))/((WallStencilWidth*2)*(WallStencilWidth*2)*(WallStencilWidth*2));
														NODE_VEL_N[0][caps][node] += weightsum*(VEL_SCALE*Uvoid_1[ii][jj][kk]);
														NODE_VEL_N[1][caps][node] += weightsum*(VEL_SCALE*Vvoid_1[ii][jj][kk]);
														NODE_VEL_N[2][caps][node] += weightsum*(VEL_SCALE*Wvoid_1[ii][jj][kk]);
													}
												}
											}
											else if(k > KM-3 && REGRESSFLAG[caps] >= -329 && REGRESSFLAG[caps] <= -189){
												ii = (int)((xc-X0_2)/DX + 1.01);
												jj = (int)((yc-Y0_2)/DX + 1.01);
												kk = (int)((zc-Z0_2)/DX + 1.01);
												
												o = DOMID_2[ii][jj][kk];
												
												if(B_index_2[ii][jj][kk] > 0){
													distance = sqrtf((xc-xp)*(xc-xp) + (yc-yp)*(yc-yp) + (zc-zp)*(zc-zp));
													if(distance <= StencilWidth*DX){
														weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
														NODE_VEL_N[0][caps][node] += weightsum*(VEL_SCALE*U_2[o]);
														NODE_VEL_N[1][caps][node] += weightsum*(VEL_SCALE*V_2[o]);
														NODE_VEL_N[2][caps][node] += weightsum*(VEL_SCALE*W_2[o]);
													}
												}
												else if(W_index_2[ii][jj][kk] > 0){
													distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
													if(distance <= WallStencilWidth*DX){
														weightsum = (1.+cosf(PI*fabsf(xc-xp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(WallStencilWidth*DX)))/((WallStencilWidth*2)*(WallStencilWidth*2)*(WallStencilWidth*2));
														NODE_VEL_N[0][caps][node] += weightsum*(VEL_SCALE*Uvoid_2[ii][jj][kk]);
														NODE_VEL_N[1][caps][node] += weightsum*(VEL_SCALE*Vvoid_2[ii][jj][kk]);
														NODE_VEL_N[2][caps][node] += weightsum*(VEL_SCALE*Wvoid_2[ii][jj][kk]);
													}
												}
											}
											else if(k > KM-3 && REGRESSFLAG[caps] >= -344 && REGRESSFLAG[caps] <= -330){
												ii = (int)((xc-X0_3)/DX + 1.01);
												jj = (int)((yc-Y0_3)/DX + 1.01);
												kk = (int)((zc-Z0_3)/DX + 1.01);
												
												o = DOMID_3[ii][jj][kk];
												
												if(B_index_3[ii][jj][kk] > 0){
													distance = sqrtf((xc-xp)*(xc-xp) + (yc-yp)*(yc-yp) + (zc-zp)*(zc-zp));
													if(distance <= StencilWidth*DX){
														weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
														NODE_VEL_N[0][caps][node] += weightsum*(VEL_SCALE*U_3[o]);
														NODE_VEL_N[1][caps][node] += weightsum*(VEL_SCALE*V_3[o]);
														NODE_VEL_N[2][caps][node] += weightsum*(VEL_SCALE*W_3[o]);
													}
												}
												else if(W_index_3[ii][jj][kk] > 0){
													distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
													if(distance <= WallStencilWidth*DX){
														weightsum = (1.+cosf(PI*fabsf(xc-xp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(WallStencilWidth*DX)))/((WallStencilWidth*2)*(WallStencilWidth*2)*(WallStencilWidth*2));
														NODE_VEL_N[0][caps][node] += weightsum*(VEL_SCALE*Uvoid_3[ii][jj][kk]);
														NODE_VEL_N[1][caps][node] += weightsum*(VEL_SCALE*Vvoid_3[ii][jj][kk]);
														NODE_VEL_N[2][caps][node] += weightsum*(VEL_SCALE*Wvoid_3[ii][jj][kk]);
													}
												}
											}
											
											if(k >= 3 && k <= KM-3){
												ii = i;
												jj = j;
												kk = k;
												
												o = DOMID[ii][jj][kk];
												
												if(B_index[ii][jj][kk] > 0 && B_index[ii][jj][kk] < 4998){
													distance = sqrtf((xc-xp)*(xc-xp) + (yc-yp)*(yc-yp) + (zc-zp)*(zc-zp));
													if(distance <= StencilWidth*DX){
														weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
														NODE_VEL_N[0][caps][node] += weightsum*(VEL_SCALE*U[o]);
														NODE_VEL_N[1][caps][node] += weightsum*(VEL_SCALE*V[o]);
														NODE_VEL_N[2][caps][node] += weightsum*(VEL_SCALE*W[o]);
													}
												}
												else if(B_index[ii][jj][kk] == 4998){//plugged cell location - this is deprecated - it is no longer used 
													distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
													if(distance <= StencilWidth*DX){
														weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
														NODE_VEL_N[0][caps][node] += weightsum*(VEL_SCALE*Uvoid[ii][jj][kk]);
														NODE_VEL_N[1][caps][node] += weightsum*(VEL_SCALE*Vvoid[ii][jj][kk]);
														NODE_VEL_N[2][caps][node] += weightsum*(VEL_SCALE*Wvoid[ii][jj][kk]);
													}
												}
												else if(W_index[ii][jj][kk] > 0){
													distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
													if(distance <= WallStencilWidth*DX){
														weightsum = (1.+cosf(PI*fabsf(xc-xp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(WallStencilWidth*DX)))/((WallStencilWidth*2)*(WallStencilWidth*2)*(WallStencilWidth*2));
														NODE_VEL_N[0][caps][node] += weightsum*(VEL_SCALE*Uvoid[ii][jj][kk]);
														NODE_VEL_N[1][caps][node] += weightsum*(VEL_SCALE*Vvoid[ii][jj][kk]);
														NODE_VEL_N[2][caps][node] += weightsum*(VEL_SCALE*Wvoid[ii][jj][kk]);
													}
												}
											}
										}
									}
								}
							}	
						}
					}///Slave cells still following parent velocity
					else if(CellStatus[caps] > -4998){//cloned slave cells
						ParentID = -CellStatus[caps];
						for(node=1;node<=NodeM;node++){
							
							xp = NODECOORD[0][caps][node]; yp = NODECOORD[1][caps][node]; zp = NODECOORD[2][caps][node];
							
							check = 0;
							if(ParentID <= 188 && zp > ZM_1+5.*DX) check = 1;//region with independent movement
							if(ParentID >= 189 && ParentID <= 329 && zp < Z0_2-5.*DX) check = 1;//region with independent movement
							if(ParentID >= 330 && ParentID <= 344 && zp < Z0_3-5.*DX) check = 1;//region with independent movement
							
							i = (int)((xp-X0)/DX + 1.01);
							j = (int)((yp-Y0)/DX + 1.01);
							k = (int)((zp-Z0)/DX + 1.01);
							
							if(check == 0){
								NODE_VEL[0][caps][node] = NODE_VEL[0][ParentID][node];
								NODE_VEL[1][caps][node] = NODE_VEL[1][ParentID][node];
								NODE_VEL[2][caps][node] = NODE_VEL[2][ParentID][node];
							}
							else{
								weightsum = 0; bvalue = 0; mindis = 1.5; region = 5000; triwall = 0;
								i1 = i-4; i2 = i+4;
								j1 = j-4; j2 = j+4;
								k1 = k-4; k2 = k+4;
								
								
								if(i1 < 1) i1 = 1; if(i2> IM) i2 = IM;
								if(j1 < 1) j1 = 1; if(j2> JM) j2 = JM;
								if(k1 < 3) k1 = 3; if(k2> KM-3) k2 = KM-3;
								
								for(i=i1;i<=i2;i++){
									for(j=j1;j<=j2;j++){
										for(k=k1;k<=k2;k++){
											xc = X0 + (i-1)*DX;
											yc = Y0 + (j-1)*DX;
											zc = Z0 + (k-1)*DX;
																						
											if(k >= 3 && k <= KM-3){
												ii = i;
												jj = j;
												kk = k;
												
												o = DOMID[ii][jj][kk];
												
												if(B_index[ii][jj][kk] > 0){
													distance = sqrtf((xc-xp)*(xc-xp) + (yc-yp)*(yc-yp) + (zc-zp)*(zc-zp));
													if(distance <= StencilWidth*DX){
														weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
														NODE_VEL[0][caps][node] += weightsum*(VEL_SCALE*U[o]);
														NODE_VEL[1][caps][node] += weightsum*(VEL_SCALE*V[o]);
														NODE_VEL[2][caps][node] += weightsum*(VEL_SCALE*W[o]);
													}
												}
												else if(W_index[ii][jj][kk] > 0){
													distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
													if(distance <= WallStencilWidth*DX){
														weightsum = (1.+cosf(PI*fabsf(xc-xp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(WallStencilWidth*DX)))/((WallStencilWidth*2)*(WallStencilWidth*2)*(WallStencilWidth*2));
														NODE_VEL[0][caps][node] += weightsum*(VEL_SCALE*Uvoid[ii][jj][kk]);
														NODE_VEL[1][caps][node] += weightsum*(VEL_SCALE*Vvoid[ii][jj][kk]);
														NODE_VEL[2][caps][node] += weightsum*(VEL_SCALE*Wvoid[ii][jj][kk]);
													}
												}
												
												if(W_index[ii][jj][kk] > 0){
													nodewall = W_index[ii][jj][kk];
													
													neighmax = WALLNODETOCELLNEIGHBORLIST[0][nodewall];
															
													for(l=1;l<=neighmax;l++){
														tri = WALLNODETOCELLNEIGHBORLIST[l][nodewall];
														n1 = WALLCELLNODE[0][tri]; n2 = WALLCELLNODE[1][tri]; n3 = WALLCELLNODE[2][tri];
														xp1 = WALLNODECOORD[0][n1];
														yp1 = WALLNODECOORD[1][n1];
														zp1 = WALLNODECOORD[2][n1];
														distance = sqrtf((xp-xp1)*(xp-xp1) + (yp-yp1)*(yp-yp1) + (zp-zp1)*(zp-zp1));
														if(distance < mindis){
															region = 0;
															mindis = distance;
															triwall = tri;
														}
														
														xp2 = WALLNODECOORD[0][n2];
														yp2 = WALLNODECOORD[1][n2];
														zp2 = WALLNODECOORD[2][n2];
														distance = sqrtf((xp-xp2)*(xp-xp2) + (yp-yp2)*(yp-yp2) + (zp-zp2)*(zp-zp2));
														if(distance < mindis){
															region = 0;
															mindis = distance;
															triwall = tri;
														}
														
														xp3 = WALLNODECOORD[0][n3];
														yp3 = WALLNODECOORD[1][n3];
														zp3 = WALLNODECOORD[2][n3];
														distance = sqrtf((xp-xp3)*(xp-xp3) + (yp-yp3)*(yp-yp3) + (zp-zp3)*(zp-zp3));
														if(distance < mindis){
															region = 0;
															mindis = distance;
															triwall = tri;
														}
														
														xp4 = WALL_CENTROID[0][tri];
														yp4 = WALL_CENTROID[1][tri];
														zp4 = WALL_CENTROID[2][tri];
														distance = sqrtf((xp-xp4)*(xp-xp4) + (yp-yp4)*(yp-yp4) + (zp-zp4)*(zp-zp4));
														if(distance < mindis){
															region = 0;
															mindis = distance;
															triwall = tri;
														}
													}
												}
											}
										}
									}
								}
								
								///Calculate wall repulsion from nearby wall triangle elements
								if(triwall > 0){
									
									if(region == 0){
										NODE_WALLNEIGH[caps][node] = triwall;
										n1 = WALLCELLNODE[0][triwall]; n2 = WALLCELLNODE[1][triwall]; n3 = WALLCELLNODE[2][triwall];
										xp1 = WALLNODECOORD[0][n1];
										yp1 = WALLNODECOORD[1][n1];
										zp1 = WALLNODECOORD[2][n1];
										
										distance = (xp-xp1)*WALLCELLNORMALS[0][triwall] + (yp-yp1)*WALLCELLNORMALS[1][triwall] + (zp-zp1)*WALLCELLNORMALS[2][triwall];
										
										if(distance < 0.25){
											if(fabsf(distance) < 1.0){
												ReflectProj = 0;
												repulf = (0.25 - distance)/1.0;
												if(distance < -0.75) repulf = 1.;
											
												for(v=0;v<=2;v++){
													NODE_WALLREPUL[v][caps][node] += repulf*(1.e-04)*WALLCELLNORMALS[v][triwall];
												}
											}
										}
									}
								}
							}	
						}		
						for(node=1;node<=NodeM_N;node++){
								
							xp = NODECOORD_N[0][caps][node]; yp = NODECOORD_N[1][caps][node]; zp = NODECOORD_N[2][caps][node];
							i = (int)((xp-X0)/DX + 1.01);
							j = (int)((yp-Y0)/DX + 1.01);
							k = (int)((zp-Z0)/DX + 1.01);
								
							check = 0;
							if(ParentID <= 188 && zp > ZM_1+5.*DX) check = 1;//region with independent movement
							if(ParentID >= 189 && ParentID <= 329 && zp < Z0_2-5.*DX) check = 1;//region with independent movement
							if(ParentID >= 330 && ParentID <= 344 && zp < Z0_3-5.*DX) check = 1;//region with independent movement
								
							if(check == 0){									
								for(i=0;i<=2;i++){
									NODE_VEL_N[i][caps][node] = NODE_VEL_N[i][ParentID][node];
								}
							}
							else{
								weightsum = 0; bvalue = 0; mindis = 1.5; region = 5000; triwall = 0;
								i1 = i-4; i2 = i+4;
								j1 = j-4; j2 = j+4;
								k1 = k-4; k2 = k+4;
								
								
								if(i1 < 1) i1 = 1; if(i2> IM) i2 = IM;
								if(j1 < 1) j1 = 1; if(j2> JM) j2 = JM;
								if(k1 < 3) k1 = 3; if(k2> KM-3) k2 = KM-3;
								
								for(i=i1;i<=i2;i++){
									for(j=j1;j<=j2;j++){
										for(k=k1;k<=k2;k++){
											xc = X0 + (i-1)*DX;
											yc = Y0 + (j-1)*DX;
											zc = Z0 + (k-1)*DX;
																						
											if(k >= 3 && k <= KM-3){
												ii = i;
												jj = j;
												kk = k;
												
												o = DOMID[ii][jj][kk];
												
												if(B_index[ii][jj][kk] > 0){
													distance = sqrtf((xc-xp)*(xc-xp) + (yc-yp)*(yc-yp) + (zc-zp)*(zc-zp));
													if(distance <= StencilWidth*DX){
														weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
														NODE_VEL_N[0][caps][node] += weightsum*(VEL_SCALE*U[o]);
														NODE_VEL_N[1][caps][node] += weightsum*(VEL_SCALE*V[o]);
														NODE_VEL_N[2][caps][node] += weightsum*(VEL_SCALE*W[o]);
													}
												}
												else if(W_index[ii][jj][kk] > 0){
													distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
													if(distance <= WallStencilWidth*DX){
														weightsum = (1.+cosf(PI*fabsf(xc-xp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(WallStencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(WallStencilWidth*DX)))/((WallStencilWidth*2)*(WallStencilWidth*2)*(WallStencilWidth*2));
														NODE_VEL_N[0][caps][node] += weightsum*(VEL_SCALE*Uvoid[ii][jj][kk]);
														NODE_VEL_N[1][caps][node] += weightsum*(VEL_SCALE*Vvoid[ii][jj][kk]);
														NODE_VEL_N[2][caps][node] += weightsum*(VEL_SCALE*Wvoid[ii][jj][kk]);
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
			
		}
		
		///RBC velocity interpolation and Wall repulsion velocity increments END
		
		
				
		///Membrane and nucelus interaction velocity increments START////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		///Membrane and nucelus interaction velocity increments START////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		if(RBCSOLVERFLAG == 1){
			#pragma omp parallel private(nzmax,nzmin,ParentID,check,e2v,dotpdt,Vneigh,ReflectProj,caps,mindis,node,xc,yc,zc,m,xp,yp,zp,xp2,yp2,zp2,xp3,yp3,zp3,distance,distance2,distance3,neighmax,l,n1,n2,n3,n1coord,n2coord,n3coord,repulf,e1v,tri,ep,tp,v,trimax,i)
			{
				#pragma omp for schedule(static,50)  nowait
				for(caps=1;caps<=CapsM;caps++){/// caps for-loop
					if(CellStatus[caps] > 0){///if conditional for updating clone parents + independent clone cells	
						nzmax = CAPS_MAX[2][caps];
						nzmin = CAPS_MIN[2][caps];
						if(iter%1000 == 0) NucMembUpdate[caps] = 0;//this flag is also zero when cells appear in the simulation for the first time
						///verify closest membrane neighbor node to nucleus node
						if(NucMembUpdate[caps] == 0){
							for(node=1;node<=NodeM_N;node++){
								xc = NODECOORD_N[0][caps][node];
								yc = NODECOORD_N[1][caps][node];
								zc = NODECOORD_N[2][caps][node];
								
								for(l=0;l<=NucMembMAX;l++) NUCLNEIGH[l][caps][node] = 0;
								
								if(CellStatus[caps] > 3){
									if((nzmax - nzmin) > 30.){//this cell is crossing periodic BC
										if(caps <= 188){
											if(zc < 0.5*(Z0_1+ZM_1)) zc += (ZM_1-Z0_1);
										}
										else if(caps >= 189 && caps <= 329){
											if(zc > 0.5*(Z0_2+ZM_2)) zc += -(ZM_2-Z0_2);
										}
										else if(caps >= 330 && caps <= 344){
											if(zc > 0.5*(Z0_3+ZM_3)) zc += -(ZM_3-Z0_3);
										}
									}
								}
								
								l = 0;
								for(m=1; m<=TriM; m++){
									xp = TRI_CENTROID[0][caps][m];
									yp = TRI_CENTROID[1][caps][m];
									zp = TRI_CENTROID[2][caps][m];
									
										if(caps <= 188){
											zp2 = zp - (ZM_1-Z0_1);
											zp3 = zp + (ZM_1-Z0_1);
										}
										else if(caps >= 189 && caps <= 329){
											zp2 = zp + (ZM_2-Z0_2);
											zp3 = zp - (ZM_2-Z0_2);
										}
										else if(caps >= 330 && caps <= 344){
											zp2 = zp + (ZM_3-Z0_3);
											zp3 = zp - (ZM_3-Z0_3);
										}
									
									
									distance = sqrtf((xc-xp)*(xc-xp) + (yc-yp)*(yc-yp) + (zc-zp)*(zc-zp));
									if(caps <= 344){
										distance2 = sqrtf((xc-xp)*(xc-xp) + (yc-yp)*(yc-yp) + (zc-zp2)*(zc-zp2));
										distance3 = sqrtf((xc-xp)*(xc-xp) + (yc-yp)*(yc-yp) + (zc-zp3)*(zc-zp3));
									}
									else{
										distance2 = 5000.;
										distance3 = 5000.;
									}
									
									if(distance < 0.5 || distance2 < 0.5 || distance3 < 0.5){
										l += 1;
										
										if(l >= NucMembMAX) printf("Nucleus neighbor list is maxed out!\n");
										NUCLNEIGH[l][caps][node] = m;
										NUCLNEIGH[0][caps][node] = l;
									}
								}
							}
							NucMembUpdate[caps] = 1;
						}
						
						
						/// neighbor node verification complete
						for(node=1;node<=NodeM_N;node++){//Node_N for loop
							mindis = 5000.;
							
							xc = NODECOORD_N[0][caps][node];
							yc = NODECOORD_N[1][caps][node];
							zc = NODECOORD_N[2][caps][node];
							
							if(CellStatus[caps] > 3){
								if((nzmax - nzmin) > 30.){//this cell is crossing periodic BC
									if(caps <= 188){
										if(zc < 0.5*(Z0_1+ZM_1)) zc += (ZM_1-Z0_1);
									}
									else if(caps >= 189 && caps <= 329){
										if(zc > 0.5*(Z0_2+ZM_2)) zc += -(ZM_2-Z0_2);
									}
									else if(caps >= 330 && caps <= 344){
										if(zc > 0.5*(Z0_3+ZM_3)) zc += -(ZM_3-Z0_3);
									}
								}
							}
							
							trimax = NUCLNEIGH[0][caps][node];
							
							if(trimax == 0) continue;
							
								for(l=1;l<=NucMembMAX;l++){///search query for-loop
									tri = NUCLNEIGH[l][caps][node];
									
									xp = TRI_CENTROID[0][caps][tri];
									yp = TRI_CENTROID[1][caps][tri];
									zp = TRI_CENTROID[2][caps][tri];
									
									
									distance = sqrtf((xc-xp)*(xc-xp) + (yc-yp)*(yc-yp) + (zc-zp)*(zc-zp));
									
									
									if(distance < 1.){
										check = 0;
										n1 = TRINODE1[tri]; n2 = TRINODE2[tri]; n3 = TRINODE3[tri];
										e1v[0] = TRI_NORMAL[0][caps][tri]; e1v[1] = TRI_NORMAL[1][caps][tri]; e1v[2] = TRI_NORMAL[2][caps][tri];
										n1coord[0] = NODECOORD[0][caps][n1];
										n1coord[1] = NODECOORD[1][caps][n1];
										n1coord[2] = NODECOORD[2][caps][n1];
										
										n2coord[0] = NODECOORD[0][caps][n2];
										n2coord[1] = NODECOORD[1][caps][n2];
										n2coord[2] = NODECOORD[2][caps][n2];
										
										n3coord[0] = NODECOORD[0][caps][n3];
										n3coord[1] = NODECOORD[1][caps][n3];
										n3coord[2] = NODECOORD[2][caps][n3];	
										
										if(CellStatus[caps] > 3){
											if((nzmax - nzmin) > 30.){//this cell is crossing periodic BC
												if(caps <= 188){
													if(n1coord[2] < 0.5*(Z0_1+ZM_1)) n1coord[2] += (ZM_1-Z0_1);
													if(n2coord[2] < 0.5*(Z0_1+ZM_1)) n2coord[2] += (ZM_1-Z0_1);
													if(n3coord[2] < 0.5*(Z0_1+ZM_1)) n3coord[2] += (ZM_1-Z0_1);
													
												}
												else if(caps >= 189 && caps <= 329){
													if(n1coord[2] > 0.5*(Z0_2+ZM_2)) n1coord[2] += -(ZM_2-Z0_2);
													if(n2coord[2] > 0.5*(Z0_2+ZM_2)) n2coord[2] += -(ZM_2-Z0_2);
													if(n3coord[2] > 0.5*(Z0_2+ZM_2)) n3coord[2] += -(ZM_2-Z0_2);
												}
												else if(caps >= 330 && caps <= 344){
													if(n1coord[2] > 0.5*(Z0_3+ZM_3)) n1coord[2] += -(ZM_3-Z0_3);
													if(n2coord[2] > 0.5*(Z0_3+ZM_3)) n2coord[2] += -(ZM_3-Z0_3);
													if(n3coord[2] > 0.5*(Z0_3+ZM_3)) n3coord[2] += -(ZM_3-Z0_3);
												}
											}
										}
											
										for(v=0;v<=2;v++){
											ep[1][v] = (n1coord[v]+n2coord[v])/2.;
											ep[2][v] = (n2coord[v]+n3coord[v])/2.;
											ep[3][v] = (n3coord[v]+n1coord[v])/2.;
											ep[4][v] = (n1coord[v]+ep[1][v])/2.;
											ep[5][v] = (ep[1][v]+n2coord[v])/2.;
											ep[6][v] = (n2coord[v]+ep[2][v])/2.;
											ep[7][v] = (ep[2][v]+n3coord[v])/2.;
											ep[8][v] = (n3coord[v]+ep[3][v])/2.;
											ep[9][v] = (ep[3][v]+n1coord[v])/2.;
											
											tp[1][v] = (n1coord[v]+n2coord[v]+n3coord[v])/3.;
											tp[2][v] = (ep[3][v]+n1coord[v]+ep[1][v])/3.;
											tp[3][v] = (ep[1][v]+n2coord[v]+ep[2][v])/3.;
											tp[4][v] = (ep[2][v]+n3coord[v]+ep[3][v])/3.;
											
											tp[5][v] = (ep[1][v]+ep[3][v])/2.;
											tp[6][v] = (ep[1][v]+ep[2][v])/2.;
											tp[7][v] = (ep[2][v]+ep[3][v])/2.;
											
											tp[8][v]  = (n1coord[v]+tp[2][v])/2.;
											tp[9][v]  = (ep[1][v]+tp[2][v])/2.;
											tp[10][v] = (ep[3][v]+tp[2][v])/2.;
											
											tp[11][v] = (ep[1][v]+tp[3][v])/2.;
											tp[12][v] = (n2coord[v]+tp[3][v])/2.;
											tp[13][v] = (ep[2][v]+tp[3][v])/2.;
											
											tp[14][v] = (ep[3][v]+tp[4][v])/2.;
											tp[15][v] = (ep[2][v]+tp[4][v])/2.;
											tp[16][v] = (n3coord[v]+tp[4][v])/2.;
											
											tp[17][v] = (ep[1][v]+tp[1][v])/2.;
											tp[18][v] = (ep[2][v]+tp[1][v])/2.;
											tp[19][v] = (ep[3][v]+tp[1][v])/2.;
										}
											
										for(v=1;v<=19;v++){
											distance = sqrtf((tp[v][0]-xc)*(tp[v][0]-xc) + (tp[v][1]-yc)*(tp[v][1]-yc) + (tp[v][2]-zc)*(tp[v][2]-zc));
											
											if(distance <= 0.2){
												repulf = (0.2 - distance)/0.2;
											
												repulf = repulf*(4.1666667e-07);
												
													if(check == 0){
														for(i=0;i<=2;i++){
															
															
																
															NODE_WALLREPUL[i][caps][n1] += repulf*e1v[i]/3.;
															NODE_WALLREPUL[i][caps][n2] += repulf*e1v[i]/3.;
															NODE_WALLREPUL[i][caps][n3] += repulf*e1v[i]/3.;
														}
													}
												
											}
										}
									}
								}///end of search query for-loop
								
								
						}///end of Node_N for-loop
						
					}///end of if conditional for updating clone parents + independent clone cells
					else if(CellStatus[caps] > -4998){
						ParentID = -CellStatus[caps];
						nzmax = CAPS_MAX[2][caps];
						nzmin = CAPS_MIN[2][caps];
						
						if(iter%1000 == 0) NucMembUpdate[caps] = 0;//this flag is also zero when cells appear in the simulation for the first time
						///verify closest membrane neighbor node to nucleus node
						if(NucMembUpdate[caps] == 0){
							for(node=1;node<=NodeM_N;node++){
								xc = NODECOORD_N[0][caps][node];
								yc = NODECOORD_N[1][caps][node];
								zc = NODECOORD_N[2][caps][node];
								
								for(l=0;l<=NucMembMAX;l++) NUCLNEIGH[l][caps][node] = 0;
								
								if(CellStatus[caps] > 3){
									if(caps <= 188){
										if(zc < 0.5*(Z0_1+ZM_1)) zc += (ZM_1-Z0_1);
									}
									else if(caps >= 189 && caps <= 329){
										if(zc > 0.5*(Z0_2+ZM_2)) zc += -(ZM_2-Z0_2);
									}
									else if(caps >= 330 && caps <= 344){
										if(zc > 0.5*(Z0_3+ZM_3)) zc += -(ZM_3-Z0_3);
									}
								}
								
								l = 0;
								for(m=1; m<=TriM; m++){
									xp = TRI_CENTROID[0][caps][m];
									yp = TRI_CENTROID[1][caps][m];
									zp = TRI_CENTROID[2][caps][m];
									
									
									distance = sqrtf((xc-xp)*(xc-xp) + (yc-yp)*(yc-yp) + (zc-zp)*(zc-zp));
									
									if(distance < 0.5){
										l += 1;
										
										if(l >= NucMembMAX) printf("Nucleus neighbor list is maxed out!\n");
										NUCLNEIGH[l][caps][node] = m;
										NUCLNEIGH[0][caps][node] = l;
									}
								}
							}
							NucMembUpdate[caps] = 1;
						}
						
						
						/// neighbor node verification complete
						for(node=1;node<=NodeM_N;node++){//Node_N for loop
							mindis = 5000.;
							
							xc = NODECOORD_N[0][caps][node];
							yc = NODECOORD_N[1][caps][node];
							zc = NODECOORD_N[2][caps][node];
							
							
							check = 0;
							if(ParentID <= 188 && zc > ZM_1+5.*DX) check = 1;//region with independent movement
							if(ParentID >= 189 && ParentID <= 329 && zc < Z0_2-5.*DX) check = 1;//region with independent movement
							if(ParentID >= 330 && ParentID <= 344 && zc < Z0_3-5.*DX) check = 1;//region with independent movement
								
								
							if(check == 0){
								trimax = NUCLNEIGH[0][caps][node];
								
								if(trimax == 0) continue;
								
								
									for(l=1;l<=NucMembMAX;l++){///search query for-loop
										tri = NUCLNEIGH[l][caps][node];
										
										xp = TRI_CENTROID[0][caps][tri];
										yp = TRI_CENTROID[1][caps][tri];
										zp = TRI_CENTROID[2][caps][tri];
										
										
										distance = sqrtf((xc-xp)*(xc-xp) + (yc-yp)*(yc-yp) + (zc-zp)*(zc-zp));
										
										if(distance < 1.){
											check = 0;
											if(AREA_UPDATE[caps][tri] == -1) check = 1;
											n1 = TRINODE1[tri]; n2 = TRINODE2[tri]; n3 = TRINODE3[tri];
											e1v[0] = TRI_NORMAL[0][caps][tri]; e1v[1] = TRI_NORMAL[1][caps][tri]; e1v[2] = TRI_NORMAL[2][caps][tri];
											n1coord[0] = NODECOORD[0][caps][n1];
											n1coord[1] = NODECOORD[1][caps][n1];
											n1coord[2] = NODECOORD[2][caps][n1];
											
											n2coord[0] = NODECOORD[0][caps][n2];
											n2coord[1] = NODECOORD[1][caps][n2];
											n2coord[2] = NODECOORD[2][caps][n2];
											
											n3coord[0] = NODECOORD[0][caps][n3];
											n3coord[1] = NODECOORD[1][caps][n3];
											n3coord[2] = NODECOORD[2][caps][n3];	
											
											
												
											for(v=0;v<=2;v++){
												ep[1][v] = (n1coord[v]+n2coord[v])/2.;
												ep[2][v] = (n2coord[v]+n3coord[v])/2.;
												ep[3][v] = (n3coord[v]+n1coord[v])/2.;
												ep[4][v] = (n1coord[v]+ep[1][v])/2.;
												ep[5][v] = (ep[1][v]+n2coord[v])/2.;
												ep[6][v] = (n2coord[v]+ep[2][v])/2.;
												ep[7][v] = (ep[2][v]+n3coord[v])/2.;
												ep[8][v] = (n3coord[v]+ep[3][v])/2.;
												ep[9][v] = (ep[3][v]+n1coord[v])/2.;
												
												tp[1][v] = (n1coord[v]+n2coord[v]+n3coord[v])/3.;
												tp[2][v] = (ep[3][v]+n1coord[v]+ep[1][v])/3.;
												tp[3][v] = (ep[1][v]+n2coord[v]+ep[2][v])/3.;
												tp[4][v] = (ep[2][v]+n3coord[v]+ep[3][v])/3.;
												
												tp[5][v] = (ep[1][v]+ep[3][v])/2.;
												tp[6][v] = (ep[1][v]+ep[2][v])/2.;
												tp[7][v] = (ep[2][v]+ep[3][v])/2.;
												
												tp[8][v]  = (n1coord[v]+tp[2][v])/2.;
												tp[9][v]  = (ep[1][v]+tp[2][v])/2.;
												tp[10][v] = (ep[3][v]+tp[2][v])/2.;
												
												tp[11][v] = (ep[1][v]+tp[3][v])/2.;
												tp[12][v] = (n2coord[v]+tp[3][v])/2.;
												tp[13][v] = (ep[2][v]+tp[3][v])/2.;
												
												tp[14][v] = (ep[3][v]+tp[4][v])/2.;
												tp[15][v] = (ep[2][v]+tp[4][v])/2.;
												tp[16][v] = (n3coord[v]+tp[4][v])/2.;
												
												tp[17][v] = (ep[1][v]+tp[1][v])/2.;
												tp[18][v] = (ep[2][v]+tp[1][v])/2.;
												tp[19][v] = (ep[3][v]+tp[1][v])/2.;
											}
												
											for(v=1;v<=19;v++){
												distance = sqrtf((tp[v][0]-xc)*(tp[v][0]-xc) + (tp[v][1]-yc)*(tp[v][1]-yc) + (tp[v][2]-zc)*(tp[v][2]-zc));
												
												if(distance < 0.2){
													repulf = (0.2 - distance)/0.2;
											
													repulf = repulf*(4.1666667e-07);
													
													
													
													if(check == 0){
														for(i=0;i<=2;i++){
															
															
															NODE_WALLREPUL_N[i][caps][node] += -repulf*e1v[i];
																
															NODE_WALLREPUL[i][caps][n1] += repulf*e1v[i]/3.;
															NODE_WALLREPUL[i][caps][n2] += repulf*e1v[i]/3.;
															NODE_WALLREPUL[i][caps][n3] += repulf*e1v[i]/3.;
														}
													}
													
													mindis = distance;
													m = tri;
												}
											}
										}
									}///end of search query for-loop
									
									
							}
						}///end of Node_N for-loop
					}
				}///end of caps for-loop
			}
		}
		///Membrane and nucelus interaction module END////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		
		///Membrane and nucelus interaction module END////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		
		/// RBC to RBC interaction module START //////////////////////////////////////////////////
		
		if(RBCSOLVERFLAG == 1){
			/// ROI
			#pragma omp parallel private(ParentID,check,checkn,Vneigh,ReflectProj,caps,node,mindisIn,mindisOut,xc,yc,zc,e1v,e2v,e3v,i,j,k,i1,j1,k1,i2,j2,k2,o,neighmax,l,check2,tri,m,n1,n2,n3,xp,yp,zp,distance,v,c,n1coord,n2coord,n3coord,ep,tp,tc,dotpdt,dotpdt2,Outneigh,Inneigh,repulf)
			{
				#pragma omp for schedule(static,50) nowait
				
				for(caps=345;caps<=CapsM;caps++){
					if(CellStatus[caps] > -4998){///if-condition filtering only cell IDs in simulation domain
						if(CellStatus[caps] > 0){///if-condition filtering only independent clone cells
							for(node=1;node<=NodeM;node++){///NodeM for-loop
								xc = NODECOORD[0][caps][node];
								yc = NODECOORD[1][caps][node];
								zc = NODECOORD[2][caps][node];
								
								e1v[0] = NODE_NORMAL[0][caps][node];
								e1v[1] = NODE_NORMAL[1][caps][node];
								e1v[2] = NODE_NORMAL[2][caps][node];
									
								mindisIn = 5000.;
								mindisOut = 5000.;
								
								i = (int)((xc-X0)/DX + 1.01);    j = (int)((yc-Y0)/DX + 1.01);    k = (int)((zc-Z0)/DX + 1.01);
								i1 = i-4; i2 = i+4;		j1 = j-4; j2 = j+4;		k1 = k-4; k2 = k+4;
								
								if(i1 < 1) i1 = 1; if(i2 > IM) i2 = IM;
								if(j1 < 1) j1 = 1; if(j2 > JM) j2 = JM;
								if(k1 < 3) k1 = 3; if(k2 > KM-3) k2 = KM-3;
								
								for(i=i1;i<=i2;i++){///start of ijk search stencil for eulerian list of neighbors
									for(j=j1;j<=j2;j++){
										for(k=k1;k<=k2;k++){
											o = DOMID[i][j][k];
											
											neighmax = EulNeigh[0][o];
												
											for(l=1;l<=neighmax;l++){
												check2 = 0;
												tri = EulNeigh[l][o];
												
												
												m = (int)((tri-1)/TriM);
														
												if(m < 345) continue;
														
												tri += -m*TriM;
														
												if(m == caps){
													n1 = TRINODE1[tri]; n2 = TRINODE2[tri]; n3 = TRINODE3[tri];
													if(n1 == node || n2 == node || n3 == node) check2 = 1;// ignore the directly adjacent tri neighbors around current node
												}
													
												if(check2 == 0){
													xp = TRI_CENTROID[0][m][tri];
													yp = TRI_CENTROID[1][m][tri];
													zp = TRI_CENTROID[2][m][tri];
													
													e2v[0] = TRI_NORMAL[0][m][tri];
													e2v[1] = TRI_NORMAL[1][m][tri];
													e2v[2] = TRI_NORMAL[2][m][tri];
														
													dotpdt = e1v[0]*e2v[0]+e1v[1]*e2v[1]+e1v[2]*e2v[2];
													
													if(dotpdt < 0){
														distance = 	sqrtf((xc-xp)*(xc-xp) + (yc-yp)*(yc-yp) + (zc-zp)*(zc-zp));
															
														if(distance < 2.){//this won't work if a particular region is sheared 3 times its orginal length
															n1 = TRINODE1[tri]; n2 = TRINODE2[tri]; n3 = TRINODE3[tri];
																
															for(v=0;v<=2;v++){
																n1coord[v] = NODECOORD[v][m][n1];
																n2coord[v] = NODECOORD[v][m][n2];
																n3coord[v] = NODECOORD[v][m][n3];
															}
															
															for(v=0;v<=2;v++){
																tp[1][v] = n1coord[v];
																tp[2][v] = n2coord[v];
																tp[3][v] = n3coord[v];
																
																tp[4][v] = 0.5*(tp[1][v]+tp[2][v]);
																tp[5][v] = 0.5*(tp[2][v]+tp[3][v]);
																tp[6][v] = 0.5*(tp[1][v]+tp[3][v]);
																
																tp[7][v] = 0.5*(tp[1][v]+tp[4][v]);
																tp[8][v] = 0.5*(tp[2][v]+tp[4][v]);
																tp[9][v] = 0.5*(tp[2][v]+tp[5][v]);
																tp[10][v] = 0.5*(tp[5][v]+tp[3][v]);
																tp[11][v] = 0.5*(tp[3][v]+tp[6][v]);
																tp[12][v] = 0.5*(tp[6][v]+tp[1][v]);
																
																tp[13][v] = 0.5*(tp[4][v]+tp[6][v]);
																tp[14][v] = 0.5*(tp[4][v]+tp[5][v]);
																tp[15][v] = 0.5*(tp[5][v]+tp[6][v]);
																
																
																tc[1][v] = (tp[1][v]+tp[7][v]+tp[12][v])/3.;
																tc[2][v] = (tp[7][v]+tp[13][v]+tp[12][v])/3.;
																tc[3][v] = (tp[7][v]+tp[4][v]+tp[13][v])/3.;
																tc[4][v] = (tp[4][v]+tp[14][v]+tp[13][v])/3.;
																tc[5][v] = (tp[4][v]+tp[8][v]+tp[14][v])/3.;
																tc[6][v] = (tp[8][v]+tp[9][v]+tp[14][v])/3.;
																tc[7][v] = (tp[8][v]+tp[2][v]+tp[9][v])/3.;
																
																tc[8][v] = (tp[12][v]+tp[13][v]+tp[6][v])/3.;
																tc[9][v] = (tp[13][v]+tp[15][v]+tp[6][v])/3.;
																tc[10][v] = (tp[13][v]+tp[14][v]+tp[15][v])/3.;
																tc[11][v] = (tp[14][v]+tp[5][v]+tp[15][v])/3.;
																tc[12][v] = (tp[14][v]+tp[9][v]+tp[5][v])/3.;
																
																tc[13][v] = (tp[6][v]+tp[15][v]+tp[11][v])/3.;
																tc[14][v] = (tp[15][v]+tp[10][v]+tp[11][v])/3.;
																tc[15][v] = (tp[15][v]+tp[5][v]+tp[10][v])/3.;
																tc[16][v] = (tp[11][v]+tp[10][v]+tp[3][v])/3.;
															}
															
															for(v=1;v<=16;v++){
																distance = sqrtf((tc[v][0]-xc)*(tc[v][0]-xc) + (tc[v][1]-yc)*(tc[v][1]-yc) + (tc[v][2]-zc)*(tc[v][2]-zc));
																
																if(distance < 1.0){
																	e3v[0] = (tc[v][0]-xc)/distance;
																	e3v[1] = (tc[v][1]-yc)/distance;
																	e3v[2] = (tc[v][2]-zc)/distance;
																	
																	repulf = -8.333e-07*(1.0 - distance)/1.0;
																	
																			
																	for(c=0;c<=2;c++){
																		NODE_WALLREPUL[c][caps][node] += repulf*e3v[c];
																		
																	}
																	
																}
															}
														}
													}
												}
											}
										}
									}
								}///end of ijk search stencil for eulerian list of neighbors
							}///end of NodeM for-loop for independent clone nodes
						}///end of if-condition filtering only independent clone cells
						else{///if-condition filtering only slave clone cells
							for(node=1;node<=NodeM;node++){///NodeM for-loop
								xc = NODECOORD[0][caps][node];
								yc = NODECOORD[1][caps][node];
								zc = NODECOORD[2][caps][node];
								
								e1v[0] = NODE_NORMAL[0][caps][node];
								e1v[1] = NODE_NORMAL[1][caps][node];
								e1v[2] = NODE_NORMAL[2][caps][node];
								
								check = 0;
								if(ParentID <= 188 && zp > ZM_1+5.*DX) check = 1;//region with independent movement
								if(ParentID >= 189 && ParentID <= 329 && zp < Z0_2-5.*DX) check = 1;//region with independent movement
								if(ParentID >= 330 && ParentID <= 344 && zp < Z0_3-5.*DX) check = 1;//region with independent movement
								
								
								if(check == 1){///start of if-condition filtering for updating non-penetrating nodes
									mindisIn = 5000.;
									mindisOut = 5000.;
									
									i = (int)((xc-X0)/DX + 1.01);    j = (int)((yc-Y0)/DX + 1.01);    k = (int)((zc-Z0)/DX + 1.01);
									i1 = i-4; i2 = i+4;		j1 = j-4; j2 = j+4;		k1 = k-4; k2 = k+4;
									
									if(i1 < 1) i1 = 1; if(i2 > IM) i2 = IM;
									if(j1 < 1) j1 = 1; if(j2 > JM) j2 = JM;
									if(k1 < 3) k1 = 3; if(k2 > KM-3) k2 = KM-3;
										
									for(i=i1;i<=i2;i++){///start of ijk search stencil for eulerian list of neighbors
										for(j=j1;j<=j2;j++){
											for(k=k1;k<=k2;k++){
												o = DOMID[i][j][k];
												
												neighmax = EulNeigh[0][o];
													
												for(l=1;l<=neighmax;l++){
													check2 = 0;
													tri = EulNeigh[l][o];
											
													m = (int)((tri-1)/TriM);
													
													if(m < 345) continue;
													
													tri += -m*TriM;
													
													if(m == caps){
														n1 = TRINODE1[tri]; n2 = TRINODE2[tri]; n3 = TRINODE3[tri];
														if(n1 == node || n2 == node || n3 == node) check2 = 1;// ignore the directly adjacent tri neighbors around current node
													}
														
													if(check2 == 0){
														xp = TRI_CENTROID[0][m][tri];
														yp = TRI_CENTROID[1][m][tri];
														zp = TRI_CENTROID[2][m][tri];
														
														e2v[0] = TRI_NORMAL[0][m][tri];
														e2v[1] = TRI_NORMAL[1][m][tri];
														e2v[2] = TRI_NORMAL[2][m][tri];
														
														dotpdt = e1v[0]*e2v[0]+e1v[1]*e2v[1]+e1v[2]*e2v[2];
														
														if(dotpdt < 0){
															distance = 	sqrtf((xc-xp)*(xc-xp) + (yc-yp)*(yc-yp) + (zc-zp)*(zc-zp));
																
															if(distance < 2.){//this won't work if a particular region is sheared 3 times its orginal length
																n1 = TRINODE1[tri]; n2 = TRINODE2[tri]; n3 = TRINODE3[tri];
																	
																for(v=0;v<=2;v++){
																	n1coord[v] = NODECOORD[v][m][n1];
																	n2coord[v] = NODECOORD[v][m][n2];
																	n3coord[v] = NODECOORD[v][m][n3];
																}
																
																for(v=0;v<=2;v++){
																	tp[1][v] = n1coord[v];
																	tp[2][v] = n2coord[v];
																	tp[3][v] = n3coord[v];
																	
																	tp[4][v] = 0.5*(tp[1][v]+tp[2][v]);
																	tp[5][v] = 0.5*(tp[2][v]+tp[3][v]);
																	tp[6][v] = 0.5*(tp[1][v]+tp[3][v]);
																	
																	tp[7][v] = 0.5*(tp[1][v]+tp[4][v]);
																	tp[8][v] = 0.5*(tp[2][v]+tp[4][v]);
																	tp[9][v] = 0.5*(tp[2][v]+tp[5][v]);
																	tp[10][v] = 0.5*(tp[5][v]+tp[3][v]);
																	tp[11][v] = 0.5*(tp[3][v]+tp[6][v]);
																	tp[12][v] = 0.5*(tp[6][v]+tp[1][v]);
																	
																	tp[13][v] = 0.5*(tp[4][v]+tp[6][v]);
																	tp[14][v] = 0.5*(tp[4][v]+tp[5][v]);
																	tp[15][v] = 0.5*(tp[5][v]+tp[6][v]);
																	
																	
																	tc[1][v] = (tp[1][v]+tp[7][v]+tp[12][v])/3.;
																	tc[2][v] = (tp[7][v]+tp[13][v]+tp[12][v])/3.;
																	tc[3][v] = (tp[7][v]+tp[4][v]+tp[13][v])/3.;
																	tc[4][v] = (tp[4][v]+tp[14][v]+tp[13][v])/3.;
																	tc[5][v] = (tp[4][v]+tp[8][v]+tp[14][v])/3.;
																	tc[6][v] = (tp[8][v]+tp[9][v]+tp[14][v])/3.;
																	tc[7][v] = (tp[8][v]+tp[2][v]+tp[9][v])/3.;
																	
																	tc[8][v] = (tp[12][v]+tp[13][v]+tp[6][v])/3.;
																	tc[9][v] = (tp[13][v]+tp[15][v]+tp[6][v])/3.;
																	tc[10][v] = (tp[13][v]+tp[14][v]+tp[15][v])/3.;
																	tc[11][v] = (tp[14][v]+tp[5][v]+tp[15][v])/3.;
																	tc[12][v] = (tp[14][v]+tp[9][v]+tp[5][v])/3.;
																	
																	tc[13][v] = (tp[6][v]+tp[15][v]+tp[11][v])/3.;
																	tc[14][v] = (tp[15][v]+tp[10][v]+tp[11][v])/3.;
																	tc[15][v] = (tp[15][v]+tp[5][v]+tp[10][v])/3.;
																	tc[16][v] = (tp[11][v]+tp[10][v]+tp[3][v])/3.;
																}
																
																for(v=1;v<=16;v++){
																	distance = sqrtf((tc[v][0]-xc)*(tc[v][0]-xc) + (tc[v][1]-yc)*(tc[v][1]-yc) + (tc[v][2]-zc)*(tc[v][2]-zc));
																	
																	if(distance < 1.0){
																		e3v[0] = (tc[v][0]-xc)/distance;
																		e3v[1] = (tc[v][1]-yc)/distance;
																		e3v[2] = (tc[v][2]-zc)/distance;
																		
																		repulf = -8.333e-07*(1.0 - distance)/1.0;
																		
																				
																		for(c=0;c<=2;c++){
																			NODE_WALLREPUL[c][caps][node] += repulf*e3v[c];
																			
																		}
																	}
																}
															}
														}
													}
												}
											}
										}
									}///end of ijk search stencil for eulerian list of neighbors
								}///end of if-condition for free regions of slave cells	
							}///end of NodeM for-loop for slave clone nodes
						}///end of if-condition filtering only slave clone cells
					}///end of if-condition for filtering only cell IDs in simulation domain
				}
			}
			
			/// DA parent
			#pragma omp parallel private(checkn,Vneigh,ReflectProj,caps,node,mindisIn,mindisOut,xc,yc,zc,e1v,e2v,e3v,ic,jc,kc,i1,j1,k1,i2,j2,k2,i,j,k,check,kk,o,neighmax,l,check2,tri,m,n1,n2,n3,xp,yp,zp,distance,v,c,n1coord,n2coord,n3coord,ep,tp,tc,dotpdt,dotpdt2,Outneigh,Inneigh,repulf)
			{
				#pragma omp for schedule(static,3) nowait
				
				for(caps=1;caps<=188;caps++){///CapsM for-loop for DA Cells
					if(CellStatus[caps] > -5000){///if-condition for filtering only cell IDs in simulation domain
						
						for(node=1;node<=NodeM;node++){///NodeM for-loop for DA Cell nodes
							///start of penetration condition check and updates
					
							xc = NODECOORD[0][caps][node];
							yc = NODECOORD[1][caps][node];
							zc = NODECOORD[2][caps][node];
							mindisIn = 5000.;
							mindisOut = 5000.;
							
							e1v[0] = NODE_NORMAL[0][caps][node];
							e1v[1] = NODE_NORMAL[1][caps][node];
							e1v[2] = NODE_NORMAL[2][caps][node];
							
							ic = (int)((xc-X0_1)/DX + 1.01);    jc = (int)((yc-Y0_1)/DX + 1.01);    kc = (int)((zc-Z0_1)/DX + 1.01);
							i1 = ic-4; i2 = ic+4;		j1 = jc-4; j2 = jc+4;		k1 = kc-4; k2 = kc+4;
							
							if(i1 < 1) i1 = 1; if(i2 > IM_1) i2 = IM_1;
							if(j1 < 1) j1 = 1; if(j2 > JM_1) j2 = JM_1;
					
							
						
							for(i=i1;i<=i2;i++){///start of ijk search stencil for eulerian list of neighbors
								for(j=j1;j<=j2;j++){
									for(k=k1;k<=k2;k++){
										check = 0;
										kk = k;
										if(k >= KM_1){// this means control point is near right bc hence querried neighbor exceeds right of right bc
											kk = k-KM_1+1;
											check = 1;
										}
										else if(k < 1){
											kk = k+KM_1-1;
											check = 2;//this means control point is near left bc hence querried neighbor exceeds left of left bc
										}
										else if(zc >= Z0_1+0.75*(ZM_1-Z0_1)) check = 3;//control point on the right quarter
										else if(zc < Z0_1+0.25*(ZM_1-Z0_1)) check = 4;//control point on the left quarter
										
										o = DOMID_1[i][j][kk];
										
										neighmax = EulNeigh_1[0][o];
										
										for(l=1;l<=neighmax;l++){
											check2 = 0;
											tri = EulNeigh_1[l][o];
											
											m = (int)((tri-1)/TriM);
											if(m > 188) continue;
											
											if(CellStatus[m] == -5000) continue;
											
											tri += -m*TriM;
											
											if(m == caps){
												n1 = TRINODE1[tri]; n2 = TRINODE2[tri]; n3 = TRINODE3[tri];
												if(n1 == node || n2 == node || n3 == node) check2 = 1;
											}
												
											if(check2 == 0){
												xp = TRI_CENTROID[0][m][tri];
												yp = TRI_CENTROID[1][m][tri];
												zp = TRI_CENTROID[2][m][tri];
												
												e2v[0] = TRI_NORMAL[0][m][tri];
												e2v[1] = TRI_NORMAL[1][m][tri];
												e2v[2] = TRI_NORMAL[2][m][tri];
												
												if(check == 1 && zp < 0.5*(ZM_1+Z0_1)){//control point is near right bc but neighbor querried is on the left bc
													zp = zp + (ZM_1-Z0_1);
												}
												if(check == 2 && zp > 0.5*(ZM_1+Z0_1)){//control point is near left bc but neighbor querried is on the right bc
													zp = zp - (ZM_1-Z0_1);
												}
												
												dotpdt = e1v[0]*e2v[0]+e1v[1]*e2v[1]+e1v[2]*e2v[2];
												
												if(dotpdt < 0){
												
													distance = 	sqrtf((xc-xp)*(xc-xp) + (yc-yp)*(yc-yp) + (zc-zp)*(zc-zp));
														
													if(distance < 2.){//this won't work if a particular region is sheared 3 times its orginal length
														n1 = TRINODE1[tri]; n2 = TRINODE2[tri]; n3 = TRINODE3[tri];
														
														for(v=0;v<=2;v++){
															n1coord[v] = NODECOORD[v][m][n1];
															n2coord[v] = NODECOORD[v][m][n2];
															n3coord[v] = NODECOORD[v][m][n3];
														}
														
														if(check == 1 || check == 3){//control point is near right bc, neighbors on the right could be virtually beyond right bc hence physically appearing on the right of left bc
															if(n1coord[2] < Z0_1+0.25*(ZM_1-Z0_1)) n1coord[2] += (ZM_1-Z0_1);
															if(n2coord[2] < Z0_1+0.25*(ZM_1-Z0_1)) n2coord[2] += (ZM_1-Z0_1);
															if(n3coord[2] < Z0_1+0.25*(ZM_1-Z0_1)) n3coord[2] += (ZM_1-Z0_1);
														}
														else if(check == 2 || check == 4){
															if(n1coord[2] >= Z0_1+0.75*(ZM_1-Z0_1)) n1coord[2] += -(ZM_1-Z0_1);
															if(n2coord[2] >= Z0_1+0.75*(ZM_1-Z0_1)) n2coord[2] += -(ZM_1-Z0_1);
															if(n3coord[2] >= Z0_1+0.75*(ZM_1-Z0_1)) n3coord[2] += -(ZM_1-Z0_1);
														}
														
														
														for(v=0;v<=2;v++){
															tp[1][v] = n1coord[v];
															tp[2][v] = n2coord[v];
															tp[3][v] = n3coord[v];
															
															tp[4][v] = 0.5*(tp[1][v]+tp[2][v]);
															tp[5][v] = 0.5*(tp[2][v]+tp[3][v]);
															tp[6][v] = 0.5*(tp[1][v]+tp[3][v]);
															
															tp[7][v] = 0.5*(tp[1][v]+tp[4][v]);
															tp[8][v] = 0.5*(tp[2][v]+tp[4][v]);
															tp[9][v] = 0.5*(tp[2][v]+tp[5][v]);
															tp[10][v] = 0.5*(tp[5][v]+tp[3][v]);
															tp[11][v] = 0.5*(tp[3][v]+tp[6][v]);
															tp[12][v] = 0.5*(tp[6][v]+tp[1][v]);
															
															tp[13][v] = 0.5*(tp[4][v]+tp[6][v]);
															tp[14][v] = 0.5*(tp[4][v]+tp[5][v]);
															tp[15][v] = 0.5*(tp[5][v]+tp[6][v]);
															
															
															tc[1][v] = (tp[1][v]+tp[7][v]+tp[12][v])/3.;
															tc[2][v] = (tp[7][v]+tp[13][v]+tp[12][v])/3.;
															tc[3][v] = (tp[7][v]+tp[4][v]+tp[13][v])/3.;
															tc[4][v] = (tp[4][v]+tp[14][v]+tp[13][v])/3.;
															tc[5][v] = (tp[4][v]+tp[8][v]+tp[14][v])/3.;
															tc[6][v] = (tp[8][v]+tp[9][v]+tp[14][v])/3.;
															tc[7][v] = (tp[8][v]+tp[2][v]+tp[9][v])/3.;
															
															tc[8][v] = (tp[12][v]+tp[13][v]+tp[6][v])/3.;
															tc[9][v] = (tp[13][v]+tp[15][v]+tp[6][v])/3.;
															tc[10][v] = (tp[13][v]+tp[14][v]+tp[15][v])/3.;
															tc[11][v] = (tp[14][v]+tp[5][v]+tp[15][v])/3.;
															tc[12][v] = (tp[14][v]+tp[9][v]+tp[5][v])/3.;
															
															tc[13][v] = (tp[6][v]+tp[15][v]+tp[11][v])/3.;
															tc[14][v] = (tp[15][v]+tp[10][v]+tp[11][v])/3.;
															tc[15][v] = (tp[15][v]+tp[5][v]+tp[10][v])/3.;
															tc[16][v] = (tp[11][v]+tp[10][v]+tp[3][v])/3.;
														}
															
														for(v=1;v<=16;v++){
															distance = sqrtf((tc[v][0]-xc)*(tc[v][0]-xc) + (tc[v][1]-yc)*(tc[v][1]-yc) + (tc[v][2]-zc)*(tc[v][2]-zc));
																
															if(distance < 1.0){
																e3v[0] = (tc[v][0]-xc)/distance;
																e3v[1] = (tc[v][1]-yc)/distance;
																e3v[2] = (tc[v][2]-zc)/distance;
																	
																repulf = -8.333e-07*(1.0 - distance)/1.0;
																	
																		
																for(c=0;c<=2;c++){
																	NODE_WALLREPUL[c][caps][node] += repulf*e3v[c];
																	
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}///end of ijk search stencil for eulerian list of neighbors
						}///end of NodeM for-loop for DA Cell nodes	
					}///if-condition for filtering only cell IDs in simulation domain
				}
			}
			
			
			/// PCV1 parent
			#pragma omp parallel private(checkn,Vneigh,ReflectProj,caps,node,mindisIn,mindisOut,xc,yc,zc,e1v,e2v,e3v,ic,jc,kc,i1,j1,k1,i2,j2,k2,i,j,k,check,kk,o,neighmax,l,check2,tri,m,n1,n2,n3,xp,yp,zp,distance,v,c,n1coord,n2coord,n3coord,ep,tp,tc,dotpdt,dotpdt2,Outneigh,Inneigh,repulf)
			{////$$$$$$$$$//////
				#pragma omp for schedule(static,3) nowait
				
				for(caps=189;caps<=329;caps++){
					if(CellStatus[caps] > -5000){
						
						for(node=1;node<=NodeM;node++){///NodeM for-loop for PCV2 Cell nodes
							///start of penetration condition check and updates
								
							
							xc = NODECOORD[0][caps][node];
							yc = NODECOORD[1][caps][node];
							zc = NODECOORD[2][caps][node];
							mindisIn = 5000.;
							mindisOut = 5000.;
								
							e1v[0] = NODE_NORMAL[0][caps][node];
							e1v[1] = NODE_NORMAL[1][caps][node];
							e1v[2] = NODE_NORMAL[2][caps][node];
								
							ic = (int)((xc-X0_2)/DX + 1.01);    jc = (int)((yc-Y0_2)/DX + 1.01);    kc = (int)((zc-Z0_2)/DX + 1.01);
							i1 = ic-4; i2 = ic+4;		j1 = jc-4; j2 = jc+4;		k1 = kc-4; k2 = kc+4;
							
							if(i1 < 1) i1 = 1; if(i2 > IM_2) i2 = IM_2;
							if(j1 < 1) j1 = 1; if(j2 > JM_2) j2 = JM_2;
							
								
								
							for(i=i1;i<=i2;i++){///start of ijk search stencil for eulerian list of neighbors
								for(j=j1;j<=j2;j++){
									for(k=k1;k<=k2;k++){
										check = 0;
										kk = k;
										if(k >= KM_2){// this means control point is near right bc hence querried neighbor exceeds right of right bc
											kk = k-KM_2+1;
											check = 1;
										}
										else if(k < 1){// this means control point is near left bc hence querried neighbor exceeds left of left bc
											kk = k+KM_2-1;
											check = 2;
										}
										else if(zc >= Z0_2+0.75*(ZM_2-Z0_2)) check = 3;
										else if(zc < Z0_2+0.25*(ZM_2-Z0_2)) check = 4;
											
										o = DOMID_2[i][j][kk];
											
										neighmax = EulNeigh_2[0][o];
											
										for(l=1;l<=neighmax;l++){
											check2 = 0;
											tri = EulNeigh_2[l][o];
												
											m = (int)((tri-1)/TriM);
											if(m < 189 || m > 329) continue;
												
											if(CellStatus[m] == -5000) continue;
												
											tri += -m*TriM;
												
											if(m == caps){
												n1 = TRINODE1[tri]; n2 = TRINODE2[tri]; n3 = TRINODE3[tri];
												if(n1 == node || n2 == node || n3 == node) check2 = 1;
											}
											
											if(check2 == 0){
												xp = TRI_CENTROID[0][m][tri];
												yp = TRI_CENTROID[1][m][tri];
												zp = TRI_CENTROID[2][m][tri];
													
												e2v[0] = TRI_NORMAL[0][m][tri];
												e2v[1] = TRI_NORMAL[1][m][tri];
												e2v[2] = TRI_NORMAL[2][m][tri];
													
													
												if(check == 1 && zp < 0.5*(ZM_2+Z0_2)){//control point is near right bc but neighbor querried is on the left bc
													zp = zp + (ZM_2-Z0_2);
												}
												if(check == 2 && zp > 0.5*(ZM_2+Z0_2)){//control point is near left bc but neighbor querried is on the right bc
													zp = zp - (ZM_2-Z0_2);
												}
													
												dotpdt = e1v[0]*e2v[0]+e1v[1]*e2v[1]+e1v[2]*e2v[2];
													
												if(dotpdt < 0){
													distance = 	sqrtf((xc-xp)*(xc-xp) + (yc-yp)*(yc-yp) + (zc-zp)*(zc-zp));
															
													if(distance < 2.){//this won't work if a particular region is sheared 3 times its orginal length
														n1 = TRINODE1[tri]; n2 = TRINODE2[tri]; n3 = TRINODE3[tri];
															
														for(v=0;v<=2;v++){
															n1coord[v] = NODECOORD[v][m][n1];
															n2coord[v] = NODECOORD[v][m][n2];
															n3coord[v] = NODECOORD[v][m][n3];
														}
														
																
														if(check == 1 || check == 3){//forward translate back half
															if(n1coord[2] < Z0_2+0.25*(ZM_2-Z0_2)) n1coord[2] += (ZM_2-Z0_2);
															if(n2coord[2] < Z0_2+0.25*(ZM_2-Z0_2)) n2coord[2] += (ZM_2-Z0_2);
															if(n3coord[2] < Z0_2+0.25*(ZM_2-Z0_2)) n3coord[2] += (ZM_2-Z0_2);
														}
														
														else if(check == 2 || check == 4){//backward translate front half
															if(n1coord[2] >= Z0_2+0.75*(ZM_2-Z0_2)) n1coord[2] += -(ZM_2-Z0_2);
															if(n2coord[2] >= Z0_2+0.75*(ZM_2-Z0_2)) n2coord[2] += -(ZM_2-Z0_2);
															if(n3coord[2] >= Z0_2+0.75*(ZM_2-Z0_2)) n3coord[2] += -(ZM_2-Z0_2);
														}
														
														for(v=0;v<=2;v++){
															tp[1][v] = n1coord[v];
															tp[2][v] = n2coord[v];
															tp[3][v] = n3coord[v];
															
															tp[4][v] = 0.5*(tp[1][v]+tp[2][v]);
															tp[5][v] = 0.5*(tp[2][v]+tp[3][v]);
															tp[6][v] = 0.5*(tp[1][v]+tp[3][v]);
															
															tp[7][v] = 0.5*(tp[1][v]+tp[4][v]);
															tp[8][v] = 0.5*(tp[2][v]+tp[4][v]);
															tp[9][v] = 0.5*(tp[2][v]+tp[5][v]);
															tp[10][v] = 0.5*(tp[5][v]+tp[3][v]);
															tp[11][v] = 0.5*(tp[3][v]+tp[6][v]);
															tp[12][v] = 0.5*(tp[6][v]+tp[1][v]);
															
															tp[13][v] = 0.5*(tp[4][v]+tp[6][v]);
															tp[14][v] = 0.5*(tp[4][v]+tp[5][v]);
															tp[15][v] = 0.5*(tp[5][v]+tp[6][v]);
															
															
															tc[1][v] = (tp[1][v]+tp[7][v]+tp[12][v])/3.;
															tc[2][v] = (tp[7][v]+tp[13][v]+tp[12][v])/3.;
															tc[3][v] = (tp[7][v]+tp[4][v]+tp[13][v])/3.;
															tc[4][v] = (tp[4][v]+tp[14][v]+tp[13][v])/3.;
															tc[5][v] = (tp[4][v]+tp[8][v]+tp[14][v])/3.;
															tc[6][v] = (tp[8][v]+tp[9][v]+tp[14][v])/3.;
															tc[7][v] = (tp[8][v]+tp[2][v]+tp[9][v])/3.;
															
															tc[8][v] = (tp[12][v]+tp[13][v]+tp[6][v])/3.;
															tc[9][v] = (tp[13][v]+tp[15][v]+tp[6][v])/3.;
															tc[10][v] = (tp[13][v]+tp[14][v]+tp[15][v])/3.;
															tc[11][v] = (tp[14][v]+tp[5][v]+tp[15][v])/3.;
															tc[12][v] = (tp[14][v]+tp[9][v]+tp[5][v])/3.;
															
															tc[13][v] = (tp[6][v]+tp[15][v]+tp[11][v])/3.;
															tc[14][v] = (tp[15][v]+tp[10][v]+tp[11][v])/3.;
															tc[15][v] = (tp[15][v]+tp[5][v]+tp[10][v])/3.;
															tc[16][v] = (tp[11][v]+tp[10][v]+tp[3][v])/3.;
														}
															
														for(v=1;v<=16;v++){
															distance = sqrtf((tc[v][0]-xc)*(tc[v][0]-xc) + (tc[v][1]-yc)*(tc[v][1]-yc) + (tc[v][2]-zc)*(tc[v][2]-zc));
																
															if(distance < 1.0){
																e3v[0] = (tc[v][0]-xc)/distance;
																e3v[1] = (tc[v][1]-yc)/distance;
																e3v[2] = (tc[v][2]-zc)/distance;
																	
																repulf = -8.333e-07*(1.0 - distance)/1.0;
																			
																for(c=0;c<=2;c++){
																	NODE_WALLREPUL[c][caps][node] += repulf*e3v[c];
																		
																}
															}
														}
														
													}
												}
											}
										}
									}
								}
							}///end of ijk search stencil for eulerian list of neighbors
						}///end of NodeM for-loop for PCV2 Cell nodes	
					}///if-condition filter for cells in simulation domain
				}///end of CapsM for-loop for PCV2 cells
			}
			////$$$$$$$$$//////
			
			/// PCV2 parent
			#pragma omp parallel private(checkn,Vneigh,ReflectProj,caps,node,mindisIn,mindisOut,xc,yc,zc,e1v,e2v,e3v,ic,jc,kc,i1,j1,k1,i2,j2,k2,i,j,k,check,kk,o,neighmax,l,check2,tri,m,n1,n2,n3,xp,yp,zp,distance,v,c,n1coord,n2coord,n3coord,ep,tp,tc,dotpdt,dotpdt2,Outneigh,Inneigh,repulf)
			{////$$$$$$$$$//////
				#pragma omp for schedule(static,1) nowait
				for(caps=330;caps<=344;caps++){
					if(CellStatus[caps] > -5000){
						for(node=1;node<=NodeM;node++){///NodeM for-loop for PCV2 Cell nodes
							///start of penetration condition check and updates
							
							xc = NODECOORD[0][caps][node];
							yc = NODECOORD[1][caps][node];
							zc = NODECOORD[2][caps][node];
							mindisIn = 5000.;
							mindisOut = 5000.;
								
							e1v[0] = NODE_NORMAL[0][caps][node];
							e1v[1] = NODE_NORMAL[1][caps][node];
							e1v[2] = NODE_NORMAL[2][caps][node];
								
							ic = (int)((xc-X0_3)/DX + 1.01);    jc = (int)((yc-Y0_3)/DX + 1.01);    kc = (int)((zc-Z0_3)/DX + 1.01);
							i1 = ic-4; i2 = ic+4;		j1 = jc-4; j2 = jc+4;		k1 = kc-4; k2 = kc+4;
							
							if(i1 < 1) i1 = 1; if(i2 > IM_3) i2 = IM_3;
							if(j1 < 1) j1 = 1; if(j2 > JM_3) j2 = JM_3;
								
								
							for(i=i1;i<=i2;i++){///start of ijk search stencil for eulerian list of neighbors
								for(j=j1;j<=j2;j++){
									for(k=k1;k<=k2;k++){
										check = 0;
										kk = k;
										if(k >= KM_3){// this means control point is near right bc hence querried neighbor exceeds right of right bc
											kk = k-KM_3+1;
											check = 1;
										}
										else if(k < 1){// this means control point is near left bc hence querried neighbor exceeds left of left bc
											kk = k+KM_3-1;
											check = 2;
										}
										else if(zc >= Z0_3+0.75*(ZM_3-Z0_3)) check = 3;
										else if(zc < Z0_3+0.25*(ZM_3-Z0_3)) check = 4;
										
										o = DOMID_3[i][j][kk];
										
										neighmax = EulNeigh_3[0][o];
										
										for(l=1;l<=neighmax;l++){
											check2 = 0;
											tri = EulNeigh_3[l][o];
											
											m = (int)((tri-1)/TriM);
											if(m < 330 || m > 344) continue;
											
											if(CellStatus[m] == -5000) continue;
											
											tri += -m*TriM;
											
											if(m == caps){
												n1 = TRINODE1[tri]; n2 = TRINODE2[tri]; n3 = TRINODE3[tri];
												if(n1 == node || n2 == node || n3 == node) check2 = 1;
											}
											
											if(check2 == 0){
												xp = TRI_CENTROID[0][m][tri];
												yp = TRI_CENTROID[1][m][tri];
												zp = TRI_CENTROID[2][m][tri];
												
												e2v[0] = TRI_NORMAL[0][m][tri];
												e2v[1] = TRI_NORMAL[1][m][tri];
												e2v[2] = TRI_NORMAL[2][m][tri];
												
												
												if(check == 1 && zp < 0.5*(ZM_3+Z0_3)){//control point is near right bc but neighbor querried is on the left bc
													zp = zp + (ZM_3-Z0_3);
												}
												if(check == 2 && zp > 0.5*(ZM_3+Z0_3)){//control point is near left bc but neighbor querried is on the right bc
													zp = zp - (ZM_3-Z0_3);
												}
												
												dotpdt = e1v[0]*e2v[0]+e1v[1]*e2v[1]+e1v[2]*e2v[2];
												
												if(dotpdt < 0){
													distance = 	sqrtf((xc-xp)*(xc-xp) + (yc-yp)*(yc-yp) + (zc-zp)*(zc-zp));
														
													if(distance < 2.){//this won't work if a particular region is sheared 3 times its orginal length
														n1 = TRINODE1[tri]; n2 = TRINODE2[tri]; n3 = TRINODE3[tri];
														
														for(v=0;v<=2;v++){
															n1coord[v] = NODECOORD[v][m][n1];
															n2coord[v] = NODECOORD[v][m][n2];
															n3coord[v] = NODECOORD[v][m][n3];
														}
														
																
														if(check == 1 || check == 3){//forward translate back half
															if(n1coord[2] < Z0_3+0.25*(ZM_3-Z0_3)) n1coord[2] += (ZM_3-Z0_3);
															if(n2coord[2] < Z0_3+0.25*(ZM_3-Z0_3)) n2coord[2] += (ZM_3-Z0_3);
															if(n3coord[2] < Z0_3+0.25*(ZM_3-Z0_3)) n3coord[2] += (ZM_3-Z0_3);
														}
														
														else if(check == 2 || check == 4){//backward translate front half
															if(n1coord[2] >= Z0_3+0.75*(ZM_3-Z0_3)) n1coord[2] += -(ZM_3-Z0_3);
															if(n2coord[2] >= Z0_3+0.75*(ZM_3-Z0_3)) n2coord[2] += -(ZM_3-Z0_3);
															if(n3coord[2] >= Z0_3+0.75*(ZM_3-Z0_3)) n3coord[2] += -(ZM_3-Z0_3);
														}
														
														for(v=0;v<=2;v++){
															tp[1][v] = n1coord[v];
															tp[2][v] = n2coord[v];
															tp[3][v] = n3coord[v];
															
															tp[4][v] = 0.5*(tp[1][v]+tp[2][v]);
															tp[5][v] = 0.5*(tp[2][v]+tp[3][v]);
															tp[6][v] = 0.5*(tp[1][v]+tp[3][v]);
															
															tp[7][v] = 0.5*(tp[1][v]+tp[4][v]);
															tp[8][v] = 0.5*(tp[2][v]+tp[4][v]);
															tp[9][v] = 0.5*(tp[2][v]+tp[5][v]);
															tp[10][v] = 0.5*(tp[5][v]+tp[3][v]);
															tp[11][v] = 0.5*(tp[3][v]+tp[6][v]);
															tp[12][v] = 0.5*(tp[6][v]+tp[1][v]);
															
															tp[13][v] = 0.5*(tp[4][v]+tp[6][v]);
															tp[14][v] = 0.5*(tp[4][v]+tp[5][v]);
															tp[15][v] = 0.5*(tp[5][v]+tp[6][v]);
															
															
															tc[1][v] = (tp[1][v]+tp[7][v]+tp[12][v])/3.;
															tc[2][v] = (tp[7][v]+tp[13][v]+tp[12][v])/3.;
															tc[3][v] = (tp[7][v]+tp[4][v]+tp[13][v])/3.;
															tc[4][v] = (tp[4][v]+tp[14][v]+tp[13][v])/3.;
															tc[5][v] = (tp[4][v]+tp[8][v]+tp[14][v])/3.;
															tc[6][v] = (tp[8][v]+tp[9][v]+tp[14][v])/3.;
															tc[7][v] = (tp[8][v]+tp[2][v]+tp[9][v])/3.;
															
															tc[8][v] = (tp[12][v]+tp[13][v]+tp[6][v])/3.;
															tc[9][v] = (tp[13][v]+tp[15][v]+tp[6][v])/3.;
															tc[10][v] = (tp[13][v]+tp[14][v]+tp[15][v])/3.;
															tc[11][v] = (tp[14][v]+tp[5][v]+tp[15][v])/3.;
															tc[12][v] = (tp[14][v]+tp[9][v]+tp[5][v])/3.;
															
															tc[13][v] = (tp[6][v]+tp[15][v]+tp[11][v])/3.;
															tc[14][v] = (tp[15][v]+tp[10][v]+tp[11][v])/3.;
															tc[15][v] = (tp[15][v]+tp[5][v]+tp[10][v])/3.;
															tc[16][v] = (tp[11][v]+tp[10][v]+tp[3][v])/3.;
														}
														
														for(v=1;v<=16;v++){
															distance = sqrtf((tc[v][0]-xc)*(tc[v][0]-xc) + (tc[v][1]-yc)*(tc[v][1]-yc) + (tc[v][2]-zc)*(tc[v][2]-zc));
															
															if(distance < 1.0){
																e3v[0] = (tc[v][0]-xc)/distance;
																e3v[1] = (tc[v][1]-yc)/distance;
																e3v[2] = (tc[v][2]-zc)/distance;
																
																repulf = -8.333e-07*(1.0 - distance)/1.0;
																			
																for(c=0;c<=2;c++){
																	NODE_WALLREPUL[c][caps][node] += repulf*e3v[c];
																		
																}
															}
														}
														
													}
												}
											}
										}
									}
								}
							}///end of ijk search stencil for eulerian list of neighbors
						}///end of NodeM for-loop for PCV2 Cell nodes	
					}///if-condition filter for cells in simulation domain
				}///end of CapsM for-loop for PCV2 cells
			}
			////$$$$$$$$$//////
		}	
		
				/// IBM_RBC() viscous force on LBM from velocity difference between membrane and immersed fluid (ROI and Slave Cells) /////////////////////////////////////////////////////////
		
		if(RBCSOLVERFLAG == 1){
			#pragma omp parallel private(o,wneighmax,nodewall,triwall,region,n1,n2,n3,xp1,yp1,zp1,xp2,yp2,zp2,xp3,yp3,zp3,xp4,yp4,zp4,xi,yi,zi,r23,ri3,r21,ri1,r31,area1,area2,area3,aratio,v,l,neighmax,tri,mindis2,xc2,yc2,zc2,xc1,yc1,zc1,iii,jjj,kkk,mindis,bvalue,dotpdt,m,caps,node,xp,yp,zp,i,j,k,i1,j1,k1,i2,j2,k2,weightsum,validinter,xc,yc,zc,ii,jj,kk,distance,interweight,repulf,ParentID,check,Ffluid) 
			{
				#pragma omp for schedule(static,50)  nowait
				for(caps=345;caps<=CapsM;caps++){
					if(CellStatus[caps] > 0){//Cloned ROI cells released for independent movement
						for(node=1;node<=NodeM;node++){
							xp = NODECOORD[0][caps][node]; yp = NODECOORD[1][caps][node]; zp = NODECOORD[2][caps][node];
							
							if(CellStatus[caps] == 1 || (zp > ZM_1 && zp < Z0_2)){
								NODE_VEL[0][caps][node] += NODE_WALLREPUL[0][caps][node];
								NODE_VEL[1][caps][node] += NODE_WALLREPUL[1][caps][node];
								NODE_VEL[2][caps][node] += NODE_WALLREPUL[2][caps][node];
								Ffluid[0] = (gamma_p*NODE_WALLREPUL[0][caps][node]+F_Total[0][caps][node]);//body force
								Ffluid[1] = (gamma_p*NODE_WALLREPUL[1][caps][node]+F_Total[1][caps][node]);//body force
								Ffluid[2] = (gamma_p*NODE_WALLREPUL[2][caps][node]+F_Total[2][caps][node]);//body force
								
								///interpolate velocities on membrane from surrounding fluid
								
								
								i = (int)((xp-X0)/DX + 1.01);
								j = (int)((yp-Y0)/DX + 1.01);
								k = (int)((zp-Z0)/DX + 1.01);
								
								i1 = i-4; i2 = i+4;
								j1 = j-4; j2 = j+4;
								k1 = k-4; k2 = k+4;
								
								if(i1 < 1) i1 = 1; if(i2> IM) i2 = IM;
								if(j1 < 1) j1 = 1; if(j2> JM) j2 = JM;
								
								weightsum = 0; bvalue = 0; mindis = 1.5; region = 5000; triwall = 0;
								for(i=i1;i<=i2;i++){
									for(j=j1;j<=j2;j++){
										for(k=k1;k<=k2;k++){
											
											xc = X0 + (i-1)*DX;
											yc = Y0 + (j-1)*DX;
											zc = Z0 + (k-1)*DX;
											
											if(k < 3 && REGRESSFLAG[caps] >= -188 && REGRESSFLAG[caps] < 0){
												ii = (int)((xc-X0_1)/DX + 1.01);
												jj = (int)((yc-Y0_1)/DX + 1.01);
												kk = (int)((zc-Z0_1)/DX + 1.01);
												
												o = DOMID_1[ii][jj][kk];
												
												if(B_index_1[ii][jj][kk] >= 0){
													distance = sqrtf((xc-xp)*(xc-xp) + (yc-yp)*(yc-yp) + (zc-zp)*(zc-zp));
													if(distance <= StencilWidth*DX){
														weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
														#pragma omp atomic update
														BODYFORCE_1[0][o] += weightsum*Ffluid[0]/UNITMASS/LBM_BFORCE_SCALE;
														#pragma omp atomic update
														BODYFORCE_1[1][o] += weightsum*Ffluid[1]/UNITMASS/LBM_BFORCE_SCALE;
														#pragma omp atomic update
														BODYFORCE_1[2][o] += weightsum*Ffluid[2]/UNITMASS/LBM_BFORCE_SCALE;
														
													}
												}
												
											}
											else if(k > KM-3 && REGRESSFLAG[caps] >= -329 && REGRESSFLAG[caps] <= -189){
												ii = (int)((xc-X0_2)/DX + 1.01);
												jj = (int)((yc-Y0_2)/DX + 1.01);
												kk = (int)((zc-Z0_2)/DX + 1.01);
												
												o = DOMID_2[ii][jj][kk];
												
												if(B_index_2[ii][jj][kk] >= 0){
													distance = sqrtf((xc-xp)*(xc-xp) + (yc-yp)*(yc-yp) + (zc-zp)*(zc-zp));
													if(distance <= StencilWidth*DX){
														weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
														#pragma omp atomic update
														BODYFORCE_2[0][o] += weightsum*Ffluid[0]/UNITMASS/LBM_BFORCE_SCALE;
														#pragma omp atomic update
														BODYFORCE_2[1][o] += weightsum*Ffluid[1]/UNITMASS/LBM_BFORCE_SCALE;
														#pragma omp atomic update
														BODYFORCE_2[2][o] += weightsum*Ffluid[2]/UNITMASS/LBM_BFORCE_SCALE;
														
													}
												}
												
											}
											else if(k > KM-3 && REGRESSFLAG[caps] >= -344 && REGRESSFLAG[caps] <= -330){
												ii = (int)((xc-X0_3)/DX + 1.01);
												jj = (int)((yc-Y0_3)/DX + 1.01);
												kk = (int)((zc-Z0_3)/DX + 1.01);
												
												o = DOMID_3[ii][jj][kk];
												
												if(B_index_3[ii][jj][kk] >= 0){
													distance = sqrtf((xc-xp)*(xc-xp) + (yc-yp)*(yc-yp) + (zc-zp)*(zc-zp));
													if(distance <= StencilWidth*DX){
														weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
														#pragma omp atomic update
														BODYFORCE_3[0][o] += weightsum*Ffluid[0]/UNITMASS/LBM_BFORCE_SCALE;
														#pragma omp atomic update
														BODYFORCE_3[1][o] += weightsum*Ffluid[1]/UNITMASS/LBM_BFORCE_SCALE;
														#pragma omp atomic update
														BODYFORCE_3[2][o] += weightsum*Ffluid[2]/UNITMASS/LBM_BFORCE_SCALE;
														
													}
												}
												
											}
											
											if(k >= 3 && k <= KM-3){
												ii = i;
												jj = j;
												kk = k;
												
												o = DOMID[ii][jj][kk];
												
												if(B_index[ii][jj][kk] >= 0){
													distance = sqrtf((xc-xp)*(xc-xp) + (yc-yp)*(yc-yp) + (zc-zp)*(zc-zp));
													if(distance <= StencilWidth*DX){
														weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
														#pragma omp atomic update
														BODYFORCE[0][o] += weightsum*Ffluid[0]/UNITMASS/LBM_BFORCE_SCALE;
														#pragma omp atomic update
														BODYFORCE[1][o] += weightsum*Ffluid[1]/UNITMASS/LBM_BFORCE_SCALE;
														#pragma omp atomic update
														BODYFORCE[2][o] += weightsum*Ffluid[2]/UNITMASS/LBM_BFORCE_SCALE;
														
													}
												}
												
											}
										}
									}
								}
							}
							
						}
						
						for(node=1;node<=NodeM_N;node++){
							xp = NODECOORD_N[0][caps][node]; yp = NODECOORD_N[1][caps][node]; zp = NODECOORD_N[2][caps][node];
							
							if(CellStatus[caps] == 1 || (zp > ZM_1 && zp < Z0_2)){
								NODE_VEL_N[0][caps][node] += NODE_WALLREPUL_N[0][caps][node];
								NODE_VEL_N[1][caps][node] += NODE_WALLREPUL_N[1][caps][node];
								NODE_VEL_N[2][caps][node] += NODE_WALLREPUL_N[2][caps][node];
								Ffluid[0] = (gamma_p*NODE_WALLREPUL_N[0][caps][node]+F_Total_N[0][caps][node]);//body force
								Ffluid[1] = (gamma_p*NODE_WALLREPUL_N[1][caps][node]+F_Total_N[1][caps][node]);//body force
								Ffluid[2] = (gamma_p*NODE_WALLREPUL_N[2][caps][node]+F_Total_N[2][caps][node]);//body force
								///interpolate velocities on membrane from surrounding fluid
								
								
								i = (int)((xp-X0)/DX + 1.01);
								j = (int)((yp-Y0)/DX + 1.01);
								k = (int)((zp-Z0)/DX + 1.01);
								
								
								i1 = i-4; i2 = i+4;
								j1 = j-4; j2 = j+4;
								k1 = k-4; k2 = k+4;
								
								if(i1 < 1) i1 = 1; if(i2> IM) i2 = IM;
								if(j1 < 1) j1 = 1; if(j2> JM) j2 = JM;
								
								weightsum = 0; bvalue = 0; 
								for(i=i1;i<=i2;i++){
									for(j=j1;j<=j2;j++){
										for(k=k1;k<=k2;k++){
											
											xc = X0 + (i-1)*DX;
											yc = Y0 + (j-1)*DX;
											zc = Z0 + (k-1)*DX;
											
											if(k < 3 && REGRESSFLAG[caps] >= -188 && REGRESSFLAG[caps] < 0){
												ii = (int)((xc-X0_1)/DX + 1.01);
												jj = (int)((yc-Y0_1)/DX + 1.01);
												kk = (int)((zc-Z0_1)/DX + 1.01);
												
												o = DOMID_1[ii][jj][kk];
												
												if(B_index_1[ii][jj][kk] >= 0){
													distance = sqrtf((xc-xp)*(xc-xp) + (yc-yp)*(yc-yp) + (zc-zp)*(zc-zp));
													if(distance <= StencilWidth*DX){
														weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
														#pragma omp atomic update
														BODYFORCE_1[0][o] += weightsum*Ffluid[0]/UNITMASS/LBM_BFORCE_SCALE;
														#pragma omp atomic update
														BODYFORCE_1[1][o] += weightsum*Ffluid[1]/UNITMASS/LBM_BFORCE_SCALE;
														#pragma omp atomic update
														BODYFORCE_1[2][o] += weightsum*Ffluid[2]/UNITMASS/LBM_BFORCE_SCALE;
														
													}
												}
											}
											else if(k > KM-3 && REGRESSFLAG[caps] >= -329 && REGRESSFLAG[caps] <= -189){
												ii = (int)((xc-X0_2)/DX + 1.01);
												jj = (int)((yc-Y0_2)/DX + 1.01);
												kk = (int)((zc-Z0_2)/DX + 1.01);
												
												o = DOMID_2[ii][jj][kk];
												
												if(B_index_2[ii][jj][kk] >= 0){
													distance = sqrtf((xc-xp)*(xc-xp) + (yc-yp)*(yc-yp) + (zc-zp)*(zc-zp));
													if(distance <= StencilWidth*DX){
														weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
														#pragma omp atomic update
														BODYFORCE_2[0][o] += weightsum*Ffluid[0]/UNITMASS/LBM_BFORCE_SCALE;
														#pragma omp atomic update
														BODYFORCE_2[1][o] += weightsum*Ffluid[1]/UNITMASS/LBM_BFORCE_SCALE;
														#pragma omp atomic update
														BODYFORCE_2[2][o] += weightsum*Ffluid[2]/UNITMASS/LBM_BFORCE_SCALE;
														
													}
												}
											}
											else if(k > KM-3 && REGRESSFLAG[caps] >= -344 && REGRESSFLAG[caps] <= -330){
												ii = (int)((xc-X0_3)/DX + 1.01);
												jj = (int)((yc-Y0_3)/DX + 1.01);
												kk = (int)((zc-Z0_3)/DX + 1.01);
												
												o = DOMID_3[ii][jj][kk];
												
												if(B_index_3[ii][jj][kk] >= 0){
													distance = sqrtf((xc-xp)*(xc-xp) + (yc-yp)*(yc-yp) + (zc-zp)*(zc-zp));
													if(distance <= StencilWidth*DX){
														weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
														#pragma omp atomic update
														BODYFORCE_3[0][o] += weightsum*Ffluid[0]/UNITMASS/LBM_BFORCE_SCALE;
														#pragma omp atomic update
														BODYFORCE_3[1][o] += weightsum*Ffluid[1]/UNITMASS/LBM_BFORCE_SCALE;
														#pragma omp atomic update
														BODYFORCE_3[2][o] += weightsum*Ffluid[2]/UNITMASS/LBM_BFORCE_SCALE;
														
													}
												}
											}
											
											if(k >= 3 && k <= KM-3){
												ii = i;
												jj = j;
												kk = k;
												
												o = DOMID[ii][jj][kk];
												
												if(B_index[ii][jj][kk] >= 0){
													distance = sqrtf((xc-xp)*(xc-xp) + (yc-yp)*(yc-yp) + (zc-zp)*(zc-zp));
													if(distance <= StencilWidth*DX){
														weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
														#pragma omp atomic update
														BODYFORCE[0][o] += weightsum*Ffluid[0]/UNITMASS/LBM_BFORCE_SCALE;
														#pragma omp atomic update
														BODYFORCE[1][o] += weightsum*Ffluid[1]/UNITMASS/LBM_BFORCE_SCALE;
														#pragma omp atomic update
														BODYFORCE[2][o] += weightsum*Ffluid[2]/UNITMASS/LBM_BFORCE_SCALE;
														
													}
												}
											}
										}
									}
								}
							}		
						}
						
						//}
					}///Slave cells still following parent velocity
					else if(CellStatus[caps] > -4998){//cloned slave cells
						ParentID = -CellStatus[caps];
						for(node=1;node<=NodeM;node++){
							xp = NODECOORD[0][caps][node]; yp = NODECOORD[1][caps][node]; zp = NODECOORD[2][caps][node];
							
							check = 0;
							if(ParentID <= 188 && zp > ZM_1+5.*DX) check = 1;//region with independent movement
							if(ParentID >= 189 && ParentID <= 329 && zp < Z0_2-5.*DX) check = 1;//region with independent movement
							if(ParentID >= 330 && ParentID <= 344 && zp < Z0_3-5.*DX) check = 1;//region with independent movement
							
								
									/// cloned slaves and parents should follow same velocity - update parent repulsion vel from clone data then copy parent value back to clone for FSI update
									
									if(check == 0) NODE_WALLREPUL[0][caps][node] = NODE_WALLREPUL[0][ParentID][node];
									if(check == 0) NODE_WALLREPUL[1][caps][node] = NODE_WALLREPUL[1][ParentID][node];
									if(check == 0) NODE_WALLREPUL[2][caps][node] = NODE_WALLREPUL[2][ParentID][node];
									
									Ffluid[0] = (gamma_p*NODE_WALLREPUL[0][caps][node]+F_Total[0][caps][node]);//body force
									Ffluid[1] = (gamma_p*NODE_WALLREPUL[1][caps][node]+F_Total[1][caps][node]);//body force
									Ffluid[2] = (gamma_p*NODE_WALLREPUL[2][caps][node]+F_Total[2][caps][node]);//body force		
									
									if(check == 1){	
										i = (int)((xp-X0)/DX + 1.01);
										j = (int)((yp-Y0)/DX + 1.01);
										k = (int)((zp-Z0)/DX + 1.01);
										
										i1 = i-4; i2 = i+4;
										j1 = j-4; j2 = j+4;
										k1 = k-4; k2 = k+4;
										
										if(i1 < 1) i1 = 1; if(i2> IM) i2 = IM;
										if(j1 < 1) j1 = 1; if(j2> JM) j2 = JM;
										if(k1 < 3) k1 = 3; if(k2> KM-3) k2 = KM-3;
										
										
										weightsum = 0; bvalue = 0; mindis = 1.5; region = 5000; triwall = 0;
										for(i=i1;i<=i2;i++){
											for(j=j1;j<=j2;j++){
												for(k=k1;k<=k2;k++){
													
													xc = X0 + (i-1)*DX;
													yc = Y0 + (j-1)*DX;
													zc = Z0 + (k-1)*DX;
													
													
													if(k >= 3 && k <= KM-3){
														
														ii = i;
														jj = j;
														kk = k;
														
														o = DOMID[ii][jj][kk];
														
														if(B_index[ii][jj][kk] >= 0){
															distance = sqrtf((xc-xp)*(xc-xp) + (yc-yp)*(yc-yp) + (zc-zp)*(zc-zp));
															if(distance <= StencilWidth*DX){
																weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
																#pragma omp atomic update
																BODYFORCE[0][o] += weightsum*Ffluid[0]/UNITMASS/LBM_BFORCE_SCALE;
																#pragma omp atomic update
																BODYFORCE[1][o] += weightsum*Ffluid[1]/UNITMASS/LBM_BFORCE_SCALE;
																#pragma omp atomic update
																BODYFORCE[2][o] += weightsum*Ffluid[2]/UNITMASS/LBM_BFORCE_SCALE;
																
															}
														}
														
													}
												}
											}
										}
									}
								//}
							
								
						}
						
						
								
							for(node=1;node<=NodeM_N;node++){
								if(CellStatus[caps] == -5000){
									for(node=1;node<=NodeM;node++){
										for(i=0;i<=2;i++){
											F_Total[i][caps][node] = 0;
												NODECOORD[i][caps][node] = -5000.;
										}
									}
									for(node=1;node<=NodeM_N;node++){
										for(i=0;i<=2;i++){
											F_Total_N[i][caps][node] = 0;
											NODECOORD_N[i][caps][node] = -5000.;
										}
									}
									break;
								}
								xp = NODECOORD_N[0][caps][node]; yp = NODECOORD_N[1][caps][node]; zp = NODECOORD_N[2][caps][node];
								
								
								check = 0;
								if(ParentID <= 188 && zp > ZM_1+5.*DX) check = 1;
								if(ParentID >= 189 && ParentID <= 329 && zp < Z0_2-5.*DX) check = 1;
								if(ParentID >= 330 && ParentID <= 344 && zp < Z0_3-5.*DX) check = 1;
								
										if(check == 0) NODE_WALLREPUL_N[0][caps][node] = NODE_WALLREPUL_N[0][ParentID][node];
										if(check == 0) NODE_WALLREPUL_N[1][caps][node] = NODE_WALLREPUL_N[1][ParentID][node];
										if(check == 0) NODE_WALLREPUL_N[2][caps][node] = NODE_WALLREPUL_N[2][ParentID][node];
										
										Ffluid[0] = (gamma_p*NODE_WALLREPUL_N[0][caps][node]+F_Total_N[0][caps][node]);//body force
										Ffluid[1] = (gamma_p*NODE_WALLREPUL_N[1][caps][node]+F_Total_N[1][caps][node]);//body force
										Ffluid[2] = (gamma_p*NODE_WALLREPUL_N[2][caps][node]+F_Total_N[2][caps][node]);//body force
										
										if(check == 1){
											i = (int)((xp-X0)/DX + 1.01);
											j = (int)((yp-Y0)/DX + 1.01);
											k = (int)((zp-Z0)/DX + 1.01);
											
											i1 = i-4; i2 = i+4;
											j1 = j-4; j2 = j+4;
											k1 = k-4; k2 = k+4;
											
											if(i1 < 1) i1 = 1; if(i2> IM) i2 = IM;
											if(j1 < 1) j1 = 1; if(j2> JM) j2 = JM;
											if(k1 < 3) k1 = 3; if(k2> KM-3) k2 = KM-3;
										
											weightsum = 0; bvalue = 0;
											for(i=i1;i<=i2;i++){
												for(j=j1;j<=j2;j++){
													for(k=k1;k<=k2;k++){
														
														xc = X0 + (i-1)*DX;
														yc = Y0 + (j-1)*DX;
														zc = Z0 + (k-1)*DX;
														
														if(k >= 3 && k <= KM-3){	
															ii = i;
															jj = j;
															kk = k;
															
															o = DOMID[ii][jj][kk];
															
															if(B_index[ii][jj][kk] >= 0){
																distance = sqrtf((xc-xp)*(xc-xp) + (yc-yp)*(yc-yp) + (zc-zp)*(zc-zp));
																if(distance <= StencilWidth*DX){
																	weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
																	#pragma omp atomic update
																	BODYFORCE[0][o] += weightsum*Ffluid[0]/UNITMASS/LBM_BFORCE_SCALE;
																	#pragma omp atomic update
																	BODYFORCE[1][o] += weightsum*Ffluid[1]/UNITMASS/LBM_BFORCE_SCALE;
																	#pragma omp atomic update
																	BODYFORCE[2][o] += weightsum*Ffluid[2]/UNITMASS/LBM_BFORCE_SCALE;
																	
																}
															}
														}											
													}
												}
											}
										}
										///	
									//}
								//}
								
							}	
							
						//}
					}
					
				
				}
			}
			
		}
		
		
		
		///IBM_RBC() viscous force on LBM from velocity difference between membrane and immersed fluid (Parent Cells)
		
		if(RBCSOLVERFLAG == 1){
			#pragma omp parallel private(o,Ffluid,caps,node,i,xp,yp,zp,j,k,i1,j1,k1,i2,j2,k2,weightsum,weightcount,bvalue,mindis,triwall,xc,yc,zc,kmap,validinter,distance,interweight,wneighmax,m,nodewall,neighmax,l,tri,n1,n2,n3,xp1,yp1,zp1,repulf,v,xp2,yp2,zp2,xp3,yp3,zp3,xp4,yp4,zp4,zpmap)
			{
				#pragma omp for schedule(static,5)  nowait
				for(caps=1;caps<=344;caps++){
					if(CellStatus[caps] == 2){//Parent cells in parent interior domains
						//DAparent
						if(caps <= 188){
							for(node=1;node<=NodeM;node++){
								NODE_VEL[0][caps][node] += NODE_WALLREPUL[0][caps][node];
								NODE_VEL[1][caps][node] += NODE_WALLREPUL[1][caps][node];
								NODE_VEL[2][caps][node] += NODE_WALLREPUL[2][caps][node];
								
								Ffluid[0] = (gamma_p*NODE_WALLREPUL[0][caps][node]+F_Total[0][caps][node]);//body force
								Ffluid[1] = (gamma_p*NODE_WALLREPUL[1][caps][node]+F_Total[1][caps][node]);//body force
								Ffluid[2] = (gamma_p*NODE_WALLREPUL[2][caps][node]+F_Total[2][caps][node]);//body force
								
								///interpolate velocities on membrane from surrounding fluid
								xp = NODECOORD[0][caps][node]; yp = NODECOORD[1][caps][node]; zp = NODECOORD[2][caps][node];
								
								i = (int)((xp-X0_1)/DX + 1.01);
								j = (int)((yp-Y0_1)/DX + 1.01);
								k = (int)((zp-Z0_1)/DX + 1.01);
								
								i1 = i-4; i2 = i+4;
								j1 = j-4; j2 = j+4;
								k1 = k-4; k2 = k+4;
								
								if(i1 < 1) i1 = 1; if(i2> IM_1) i2 = IM_1;
								if(j1 < 1) j1 = 1; if(j2> JM_1) j2 = JM_1;
								
								
								weightsum = 0; weightcount = 0; bvalue = 0; mindis = 1.5; triwall = 0;
								for(i=i1;i<=i2;i++){
									for(j=j1;j<=j2;j++){
										for(k=k1;k<=k2;k++){
											
											xc = X0_1 + (i-1)*DX;
											yc = Y0_1 + (j-1)*DX;
											zc = Z0_1 + (k-1)*DX;
											
											if(k < 1) kmap = k + (KM_1-1);
											else if(k > KM_1) kmap = k - (KM_1-1);
											else kmap = k;
											
											o = DOMID_1[i][j][kmap];
											
											if(B_index_1[i][j][kmap] >= 0){
												
												distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
												if(distance <= StencilWidth*DX){
													
													weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
													#pragma omp atomic update
													BODYFORCE_1[0][o] += weightsum*Ffluid[0]/UNITMASS/LBM_BFORCE_SCALE;
													#pragma omp atomic update
													BODYFORCE_1[1][o] += weightsum*Ffluid[1]/UNITMASS/LBM_BFORCE_SCALE;
													#pragma omp atomic update
													BODYFORCE_1[2][o] += weightsum*Ffluid[2]/UNITMASS/LBM_BFORCE_SCALE;
													
												}
											}
											
											
										}
									}
								}
								
								
							}
						
							for(node=1;node<=NodeM_N;node++){
								NODE_VEL_N[0][caps][node] += NODE_WALLREPUL_N[0][caps][node];
								NODE_VEL_N[1][caps][node] += NODE_WALLREPUL_N[1][caps][node];
								NODE_VEL_N[2][caps][node] += NODE_WALLREPUL_N[2][caps][node];
								
								Ffluid[0] = (gamma_p*NODE_WALLREPUL_N[0][caps][node]+F_Total_N[0][caps][node]);//body force
								Ffluid[1] = (gamma_p*NODE_WALLREPUL_N[1][caps][node]+F_Total_N[1][caps][node]);//body force
								Ffluid[2] = (gamma_p*NODE_WALLREPUL_N[2][caps][node]+F_Total_N[2][caps][node]);//body force
								///interpolate velocities on membrane from surrounding fluid
								xp = NODECOORD_N[0][caps][node]; yp = NODECOORD_N[1][caps][node]; zp = NODECOORD_N[2][caps][node];
								
								i = (int)((xp-X0_1)/DX + 1.01);
								j = (int)((yp-Y0_1)/DX + 1.01);
								k = (int)((zp-Z0_1)/DX + 1.01);
								
								
								i1 = i-4; i2 = i+4;
								j1 = j-4; j2 = j+4;
								k1 = k-4; k2 = k+4;
								
								if(i1 < 1) i1 = 1; if(i2> IM_1) i2 = IM_1;
								if(j1 < 1) j1 = 1; if(j2> JM_1) j2 = JM_1;
								
								weightsum = 0; weightcount = 0; bvalue = 0;
								for(i=i1;i<=i2;i++){
									for(j=j1;j<=j2;j++){
										for(k=k1;k<=k2;k++){
											
											xc = X0_1 + (i-1)*DX;
											yc = Y0_1 + (j-1)*DX;
											zc = Z0_1 + (k-1)*DX;
											
											if(k < 1) kmap = k + (KM_1-1);
											else if(k > KM_1) kmap = k - (KM_1-1);
											else kmap = k;
											
											o = DOMID_1[i][j][kmap];
											
											if(B_index_1[i][j][kmap] >= 0){
												
												distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
												if(distance <= StencilWidth*DX){
													weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
													#pragma omp atomic update
													BODYFORCE_1[0][o] += weightsum*Ffluid[0]/UNITMASS/LBM_BFORCE_SCALE;
													#pragma omp atomic update
													BODYFORCE_1[1][o] += weightsum*Ffluid[1]/UNITMASS/LBM_BFORCE_SCALE;
													#pragma omp atomic update
													BODYFORCE_1[2][o] += weightsum*Ffluid[2]/UNITMASS/LBM_BFORCE_SCALE;
													
												}
											}
										}
									}
								}
								
								
							}
						
						}
							//PCV1 parent
						if(caps >= 189 && caps <= 329){
							for(node=1;node<=NodeM;node++){
								NODE_VEL[0][caps][node] += NODE_WALLREPUL[0][caps][node];
								NODE_VEL[1][caps][node] += NODE_WALLREPUL[1][caps][node];
								NODE_VEL[2][caps][node] += NODE_WALLREPUL[2][caps][node];
								
								Ffluid[0] = (gamma_p*NODE_WALLREPUL[0][caps][node]+F_Total[0][caps][node]);//body force
								Ffluid[1] = (gamma_p*NODE_WALLREPUL[1][caps][node]+F_Total[1][caps][node]);//body force
								Ffluid[2] = (gamma_p*NODE_WALLREPUL[2][caps][node]+F_Total[2][caps][node]);//body force
								///interpolate velocities on membrane from surrounding fluid
								xp = NODECOORD[0][caps][node]; yp = NODECOORD[1][caps][node]; zp = NODECOORD[2][caps][node];
								
								i = (int)((xp-X0_2)/DX + 1.01);
								j = (int)((yp-Y0_2)/DX + 1.01);
								k = (int)((zp-Z0_2)/DX + 1.01);
								
								
								i1 = i-4; i2 = i+4;
								j1 = j-4; j2 = j+4;
								k1 = k-4; k2 = k+4;
								
								if(i1 < 1) i1 = 1; if(i2> IM_2) i2 = IM_2;
								if(j1 < 1) j1 = 1; if(j2> JM_2) j2 = JM_2;
								
								weightsum = 0; weightcount = 0; bvalue = 0; mindis = 1.5; triwall = 0;
								for(i=i1;i<=i2;i++){
									for(j=j1;j<=j2;j++){
										for(k=k1;k<=k2;k++){
											xc = X0_2 + (i-1)*DX;
											yc = Y0_2 + (j-1)*DX;
											zc = Z0_2 + (k-1)*DX;
											
											if(k < 1) kmap = k + (KM_2-1);
											else if(k > KM_2) kmap = k - (KM_2-1);
											else kmap = k;
											
											o = DOMID_2[i][j][kmap];
											
											if(B_index_2[i][j][kmap] >= 0){
												
												distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
												if(distance <= StencilWidth*DX){
													weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
													#pragma omp atomic update
													BODYFORCE_2[0][o] += weightsum*Ffluid[0]/UNITMASS/LBM_BFORCE_SCALE;
													#pragma omp atomic update
													BODYFORCE_2[1][o] += weightsum*Ffluid[1]/UNITMASS/LBM_BFORCE_SCALE;
													#pragma omp atomic update
													BODYFORCE_2[2][o] += weightsum*Ffluid[2]/UNITMASS/LBM_BFORCE_SCALE;
													
												}
												
											}
											
										}
									}
								}
								
								
							}
						
							for(node=1;node<=NodeM_N;node++){
								NODE_VEL_N[0][caps][node] += NODE_WALLREPUL_N[0][caps][node];
								NODE_VEL_N[1][caps][node] += NODE_WALLREPUL_N[1][caps][node];
								NODE_VEL_N[2][caps][node] += NODE_WALLREPUL_N[2][caps][node];
								
								Ffluid[0] = (gamma_p*NODE_WALLREPUL_N[0][caps][node]+F_Total_N[0][caps][node]);//body force
								Ffluid[1] = (gamma_p*NODE_WALLREPUL_N[1][caps][node]+F_Total_N[1][caps][node]);//body force
								Ffluid[2] = (gamma_p*NODE_WALLREPUL_N[2][caps][node]+F_Total_N[2][caps][node]);//body force
								///interpolate velocities on membrane from surrounding fluid
								xp = NODECOORD_N[0][caps][node]; yp = NODECOORD_N[1][caps][node]; zp = NODECOORD_N[2][caps][node];
								
								i = (int)((xp-X0_2)/DX + 1.01);
								j = (int)((yp-Y0_2)/DX + 1.01);
								k = (int)((zp-Z0_2)/DX + 1.01);
								
								
								i1 = i-4; i2 = i+4;
								j1 = j-4; j2 = j+4;
								k1 = k-4; k2 = k+4;
								
								if(i1 < 1) i1 = 1; if(i2> IM_2) i2 = IM_2;
								if(j1 < 1) j1 = 1; if(j2> JM_2) j2 = JM_2;
								
								weightsum = 0; weightcount = 0; bvalue = 0;
								for(i=i1;i<=i2;i++){
									for(j=j1;j<=j2;j++){
										for(k=k1;k<=k2;k++){
											xc = X0_2 + (i-1)*DX;
											yc = Y0_2 + (j-1)*DX;
											zc = Z0_2 + (k-1)*DX;
											
											if(k < 1) kmap = k + (KM_2-1);
											else if(k > KM_2) kmap = k - (KM_2-1);
											else kmap = k;
											
											o = DOMID_2[i][j][kmap];
											
											if(B_index_2[i][j][kmap] >= 0){
												
												distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
												if(distance <= StencilWidth*DX){
													weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
													#pragma omp atomic update
													BODYFORCE_2[0][o] += weightsum*Ffluid[0]/UNITMASS/LBM_BFORCE_SCALE;
													#pragma omp atomic update
													BODYFORCE_2[1][o] += weightsum*Ffluid[1]/UNITMASS/LBM_BFORCE_SCALE;
													#pragma omp atomic update
													BODYFORCE_2[2][o] += weightsum*Ffluid[2]/UNITMASS/LBM_BFORCE_SCALE;
													
												}
											}
										}
									}
								}
								
							}
						
						}
						//PCV2 parent
						if(caps >= 330 && caps <= 344){
							for(node=1;node<=NodeM;node++){
								NODE_VEL[0][caps][node] += NODE_WALLREPUL[0][caps][node];
								NODE_VEL[1][caps][node] += NODE_WALLREPUL[1][caps][node];
								NODE_VEL[2][caps][node] += NODE_WALLREPUL[2][caps][node];
								
								Ffluid[0] = (gamma_p*NODE_WALLREPUL[0][caps][node]+F_Total[0][caps][node]);//body force
								Ffluid[1] = (gamma_p*NODE_WALLREPUL[1][caps][node]+F_Total[1][caps][node]);//body force
								Ffluid[2] = (gamma_p*NODE_WALLREPUL[2][caps][node]+F_Total[2][caps][node]);//body force
								///interpolate velocities on membrane from surrounding fluid
								xp = NODECOORD[0][caps][node]; yp = NODECOORD[1][caps][node]; zp = NODECOORD[2][caps][node];
								
								i = (int)((xp-X0_3)/DX + 1.01);
								j = (int)((yp-Y0_3)/DX + 1.01);
								k = (int)((zp-Z0_3)/DX + 1.01);
								
								
								i1 = i-4; i2 = i+4;
								j1 = j-4; j2 = j+4;
								k1 = k-4; k2 = k+4;
								
								if(i1 < 1) i1 = 1; if(i2> IM_3) i2 = IM_3;
								if(j1 < 1) j1 = 1; if(j2> JM_3) j2 = JM_3;
								
								weightsum = 0; weightcount = 0; bvalue = 0; mindis = 1.5; triwall = 0;
								for(i=i1;i<=i2;i++){
									for(j=j1;j<=j2;j++){
										for(k=k1;k<=k2;k++){
											xc = X0_3 + (i-1)*DX;
											yc = Y0_3 + (j-1)*DX;
											zc = Z0_3 + (k-1)*DX;
											
											if(k < 1) kmap = k + (KM_3-1);
											else if(k > KM_3) kmap = k - (KM_3-1);
											else kmap = k;
										
											o = DOMID_3[i][j][kmap];
											
											if(B_index_3[i][j][kmap] >= 0){
												distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
												if(distance <= StencilWidth*DX){
													weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
													#pragma omp atomic update
													BODYFORCE_3[0][o] += weightsum*Ffluid[0]/UNITMASS/LBM_BFORCE_SCALE;
													#pragma omp atomic update
													BODYFORCE_3[1][o] += weightsum*Ffluid[1]/UNITMASS/LBM_BFORCE_SCALE;
													#pragma omp atomic update
													BODYFORCE_3[2][o] += weightsum*Ffluid[2]/UNITMASS/LBM_BFORCE_SCALE;
													
												}
											}
											
										}
									}
								}
								
								
							}
						
							for(node=1;node<=NodeM_N;node++){
								NODE_VEL_N[0][caps][node] += NODE_WALLREPUL_N[0][caps][node];
								NODE_VEL_N[1][caps][node] += NODE_WALLREPUL_N[1][caps][node];
								NODE_VEL_N[2][caps][node] += NODE_WALLREPUL_N[2][caps][node];
								
								Ffluid[0] = (gamma_p*NODE_WALLREPUL_N[0][caps][node]+F_Total_N[0][caps][node]);//body force
								Ffluid[1] = (gamma_p*NODE_WALLREPUL_N[1][caps][node]+F_Total_N[1][caps][node]);//body force
								Ffluid[2] = (gamma_p*NODE_WALLREPUL_N[2][caps][node]+F_Total_N[2][caps][node]);//body force
								///interpolate velocities on membrane from surrounding fluid
								xp = NODECOORD_N[0][caps][node]; yp = NODECOORD_N[1][caps][node]; zp = NODECOORD_N[2][caps][node];
								
								i = (int)((xp-X0_3)/DX + 1.01);
								j = (int)((yp-Y0_3)/DX + 1.01);
								k = (int)((zp-Z0_3)/DX + 1.01);
								
								
								i1 = i-4; i2 = i+4;
								j1 = j-4; j2 = j+4;
								k1 = k-4; k2 = k+4;
								
								if(i1 < 1) i1 = 1; if(i2> IM_3) i2 = IM_3;
								if(j1 < 1) j1 = 1; if(j2> JM_3) j2 = JM_3;
								
								weightsum = 0; weightcount = 0; bvalue = 0;
								for(i=i1;i<=i2;i++){
									for(j=j1;j<=j2;j++){
										for(k=k1;k<=k2;k++){
											xc = X0_3 + (i-1)*DX;
											yc = Y0_3 + (j-1)*DX;
											zc = Z0_3 + (k-1)*DX;
											
											if(k < 1) kmap = k + (KM_3-1);
											else if(k > KM_3) kmap = k - (KM_3-1);
											else kmap = k;
											
											o = DOMID_3[i][j][kmap];
											
											if(B_index_3[i][j][kmap] >= 0){
												distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
												if(distance <= StencilWidth*DX){
													weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
													#pragma omp atomic update
													BODYFORCE_3[0][o] += weightsum*Ffluid[0]/UNITMASS/LBM_BFORCE_SCALE;
													#pragma omp atomic update
													BODYFORCE_3[1][o] += weightsum*Ffluid[1]/UNITMASS/LBM_BFORCE_SCALE;
													#pragma omp atomic update
													BODYFORCE_3[2][o] += weightsum*Ffluid[2]/UNITMASS/LBM_BFORCE_SCALE;
													
												}
											}
										}
									}
								}
								
							}
						}
					}
					else if(CellStatus[caps] > 2){//Parent cells crossing periodic boundaries
						if(caps <= 188){
							for(node=1;node<=NodeM;node++){
								NODE_VEL[0][caps][node] += NODE_WALLREPUL[0][caps][node];
								NODE_VEL[1][caps][node] += NODE_WALLREPUL[1][caps][node];
								NODE_VEL[2][caps][node] += NODE_WALLREPUL[2][caps][node];
								
								Ffluid[0] = (gamma_p*NODE_WALLREPUL[0][caps][node]+F_Total[0][caps][node]);//body force
								Ffluid[1] = (gamma_p*NODE_WALLREPUL[1][caps][node]+F_Total[1][caps][node]);//body force
								Ffluid[2] = (gamma_p*NODE_WALLREPUL[2][caps][node]+F_Total[2][caps][node]);//body force
								///interpolate velocities on membrane from surrounding fluid
								xp = NODECOORD[0][caps][node]; yp = NODECOORD[1][caps][node]; zp = NODECOORD[2][caps][node];
								
								i = (int)((xp-X0_1)/DX + 1.01);
								j = (int)((yp-Y0_1)/DX + 1.01);
								k = (int)((zp-Z0_1)/DX + 1.01);
											
											
								i = (int)((xp-X0_1)/DX + 1.01);
								j = (int)((yp-Y0_1)/DX + 1.01);
								k = (int)((zp-Z0_1)/DX + 1.01);
								
								i1 = i-4; i2 = i+4;
								j1 = j-4; j2 = j+4;
								k1 = k-4; k2 = k+4;
								if(i1 < 1) i1 = 1; if(i2 > IM_1) i2 = IM_1;
								if(j1 < 1) j1 = 1; if(j2 > JM_1) j2 = JM_1;
								
								weightsum = 0; weightcount = 0; bvalue = 0; mindis = 1.5; triwall = 0;
								for(i=i1;i<=i2;i++){
									for(j=j1;j<=j2;j++){
										for(k=k1;k<=k2;k++){
											xc = X0_1 + (i-1)*DX;
											yc = Y0_1 + (j-1)*DX;
											zc = Z0_1 + (k-1)*DX;
											
											
											if(k > KM_1){
												kmap = k - (KM_1-1);
											}
											else if(k < 1) kmap = k + (KM_1-1);
											else kmap = k;
											
											o = DOMID_1[i][j][kmap];
											
											if(B_index_1[i][j][kmap] >= 0){
												distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
												if(distance <= StencilWidth*DX){
													weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
													#pragma omp atomic update
													BODYFORCE_1[0][o] += weightsum*Ffluid[0]/UNITMASS/LBM_BFORCE_SCALE;
													#pragma omp atomic update
													BODYFORCE_1[1][o] += weightsum*Ffluid[1]/UNITMASS/LBM_BFORCE_SCALE;
													#pragma omp atomic update
													BODYFORCE_1[2][o] += weightsum*Ffluid[2]/UNITMASS/LBM_BFORCE_SCALE;
													
												}
											}
											
										}
									}
								}
								
							}
							
							for(node=1;node<=NodeM_N;node++){
								NODE_VEL_N[0][caps][node] += NODE_WALLREPUL_N[0][caps][node];
								NODE_VEL_N[1][caps][node] += NODE_WALLREPUL_N[1][caps][node];
								NODE_VEL_N[2][caps][node] += NODE_WALLREPUL_N[2][caps][node];
								
								Ffluid[0] = (gamma_p*NODE_WALLREPUL_N[0][caps][node]+F_Total_N[0][caps][node]);//body force
								Ffluid[1] = (gamma_p*NODE_WALLREPUL_N[1][caps][node]+F_Total_N[1][caps][node]);//body force
								Ffluid[2] = (gamma_p*NODE_WALLREPUL_N[2][caps][node]+F_Total_N[2][caps][node]);//body force
								///interpolate velocities on membrane from surrounding fluid
								xp = NODECOORD_N[0][caps][node]; yp = NODECOORD_N[1][caps][node]; zp = NODECOORD_N[2][caps][node];
								
																
								i = (int)((xp-X0_1)/DX + 1.01);
								j = (int)((yp-Y0_1)/DX + 1.01);
								k = (int)((zp-Z0_1)/DX + 1.01);
								
								i1 = i-4; i2 = i+4;
								j1 = j-4; j2 = j+4;
								k1 = k-4; k2 = k+4;
								if(i1 < 1) i1 = 1; if(i2 > IM_1) i2 = IM_1;
								if(j1 < 1) j1 = 1; if(j2 > JM_1) j2 = JM_1;
								
								weightsum = 0; weightcount = 0; bvalue = 0;
								for(i=i1;i<=i2;i++){
									for(j=j1;j<=j2;j++){
										for(k=k1;k<=k2;k++){
											xc = X0_1 + (i-1)*DX;
											yc = Y0_1 + (j-1)*DX;
											zc = Z0_1 + (k-1)*DX;
											
											
											if(k > KM_1) kmap = k - (KM_1-1);
											else if(k < 1) kmap = k + (KM_1-1);
											else kmap = k;
											
											o = DOMID_1[i][j][kmap];
											
											if(B_index_1[i][j][kmap] >= 0){
												distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
												
												if(distance <= StencilWidth*DX){
													weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
													#pragma omp atomic update
													BODYFORCE_1[0][o] += weightsum*Ffluid[0]/UNITMASS/LBM_BFORCE_SCALE;
													#pragma omp atomic update
													BODYFORCE_1[1][o] += weightsum*Ffluid[1]/UNITMASS/LBM_BFORCE_SCALE;
													#pragma omp atomic update
													BODYFORCE_1[2][o] += weightsum*Ffluid[2]/UNITMASS/LBM_BFORCE_SCALE;
													
												}
											}
										}
									}
								}
								
							}
							 
						}
						if(caps >= 189 && caps <= 329){
							for(node=1;node<=NodeM;node++){
								NODE_VEL[0][caps][node] += NODE_WALLREPUL[0][caps][node];
								NODE_VEL[1][caps][node] += NODE_WALLREPUL[1][caps][node];
								NODE_VEL[2][caps][node] += NODE_WALLREPUL[2][caps][node];
								
								Ffluid[0] = (gamma_p*NODE_WALLREPUL[0][caps][node]+F_Total[0][caps][node]);//body force
								Ffluid[1] = (gamma_p*NODE_WALLREPUL[1][caps][node]+F_Total[1][caps][node]);//body force
								Ffluid[2] = (gamma_p*NODE_WALLREPUL[2][caps][node]+F_Total[2][caps][node]);//body force
								///interpolate velocities on membrane from surrounding fluid
								xp = NODECOORD[0][caps][node]; yp = NODECOORD[1][caps][node]; zp = NODECOORD[2][caps][node];
								
								i = (int)((xp-X0_2)/DX + 1.01);
								j = (int)((yp-Y0_2)/DX + 1.01);
								k = (int)((zp-Z0_2)/DX + 1.01);
								
								
								i = (int)((xp-X0_2)/DX + 1.01);
								j = (int)((yp-Y0_2)/DX + 1.01);
								k = (int)((zp-Z0_2)/DX + 1.01);
								
								i1 = i-4; i2 = i+4;
								j1 = j-4; j2 = j+4;
								k1 = k-4; k2 = k+4;
								if(i1 < 1) i1 = 1; if(i2 > IM_2) i2 = IM_2;
								if(j1 < 1) j1 = 1; if(j2 > JM_2) j2 = JM_2;
								
								weightsum = 0; weightcount = 0; bvalue = 0; mindis = 1.5; triwall = 0;
								for(i=i1;i<=i2;i++){
									for(j=j1;j<=j2;j++){
										for(k=k1;k<=k2;k++){
											xc = X0_2 + (i-1)*DX;
											yc = Y0_2 + (j-1)*DX;
											zc = Z0_2 + (k-1)*DX;
											
											
											if(k < 1) kmap = k + (KM_2-1);
											else if(k > KM_2) kmap = k - (KM_2-1);
											else kmap = k;
											
											o = DOMID_2[i][j][kmap];
											
											if(B_index_2[i][j][kmap] >= 0){
												distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
												
												if(distance <= StencilWidth*DX){
													weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
													#pragma omp atomic update
													BODYFORCE_2[0][o] += weightsum*Ffluid[0]/UNITMASS/LBM_BFORCE_SCALE;
													#pragma omp atomic update
													BODYFORCE_2[1][o] += weightsum*Ffluid[1]/UNITMASS/LBM_BFORCE_SCALE;
													#pragma omp atomic update
													BODYFORCE_2[2][o] += weightsum*Ffluid[2]/UNITMASS/LBM_BFORCE_SCALE;
													
												}
											}
											
										}
									}
								}
								
								

							}
							
							for(node=1;node<=NodeM_N;node++){
								NODE_VEL_N[0][caps][node] += NODE_WALLREPUL_N[0][caps][node];
								NODE_VEL_N[1][caps][node] += NODE_WALLREPUL_N[1][caps][node];
								NODE_VEL_N[2][caps][node] += NODE_WALLREPUL_N[2][caps][node];
								
								Ffluid[0] = (gamma_p*NODE_WALLREPUL_N[0][caps][node]+F_Total_N[0][caps][node]);//body force
								Ffluid[1] = (gamma_p*NODE_WALLREPUL_N[1][caps][node]+F_Total_N[1][caps][node]);//body force
								Ffluid[2] = (gamma_p*NODE_WALLREPUL_N[2][caps][node]+F_Total_N[2][caps][node]);//body force
								///interpolate velocities on membrane from surrounding fluid
								xp = NODECOORD_N[0][caps][node]; yp = NODECOORD_N[1][caps][node]; zp = NODECOORD_N[2][caps][node];
								
								
								i = (int)((xp-X0_2)/DX + 1.01);
								j = (int)((yp-Y0_2)/DX + 1.01);
								k = (int)((zp-Z0_2)/DX + 1.01);
								
								i1 = i-4; i2 = i+4;
								j1 = j-4; j2 = j+4;
								k1 = k-4; k2 = k+4;
								if(i1 < 1) i1 = 1; if(i2 > IM_2) i2 = IM_2;
								if(j1 < 1) j1 = 1; if(j2 > JM_2) j2 = JM_2;
								
								weightsum = 0; weightcount = 0; bvalue = 0;
								for(i=i1;i<=i2;i++){
									for(j=j1;j<=j2;j++){
										for(k=k1;k<=k2;k++){
											xc = X0_2 + (i-1)*DX;
											yc = Y0_2 + (j-1)*DX;
											zc = Z0_2 + (k-1)*DX;
											
											
											if(k < 1) kmap = k + (KM_2-1);
											else if(k > KM_2) kmap = k - (KM_2-1);
											else kmap = k;
											
											o = DOMID_2[i][j][kmap];
											
											if(B_index_2[i][j][kmap] >= 0){
												distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
												
												if(distance <= StencilWidth*DX){
													weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
													#pragma omp atomic update
													BODYFORCE_2[0][o] += weightsum*Ffluid[0]/UNITMASS/LBM_BFORCE_SCALE;
													#pragma omp atomic update
													BODYFORCE_2[1][o] += weightsum*Ffluid[1]/UNITMASS/LBM_BFORCE_SCALE;
													#pragma omp atomic update
													BODYFORCE_2[2][o] += weightsum*Ffluid[2]/UNITMASS/LBM_BFORCE_SCALE;
													
												}
											}
										}
									}
								}
								
							}
						}
						if(caps >= 330 && caps <= 344){
							for(node=1;node<=NodeM;node++){
								NODE_VEL[0][caps][node] += NODE_WALLREPUL[0][caps][node];
								NODE_VEL[1][caps][node] += NODE_WALLREPUL[1][caps][node];
								NODE_VEL[2][caps][node] += NODE_WALLREPUL[2][caps][node];
								
								Ffluid[0] = (gamma_p*NODE_WALLREPUL[0][caps][node]+F_Total[0][caps][node]);//body force
								Ffluid[1] = (gamma_p*NODE_WALLREPUL[1][caps][node]+F_Total[1][caps][node]);//body force
								Ffluid[2] = (gamma_p*NODE_WALLREPUL[2][caps][node]+F_Total[2][caps][node]);//body force
								///interpolate velocities on membrane from surrounding fluid
								triwall = 0;
								xp = NODECOORD[0][caps][node]; yp = NODECOORD[1][caps][node]; zp = NODECOORD[2][caps][node];
								
								i = (int)((xp-X0_3)/DX + 1.01);
								j = (int)((yp-Y0_3)/DX + 1.01);
								k = (int)((zp-Z0_3)/DX + 1.01);
								
								
								i = (int)((xp-X0_3)/DX + 1.01);
								j = (int)((yp-Y0_3)/DX + 1.01);
								k = (int)((zp-Z0_3)/DX + 1.01);
								
								i1 = i-4; i2 = i+4;
								j1 = j-4; j2 = j+4;
								k1 = k-4; k2 = k+4;
								if(i1 < 1) i1 = 1; if(i2 > IM_3) i2 = IM_3;
								if(j1 < 1) j1 = 1; if(j2 > JM_3) j2 = JM_3;
								
								weightsum = 0; weightcount = 0; bvalue = 0; mindis = 5000.;
								for(i=i1;i<=i2;i++){
									for(j=j1;j<=j2;j++){
										for(k=k1;k<=k2;k++){
											xc = X0_3 + (i-1)*DX;
											yc = Y0_3 + (j-1)*DX;
											zc = Z0_3 + (k-1)*DX;
											
											
											if(k < 1) kmap = k + (KM_3-1);
											else if(k > KM_3) kmap = k - (KM_3-1);
											else kmap = k;
											
											o = DOMID_3[i][j][kmap];
											
											if(B_index_3[i][j][kmap] >= 0){
												distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
												
												if(distance <= StencilWidth*DX){
													weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
													#pragma omp atomic update
													BODYFORCE_3[0][o] += weightsum*Ffluid[0]/UNITMASS/LBM_BFORCE_SCALE;
													#pragma omp atomic update
													BODYFORCE_3[1][o] += weightsum*Ffluid[1]/UNITMASS/LBM_BFORCE_SCALE;
													#pragma omp atomic update
													BODYFORCE_3[2][o] += weightsum*Ffluid[2]/UNITMASS/LBM_BFORCE_SCALE;
													
												}
											}
											
										}
									}
								}
								
							}
														
							for(node=1;node<=NodeM_N;node++){
								NODE_VEL_N[0][caps][node] += NODE_WALLREPUL_N[0][caps][node];
								NODE_VEL_N[1][caps][node] += NODE_WALLREPUL_N[1][caps][node];
								NODE_VEL_N[2][caps][node] += NODE_WALLREPUL_N[2][caps][node];
								
								Ffluid[0] = (gamma_p*NODE_WALLREPUL_N[0][caps][node]+F_Total_N[0][caps][node]);//body force
								Ffluid[1] = (gamma_p*NODE_WALLREPUL_N[1][caps][node]+F_Total_N[1][caps][node]);//body force
								Ffluid[2] = (gamma_p*NODE_WALLREPUL_N[2][caps][node]+F_Total_N[2][caps][node]);//body force
								///interpolate velocities on membrane from surrounding fluid
								xp = NODECOORD_N[0][caps][node]; yp = NODECOORD_N[1][caps][node]; zp = NODECOORD_N[2][caps][node];
								
								
								i = (int)((xp-X0_3)/DX + 1.01);
								j = (int)((yp-Y0_3)/DX + 1.01);
								k = (int)((zp-Z0_3)/DX + 1.01);
								
								i1 = i-4; i2 = i+4;
								j1 = j-4; j2 = j+4;
								k1 = k-4; k2 = k+4;
								if(i1 < 1) i1 = 1; if(i2 > IM_3) i2 = IM_3;
								if(j1 < 1) j1 = 1; if(j2 > JM_3) j2 = JM_3;
								
								weightsum = 0; weightcount = 0; bvalue = 0;
								for(i=i1;i<=i2;i++){
									for(j=j1;j<=j2;j++){
										for(k=k1;k<=k2;k++){
											xc = X0_3 + (i-1)*DX;
											yc = Y0_3 + (j-1)*DX;
											zc = Z0_3 + (k-1)*DX;
											
											validinter[i-i1][j-j1][k-k1] = 0;
											
											if(k < 1) kmap = k + (KM_3-1);
											else if(k > KM_3) kmap = k - (KM_3-1);
											else kmap = k;
											
											o = DOMID_3[i][j][kmap];
											
											if(B_index_3[i][j][kmap] >= 0){
												distance = sqrtf((xc-xp)*(xc-xp)+(yc-yp)*(yc-yp)+(zc-zp)*(zc-zp));
												
												if(distance <= StencilWidth*DX){
													weightsum = (1.+cosf(PI*fabsf(xc-xp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(yc-yp)/(StencilWidth*DX)))*(1.+cosf(PI*fabsf(zc-zp)/(StencilWidth*DX)))/((StencilWidth*2)*(StencilWidth*2)*(StencilWidth*2));
													#pragma omp atomic update
													BODYFORCE_3[0][o] += weightsum*Ffluid[0]/UNITMASS/LBM_BFORCE_SCALE;
													#pragma omp atomic update
													BODYFORCE_3[1][o] += weightsum*Ffluid[1]/UNITMASS/LBM_BFORCE_SCALE;
													#pragma omp atomic update
													BODYFORCE_3[2][o] += weightsum*Ffluid[2]/UNITMASS/LBM_BFORCE_SCALE;
													
												}
											}
										}
									}
								}
								
							}
						}
					}
				}
			}
		}
		
		
		/// Updating node nodal positions of ROI cells//////////////////////////////////////////////////////////////////////////////////////////
		
		if(RBCSOLVERFLAG == 1){
			#pragma omp parallel private(i, caps, nzmax, nzmin, node, ParentID, check, tri, cellvel,counter,repulf)
			{
				#pragma omp for schedule(static,50)  nowait
				for(caps=345;caps<=CapsM;caps++){
					if(CellStatus[caps] < 0 && CellStatus[caps] > -4998){//these are clone cells still slave to parents
						nzmax = -5000.; nzmin = 5000.;
						ParentID = -CellStatus[caps];
						for(i=0;i<=2;i++){
							CAPS_MAX[i][caps] = -5000.;
							CAPS_MIN[i][caps] = 5000.;
						}
						
						for(node=1;node<=NodeM;node++){
							check = 0;
							if(ParentID <= 188 && NODECOORD[2][caps][node] > ZM_1+5.*DX) check = 1;
							if(ParentID >= 189 && ParentID <= 329 && NODECOORD[2][caps][node] < Z0_2-5.*DX) check = 1;
							if(ParentID >= 330 && ParentID <= 344 && NODECOORD[2][caps][node] < Z0_3-5.*DX) check = 1;
							for(i=0;i<=2;i++){
								if(check == 0) NODE_VEL[i][caps][node] = NODE_VEL[i][ParentID][node];
								NODECOORD[i][caps][node] += NODE_VEL[i][caps][node]*DT*((float)SOLVERSTAGGER);
								if(NODECOORD[i][caps][node] > CAPS_MAX[i][caps]) CAPS_MAX[i][caps] = NODECOORD[i][caps][node];
								if(NODECOORD[i][caps][node] < CAPS_MIN[i][caps]) CAPS_MIN[i][caps] = NODECOORD[i][caps][node];
							}
							if(NODECOORD[2][caps][node] > nzmax) nzmax = NODECOORD[2][caps][node];
							if(NODECOORD[2][caps][node] < nzmin) nzmin = NODECOORD[2][caps][node];
						}
						
						for(node=1;node<=NodeM_N;node++){
							check = 0;
							if(ParentID <= 188 && NODECOORD_N[2][caps][node] > ZM_1+5.*DX) check = 1;
							if(ParentID >= 189 && ParentID <= 329 && NODECOORD_N[2][caps][node] < Z0_2-5.*DX) check = 1;
							if(ParentID >= 330 && ParentID <= 344 && NODECOORD_N[2][caps][node] < Z0_3-5.*DX) check = 1;
							for(i=0;i<=2;i++){
								if(check == 0) NODE_VEL_N[i][caps][node] = NODE_VEL_N[i][ParentID][node];
								NODECOORD_N[i][caps][node] += NODE_VEL_N[i][caps][node]*DT*((float)SOLVERSTAGGER);
							}
						}
						
						if(nzmin > ZM_1 && nzmax < Z0_2){//they have fully crossed the cloning line
							CellStatus[caps] = 1;
						}
					}
					else if(CellStatus[caps] == 1){//check if they exit ROI domain, if so then delete
						for(i=0;i<=2;i++){
							CAPS_MAX[i][caps] = -5000.;
							CAPS_MIN[i][caps] = 5000.;
						}
						for(node=1;node<=NodeM_N;node++){
							for(i=0;i<=2;i++){
								NODECOORD_N[i][caps][node] += NODE_VEL_N[i][caps][node]*DT*((float)SOLVERSTAGGER);
							}
						}
						for(node=1;node<=NodeM;node++){
							repulf = sqrtf(NODE_VEL[0][caps][node]*NODE_VEL[0][caps][node] + NODE_VEL[1][caps][node]*NODE_VEL[1][caps][node] + NODE_VEL[2][caps][node]*NODE_VEL[2][caps][node]);
							
							if(repulf > 5.e-3 || isnan(repulf) == 1 || isinf(repulf) == 1){
								CellStatus[caps] = -5000;
								break;
							}
							
							for(i=0;i<=2;i++){
								
								
								NODECOORD[i][caps][node] += NODE_VEL[i][caps][node]*DT*((float)SOLVERSTAGGER);
								if(NODECOORD[i][caps][node] > CAPS_MAX[i][caps]) CAPS_MAX[i][caps] = NODECOORD[i][caps][node];
								if(NODECOORD[i][caps][node] < CAPS_MIN[i][caps]) CAPS_MIN[i][caps] = NODECOORD[i][caps][node];
							}
							if((NODECOORD[2][caps][node] < ZM_1) && REGRESSFLAG[caps] == 0){
								CellStatus[caps] = 3001;///flag cell for eventual removal from CV
							}
							if((NODECOORD[2][caps][node] > Z0_2) && REGRESSFLAG[caps] == 0){
								CellStatus[caps] = 3002;///flag cell for eventual removal from CA
							}
							
						}
						
						if(CellStatus[caps] > 2000){
							for(tri=1;tri<=TriM;tri++){
								AREA_UPDATE[caps][tri] = -1;
							}
						}
						
						
						if(CellStatus[caps] == -5000){//Let's reposition all nodes to null domain
							REGRESSFLAG[caps] = 0;
							NucMembUpdate[caps] = 0;
							for(node=1;node<=NodeM;node++){
								NODE_WALLNEIGH[caps][node] = 0;
								MembIn_NeighFlag[caps][node] = 0;
								MembOut_NeighFlag[caps][node] = 0;
								
								for(i=0;i<=2;i++){
									NODE_WALLREPUL[i][caps][node] = 0;
									NODECOORD[i][caps][node] = -5000.;
									F_Total[i][caps][node] = 0.;
									NODE_VEL[i][caps][node] = 0.;
								}
							}
							
							for(node=1;node<=NodeM_N;node++){
								NUC_NeighFlag[caps][node] = 0;
								for(i=0;i<=2;i++){
									NODECOORD_N[i][caps][node] = -5000.;
									F_Total_N[i][caps][node] = 0.;
									NODE_VEL_N[i][caps][node] = 0.;
								}
							}
						}
						
					}
					else if(CellStatus[caps] == 3001 || CellStatus[caps] == 3002){//check if they exit ROI domain, if so then delete
						for(i=0;i<=2;i++){
							CAPS_MAX[i][caps] = -5000.;
							CAPS_MIN[i][caps] = 5000.;
							cellvel[i] = 0;
						}
						counter = 0;
						if(CellStatus[caps] == 3001){///CV exit
							for(node=1;node<=NodeM;node++){
								
								if(NODECOORD[2][caps][node] > ZM_1){
									counter += 1;
									for(i=0;i<=2;i++) cellvel[i] += NODE_VEL[i][caps][node];
								}
							}
							
							if(counter > 0){
								for(i=0;i<=2;i++) cellvel[i] = cellvel[i]/counter;
								
								cellvel[2] += -10.e-6;
							}
							else{
								cellvel[0] = 0; cellvel[1] = 0; cellvel[2] = -10.e-6;
							}
							
							
							for(node=1;node<=NodeM_N;node++){
								if(NODECOORD_N[2][caps][node] > ZM_1){
									for(i=0;i<=2;i++){
										NODECOORD_N[i][caps][node] += NODE_VEL_N[i][caps][node]*DT*((float)SOLVERSTAGGER);
									}
								}
								else{
									for(i=0;i<=2;i++){
										NODECOORD_N[i][caps][node] += cellvel[i]*DT*((float)SOLVERSTAGGER);
									}
								}
							}
							for(node=1;node<=NodeM;node++){
								if(NODECOORD[2][caps][node] > ZM_1){
									for(i=0;i<=2;i++){
										NODECOORD[i][caps][node] += NODE_VEL[i][caps][node]*DT*((float)SOLVERSTAGGER);
									
										if(NODECOORD[i][caps][node] > CAPS_MAX[i][caps]) CAPS_MAX[i][caps] = NODECOORD[i][caps][node];
										if(NODECOORD[i][caps][node] < CAPS_MIN[i][caps]) CAPS_MIN[i][caps] = NODECOORD[i][caps][node];
									}
								}
								else{
									for(i=0;i<=2;i++){
										NODECOORD[i][caps][node] += cellvel[i]*DT*((float)SOLVERSTAGGER);
										
										if(NODECOORD[i][caps][node] > CAPS_MAX[i][caps]) CAPS_MAX[i][caps] = NODECOORD[i][caps][node];
										if(NODECOORD[i][caps][node] < CAPS_MIN[i][caps]) CAPS_MIN[i][caps] = NODECOORD[i][caps][node];
									}
								}
								if(CAPS_MAX[2][caps] <= ZM_1){
									CellStatus[caps] = -5000;///remove cell from simulation
									
									break;
								}
							}
							
						}
						else if(CellStatus[caps] == 3002){///CA exit
							for(node=1;node<=NodeM;node++){
								
								if(NODECOORD[2][caps][node] < Z0_2){
									counter += 1;
									for(i=0;i<=2;i++) cellvel[i] += NODE_VEL[i][caps][node];
								}
							}
							
							if(counter > 0){
								for(i=0;i<=2;i++) cellvel[i] = cellvel[i]/counter;
								
								cellvel[2] += 10.e-6;
							}
							else{
								cellvel[0] = 0; cellvel[1] = 0; cellvel[2] = 10.e-6;
							}
							
							
							for(node=1;node<=NodeM_N;node++){
								if(NODECOORD_N[2][caps][node] < Z0_2){
									for(i=0;i<=2;i++){
										NODECOORD_N[i][caps][node] += NODE_VEL_N[i][caps][node]*DT*((float)SOLVERSTAGGER);
									}
								}
								else{
									for(i=0;i<=2;i++){
										NODECOORD_N[i][caps][node] += cellvel[i]*DT*((float)SOLVERSTAGGER);
									}
								}
							}
							for(node=1;node<=NodeM;node++){
								if(NODECOORD[2][caps][node] < Z0_2){
									for(i=0;i<=2;i++){
										NODECOORD[i][caps][node] += NODE_VEL[i][caps][node]*DT*((float)SOLVERSTAGGER);
									
										if(NODECOORD[i][caps][node] > CAPS_MAX[i][caps]) CAPS_MAX[i][caps] = NODECOORD[i][caps][node];
										if(NODECOORD[i][caps][node] < CAPS_MIN[i][caps]) CAPS_MIN[i][caps] = NODECOORD[i][caps][node];
									}
								}
								else{
									for(i=0;i<=2;i++){
										NODECOORD[i][caps][node] += cellvel[i]*DT*((float)SOLVERSTAGGER);
										
										if(NODECOORD[i][caps][node] > CAPS_MAX[i][caps]) CAPS_MAX[i][caps] = NODECOORD[i][caps][node];
										if(NODECOORD[i][caps][node] < CAPS_MIN[i][caps]) CAPS_MIN[i][caps] = NODECOORD[i][caps][node];
									}
								}
								if(CAPS_MIN[2][caps] >= Z0_2){
									CellStatus[caps] = -5000;///remove cell from simulation
									
									break;
								}
							}
							
						}
						
						if(CellStatus[caps] > 2000){
							for(tri=1;tri<=TriM;tri++){
								AREA_UPDATE[caps][tri] = -1;
							}
						}
						
						if(CellStatus[caps] == -5000){//Let's reposition all nodes to null domain
							REGRESSFLAG[caps] = 0;
							NucMembUpdate[caps] = 0;
							for(node=1;node<=NodeM;node++){
								NODE_WALLNEIGH[caps][node] = 0;
								MembIn_NeighFlag[caps][node] = 0;
								MembOut_NeighFlag[caps][node] = 0;
								
								for(i=0;i<=2;i++){
									NODE_WALLREPUL[i][caps][node] = 0;
									NODECOORD[i][caps][node] = -5000.;
									F_Total[i][caps][node] = 0.;
									NODE_VEL[i][caps][node] = 0.;
								}
							}
							
							for(node=1;node<=NodeM_N;node++){
								NUC_NeighFlag[caps][node] = 0;
								for(i=0;i<=2;i++){
									NODECOORD_N[i][caps][node] = -5000.;
									F_Total_N[i][caps][node] = 0.;
									NODE_VEL_N[i][caps][node] = 0.;
								}
							}
						}
					}
				
				}
				
			}
			
		}
		
		
		/// Generation of ROI cells from parent cells at periodic BCs //////////////////////////////////////////
		/// Updating node positions of parent cells//////////////////////////////////////////////////////////////////////////////////////////
		
		if(RBCSOLVERFLAG == 1){
			#pragma omp parallel private(i, nzmax, nzmin, caps, node, m)
			{
				#pragma omp for schedule(static,50)  nowait	
				for(caps=1;caps<=344;caps++){
					CLONEFLAG[caps] = 0;
					if(CellStatus[caps] == 2){//once cells cross periodic boundaries the partial domains need to be looped to the to the other side of periodic boundary, the first crossing instance needs to be flagged for clone cell generation
						for(node=1;node<=NodeM;node++){
							for(i=0;i<=2;i++){
								NODECOORD[i][caps][node] += NODE_VEL[i][caps][node]*DT*((float)SOLVERSTAGGER);
							}
							if(caps <= 188){
								if(NODECOORD[2][caps][node] >= ZM_1){
									if(REGRESSFLAG[caps] == 0) CLONEFLAG[caps] = 1;
									NODECOORD[2][caps][node] = NODECOORD[2][caps][node] - (ZM_1-Z0_1);
									// in the next part of the code (after this for(caps=1;caps<=344;caps++) loop ends), 
									// the clone cell will be generated and at the end of the cloning the following flags
									// will be changed to prevent multiple cloning from the same parent cell crossing
									// the BC multiple times due to noise:
									// CLONEFLAG[caps] = 0 (1 previously);REGRESSFLAG[caps] = 1 (0 previously); CellStatus[caps] = 4 (2 previously);
									// only parents with CLONEFLAG[caps] == 1 & REGRESSFLAG[caps] == 0 & CellStatus[caps] == 2 are marked for cloning 
								}
								else if(NODECOORD[2][caps][node] < Z0_1){//this node moved backwards in one update frame
									NODECOORD[2][caps][node] = NODECOORD[2][caps][node] + (ZM_1-Z0_1);
									//REGRESSFLAG[caps] = 1;
								}
							}
							else if(caps >= 189 && caps <= 329){
								if(NODECOORD[2][caps][node] <= Z0_2){
									if(REGRESSFLAG[caps] == 0) CLONEFLAG[caps] = 2;
									NODECOORD[2][caps][node] = NODECOORD[2][caps][node] + (ZM_2-Z0_2);
								}
								else if(NODECOORD[2][caps][node] > ZM_2){//this node moved backwards in one update frame
									NODECOORD[2][caps][node] = NODECOORD[2][caps][node] - (ZM_2-Z0_2);
								}
							}
							else if(caps >= 330 && caps <= 344){
								if(NODECOORD[2][caps][node] <= Z0_3){
									if(REGRESSFLAG[caps] == 0) CLONEFLAG[caps] = 3;
									NODECOORD[2][caps][node] = NODECOORD[2][caps][node] + (ZM_3-Z0_3);
								}
								else if(NODECOORD[2][caps][node] > ZM_3){//this node moved backwards in one update frame
									NODECOORD[2][caps][node] = NODECOORD[2][caps][node] - (ZM_3-Z0_3);
								}
							}
						}
						for(node=1;node<=NodeM_N;node++){
							for(i=0;i<=2;i++){
								NODECOORD_N[i][caps][node] += NODE_VEL_N[i][caps][node]*DT*((float)SOLVERSTAGGER);
							}
							if(caps <= 188){
								if(NODECOORD_N[2][caps][node] >= ZM_1){
									if(REGRESSFLAG[caps] == 0) CLONEFLAG[caps] = 1;
									NODECOORD_N[2][caps][node] = NODECOORD_N[2][caps][node] - (ZM_1-Z0_1);
								}
								else if(NODECOORD_N[2][caps][node] < Z0_1){//this node moved backwards in one update frame
									NODECOORD_N[2][caps][node] = NODECOORD_N[2][caps][node] + (ZM_1-Z0_1);
								}
							}
							else if(caps >= 189 && caps <= 329){
								if(NODECOORD_N[2][caps][node] <= Z0_2){
									if(REGRESSFLAG[caps] == 0) CLONEFLAG[caps] = 2;
									NODECOORD_N[2][caps][node] = NODECOORD_N[2][caps][node] + (ZM_2-Z0_2);
								}
								else if(NODECOORD_N[2][caps][node] > ZM_2){//this node moved backwards in one update frame
									NODECOORD_N[2][caps][node] = NODECOORD_N[2][caps][node] - (ZM_2-Z0_2);
								}
							}
							else if(caps >= 330 && caps <= 344){
								if(NODECOORD_N[2][caps][node] <= Z0_3){
									if(REGRESSFLAG[caps] == 0) CLONEFLAG[caps] = 3;
									NODECOORD_N[2][caps][node] = NODECOORD_N[2][caps][node] + (ZM_3-Z0_3);
								}
								else if(NODECOORD_N[2][caps][node] > ZM_3){//this node moved backwards in one update frame
									NODECOORD_N[2][caps][node] = NODECOORD_N[2][caps][node] - (ZM_3-Z0_3);
								}
							}
						}
					}
					else if(CellStatus[caps] == 4){
						nzmax = -5000.; nzmin = 5000.;
						REGRESSFLAG[caps] = 1;
						for(node=1;node<=NodeM;node++){
							for(i=0;i<=2;i++){
								NODECOORD[i][caps][node] += NODE_VEL[i][caps][node]*DT*((float)SOLVERSTAGGER);
							}
							
							if(caps <= 188){
								if(NODECOORD[2][caps][node] >= ZM_1) NODECOORD[2][caps][node] = NODECOORD[2][caps][node] - (ZM_1-Z0_1);
								else if(NODECOORD[2][caps][node] < Z0_1){
									NODECOORD[2][caps][node] = NODECOORD[2][caps][node] + (ZM_1-Z0_1);
								}
							}
							else if(caps >= 189 && caps <= 329){
								if(NODECOORD[2][caps][node] <= Z0_2) NODECOORD[2][caps][node] = NODECOORD[2][caps][node] + (ZM_2-Z0_2);
								else if(NODECOORD[2][caps][node] > ZM_2){
									NODECOORD[2][caps][node] = NODECOORD[2][caps][node] - (ZM_2-Z0_2);
								}
								
							}
							else if(caps >= 330 && caps <= 344){
								if(NODECOORD[2][caps][node] <= Z0_3) NODECOORD[2][caps][node] = NODECOORD[2][caps][node] + (ZM_3-Z0_3);
								else if(NODECOORD[2][caps][node] > ZM_3){
									NODECOORD[2][caps][node] = NODECOORD[2][caps][node] - (ZM_3-Z0_3);
								}
							}
							
							if(NODECOORD[2][caps][node] > nzmax) nzmax = NODECOORD[2][caps][node];
							if(NODECOORD[2][caps][node] < nzmin) nzmin = NODECOORD[2][caps][node];
						}
						for(node=1;node<=NodeM_N;node++){
							for(i=0;i<=2;i++){
								NODECOORD_N[i][caps][node] += NODE_VEL_N[i][caps][node]*DT*((float)SOLVERSTAGGER);
							}
							if(caps <= 188){
								if(NODECOORD_N[2][caps][node] >= ZM_1) NODECOORD_N[2][caps][node] = NODECOORD_N[2][caps][node] - (ZM_1-Z0_1);
								else if(NODECOORD_N[2][caps][node] < Z0_1) NODECOORD_N[2][caps][node] = NODECOORD_N[2][caps][node] + (ZM_1-Z0_1);
							}
							else if(caps >= 189 && caps <= 329){
								if(NODECOORD_N[2][caps][node] <= Z0_2) NODECOORD_N[2][caps][node] = NODECOORD_N[2][caps][node] + (ZM_2-Z0_2);
								else if(NODECOORD_N[2][caps][node] > ZM_2) NODECOORD_N[2][caps][node] = NODECOORD_N[2][caps][node] - (ZM_2-Z0_2);
							}
							else if(caps >= 330 && caps <= 344){
								if(NODECOORD_N[2][caps][node] <= Z0_3) NODECOORD_N[2][caps][node] = NODECOORD_N[2][caps][node] + (ZM_3-Z0_3);
								else if(NODECOORD_N[2][caps][node] > ZM_3) NODECOORD_N[2][caps][node] = NODECOORD_N[2][caps][node] - (ZM_3-Z0_3);
							}
						}
						if(caps <= 188){
							if(nzmin > (Z0_1+5.) && nzmax < (Z0_1+0.75*(ZM_1-Z0_1))){
								CellStatus[caps] = 2;
							}
						}
						else if(caps >= 189 && caps <= 329){
							if(nzmax < (ZM_2-5.) && nzmin > (ZM_2-0.75*(ZM_2-Z0_2))){
								CellStatus[caps] = 2; 
							}
						}
						else if(caps >= 330 && caps <= 344){
							if(nzmax < (ZM_3-5.) && nzmin > (ZM_3-0.75*(ZM_3-Z0_3))){
								CellStatus[caps] = 2; 
							}
						}
					}
					
				}
			}
			
		}
		
		///cell-cloning, perform serially - cannot parallelize
		
		if(RBCSOLVERFLAG == 1){
			for(caps=1;caps<=344;caps++){
				if(CLONEFLAG[caps] > 0){
					for(m=345;m<=CapsM;m++){
						if(CellStatus[m] == -5000){
							CellStatus[m] = -caps;
							REGRESSFLAG[m] = -caps;
							
							for(node=1;node<=NodeM;node++){///start cloning the parent cell nodes to new ROI slave cell nodes
								for(i=0;i<=2;i++){
									NODECOORD[i][m][node] = NODECOORD[i][caps][node];
									NODE_VEL[i][m][node] = NODE_VEL[i][caps][node];
								}
								if(CLONEFLAG[caps] == 1 && NODECOORD[2][m][node] < 0.5*(Z0_1 + ZM_1)){
									NODECOORD[2][m][node] += (ZM_1-Z0_1);
								}
								else if(CLONEFLAG[caps] == 2 && NODECOORD[2][m][node] > 0.5*(Z0_2 + ZM_2)){
									NODECOORD[2][m][node] += -(ZM_2-Z0_2);
								}
								else if(CLONEFLAG[caps] == 3 && NODECOORD[2][m][node] > 0.5*(Z0_3 + ZM_3)){
									NODECOORD[2][m][node] += -(ZM_3-Z0_3);
								}
							}
							for(node=1;node<=NodeM_N;node++){///start cloning the parent cell nuclues nodes to new ROI slave cell nucleus nodes
								for(i=0;i<=2;i++){
									NODECOORD_N[i][m][node] = NODECOORD_N[i][caps][node];
									NODE_VEL_N[i][m][node] = NODE_VEL_N[i][caps][node];
								}
								if(CLONEFLAG[caps] == 1 && NODECOORD_N[2][m][node] < 0.5*(Z0_1 + ZM_1)){
									NODECOORD_N[2][m][node] += (ZM_1-Z0_1);
								}
								else if(CLONEFLAG[caps] == 2 && NODECOORD_N[2][m][node] > 0.5*(Z0_2 + ZM_2)){
									NODECOORD_N[2][m][node] += -(ZM_2-Z0_2);
								}
								else if(CLONEFLAG[caps] == 3 && NODECOORD_N[2][m][node] > 0.5*(Z0_3 + ZM_3)){
									NODECOORD_N[2][m][node] += -(ZM_3-Z0_3);
								}
							}
							
							for(tri=1;tri<=TriM;tri++){
								AREA_UPDATE[m][tri] = AREA_UPDATE[caps][tri];
							}
							
							
							CAPS_AREA[m] = CAPS_AREA[caps];
							CAPS_VOL[m] = CAPS_VOL[caps];
							break;
						}
					}
					CellStatus[caps] = 4;
					CLONEFLAG[caps] = 0;
					REGRESSFLAG[caps] = 1;
					
				}
			}	
			
		}
		
		
		
		
		if(iter%outputfreq == 0) OUTPUT();
		
		
		
	}while(iter < itermax);
	///end of LBM loop
	
	
	
	return 0;
}
