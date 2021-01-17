
/*
 This is the model to simulate the PC network with GCs and MLIs.
 17/01/2021
 */

#include <omp.h>
//#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
//#include <conio.h>  //  getch();   
#include <cstdlib>  // pause

using namespace std;

#include "param.h"          // file for all parameters 
#include "synapse.h"

//------------Number of ODE for each cell -------------------------------
const int N_EQ_G = 4;
const int N_EQ_M = 2;
const int N_EQ_MLI = 4;
const int N_EQ_U = 5;
const int N_PC_U = 4;

// Complete number of ODE
#define N_EQ   ( N_EQ_MLI*NMLI + N_EQ_G*NGC + N_PC_U*NPC + N_EQ_M*NMF ) 
//+++++++++++++++++++ MAIN PROGRAM +++++++++++++++++++++++++++++++++++++++++++

//----------external variables ---------------------------------------------
int no_GC[NGC][N_EQ_G], no_MLI[NMLI][N_EQ_MLI],
    no_eMF[NMF][N_EQ_M],no_PC[NPC][N_PC_U]; 

MatDoub g_eG(NGC,N_eG, 0.),    g_iG(NGC,N_iG, 0.),    g_MLIPC(NPC,N_MLIPC, 0.), 
        g_MMLI(NMLI,N_MMLI, 0.), g_GMLI(NMLI,N_GMLI, 0.),
		g_PC(NPC, N_GCPC, 0.);
MatInt pre_eG(NGC,N_eG, 0),  pre_iG(NGC,N_iG, 0), pre_MLIG(NGC,N_MLIG, 0), 
       pre_MMLI(NMLI,N_MMLI, 0),   pre_GMLI(NMLI,N_GMLI, 0),
	    pre_GCPC(NPC, N_GCPC, 0), pre_MLIPC(NPC, N_MLIPC, 0);

// learning parameters
VecDoub d_GG(NMLI, 0.), d_GCGC(NGC, 0.);

//----------external classes (beginning of initialization)------------------
 
GCBRUNEL GC[NGC];  
MLIBRUNEL MLI[NMLI];  
Purkinjeone PC[NPC];//one-compartment
int ii;
//mossy firbers
MFIAF  eMF[NMF];

// postsynapse currernt 
IMLI   *syn_MLIPC[NPC]; // MLI->GC
IMFGC *syn_iG[NGC];     // iMF -> GC
IMFGC *syn_eG[NGC];     // eMF->GC
IMF   *syn_GMLI[NMLI];  // GC ->MLI
IPC   *syn_GcPC[NPC];   //Gc->PC

//   -----external functions----------------------------------------------
void fun(double x, double *y_ini, double *f_ini);
void rkForth(int n, void fun(double, double*, double*),
		double h, double x, double* y, double* f, double* s, double* yk);
void rkSecond(int n, void fun(double, double*, double*),
		double h, double x, double* y, double* f, double* s, double* yk);

void euler(int n, void fun(double, double*, double*),
		double h, double x, double* y, double* f);

////////////////////////////////////////////////////////////////////////////////////
//++++++Main program+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
////////////////////////////////////////////////////////////////////////////////////
double EPSP=0.;
int main(int argc,char **argv) 
{
	//MPI
	int num_procs;
	double zz,ww;

	// running time counter
	clock_t start,finish;
	double totaltime;
	//start = clock();
	start = omp_get_wtime();

	FILE *f1, *f2, *f3, *f31, *f32, *f33, *f4, *f5, *f6, *f7, *f8, *f9, *f10, *f11, *f12, *f13, *f14, *f15, *f16;

	//---------allocate place for ALL variables and functions-----------
	double y_ini[N_EQ], f_ini[N_EQ];

	//---------allocate place for TWO temporal arrays using by RK.C solver 
	double y1[N_EQ], y2[N_EQ];

	//---------general parameters----------------------------------------------
	double t = 0, h = DT;
	int i, j, k, ii = 0;


   //
   //  How many processors are available?
   //
   num_procs = omp_get_num_procs ( );

   cout << "\n";
   cout << "  The number of processors available:\n";
   cout << "  OMP_GET_NUM_PROCS () = " << num_procs << "\n";
   
  int nthreads, thread_id;
#if TEST
  nthreads = 2;
#else
  nthreads = 4;
#endif
  omp_set_num_threads ( nthreads );

// Fork the threads and give them their own copy of the variables 
#pragma omp parallel private(nthreads, thread_id)
  {
    
    //* Get the thread number 
    thread_id = omp_get_thread_num();
    printf("Thread %d says: Hello World\n", thread_id);

    if (thread_id == 0)     // This is only done by the master 
	{
		nthreads = omp_get_num_threads();
		printf("Thread %d reports: the number of threads are %d\n", thread_id, nthreads);
	}
	
  }    // Now all threads join master thread and they are disbanded 


  
  /* Seed the random-number generator with current time so that
    * the numbers will be different every time we run.
    */
   //srand( (unsigned)time( NULL ) );
   srand( SEED );

   //----------arrays initialization----------------------------------------------
#pragma omp parallel private(i,j,k)
{

#pragma omp for nowait
   for(i=0; i<N_EQ; i++){
	   y_ini[i] = 0,  f_ini[i] = 0; 
	   y1[i] = 0,     y2[i] = 0;
   }


   //----------classes initialization (continue)----------------------------
   // depend on the order of defination synapse class
   #pragma omp for nowait
   for (i = 0; i < NPC; i++) {
	   syn_GcPC[i] = new IPC[N_GCPC];
	   syn_MLIPC[i] = new IMLI[N_MLIPC];
   }
#pragma omp for nowait
   for(i=0; i < NGC; i++) {
	   syn_iG[i] = new IMFGC[N_iG];
	   syn_eG[i] = new IMFGC[N_eG];

   }
//

#pragma omp for nowait
   for(i=0; i < NMLI; i++) {
	   syn_GMLI[i] = new IMF[N_GMLI];
   }

   //----------creating the integer arrays containing the addresses--------
   //----------of  ALL internal variables for ALL objects RE, TC ----------
   //----------and GB classes (e.g., no_re[i][j][k] is the address --------	
   //----------of the variable y[k] for the object re_cell[i][j]) ---------
   //----------NOTE: this is the relative addresses and you should ---------
   //----------add the REAL address of the first element of the -----------
   //----------original 1D array to use these arrays----------------------- 
#pragma omp for nowait
   for(i=0; i < NGC; i++)
	   for(k=0; k < N_EQ_G; k++)
		   no_GC[i][k] = k + i * N_EQ_G;



#pragma omp for nowait
   for(i=0; i < NMF; i++)
	   for(k=0; k < N_EQ_M; k++)
		   no_eMF[i][k] = NGC*N_EQ_G + k + i * N_EQ_M;

 #pragma omp for nowait
    for(i=0; i < NMLI; i++)
	    for(k=0; k < N_EQ_MLI; k++)
		    no_MLI[i][k] = NGC*N_EQ_G+ NMF*N_EQ_M + k + i * N_EQ_MLI;
#pragma omp for nowait
   for (i = 0; i < NPC; i++)
	   for (k = 0; k < N_PC_U; k++)
		   no_PC[i][k] = NGC*N_EQ_G+NMF*N_EQ_M +NMLI*N_EQ_MLI +k + i * N_PC_U;

#pragma omp barrier

   //---variable initialization (additional for standard constructor)----------
#pragma omp for nowait
   for(i=0; i<NGC; i++) {
	   GC[i].init(y_ini+no_GC[i][0]);
	   GC[i].reset();
   }

	 
#pragma omp for nowait
   for(i=0; i<NMF; i++) {
	   eMF[i].init(y_ini+no_eMF[i][0]);
	   eMF[i].reset();
   }
 #pragma omp for nowait
    for(i=0; i<NMLI; i++) {
	    MLI[i].init(y_ini+no_MLI[i][0]);
	    MLI[i].reset();
    }
#pragma omp for nowait
   for (i = 0; i < NPC; i++) {
	   PC[i].init(y_ini + no_PC[i][0]);
	   PC[i].reset();
   }

} //end omp

   //--------------open ALL files-------------------------------------
   if(    ( f1 = fopen("RasterGC.dat", "w") ) == NULL
	   || ( f2 = fopen("RasterMF.dat", "w") )      ==NULL
	   || ( f3 = fopen("connectionMG.dat", "w") )==NULL
	   || (f31 = fopen("connectionGMIL.dat", "w")) == NULL
	   || (f32 = fopen("connectionGPC.dat", "w")) == NULL
	   || (f33 = fopen("connectionMILPC.dat", "w")) == NULL
	   || ( f4 = fopen("wcon.dat", "w") )==NULL
	   || ( f5 = fopen("RasterPC.dat", "w") )==NULL
	   || ( f6 = fopen("PSC.dat", "w") )==NULL
	   || ( f8 = fopen("RasterMLI.dat", "w") )==NULL
	   || ( f9 = fopen("EPSC.dat", "w") )==NULL
	   || ( f7 = fopen("Voltage.dat", "w") )==NULL
	   || ( f10 = fopen("meanrate.dat", "w")) == NULL
	   || (f11 = fopen("rate.dat", "w")) == NULL
	   || (f12 = fopen("gampa.dat", "w")) == NULL
	   || (f13 = fopen("gampasl.dat", "w")) == NULL
	   || (f15 = fopen("meanMLIn.dat", "w")) == NULL
	   || (f16 = fopen("R.dat", "w")) == NULL
	   || (f14 = fopen("pcrate.dat", "w")) == NULL) {
		   printf("can't open files\n");
		   exit(0);
   }

   //fprintf(f6, "%4d\n",9999 );
   //fprintf(f7, "%4d\n",9999 );
   //fprintf(f8, "%4d %4d \n", 9999, 0 );

   //----------Connection matrix-------------------------------
   printf("\n Begin Connect Matrix");

   //-------------- Network topology ------------------------------------
   int pre, exists;
   for (i = 0; i < 1000; i++) {
	   pre = GetRandom(NMF);  
   }
   for(i=0; i<NGC; i++) {
	   // MF -> GC
	   for(j=0; j<N_eG; j++) {
		   do{
			   exists = 0;		// avoid multiple synapses
			   pre = GetRandom(NMF);
			   for (k=0;k<j;k++) if (pre_eG[i][k]==pre) exists = 1;	// synapse already exists  
		   }while (exists == 1);
		   pre_eG[i][j]= pre;
		   fprintf(f3, "%4d %4d \n", i, pre);		
	   }
	   
   }


   for(i=0; i<NMLI; i++) { 
	 

	   // GC -> MLI
	   for(j=0; j<N_GMLI; j++) {
		   do{
			   exists = 0;		// avoid multiple synapses
			   pre = GetRandom(NGC);
			   for (k=0;k<j;k++) if (pre_GMLI[i][k]==pre) exists = 1;	// synapse already exists  
		   }while (exists == 1);
		   pre_GMLI[i][j]= pre;
		   fprintf(f31, "%4d %4d \n", i, pre);
	   }
   }
		   
	for (i = 0; i < NPC; i++) {
	   for (j = 0; j < N_GCPC; j++) {
		   do {
			   exists = 0;		// avoid multiple synapses
			   pre = GetRandom(NGC);
			   for (k = 0; k < j; k++) if (pre_GCPC[i][k] == pre) exists = 1;	// synapse already exists  
		   } while (exists == 1);
		   pre_GCPC[i][j] = pre;
		   //fprintf(f5, "%4d %4d \n", i, pre);
		   fprintf(f32, "%4d %4d \n", i, pre);
	   }
	   for (j = 0; j < N_MLIPC; j++) {
		   do {
			   exists = 0;		// avoid multiple synapses
			   pre = GetRandom(NMLI);
			   for (k = 0; k < j; k++) if (pre_MLIPC[i][k] == pre) exists = 1;	// synapse already exists  
		   } while (exists == 1);
		   pre_MLIPC[i][j] = pre;
		   fprintf(f33, "%4d %4d \n", i, pre);
		   //fprintf(f5, "%4d %4d \n", i, pre);
	   }
   }

	
	  
	  printf("\n End Connect Matrix");


   //--------the initial conductances----------------------------------- 
   // normal distribution
   double weight;
   for (weight =3.5; weight <=3.5; weight = weight +0.5) {
	   //AMPArec = weight;
   cout << "  Weight  = " << weight << "\n";
   double WiG = wig;    
   Normaldev grandWiG( WiG, WiG * var, SEED);
  
   double WeG = weg;   
   Normaldev grandWeG( WeG, WeG * var, SEED);
     
   double WMLIG = weight;
   Normaldev grandWMLIG( WMLIG, WMLIG * var, SEED);////MLI-PC

   double WMMLI = wmMLI;   
   Normaldev grandWMMLI( WMMLI, WMMLI * var, SEED);

   double WGMLI = wgMLI;   
   Normaldev grandWGMLI( WGMLI, WGMLI * var, SEED);

   double WPC = wpc;
   Normaldev grandWPC(WPC, WPC * var, SEED);

   double Starttime = starttime;
   Normaldev grandstart(Starttime, 3*Starttime * var, SEED);



   for(i=0; i<NGC; i++) {
	   for(j=0; j<N_iG; j++) {
#if TEST
		   g_iG[i][j] = WiG;
#else
		   g_iG[i][j] = grandWiG.dev(); 
		   if ( g_iG[i][j] < 0. )    g_iG[i][j] = GetRandomDoub(2.*WiG);
#endif
	   }

	   for(j=0; j<N_eG; j++) {
#if TEST
		   g_eG[i][j] = WeG; 
#else
		   g_eG[i][j] = grandWeG.dev(); 
		   if ( g_eG[i][j] < 0. )    g_eG[i][j] = GetRandomDoub(2.*WeG);
#endif
	   }
   }

   

   Normaldev bigMMLI(15, 15.*var, SEED);
   for(i=0; i<NMLI; i++) {
//	  
//	   }

	   for (j = 0; j<N_GMLI; j++) {
#if TEST
		   g_MMLI[i][j] = WMMLI;
#else
		   g_GMLI[i][j] = grandWGMLI.dev();
		   if (g_GMLI[i][j] < 0.)    g_GMLI[i][j] = GetRandomDoub(2.*WGMLI);
		   //}
#endif
	   }


   }
   for (i = 0; i < NPC; i++) {
	   for (j = 0; j < N_GCPC; j++) {
		   g_PC[i][j] = grandWPC.dev();
		   if (g_PC[i][j] < 0) g_PC[i][j] = GetRandomDoub(2.*WPC);
	   }
	   for (j = 0; j < N_MLIPC; j++) {
		   g_MLIPC[i][j] = grandWMLIG.dev();
		   if (g_MLIPC[i][j] < 0.)    g_MLIPC[i][j] = GetRandomDoub(2.*WMLIG);
	   }
		  
   }

   
   fprintf(f4,"%5d \n", 9999 );
   printf("\n MLI network is done ! ");

   //**********************************************************

   //**********************************************************
   
   printf("\n GC network is done ! ");
   //**********************************************************

   
   printf("\n Stimulus is done ! ");
  
   //---changes of the parameters to get variability---------------------------
   int cy;
   //double cy;
   int len =2;
   double meanfre, UUU;
   for (UUU =0.4; UUU <= 0.4; UUU = UUU + 0.05) {
   AMPAfastU =  UUU;
   AMPAslowU =  UUU;
	   //AMPAfac = UUU;
   cout << "  UUU  = " << UUU << "\n";
   for (cy = 20; cy <= 20; cy = cy +20) {  //regular 100=100hz irregular 140=100hz  regular 10:20:190  irregular 10:40:300
	   t = 0.; ii = 0;
	   meanfre = 1000 /cy;
	  //GABAfastU = cy;
	  //GABAslowU = cy;
	   cout << "  cy  = " << cy << "\n";

//#if MFSTIM
	   for (i = 0; i < NMF; i++) {
		   eMF[i].TYPE = 1;  //randomly select TYPE I or II
		   eMF[i].k = 1.0;


		   eMF[i].slow_invl = meanfre;//meanfre;//Poisson mean ISI
		   eMF[i].start = 0;// grandstart.dev();//0; //ST*T+450;  //500;
		   eMF[i].end = TMAX; //ST*T+1500; //2550;
		   eMF[i].fast_invl = 11.1;// len;   //ST*T+1500; //2550;
		   eMF[i].burst_len = 25;// len * 5;    // spike number, 1 for single stimulus 5:2Brurst;8:3Burst
		   //1->28 5->45 10->90 15->130

		   eMF[i].noise = 0;
		   eMF[i].reset_input();
		   eMF[i].lastspk = -9999;
		   eMF[i].q = 0;
		   eMF[i].init(y_ini + no_eMF[i][0]);
		   //eMF[i].postB = 0;
	   }
	   


	   Normaldev grand1(THR_GC, fabs(THR_GC * 0.05), SEED);
	   //Normaldev NMDA1( NAR_MG, NAR_MG*var, SEED);
	   for (i = 0; i < NGC; i++) {
#if TEST
		   GC[i].Vthr = THR_GC;
#else
		   GC[i].Vthr = grand1.dev();
		   GC[i].lastspk = -9999;
		   GC[i].TimeCounter = -1.0;
		   GC[i].q = 0;
		   GC[i].init(y_ini + no_GC[i][0]);
		   GC[i].postB = 0;
		  // GC[i].U = 0.7;
		   
#endif

		   for (j = 0; j < N_eG; ++j) {
			   syn_eG[i][j].AMPA_NMDA_RATIO = NAR_MG;
			   syn_eG[i][j].AMPA_AMPAsl_RATIO = NAR_MG2;
			   syn_eG[i][j].I = 0;
		   }
		
	   }
	   printf("\n Setting up GC is done ! ");
	  Normaldev grand2(THR_PC, fabs(THR_PC * 0.02), SEED);
	   for (i = 0; i < NPC; i++) {

		   PC[i].Vthr = grand2.dev();
		   PC[i].lastspk = -9999;
		   PC[i].TimeCounter = -1.0;
		   PC[i].q = 0;
		   PC[i].postB = 0;
		   //PC[1].postB = 2;
		   //PC[2].postB = 3;
		   PC[i].init(y_ini + no_PC[i][0]);
		   for (j = 0; j < N_GCPC; ++j) {
			   syn_GcPC[i][j].AMPA_NMDA_RATIO = NAR_PC;
			   syn_GcPC[i][j].AMPA_AMPAsl_RATIO = NAR_PC;
			   syn_GcPC[i][j].I = 0;

			   // syn_eG[i][j].AMPA_AMPAsl_RATIO = NAR_MG2;
		   }
		   for (j = 0; j < N_MLIPC; ++j) {
#if TEST
			   syn_MLIPC[i][j].mix = 1; //mixed fast and slow component
#else
			   syn_MLIPC[i][j].mix = GetRandomDoub(1); //mixed fast and slow component
			   syn_MLIPC[i][j].I = 0;
#endif
		   }
	   }
	 /*  Normaldev gka(GKA, 0.3*GKA, SEED);
	   Normaldev gL(3.0, 0.05*3.0, SEED);
	   Normaldev tAHP(20.0, 0.1 * 20, SEED);*/
	   Normaldev delay1(1.0, 0.2, SEED);
	   Normaldev grand3(THR_MLI, fabs(THR_MLI * 0.05), SEED);
	   //Normaldev NMDA2( NAR_MMLI, NAR_MMLI*var, SEED);
	   //Normaldev NMDA3( NAR_GMLI, NAR_MG*var, SEED);
	   for (i = 0; i < NMLI; i++) {
#if TEST
		   MLI[i].delayA = 1.0;
		   MLI[i].delayB = 1.0;
		   MLI[i].para();
		   MLI[i].G_ka = GKA;
		   MLI[i].G_l = 3;
		   MLI[i].Vthr = THR_MLI;
		   MLI[i].tauAHP = 20;
#else
		   MLI[i].delayA = delay1.dev();
		   MLI[i].delayB = delay1.dev();
		   MLI[i].para();
		   MLI[i].G_ka = 0;//gka.dev(); 
		   MLI[i].Cm = 20;// GetRandomDoub(40.) + 20.; //20-60
		   MLI[i].G_l = MLI[i].Cm / MLI[i].taum;
		   MLI[i].tauAHP = 20; //tAHP.dev(); 
		   MLI[i].Vthr = grand3.dev();
		   MLI[i].E_l = -1.0*EL;

		   MLI[i].lastspk = -9999;
		   MLI[i].TimeCounter = -1.0;
		   MLI[i].q = 0;
		   MLI[i].postB = 0;
		   MLI[i].init(y_ini + no_MLI[i][0]);
		   MLI[i].reset();
#endif
		
		   for (j = 0; j < N_GMLI; ++j) {
			   syn_GMLI[i][j].AMPA_NMDA_RATIO = NAR_GMLI; //mixed fast and slow component
			   syn_GMLI[i][j].I = 0;
		   }
	   }
	   printf("\n Setting up MLI is done ! ");
	   //--------------end variability------------------------------------



	   //---------- end: changing variables----------------------
	   int it = 0, trial = 0;
	   double Vstim = 0.0, PSC = 0.0, GCFR = 0.0, GCFRmean = 0.0, GCfired = 0.0, Gcontrol = 0.0;
	   double PCFRmean = 0.0, PCfired = 0.0, MFFRmean = 0.0, MFfired = 0.0, UBCFRmean = 0.0, UBCfired = 0.0, MLIFRmean = 0.0, MLIfired = 0.0;
	   double PCrate = 0.0, MFrate = 0.0, UBCrate = 0.0, GCrate = 0.0, MLIrate = 0.0;
	   double pcrate1 = 0.0;// single cell rate
	   double gamamean = 0.0, gamaslmean = 0.0, meanMLIn=0.0, meanMLInR=0.0;
	   //----------------CALCULATION----------------------------------------
	   printf("\n CALCULATION IN PROGICSS!!!: TMAX= %lf (min)", TMAX / (1000.*60.));


	   for (i = 0; i < NMF; i++) {
		   zz = GetRandomDoub(1.0);
		   
		   if (zz <= 1) {
			   eMF[i].Inoise = -200;
			 
		   }
	   }

	   while (t <= TMAX) {

		   //add burst in Poisson Burst 1:3ms;3:8ms;7:20ms;10;28ms
		   //20-100;1-30. 28-130.10 spike
		   //20-20-1��20-30-3��20-70-7
		 if (t >=1200 && t <= (1200+25)) {
			   for (i = 0; i < NMF; i++) {
				   eMF[i].Startburst = 1;
				   if (t >= 1200 && t <= 1201) {
					   eMF[i].TimeCounter = -3;
				   }
			   }
		   }
		   else {
			   for (i = 0; i < NMF; i++) {
				   eMF[i].Startburst = 0;
				  
			   }
		   }
		   //len = len +4;

		   //printf("time= %lf (ms)", t);
#if HOLD
		   y_ini[no_GC[0][0]] = -70;
		   y_ini[no_UBC[1][0]] = -70;
		   y_ini[no_MLI[0][0]] = -70;
		   y_ini[no_PC[0][0]] = -70;
#endif
		   if (t >= ST*SHORTTIME) {
			   Vstim = A*1e-3 * sin(3 * 2*PI * t / T);
			   //Vstim = MFRATE*1e-3 * ( A* 1./T*1e3 * sin( 2*PI * t/T) );	
		   }

#if MFSTIM
		   //---------- modulate stimulating patterns ----------------------
		   for (i = 0; i < NMF; i++) {
			   eMF[i].mu = eMF[i].k * fabs(Vstim);
		   }
		   // bidirectional varying stimuli
		   if (Vstim >= THR && t >= ST*SHORTTIME) {

			   for (i = 0; i < NMF; i++) {

				   zz = GetRandomDoub(1.0);
				   //printf("%lf\n", zz);
				   ww = DT * (MFRATE*1e-3 + eMF[i].mu);
				   if ((eMF[i].TYPE == 1) && (zz - ww < 0)) {
					   eMF[i].Inoise = -200;
					   //fprintf(f2, "%5f %5d \n", t,i);
				   }
				   else if ((eMF[i].TYPE == 2) && (GetRandomDoub(1.0) < DT * (MFRATE*1e-3 - eMF[i].mu))) {
					   eMF[i].Inoise = -200;
				   }
				   else {
					   eMF[i].Inoise = GetRandomDoub(1.0) - 0.5;
				   }
			   }

		   }
		   else if (Vstim < THR && t >= ST*SHORTTIME) {
			   for (i = 0; i < NMF; i++) {
				   if ((eMF[i].TYPE == 2) && (GetRandomDoub(1.0) < DT * (MFRATE*1e-3 + eMF[i].mu))) {
					   eMF[i].Inoise = -200;
				   }
				   else if ((eMF[i].TYPE == 1) && (GetRandomDoub(1.0) < DT * (MFRATE*1e-3 - eMF[i].mu))) {
					   eMF[i].Inoise = -200;
					   //fprintf(f2, "%5f %5d \n", t, i);
				   }
				   else {
					   eMF[i].Inoise = GetRandomDoub(1.0) - 0.5;
				   }
			   }
		   }
#endif
#if BURST
		   //---------- modulate stimulating patterns ----------------------
		   for (i = 0; i < NMF; i++) {
			   eMF[i].mu = eMF[i].k * fabs(Vstim);
		   }
		   // bidirectional varying stimuli
		   if (Vstim >= THR && t >= ST*SHORTTIME) {

			   for (i = 0; i < NMF; i++) {

				   zz = GetRandomDoub(1.0);
				   //printf("%lf\n", zz);
				   ww = DT * (MFRATE*1e-3 + eMF[i].mu);
				   if ((eMF[i].TYPE == 1) && (zz - ww < 0)) {
					   eMF[i].Inoise = -200;
					   //fprintf(f2, "%5f %5d \n", t,i);
				   }
				   else if ((eMF[i].TYPE == 2) && (GetRandomDoub(1.0) < DT * (MFRATE*1e-3 - eMF[i].mu))) {
					   eMF[i].Inoise = -200;
				   }
				   else {
					   eMF[i].Inoise = GetRandomDoub(1.0) - 0.5;
				   }
			   }

		   }
		   else if (Vstim < THR && t >= ST*SHORTTIME) {
			   for (i = 0; i < NMF; i++) {
				   if ((eMF[i].TYPE == 2) && (GetRandomDoub(1.0) < DT * (MFRATE*1e-3 + eMF[i].mu))) {
					   eMF[i].Inoise = -200;
				   }
				   else if ((eMF[i].TYPE == 1) && (GetRandomDoub(1.0) < DT * (MFRATE*1e-3 - eMF[i].mu))) {
					   eMF[i].Inoise = -200;
					   //fprintf(f2, "%5f %5d \n", t, i);
				   }
				   else {
					   eMF[i].Inoise = GetRandomDoub(1.0) - 0.5;
				   }
			   }
		   }
#endif

		   /*y_ini[no_PC[0][0]] = -70;
		   EPSP = 0.;*/

		   //rkForth(N_EQ, fun, h, t, y_ini, f_ini, y1, y2);
		   //rkSecond(N_EQ, fun, h, t, y_ini, f_ini, y1, y2);
		   euler(N_EQ, fun, h, t, y_ini, f_ini);
		   t += h;
		   ii++;

#if SAVE_VOLTAGE 
		   //------save voltage data -----------------------------------------------------------------
		   if (ii % int(1) == 0) {
			
			   if (t > 1000) {
				   //fprintf(f7, "%lf %lf %lf %lf\n", t, y_ini[no_PC[0][0]], y_ini[no_PC[1][0]], EPSP);
			   }
			   fprintf(f7, "%lf %lf\n", t, y_ini[no_PC[0][0]]);

			 //fprintf(f7, "%lf %lf  %lf\n", t, GC[0].AMPAfR, GC[0].AMPAfu);
			// fprintf(f7, "%lf %lf  %lf\n", t, GC[0].AMPAsR, GC[0].AMPAsu);
			   gamamean = 0.0; gamaslmean = 0.0, meanMLIn = 0.0,meanMLInR = 0.0;
			   if (t >= 000) {
				/*   fprintf(f12, "%lf ", t);
				   fprintf(f13, "%lf ", t);
				   fprintf(f15, "%lf ", t);*/
				   //N_GCPC
				   for (j = 0; j < NPC; j++) {
					   for (i = 0; i < N_GCPC; i++) {
						   //fprintf(f12, "%lf ", GC[pre_GCPC[0][i]].ampaGP);
						   //fprintf(f13, "%lf ", GC[pre_GCPC[0][i]].ampaMGsl*0.3);
						   gamamean += GC[pre_GCPC[j][i]].ampaGP;
						   gamaslmean += GC[pre_GCPC[j][i]].ampaslGP*1.4;
						   meanMLInR += GC[pre_GCPC[j][i]].AMPAfR;
						   meanMLInR += GC[pre_GCPC[j][i]].AMPAsR;
					   }
					   meanMLIn = gamaslmean + gamamean;
					   meanMLIn = meanMLIn / N_GCPC;
					   meanMLInR = meanMLInR / N_GCPC;
					   gamamean = gamamean / N_GCPC;
					   gamaslmean = gamaslmean / N_GCPC;

					   fprintf(f12, "%lf ", gamamean);
					   fprintf(f13, "%lf ", gamaslmean);
					   fprintf(f15, "%lf ", meanMLIn);
					   fprintf(f16, "%lf ", meanMLInR/2);

								
				   }
				   fprintf(f12, "\n ");
				   fprintf(f13, "\n ");
				   fprintf(f15, "\n ");
				   fprintf(f16, "\n ");
			   }
			
			  //fprintf(f7, "%lf %lf %lf %lf %lf\n", t, eMF[1].ampaMG, eMF[1].nmdaMG, MLI[1].MLIBR, MLI[1].MLIBu);
			  // fprintf(f3, "%lf %lf %lf %lf %lf\n", t, eMF[0].ampaMG, eMF[0].nmdaMG, MLI[0].gabaa, MLI[0].gabab);
			 // fprintf(f3, "%lf %lf\n", t, y_ini[no_PC[2][0]]);
		   }
#endif

#if SAVE_EPSC
		   if (ii % int(1 / h) == 0) {

			  
			   PSC = 0.0;
			   //--------GC-PC-----------
			   if (t >= 0000) {
				   //for (i = 0; i < NPC; ++i) {
					   for (j = 0; j < N_GCPC; ++j) {
						   PSC += syn_GcPC[0][j].I;
					   }
					   fprintf(f9, "%lf %lf\n ",t, PSC);
					   
				  // }
					   PSC = 0.0;
					   for (j = 0; j < N_MLIPC; ++j) {
						   PSC += syn_MLIPC[0][j].I;
					   }
					   fprintf(f6, "%lf %lf\n ",t, PSC);

			

			   }
			   			
		   }
#endif
#if SAVE_RATE
		   if (ii % int(1/ h) == 0) {
			   PCrate = 0.0, MFrate = 0.0, UBCrate = 0.0, GCrate = 0.0,MLIrate=0.;
			   //GC rate in bin=1ms
			   for (i = 0; i < NGC; i++) {
				   for (j = 0; j < GC[i].Ca; j++) {
					   if (GC[i].spiketimes[j]>(t -1)) {
						   GCrate += 1;
					   }
				   }
			   }
			   
			   ////////////////////////////////////
			
			   for (i = 0; i < NPC; i++) {
				   for (j = 0; j < PC[i].Ca; j++) {
					   if (PC[i].spiketimes[j]>(t - 1)) {
						   PCrate += 1;
					   }
				   }
			   }
			
			   //MF rate in bin=1ms
			   for (i = 0; i < NMF; i++) {
				   for (j = 0; j < eMF[i].Ca; j++) {
					   if (eMF[i].spiketimes[j]>(t - 1)) {
						   MFrate += 1;
					   }
				   }
			   }
			   //MLI rate in bin=1ms
			   for (i = 0; i < NMLI; i++) {
				   for (j = 0; j < MLI[i].Ca; j++) {
					   if (MLI[i].spiketimes[j]>(t - 1)) {
						   MLIrate += 1;
					   }
				   }
			   }
			  fprintf(f11, "%lf %lf %lf %lf %lf\n", t, MFrate, MLIrate, GCrate, PCrate);

		   }
#endif

		   //--------  at the end of one trial -----------------------------------
		   if (ii % int(SHORTTIME / h) == 0) {
			   it++;
			   trial = (int)ceil((double)it / (double)NSTIM);

			   printf("\n Trial = %d ", (int)ceil((double)it / (double)NSTIM));
			   //if ( it ==40 ) system("PAUSE"); 		  

			   //fprintf(f7, "%4d\n",9999 );
			   //fprintf(f8, "%4d %4d \n",9999, 0 );

			   //------save spike data -----------------------------------------------------------------
			   fprintf(f1, " %4d %4d %4d\n", 9999, it, 1);
			   GCFRmean = 0.0;
			   GCfired = 0.0;
			   PCFRmean = 0.0;
			   PCfired = 0.0;
			   MFFRmean = 0.0;
			   MFfired = 0.0;
			   UBCFRmean = 0.0;
			   UBCfired = 0.0;
			   MLIFRmean = 0.0, MLIfired = 0.0;
			   pcrate1 = 0.0;
		


			   for (i = 0; i < NGC; i++) {
				   for (j = 0; j < GC[i].Ca; j++) {
					   //fprintf(f1,"%4d %4d %lf \n", trial, i, GC[i].spiketimes[j]- (it-1)*SHORTTIME );
					   fprintf(f1, "%4d %4d %lf \n", trial, i, GC[i].spiketimes[j]);

					   GCfired += 1.0;
					   GCFRmean += 1;

				   }
				
			   }
			   if (GCfired > 0.0) {
				   GCFRmean /= NGC;
			   }
			   else {
				   GCFRmean = 0.0;
			   }
			   Gcontrol += ALPHA * (GCFRmean - GCRATE);
			   for (i = 0; i < NGC; i++) {
				   GC[i].Ginh = Gcontrol;
			   }

			   
			   //save single PC rate
			   if (t > 1999) {
				   for (i = 0; i < NPC; i++) {
					   for (j = 0; j < PC[i].Ca; j++) {
						   pcrate1 += 1;
					   }
					   fprintf(f14, " %1f", pcrate1);//
					   pcrate1 = 0.0;
				   }
			   }
			   fprintf(f14, "\n");
			   //PC mean rate
			   for (i = 0; i < NPC; i++) {
				   for (j = 0; j < PC[i].Ca; j++) {
					   fprintf(f5, "%4d %4d %lf \n", trial, i, PC[i].spiketimes[j]);//
					   //fprintf(f5, "%4d %4d %lf \n", trial, i, PC[i].spiketimes[j] - (it - 1)*SHORTTIME);//
					   //fprintf(f1, "%4d %4d %lf \n", trial, i, PC[i].spiketimes[j] - (it - 1)*SHORTTIME);
					   PCFRmean += 1;
					   PCfired += 1.0;

				   }
				  
			   }
			   if (PCfired > 0.0) {
				   PCFRmean /= NPC;
			   }
			   else {
				   PCFRmean = 0.0;
			   }
			   //MF mean rate
			   for (i = 0; i < NMF; i++) {
				   for (j = 0; j < eMF[i].Ca; j++) {
					   //fprintf(f2, "%4d %4d %lf \n", trial, i, eMF[i].spiketimes[j] - (it - 1)*SHORTTIME);
					   fprintf(f2, "%4d %4d %lf \n", trial, i, eMF[i].spiketimes[j]);

					   MFFRmean += 1;
					   MFfired += 1.0;

				   }
			   }
			   if (MFfired > 0.0) {
				   MFFRmean /= NMF;
			   }
			   else {
				   MFFRmean = 0.0;
			   }
				//MLI mean rate
			   for (i = 0; i < NMLI; i++) {
				   for (j = 0; j < MLI[i].Ca; j++) {
					   //fprintf(f2, "%4d %4d %lf \n", trial, i, eMF[i].spiketimes[j] - (it - 1)*SHORTTIME);	
					   fprintf(f8, "%4d %4d %lf \n", trial, i, MLI[i].spiketimes[j]);

					   MLIFRmean += 1;
					   MLIfired += 1.0;

				   }
			   }
			   if (MLIfired > 0.0) {
				   MLIFRmean /= NMLI;
			   }
			   else {
				   MLIFRmean = 0.0;
			   }
			   printf("Trial=%4d, G=%lf, MFrate=%lf \n", trial, Gcontrol, MFFRmean);
			   printf("Trial=%4d, G=%lf, GCrate=%lf \n", trial, Gcontrol, GCFRmean);
			  // printf("Trial=%4d, G=%lf, UBCrate=%lf \n", trial, Gcontrol, UBCFRmean);
			   printf("Trial=%4d, G=%lf, MLIrate=%lf \n", trial, Gcontrol, MLIFRmean);
			   printf("Trial=%4d, G=%lf, PCrate=%lf \n", trial, Gcontrol, PCFRmean);
			   if (t > 1999) {
				   fprintf(f10, "%lf %lf %lf %lf\n", MFFRmean, MLIFRmean, GCFRmean, PCFRmean);
			   }


			   

			   // reset variables
			   for (i = 0; i < NMF; i++)    eMF[i].reset();
			   for (i = 0; i < NGC; i++)    GC[i].reset();
			   //for (i = 0; i < NUBC; i++)   UBC[i].reset();
			   for (i = 0; i < NMLI; i++)   MLI[i].reset();
			   for (i = 0; i < NPC; i++)PC[i].reset();
			   //fprintf(f12, "%lf %lf\n", 999.0, 0.0);
			   //fprintf(f13, "%lf %lf\n", 999.0, 0.0);

		   }// end of one trial 
		  		 
	   }
	  /* for (i = 0; i < NPC+1; i++) {
		  
		   fprintf(f11, "%4d ", -100);
	
	   }

	   fprintf(f11, "%lf\n");*/
	   //fprintf(f11, "%ld %ld %ld %ld  %ld\n", -1, -1, -1, -1,-1);
	   //fprintf(f11, "\n");
	   fprintf(f11, "%ld %ld %ld %ld  %ld\n", -1, -1, -1, -1,-1);
	   fprintf(f12, "%lf %lf\n", -999.0, -999.0);
	   fprintf(f13, "%lf %lf\n", -999.0, -999.0);
	  // fprintf(f15, "%lf %lf\n", -999.0, -999.0);
	   fprintf(f5, "%4d %4d  %4d \n", -99, -99, -99);//

   }
   //U 
   fprintf(f12, "%lf %lf\n", -9999.0, -9999.0);
   fprintf(f13, "%lf %lf\n", -9999.0, -9999.0);
  // fprintf(f15, "%lf %lf\n", -9999.0, -9999.0);
  // fprintf(f10, "%ld %ld %ld %ld\n", -1, -1, -1, -1);
   //fprintf(f7, "%ld %ld %ld %ld\n", -1, -1, -1, -1);
   fprintf(f5, "%4d %4d %4d \n", -999, -999, -999);//
   fprintf(f11, "%ld %ld %ld %ld  %ld\n", -11, -11, -11, -11, -11);
   }
   //Inh
   fprintf(f12, "%lf %lf\n", -99999.0, -99999.0);
   fprintf(f13, "%lf %lf\n", -99999.0, -99999.0);
   //fprintf(f15, "%lf %lf\n", -99999.0, -99999.0);
   //fprintf(f10, "%ld %ld %ld %ld\n", -2, -2, -2, -2);
   fprintf(f5, "%4d %4d  %4d \n", -9999, -9999, -9999);//
   fprintf(f11, "%ld %ld %ld %ld  %ld\n", -111, -111, -111, -111, -111);
   }//--------------------END CALCULATION-------------------------------
      
  //-----------------close ALL files-----------------------------------
  fprintf(f1, " %4d %4d %4d\n",9999, NTRIAL, NSTIM );
  fclose(f1);
  fclose(f2);
  fclose(f3);
  fclose(f31);
  fclose(f32);
  fclose(f33);
  fclose(f4);
  fclose(f5);
  fclose(f6);
  fclose(f7);
  fclose(f8);
  fclose(f9);
  fclose(f10);
  fclose(f11);
  fclose(f12);
  fclose(f13);
  fclose(f14);
  fclose(f15);
  fclose(f16);
  // free memory


   //free memory
  for(i=0; i < NGC; i++) {
	  delete [] syn_iG[i];
	  delete [] syn_eG[i];

  }


  for(i=0; i < NMLI; i++) {

	  delete [] syn_GMLI[i];

  }
  for (i = 0; i < NPC; i++) {
	  delete[] syn_GcPC[i];
	  delete[] syn_MLIPC[i];


  }


 

   finish = omp_get_wtime();
   //totaltime = (double)(finish-start)/CLOCKS_PER_SEC;
   totaltime = (double)(finish-start);
   cout<<"\n running time is "<<totaltime<<" sec"<<endl;

   system("PAUSE"); 
   return 0;
}



//----------external functions----------------------------------------------
void fun(double x, double *y_ini, double *f_ini){
	int i, j; 
	


		//========here the MAIN loop to calculate intrinsic conductances===========
	//--------(f_ini IS changed, y_ini IS NOT changed)-------------------------
#pragma omp parallel private(i,j)
{
#pragma omp for nowait
	for(i=0; i<NGC; ++i)    GC[i].calc( x, y_ini+no_GC[i][0],  f_ini+no_GC[i][0]); 

#pragma omp for nowait
	
	for (i = 0; i < NMF; ++i) { 
		
		eMF[i].calc(x, y_ini + no_eMF[i][0], f_ini + no_eMF[i][0]);
	}
#pragma omp for nowait
	for(i=0; i<NMLI; ++i) 	MLI[i].calc(x, y_ini+no_MLI[i][0], f_ini+no_MLI[i][0]); 
#pragma omp for nowait
	for (i = 0; i < NPC; ++i)  PC[i].calc(x, y_ini + no_PC[i][0], f_ini + no_PC[i][0]);

#pragma omp barrier

		//========here the MAIN loop to calculate synaptic conductances=============
		//--------(f_ini IS changed, y_ini IS NOT changed) -------------------------
	    // i-post; j-pre
		// update GC post synapse
#pragma omp for nowait
		for(i = 0; i < NGC; ++i)	{						

		

			//-----------AMPA and NMDA from eMF to GC cells--------------------------------------
			for (j = 0; j < N_eG; ++j) {
				syn_eG[i][j].calc(g_eG[i][j], eMF[pre_eG[i][j]].ampaMG, eMF[pre_eG[i][j]].nmdaMG, eMF[pre_eG[i][j]].ampaMGsl,
					y_ini[no_GC[i][0]], GC[i].postB);
				f_ini[no_GC[i][0]] -= syn_eG[i][j].I / GC[i].Cm;
			}
		
		}


#pragma omp for nowait
		// update MLI post synapse
		for(i = 0; i < NMLI; ++i)	{						

			//-----------AMPA and NMDA from GC to MLI cells--------------------------------------
			for(j = 0; j < N_GMLI; ++j){
				syn_GMLI[i][j].calc( g_GMLI[i][j], GC[ pre_GMLI[i][j] ].ampaGMLI, GC[ pre_GMLI[i][j] ].nmdaGMLI, y_ini[no_MLI[i][0]], MLI[i].postB );
				f_ini[no_MLI[i][0]] -= syn_GMLI[i][j].I / MLI[i].Cm;	
				
			}
		
		}
#pragma omp for nowait
		// update pc post synapse
		for (i = 0; i < NPC; ++i) {
		
			//-----------ampa and nmda from GC to PC cells--------------------------------------
			for (j = 0; j < N_GCPC; ++j) {
				syn_GcPC[i][j].calc(g_PC[i][j], GC[pre_GCPC[i][j]].ampaGP, GC[pre_GCPC[i][j]].ampaslGP, y_ini[no_PC[i][0]], PC[i].postB);
				//f_ini[no_PC[i][0]] -= syn_GcPC[i][j].I / PC[i].Cms;
				f_ini[no_PC[i][0]] -= syn_GcPC[i][j].I / PC[i].Cm;//one-compartment
				EPSP += -syn_GcPC[i][j].I / PC[i].Cm;
			
			}
			//-----------GABAa and GABAb from Intermediate neuron to PC cells--------------------------------------	 
			for (j = 0; j < N_MLIPC; ++j) {
				//cout<<" pre_MLI is "<<pre_MLIG[i][j] <<endl;
				syn_MLIPC[i][j].calc(g_MLIPC[i][j], MLI[pre_MLIPC[i][j]].gabaa, MLI[pre_MLIPC[i][j]].gabab, y_ini[no_PC[i][0]]);//
				f_ini[no_PC[i][0]] -= syn_MLIPC[i][j].I / PC[i].Cm;
			}

		}
		//--------------

}
		//=============END of MAIN loop==============================================

}


// 1st order Rounge-Kutta ((forward) Euler method) solver for ODE -------------------------------------------------
//   rk(N_EQ, fun, h, t, y_ini, f_ini, y1, y2);
void euler(int n, void fun(double, double*, double*), 
        double h, double x, double* y, double* f)
{
	int i;


#pragma omp parallel for
	for(i = 0; i < n; ++i) 	{
		y[i] += h * f[i]; 
	}

	//k1
	fun(x, y, f);
}



// 2nd order Rounge-Kutta solver for ODE -------------------------------------------------
//   rk(N_EQ, fun, h, t, y_ini, f_ini, y1, y2);
void rkSecond(int n, void fun(double, double*, double*), 
        double h, double x, double* y, double* f, double* s, double* yk)
{
	int i;
	double xk;

	//k1
	fun(x, y, f);    

	for(i = 0; i < n; i++)
	{
		yk[i] = y[i] + (h/2.)*f[i]; 
	}
	xk = x + h/2.;    
	//k2
	fun(xk, yk, f);

	for(i = 0; i < n; ++i) y[i] += h * f[i];
}

