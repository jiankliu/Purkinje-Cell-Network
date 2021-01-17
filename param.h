
#include <math.h>

#ifndef PARAM_H
#define PARAM_H

/************** Defines some global paramaters -- *************/
	
//  units:
//  (nS) = 
//  (nA) = (nanoamp)
//  (mV) = (millivolt)
//  (umho) = (micromho)
//  (mM) = (milli/liter)
//  t = msec

#define WINVC 1
#define UNIX  0

//==============================================================================
//     [simulation parameters and settings]
//==============================================================================

//------------------------------------------------------------------------------
// general settings
//------------------------------------------------------------------------------
#define TEST   0
#define MFSTIM   0
#define POISSON 1
#define BURST 0
#define SAVE_rate 0
#define SAVE_VOLTAGE  1
#define SAVE_EPSC  1
#define SAVE_RATE 1
#define CURRENT_CLAMP 0
#define V_CLAMP 0
#define NOUBC  0
#define HOLD   0
#define GABAA  0       // 1: only GABAA from MLI->GC

const int NSTIM = 1;
const int NTRIAL =2;
const int ST = 0;        // time without stimuli
const int NT = 1;
const double T = 1000;  // [5000 2500 1667 1250 1000 833 714 625 555 500] (ms) period of stimulus 
const double SHORTTIME = T;
const double TMAX =  NTRIAL * SHORTTIME + ST*T;       //simulation time (m sec)

const int SEED = 2947;  // [43 433 2947] randome seed 
 int SEED1 = 1112; 
const int RANDIN = 1;
const double GK = 2;

//MLI
const double GKA = 0; //1.2;
const double EL = 70;  // 
const double A = 30;// 10;//10;      // ( o/s) amplitude of stimulus velocity

//
const double TONICI = 900;  // 
const double TSTIM  = 100;  // ty
const double GCRATE  = 5;     // average GC firing rate
const double ALPHA = 0.0;                 // time constant of controling

//------------------------------------------------------------------------------
//  network simulation parameter
//------------------------------------------------------------------------------
const int NMF = 500;// 500;//500; //56;
const int NGC = 1000;//1500;//4900; //900;
const int NMLI = 500;// 500; //144;
const int NPC =50 ;// 50;//if add one compartment is 2 (Enrico Mugnaini,2012) the ratios of UBCs/MLIlgi cells and UBCs/Purkinje cells increase to 7.3 and 2,respectively

const double stimp = 10;
const int NINPUT = (int)ceil(NMLI * stimp/100);


//------------------------------------------------------------------------------
//  stimulus parameter
//------------------------------------------------------------------------------
const double THR = 0;   //  threshold to cutoff the amplitude of velocity

// Number of input synapses for post cell 
const int N_iG  = 4; //4;  // iMF -> GC
const int N_eG  = 1; //4;  // eMF -> GC
const int N_MMLI = 5; //10;   //  MF -> MLI
const int N_GMLI = 4;// (int)(NGC / 16.); // 8=0.6; 16=0.3;   //  GC -> MLI
const int N_MLIG = (int)(NMLI/4.8); //4.8;   //  MLI -> GC
const int N_MLIU = 1; //6;   //  MLI -> UBC
const int N_GCPC = 100;// 100;// 1000;   //GC->PC
const int N_MLIPC = 8; //interneuron -> PC 
const double DIMGG = 0.4;


//------------------------------------------------------------------------------
//  synaptic simulation parameter
//------------------------------------------------------------------------------
// synapse weight in units of nS
const double AMPA_MAXiG = 30;     // 1.6  for single MF-GC
const double AMPA_MAXeG = 30;     // 1.6  for single MF-GC
const double AMPA_MAXMMLI = 50;     // 2.5 
const double AMPA_MAXUMLI = 50;     // 2.5 
const double AMPA_MAXGMLI = 50;     // 2.5 
const double AMPA_MAXMLIG = 400;    // <20 
const double AMPA_MAXGG  = 600;    // <30 
const double weg = 1000;    
const double wig =1.6;
const double wmMLI = 3;
const double wgMLI =3.2;
const double wMLIg =  3.5;                       

const double wpc = 0.5;
const double NAR_MG   = 2.4;  // MF-GC   = NMDA  / AMPA 
const double NAR_PC = 1.4;      //GC-PC  = NMDA  / AMPA 
const double NAR_MG2  = 2.0;       // MF-GC   = AMPAslow  / AMPA  
const double NAR_MMLI = 2;// 0.0;        // GC-MLI  =  NMDA / AMPA  
const double NAR_MMLI2  = 2;        // GC-MLI  = AMPAslow  / AMPA  

const double NAR_GMLI = 3;// 0.5;// 0.0;        // GC-MLI  = NMDA  / AMPA  
const double NAR_GABA = 1.5;//0.15;     // MLI-GC(UBC)

// for gaussian
const double var = 0.1;  //this is sacle variance for weights with var = avgW*SDW

const double TRECNMDA = 12;
const double TRECAMPA = 12;

//------------------------------------------------------------------------------
//  intfire cell properties [2]
//------------------------------------------------------------------------------
const double THR_GC = -49;   
const double THR_MLI = -45;
const int SPKTIME = 1000;
const double THR_PC = -50;

const double DT = 0.1;    //IAF
const double TERR = DT; 

//------------------------------------------------------------------------------
const double starttime = 200;
 double  AMPAfastU = 0.4;
 double  AMPAslowU = 0.4;
 double  AMPArec = 50; //0.01;//
 double AMPAfac = 400;
 double GABAfastU =  0.1;
 double GABAslowU =  0.05;
#endif
