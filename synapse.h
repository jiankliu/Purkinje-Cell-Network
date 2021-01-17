

//use nr MLIdes
#include "../nr3/nr3.h"
#include "../nr3/ran.h"
#include "../nr3/gamma.h"
#include "../nr3/deviates.h"
#include "../nr3/distributions.h"
#pragma once
//#include <MLInio.h>  //  getch();   
#include <cstdlib>  // pause

#include "param.h"          // file for all parameters 

//---------------constANTs initialization----------------------------
const double  PI = 3.1415926;

#define MIN(X,Y)  ( X<Y ? X:Y )
#define MAX(X,Y)  ((X) > (Y) ?  (X) : (Y))
#define GetRandom(max) ((rand()%(int)((max)))) // random integer between 0 and max-1
//#define getrandom( min, max ) ((rand() % (int)(((max) + 1) - (min))) + (min)) // random number between min and max,
#define GetRandomDoub(max) ((double)(rand()/(double)(RAND_MAX)*((max)))) // random double between 0 and max


class cls_random{
	public:
	double pV = 0.0;
	int lambda;
	cls_random(){};
	void randomExponential();
};
void cls_random::randomExponential()
{
	
	pV = (double)rand() / (double)RAND_MAX;
	pV = (-1.0 / 1)*log(1 - pV);
}

//---- simple model Wang XJ. for GABA-A synapse with the short term plasticity (STP)---------
class Gaba_A2 {
	static double Cdur, Cmax, Deadtime;  
	double alpha, lastrelease, releaseat, TimeMLIunt;
	double C, q, S, tauR, tauD;
    
	//SHORT-TERM PLASTICITY
    double  R, u; 	

public:
	double gabaa, Delay;
    double  U, trec, tfac;
	double  RR, uu;
	Gaba_A2() { // set all parameters
		S = 0;
		tauR = 1;// 0.5;// 1;
		tauD = 10;// 1.2;//5;
		alpha = 1./tauR * pow(tauD/tauR, tauR/(tauD-tauR));
		q = 0;
		gabaa = 0, C = 0;
		Delay = 1.4;//0.6;
		releaseat = 0;
		lastrelease = -9999;
		TimeMLIunt = -1.0;	//(ms)		: time MLIunter

		// Values based on Basket Cell/type F2 from Gupta, Wang, Markram (2000)
		tfac = 20;    // (ms)	this value should be close to 0 for no facilitation 
		trec = 700;   // (ms) 	reMLIvery from depression time constant 
		U    = 0.2;  //	    percent of transmitter released on first pulse
		R    = 1;     //  Releasable pool 
	}
	void reset_gaba_U() {
		u    = U;     //  for running value of U 	
		q = 0;
		R = 1;
		C = 0;
		gabaa = 0;
		S = 0;
		releaseat = 0;
		TimeMLIunt = -1.0;
		lastrelease = -9999;
	}
	void calc(double lastspike, double x);
};

double Gaba_A2::Cdur = 1., Gaba_A2::Cmax = 1., Gaba_A2::Deadtime = 1;

void Gaba_A2::calc(double lastspike, double x) {

	releaseat = lastspike + Delay; 
	q = (x - lastrelease - Cdur);  //time since last release ended
    TimeMLIunt -= DT;                  //time since last release ended

	if (TimeMLIunt < -Deadtime+TERR) {
	//if ( q > Deadtime-TERR) {
		if ( fabs(x - releaseat) < TERR) {
			lastrelease = x;
			TimeMLIunt = Cdur;

			// Synaptic MLInnections displaying depression are characterized by
			// negligible values of tau facil and hence 
			
			u = U + (1 - U) * u * exp(-q / tfac);
			R = 1 + (R*(1 - u) - 1)*exp(-q / trec);
			
			//u= U + (1-U) * u * exp(-q/tfac);
			RR = R;
			uu = u;
			C = R*u;	//			: start new release, turn on 
			S += C;
		}
	//} else if (q < 0) { 
    } else if (TimeMLIunt > 0) {		//: still releasing?
		// do nothing
	} else if (C > 0) {                  
		C = 0.;
	}
	gabaa += DT * (- gabaa /tauD + 3 * S * (1-gabaa) );
	//gabaa += DT * (- gabaa /tauD + alpha*S * (1-gabaa));
	S += -DT * S/tauR;
	
}
//--- end gaba_a stp --



//---- first order kinet model for GABA-B synapse with the short term plasticity (STP)---------
class Gaba_B2 {
	static double Cdur, Cmax, Deadtime;  
	double alpha, lastrelease, releaseat, TimeMLIunt;
	double C, q, S, tauR, tauD;
    
	//SHORT-TERM PLASTICITY
    double  R, u; 	

public:
	double gabab;
	double Delay;
    double  U, trec, tfac;
	double  RR, uu;
	Gaba_B2() { // set all parameters
		S = 0;
		tauR = 5;// 8;
		tauD = 100;// 30;// 100;
		alpha = 1./tauR * pow(tauD/tauR, tauR/(tauD-tauR));
		Delay = 1.4;// 0.6;
		q = 0;
		C = 0;
		releaseat = 0;
		lastrelease = -9999;
		TimeMLIunt = -1.0;	//(ms)		: time MLIunter

		// Values based on Basket Cell/type F2 from Gupta, Wang, Markram (2000)
		tfac = 20;    // (ms)	this value should be close to 0 for no facilitation 
		trec = 700;   // (ms) 	reMLIvery from depression time constant 
		U    = 0.2;  //	    percent of transmitter released on first pulse
		R    = 1;     //  Releasable pool 
	}
	void reset_gaba_U() {
		u    = U;     //  for running value of U 

		q = 0;
		R = 1;
		C = 0;
		gabab = 0;
		S = 0;
		releaseat = 0;
		TimeMLIunt = -1.0;
		lastrelease = -9999;
	}
	void calc(double lastspike, double x);
};

double Gaba_B2::Cdur = 1., Gaba_B2::Cmax = 1., Gaba_B2::Deadtime = 1;

void Gaba_B2::calc(double lastspike, double x) {

	releaseat = lastspike + Delay; 
	q = (x - lastrelease - Cdur);  //time since last release ended
    TimeMLIunt -= DT;                  //time since last release ended

	if (TimeMLIunt < -Deadtime+TERR) {
	//if ( q > Deadtime-TERR) {
		if ( fabs(x - releaseat) < TERR) {
			lastrelease = x;
			TimeMLIunt = Cdur;
			
			// Synaptic MLInnections displaying depression are characterized by
			// negligible values of tau facil and hence 
		
			u = U + (1 - U) * u * exp(-q / tfac);
			R = 1 + (R*(1 - u) - 1)*exp(-q / trec);
		
			RR = R;
			uu = u;
			
			C = R*u;	//			: start new release, turn on 
			S += C;
		}
	//} else if (q < 0) { 
    } else if (TimeMLIunt > 0) {		//: still releasing?
		// do nothing
	} else if (C > 0) {                  
		C = 0.;
	}
	gabab += DT * (- gabab /tauD + 0.3 * S * (1-gabab) );
	//gabab += DT * (- gabab /tauD + alpha*S * (1-gabab));
	S += -DT * S/tauR;	
}
//--- end gaba_a stp --


//------------first order kiner model for slow AMPA synapse---------------------
class AMPAsl {
	static double Cdur, Cmax, Deadtime; 
	double alpha, C, S_nmda;
	double lastrelease, releaseat, TimeMLIunt;

	//SHORT-TERM PLASTICITY
    double  q, R, u; 

public:
	double nmda, tauR, tauD;
    double  U, trec, tfac;
    double  Delay;
	double  RR, uu;
	AMPAsl() {
		Delay = 1.4;// 1.0;
		q = 0;

		C = 0;
		nmda = 0;
		S_nmda = 0;
		tauR = 163;
		tauD = 504;
		alpha = 1./tauR * pow(tauD/tauR, tauR/(tauD-tauR));

		lastrelease = -9999;
		releaseat = 0;
		TimeMLIunt = -1.0;	//(ms)		: time MLIunter

		// Values based on Basket Cell/type F2 from Gupta, Wang, Markram (2000)
		tfac = 0;    // (ms)	this value should be close to 0 for no facilitation 
		trec = 500;   // (ms) 	reMLIvery from depression time constant 
		U    = 0.5;  //	    percent of transmitter released on first pulse
		R    = 1;     //  Releasable pool 
	}

	void reset_exsyn_U() {
		u    = U;     //  for running value of U	
		q = 0;
		R = 1;
		C = 0;
		nmda = 0;
		S_nmda = 0;
		releaseat = 0;
		TimeMLIunt = -1.0;
		lastrelease = -9999;
	}

	void calc(double lastspike, double x);
};

double AMPAsl::Cdur = 1, AMPAsl::Cmax = 1., AMPAsl::Deadtime = 1;

void AMPAsl::calc(double lastspike, double x) {
	
    // FIND OUT THERE WAS A SPIKE 
	releaseat = lastspike + Delay; 

	q = ((x - lastrelease) - Cdur);  //time since last receptor replease
    TimeMLIunt -= DT;                  //time since last release ended

	if (TimeMLIunt < -Deadtime+TERR) {
	//if ( q > Deadtime-TERR) {
		if ( fabs(x - releaseat) < TERR) {
			lastrelease = x;
			TimeMLIunt = Cdur;
			
			// Synaptic MLInnections displaying depression are characterized by
			// negligible values of tfacil and hence 
			
			u = U + (1 - U) * u * exp(-q / tfac);
			R = 1 + (R*(1 - u) - 1)*exp(-q / trec);
			//printf("R=%1f,U=%1f\n", R, u);
		
			RR = R;
			uu = u;
			
			C = R*u;	//			: start new release, turn on 
			S_nmda += C;
		}
	//} else if (q < 0) { 
    } else if (TimeMLIunt > 0) {		//: still releasing?
		// do nothing
	} else if (C > 0) {                  
		C = 0.;
	}
	
	nmda += DT * (- nmda /tauD + 0.3 * S_nmda * (1-nmda) );
	//nmda +=2DT * (- nmda /tauD + alpha*S_nmda * (1-nmda));
	S_nmda += -DT * S_nmda/tauR;
}

//------------first order kiner model for AMPA and NMDA synapse---------------------
class AMPA_NMDA2 {
	static double Cdur, Cmax, Deadtime; 
	double Alpha1, Alpha2, C, S_ampa, S_nmda;
	double lastrelease, releaseat, TimeMLIunt;

	//SHORT-TERM PLASTICITY
    double  q, R, u, R2, u2; 

public:
	double ampa, nmda, tauR1, tauD1, tauR2, tauD2;
    double  U, trec, tfac, U2, trec2, tfac2;
    double  Delay;
	double  RR, uu, RR1, uu1;
	AMPA_NMDA2() {
		Delay = 1.4;// 1.0;
		q = 0;

		ampa = 0, C = 0;
		nmda = 0;
		S_ampa = 0;
		S_nmda = 0;
		tauR1 = 0.4;
		tauD1 = 2;
		Alpha1 = 1./tauR1 * pow(tauD1/tauR1, tauR1/(tauD1-tauR1));
		tauR2 = 5;
		tauD2 = 100;
		Alpha2 = 1./tauR2 * pow(tauD2/tauR2, tauR2/(tauD2-tauR2));

		lastrelease = -9999;
		releaseat = 0;
		TimeMLIunt = -1.0;	//(ms)		: time MLIunter

		// Values based on Basket Cell/type F2 from Gupta, Wang, Markram (2000)
		tfac = 0;    // (ms)	this value should be close to 0 for no facilitation 
		trec = 500;   // (ms) 	reMLIvery from depression time constant 
		U    = 0.5;  //	    percent of transmitter released on first pulse
		R    = 1;     //  Releasable pool 

		// Values based on Basket Cell/type F2 from Gupta, Wang, Markram (2000)
		tfac2 = 0.01;     // (ms)	this value should be close to 0 for no facilitation 
		trec2 = 0.01;   // (ms) 	reMLIvery from depression time constant 
		U2    = 0.05;  //	    percent of transmitter released on first pulse
		R2    = 1;     //  Releasable pool 
	}

	void reset_exsyn_U() {
		u    = U;     //  for running value of U	
		u2   = U2;     //  for running value of U	
		q = 0;

		ampa = 0, C = 0;
		nmda = 0;
		S_ampa = 0;
		S_nmda = 0;
		//U = 0.5;  //	    percent of transmitter released on first pulse
		R = 1;
		//U2 = 0.05;  //	    percent of transmitter released on first pulse
		R2 = 1;
		lastrelease = -9999;
	}

	void calc(double lastspike, double x);
};

double AMPA_NMDA2::Cdur = 1, AMPA_NMDA2::Cmax = 1., AMPA_NMDA2::Deadtime = 1;

void AMPA_NMDA2::calc(double lastspike, double x) {
	
    // FIND OUT THERE WAS A SPIKE 
	releaseat = lastspike + Delay; 

	q = ((x - lastrelease) - Cdur);  //time since last receptor replease
    TimeMLIunt -= DT;                  //time since last release ended

	if ( TimeMLIunt < -Deadtime+TERR) {
	//if ( q > Deadtime-TERR) {
		if ( fabs(x - releaseat) < TERR) {
			lastrelease = x;
			TimeMLIunt = Cdur;
			
			// Synaptic MLInnections displaying depression are characterized by
			// negligible values of tfacil and hence 
			
			u = U + (1 - U) * u * exp(-q / tfac);
			R = 1 + (R*(1 - u) - 1)*exp(-q / trec);
		
			
			RR = R;
			uu = u;
			C = R*u;	//			: start new release, turn on 
			S_ampa += R*u;

			
			u2 = U2 + (1-U2) * u2 * exp(-q/tfac2);		
			R2 = 1 + (R2*(1 - u2) - 1) * exp(-q / trec2);
			S_nmda += R2*u2;
			RR1 = R2;
			uu1 = u2;
		}
	//} else if (q < 0) { 
    } else if (TimeMLIunt > 0) {		//: still releasing?
		// do nothing
	} else if (C > 0) {                  
		C = 0.;
	}
	
	ampa += DT * (- ampa /tauD1 + 3. * S_ampa * (1-ampa) );
	//ampa += DT * (- ampa /tauD1 + Alpha1 * S_ampa* (1-ampa) );
	S_ampa += -DT * S_ampa/tauR1;
	nmda += DT * (- nmda /tauD2 + 0.35 * S_nmda * (1-nmda) );
	//nmda += DT * (- nmda /tauD2 + Alpha2 * S_nmda* (1-nmda));
	S_nmda += -DT * S_nmda/tauR2;
}

//------------UBC type: first order kiner model for mGluR2 Inh synapse---------------------
class MGR2 {
	static double Cdur, Cmax, Deadtime; 
	double alpha, C, lastrelease, releaseat, TimeMLIunt;

	//SHORT-TERM PLASTICITY
    double  q, R, u; 

public:
	double mGluR2;
    double  U, trec, tfac;
    double  xx, Delay, tauD, tauR;
	MGR2() {
		Delay = 1.;
		tauD = 500; 
		tauR = 150; 
		alpha = 1./tauR * pow(tauD/tauR, tauR/(tauD-tauR));
		q = 0;

		C = 0;
		mGluR2 = 0;
		xx = 0;

		lastrelease = -9999;
		releaseat = 0;
		TimeMLIunt = -1.0;	//(ms)		: time MLIunter

		// Values based on Basket Cell/type F2 from Gupta, Wang, Markram (2000)
		tfac = 0.001;    // (ms)	this value should be close to 0 for no facilitation 
		trec = 0.001;   // (ms) 	reMLIvery from depression time constant 
		U    = 0.1; //0.05;  //	    percent of transmitter released on first pulse
		R    = 1;     //  Releasable pool 
	}

	void reset_exsyn_U() {
		u    = U;     //  for running value of U	
		q = 0;

		C = 0;
		mGluR2 = 0;
		xx = 0;

		lastrelease = -9999;
		releaseat = 0;
		TimeMLIunt = -1.0;	//(ms)		: time MLIunter
		U = 0.1; //0.05;  //	    percent of transmitter released on first pulse
		R = 1;     //  Releasable pool 
	}

	void calc(double lastspike, double x);
};

double MGR2::Cdur = 1, MGR2::Cmax = 1., MGR2::Deadtime = 1;

// expilcity solving M2
void MGR2::calc(double lastspike, double x) {
	
    // FIND OUT THERE WAS A SPIKE 
	releaseat = lastspike + Delay; 

	q = ((x - lastrelease) - Cdur);  //time since last receptor replease
    TimeMLIunt -= DT;                  //time since last release ended

	if (TimeMLIunt < -Deadtime+TERR) {
		if ( fabs(x - releaseat) < TERR) {
			lastrelease = x;
			TimeMLIunt = Cdur;
			
			// Synaptic MLInnections displaying depression are characterized by
			// negligible values of tfacil and hence 
			//R=1+(R*(1-u)-1)*exp(-q/trec);
			//u= U + (1-U) * u * exp(-q/tfac);			
			
			//C = R*u;	//			: start new release, turn on 
			xx += U;// C;
		}
	//} else if (q < 0) { 
    } else if (TimeMLIunt > 0) {		//: still releasing?
		// do nothing
	} else if (C > 0) {                  
		C = 0.;
	}
	
	// R2'=-(R2/tauNMDA)
	// Update nmda
	mGluR2 += DT * (- mGluR2 /tauD + 0.02 * xx * (1-mGluR2) );
	//mGluR2 += DT * (- mGluR2 /tauD + alpha*xx * (1-mGluR2));
	xx += -DT * xx/tauR;
}

//------------UBC type: first order kiner model for mGluR2 Inh synapse---------------------
class MGR1 {
	static double Cdur, Cmax, Deadtime; 
	double C, lastrelease, releaseat, TimeMLIunt;

	//SHORT-TERM PLASTICITY
    double  q, R, u; 

public:
	double mGluR1;
    double  U, trec, tfac;
    double  xx, Delay, tauD, tauR;
	MGR1() {
		Delay = 1.;
		tauD = 1290.*2.0;
		tauR = 146.*2.0;
		q = 0;

		C = 0;
		mGluR1 = 0;
		xx = 0;

		lastrelease = -9999;
		releaseat = 0;
		TimeMLIunt = -1.0;	//(ms)		: time MLIunter

		// Values based on Basket Cell/type F2 from Gupta, Wang, Markram (2000)
		tfac = 0.001;    // (ms)	this value should be close to 0 for no facilitation 
		trec = 0.001;   // (ms) 	reMLIvery from depression time constant 
		U    = 0.03;  //	    percent of transmitter released on first pulse
		R    = 1;     //  Releasable pool 
	}

	void reset_exsyn_U() {
		u    = U;     //  for running value of U	
		q = 0;

		C = 0;
		mGluR1 = 0;
		xx = 0;

		lastrelease = -9999;
		releaseat = 0;
		TimeMLIunt = -1.0;	//(ms)		: time MLIunter
		U = 0.03;  //	    percent of transmitter released on first pulse
		R = 1;     //  Releasable pool 
	}

	void calc(double lastspike, double x);
};

double MGR1::Cdur = 1, MGR1::Cmax = 1., MGR1::Deadtime = 1;

// expilcity solving M2
void MGR1::calc(double lastspike, double x) {
	
    // FIND OUT THERE WAS A SPIKE 
	releaseat = lastspike + Delay; 

	q = ((x - lastrelease) - Cdur);  //time since last receptor replease
    TimeMLIunt -= DT;                  //time since last release ended

	if (TimeMLIunt < -Deadtime+TERR) {
		if ( fabs(x - releaseat) < TERR) {
			lastrelease = x;
			TimeMLIunt = Cdur;
			
			// Synaptic MLInnections displaying depression are characterized by
			// negligible values of tfacil and hence 
			u = U + (1 - U) * u * exp(-q / tfac);
			R=1+(R*(1-u)-1)*exp(-q/trec);
				
			
			C = R*u;	//			: start new release, turn on 
			xx += C;
		}
	//} else if (q < 0) { 
    } else if (TimeMLIunt > 0) {		//: still releasing?
		// do nothing
	} else if (C > 0) {                  
		C = 0.;
	}	
	mGluR1 += DT * (- mGluR1 /tauD + 0.004 * xx * (1-mGluR1) );
	xx += -DT * xx/tauR;
}




//------------MF-GC inputs plastic soma ---------------------
class MGsoma : public AMPA_NMDA2, public AMPAsl {

public:
	double ampaMG, nmdaMG, ampaMGsl;
	double AMPAsR, AMPAsu, AMPAfR, AMPAfu, NMDAfR, NMDAfu;

	MGsoma() {
		ampaMG = 0;
		nmdaMG = 0;
		AMPA_NMDA2::U = 0.5;// 0.5;
		AMPA_NMDA2::tfac = 0.01;//12;// 0.01;// 12;
		AMPA_NMDA2::trec = 0.01;//12;// 0.01;// 12;
		AMPA_NMDA2::tauR1 = 0.1;// 0.5;
		AMPA_NMDA2::tauD1 = 0.1;// 1.2;

		AMPA_NMDA2::U2 = 0.05;
		AMPA_NMDA2::tfac2 =  0.001;
		AMPA_NMDA2::trec2 = 0.001; //12
		AMPA_NMDA2::tauR2 = 0.1;// 8;
		AMPA_NMDA2::tauD2 = 0.1;//30;

		ampaMGsl = 0;
		AMPAsl::U = 0.1;
		AMPAsl::tfac = 0.01;//12;// 0.01;// 12;
		AMPAsl::trec = 0.01;//12;// 0.01;// 12; //12s
		AMPAsl::tauR = 0.1;// 0.5;
		AMPAsl::tauD = 0.1;// 5.0;
	}
	void reset() {
		AMPA_NMDA2::reset_exsyn_U();
		AMPAsl::reset_exsyn_U();

	}
	void calc(double lastspk, double x) {
		AMPA_NMDA2::calc(lastspk, x);
		ampaMG = AMPA_NMDA2::ampa;
		nmdaMG = AMPA_NMDA2::nmda;
		AMPAsl::calc(lastspk, x);
		ampaMGsl = AMPAsl::nmda;
		AMPAsR = AMPAsl::RR;
		AMPAsu = AMPAsl::uu;
		AMPAfR = AMPA_NMDA2::RR;
		AMPAfu = AMPA_NMDA2::uu;
		NMDAfR = AMPA_NMDA2::RR1;
		NMDAfu = AMPA_NMDA2::uu1;


	}
};

//------------ MLI-PC plastic soma ---------------------
class MLIsoma : public Gaba_A2, public Gaba_B2 {

public:
	double delayA, delayB;
	double MLIAR, MLIAu, MLIBR, MLIBu;
	MLIsoma() {
		//fast synapse depression
		Gaba_A2::U = 0.1;
		Gaba_A2::tfac =  800;
		Gaba_A2::trec =  100;


		//slow synapse facilitation
		Gaba_B2::U = 0.05;
		Gaba_B2::tfac =  100;
		Gaba_B2::trec = 800;
	}
	void reset() {
		Gaba_A2::reset_gaba_U();
		Gaba_B2::reset_gaba_U();
	}
	void para() {
		Gaba_A2::Delay = delayA;
		Gaba_B2::Delay = delayB;
	}
	void assignment() {
		Gaba_A2::U = GABAfastU;
		Gaba_B2::U = GABAslowU;//0.4
		//AMPA_NMDA2::tfac = AMPAfac;// 0.01;// 400;//0.01;
		//AMPA_NMDA2::trec = AMPArec;// 0.01;// 50;// 0.01; //12
		//AMPAsl::tfac = AMPAfac;// 0.01;// 400;
		//AMPAsl::trec = AMPArec;// 0.01;// 50;//0.01;400;
	}
	void calc(double lastspike, double x) {
		Gaba_A2::calc(lastspike, x);
		Gaba_B2::calc(lastspike, x);
		MLIAR = Gaba_A2::RR;
		MLIAu = Gaba_A2::uu;
		MLIBR = Gaba_B2::RR;
		MLIBu = Gaba_B2::uu;


	}
};

//------------XXXXXXXXXX GC-PC inputs plastic somaXXXXXXXXXXXX ---------------------
class GCsoma: public AMPA_NMDA2,public AMPAsl {

public:
	double ampaGP, nmdaGP, ampaslGP;
	double AMPAsR, AMPAsu, AMPAfR, AMPAfu;
	GCsoma() {
		ampaGP=0;
		ampaslGP = 0;
		nmdaGP=0;
		//AMPAfastU = 0.4;
		//AMPAslowU = 0.4;
		// fast AMPA
		
		AMPA_NMDA2::U = 0.4;
		AMPA_NMDA2::tfac = 0.01;// 400;// 0.01;// 400;//0.01;
		AMPA_NMDA2::trec = 0.01;//50;// 0.01;// 50;// 0.01; //12
		//NMDA
		AMPA_NMDA2::tfac2 = 30;//  500;//0.01;
		AMPA_NMDA2::trec2 = 500;// 30;// 0.01; //12

		AMPA_NMDA2::tauR1 =0.8;
		AMPA_NMDA2::tauD1 = 3;
		// if NMDA tauR2 and tauD2 decrease result in PC Firing rate decrease;
		AMPA_NMDA2::tauR2 = 0.01;// 8;//8.;
		AMPA_NMDA2::tauD2 = 0.01;//30;//30.;
		// AMPA Slow
		AMPAsl::U = 0.4;//0.4
		AMPAsl::tfac = 0.01;// 400;// 0.01;// 400;
		AMPAsl::trec = 0.01;//50;// 0.01;// 50;//0.01;400;
		AMPAsl::tauD = 5;
		AMPAsl::tauR =10;
		//assignment();
	}
	void reset() {
		AMPA_NMDA2::reset_exsyn_U();
		AMPAsl::reset_exsyn_U();
	}
	void assignment() {
		AMPA_NMDA2::U = AMPAfastU;
		AMPAsl::U = AMPAslowU;//0.4
		AMPA_NMDA2::tfac =  AMPAfac;// 0.01;// 400;//0.01;
		AMPA_NMDA2::trec =  AMPArec;// 0.01;// 50;// 0.01; //12
		AMPAsl::tfac =  AMPAfac;// 0.01;// 400;
		AMPAsl::trec =  AMPArec;// 0.01;// 50;//0.01;400;
	}

 	void calc(double lastspk, double x){
		AMPA_NMDA2::calc(lastspk, x); 
		ampaGP = AMPA_NMDA2::ampa;
		nmdaGP = AMPA_NMDA2::nmda;
		AMPAsl::calc(lastspk, x);
		ampaslGP = AMPAsl::nmda;
		AMPAsR = AMPAsl::RR;
		AMPAsu = AMPAsl::uu;
		AMPAfR = AMPA_NMDA2::RR;
		AMPAfu = AMPA_NMDA2::uu;
		//nmdaMG = 0;
	} 
};



//********GC-interneuron
class GMLIsoma : public AMPA_NMDA2 {

public:
	double ampaGMLI, nmdaGMLI;
	GMLIsoma() {
		ampaGMLI = 0;
		nmdaGMLI = 0;
		AMPA_NMDA2::U =  0.1;
		AMPA_NMDA2::tfac =  50;//0.01;
		AMPA_NMDA2::trec = 100;// 0.01; //12
		AMPA_NMDA2::tauR1 = 1;
		AMPA_NMDA2::tauD1 = 1.5;
		// if NMDA tauR2 and tauD2 decrease result in PC Firing rate decrease;
		AMPA_NMDA2::tfac2 = 100;//0.01;
		AMPA_NMDA2::trec2 =  50;// 0.01; //12
		AMPA_NMDA2::U2 = 0.07;
		AMPA_NMDA2::tauR2 = 5;//8.;
		AMPA_NMDA2::tauD2 = 20;//30.;
	}
	void reset() {
		AMPA_NMDA2::reset_exsyn_U();
	}

	void calc(double lastspk, double x) {
		AMPA_NMDA2::calc(lastspk, x);
		ampaGMLI = AMPA_NMDA2::ampa;
		nmdaGMLI = AMPA_NMDA2::nmda;
		//nmdaMG = 0;
	}
};








//------------  GC-MLI current ---------------------
class IMF {
	static double E1, E2;

public:
	double I, AMPA_NMDA_RATIO;
	IMF() {
		I = 0;
		AMPA_NMDA_RATIO = NAR_MG;
	}
	void calc(double g_AMPA, double ampa, double nmda, double y_post,  double postB);
};

double IMF::E1 = 0., IMF::E2 = 2.404;// 0.;

void IMF::calc(double g_AMPA, double ampa, double nmda, double y_post,  double postB) {
	I = g_AMPA * ampa * (y_post - E1) + g_AMPA * AMPA_NMDA_RATIO * nmda * postB* (y_post - E2);
	//I = g_AMPA * 0 * (y_post - E1) + g_AMPA * AMPA_NMDA_RATIO * nmda * postB* (y_post - E2);
}


//------------ Ex MF current ---------------------



//------------ Ex MF-UBC current ---------------------

///####################################**************PC network***************#############################

class IMFGC {
	static double E1, E2;

public:
	double I, AMPA_NMDA_RATIO, AMPA_AMPAsl_RATIO;
	double I1;
	IMFGC() {
		I = 0;
		I1 = 0;
		AMPA_NMDA_RATIO = NAR_MG;//2.4
		AMPA_AMPAsl_RATIO = NAR_MG2;//2.0
	}
	void calc(double g_AMPA, double ampa, double nmda, double ampasl, double y_post, double postB);
};

double IMFGC::E1 = 0., IMFGC::E2 = 2.404;

void IMFGC::calc(double g_AMPA, double ampa, double nmda, double ampasl, double y_post, double postB) {
	//	I = g_AMPA * ampa * (y_post - E1) + g_AMPA * AMPA_NMDA_RATIO * nmda * postB* (y_post - E2)
	//		+ g_AMPA * AMPA_AMPAsl_RATIO * ampasl * (y_post - E1);
	I = 1 * g_AMPA * ampa * (y_post - E1) + 1 * g_AMPA * AMPA_NMDA_RATIO * nmda * postB* (y_post - E2)
		+ 1 * g_AMPA * AMPA_AMPAsl_RATIO * ampasl * (y_post - E1);
	//I1 = g_AMPA * AMPA_NMDA_RATIO * nmda * postB* (y_post - E2);

}
//-----------------GC->PC------------
class IPC {
	static double E1, E2;

public:
	double I, AMPA_NMDA_RATIO, AMPA_AMPAsl_RATIO;
	IPC() {
		I = 0;
		AMPA_NMDA_RATIO = NAR_PC;
		AMPA_AMPAsl_RATIO = NAR_PC;
	}
	void calc(double g_AMPA, double ampa, double nmda, double y_post, double postB);
};

double IPC::E1 = 0., IPC::E2 = 0.;

void IPC::calc(double g_AMPA, double ampa, double ampasl, double y_post, double postB) {
	//I = g_AMPA * ampa * (y_post - E1) + g_AMPA * AMPA_NMDA_RATIO * nmda * 1* (y_post - E2);
	I = g_AMPA * ampa * (y_post - E1) + g_AMPA * AMPA_AMPAsl_RATIO * ampasl * (y_post - E2);
	//I = g_AMPA * 0 * (y_post - E1) + g_AMPA * AMPA_NMDA_RATIO * nmda * postB* (y_post - E2);
}

//------------ Inh current Inter->PC ---------------------
class IMLI {
	static double E, EB;
public:
	double I, GABA_AB_RATIO, mix;
	IMLI() {
		I = 0;
		GABA_AB_RATIO = NAR_GABA;//0.15
		//GABA_AB_RATIO = 1;//0.15
	}
	void calc(double g_GABAA, double gabaa, double gabab, double y_post);
};

double IMLI::E = -80, IMLI::EB = -90;// -75., IMLI::EB = -90.;

void IMLI::calc(double g_GABAA, double gabaa, double gabab, double y_post) {
	//I = mix * g_GABAA * gabaa * (y_post - E) + (1-mix) * g_GABAA * GABA_AB_RATIO * gabab * (y_post - EB);
	I = g_GABAA * gabaa * (y_post - E) + g_GABAA * GABA_AB_RATIO * gabab * (y_post - EB);
}
/*

modified from pregen.mod - Mainen and Destexhe 
modified to inMLIrporate arbitrary sequences of pulses with 
assigned inter-burst intervals
spkgenmode = 0 is the same as the previous version
spkgenmode = 1 permits the user to assign inter-burst intervals
in an aperiodic fashion. In this mode the on_times variable
must be initialized to +9e9 in the .oc or .hoc file.
Potential bug was fixed by adding dt/2
Dean Buonomano 5/20/99

bugs were fixed
Jian Liu 
*/


class INPUT : public Poissondev, public Expondist, public Expondev, public cls_random{
  static double spikedur, refact;
  double on_times[10];
  int spkgenmode;
  double burst;
  int    pulsecntr;

  double lcntr, scntr, bcntr;

public:
  //int SEED1=2497;
  double Istim, Spike;
  double noise, fast_invl, slow_invl, burst_len, start, end;
  double Startburst;
  INPUT() : Poissondev(1, SEED), Expondist(SEED), Expondev(1, SEED1){
	  spkgenmode = 0;
	  fast_invl  = 1;      //: time between spikes in a burst (msec)
	  slow_invl  = 50;     //: burst period (msec)
	  burst_len  = 3.;     //: burst length (# spikes)
	  start      = 50.;    //: location for first burst
	  end        = 200.;   //: time to stop bursting
	  noise = 0;// 0;
	  Spike      = 0.;  
	  pulsecntr  = -1;
	  burst      = 0;
	  Startburst = 0;
  }

  void reset_input() {
	  scntr = fast_invl;
	  bcntr = burst_len;

	  if (noise != 0) {
		  scntr	= noise * Poissondev::dev(fast_invl ) + (1 - noise) * fast_invl;
		  bcntr	= noise * Poissondev::dev(burst_len ) + (1 - noise) * burst_len;
	  }
	  lcntr = start;
  }

  void calc( double y_post, double x);
};

void INPUT::calc(double y_post, double x) {
	if (x < end) {
		
		Spike = 0.;
		scntr -= DT;//if burst
		lcntr -= DT;
#if POISSON
		if (Startburst) {
			scntr -= DT;
			lcntr += DT;
			if (burst) {
				// in a burst
				if (scntr < DT + TERR) {
					if (bcntr <= 1.0) {
						// last spike in burst?
						burst = 0;
					}
					Spike = 10;
					bcntr -= 11.1;
					if (noise == 0) {
						scntr = fast_invl + DT;
					}
					else {
						scntr = noise * Poissondev::dev(fast_invl) + (1 - noise) * fast_invl;
						//scntr = 1 * Expondev::dev()*slow_invl + (1 - 1) * slow_invl;
					}
				}

			}
			else {
				bcntr = burst_len;
				burst = 1;
			}
		
		}
		else {
			if (lcntr <= DT + TERR) {
				if (noise == 0) {
					//lcntr = slow_invl - (burst_len-1) * fast_invl;
					cls_random::randomExponential();
					lcntr = 1 * cls_random::pV*slow_invl;
				    //lcntr = slow_invl;

				}
				else {
					//lcntr = noise * Poissondev::dev(slow_invl) + (1 - noise) * slow_invl;
					lcntr = noise * Expondev::dev()*slow_invl + (1 - noise) * slow_invl;
					//printf("%4f\n", lcntr);	
				}
				Spike = 10;
			}
		}
#endif
		/****************Generate Burst spike*****************/
#if BURST
		if ( burst ) {
			// in a burst
			if ( scntr < DT + TERR ) {     
			 //if (lcntr < DT + TERR) {
				// a spike
				//if ( fabs(bcntr-spikedur) < TERR ) {
				if ( bcntr <= 1.0 ) {
					// last spike in burst?
					burst = 0;
					if ( spkgenmode ==0 ) {
						if ( noise ==0 ) {
							lcntr = slow_invl;
							//lcntr = slow_invl - (burst_len-1) * fast_invl;
							//lcntr=1 * Expondev::dev()*slow_invl;
							
						} else {
							//lcntr = noise * Poissondev::dev(slow_invl) + (1 - noise) * slow_invl;
							lcntr = noise * Poissondev::dev(slow_invl) + (1 - noise) * slow_invl;
							//printf("%4f\n", lcntr);							
						}					
					} else if ( spkgenmode ==1 ) {
						lcntr = on_times[pulsecntr] + DT;
					}
				}

				Spike = 10;
				bcntr -= 11.1;

				if (noise==0) {
					scntr = fast_invl + DT;
				} else {
					scntr = noise * Poissondev::dev(fast_invl) + (1 - noise) * fast_invl; 
					//scntr = 1 * Expondev::dev()*slow_invl + (1 - 1) * slow_invl;
				}
			}
		} else {                  
			//  between burst
			//if ( fabs(lcntr-DT) < TERR) {	//	: there is a round off problem here
			if (lcntr <= DT + TERR) {	//	: there is a round off problem here
				if (noise==0) { 
					bcntr = burst_len;
				} else {
					bcntr = noise * Poissondev::dev(burst_len) + (1 - noise) * burst_len;
					//bcntr = 1 * Expondev::dev()*slow_invl + (1 - 1) * slow_invl;
				}
				burst = 1;
				//Spike = 10;
				pulsecntr = pulsecntr + 1;	
			}
		} 
#endif

	} // x end
	Istim = Spike * ( y_post  + 20 );
}


//-----------------------------------------------------------------------
//            Define Cells
//-----------------------------------------------------------------------

//-------------------extral MF -----------------------------------------------
class MFIAF: public MGsoma, public GMLIsoma, public INPUT {
	static double V0; 

	static double VON, VOFF; 
	static double spkdur, refact; 
	double Imax,rr;

public:
	double Cm, G_l, E_l;
	int    TYPE, Ca;  //Ca is spike number
	double spiketimes[SPKTIME], postB;
	//using Expondev::Expondev;
    double Vthr, q, lastspk;
	double tauAHP, gAHPbar, eAHP;
	double k, mu, Inoise, TimeCounter;

	MFIAF(){
		TYPE = 0;
		k = 0;
		mu = 0;
		
		Cm = 1;
		E_l  = -60.;      // mV    10^-3 
        G_l  = 0.025;        // mS/cm^2  The maximum specific leakage MLInductance
		lastspk = -9999;
        TimeCounter = -1.0;
		Ca   = 0;
		postB = 0;

		Imax = 1;     //  uA/cm^2
		Inoise = 0;   //  uA/cm^2
		Vthr = -40;

  	    tauAHP   = 10;
		gAHPbar = 0.00007 * 1e3;   // mS/cm2; peak of AHP current
		eAHP = -90;

		q= 0;
		for (int i=0; i<SPKTIME; i++) spiketimes[i] = 0.;
	}

	void reset(){
		Ca = 0;
		for (int i=0; i<SPKTIME; i++) spiketimes[i] = 0.;
	}

	void init(double *y) {
        y[0] = V0;
        y[1] = 0; //gAHP;
		ampaMG = 0;
		nmdaMG = 0;
		ampaMGsl = 0;
		ampaGMLI = 0;
		nmdaGMLI = 0;		
		MGsoma::reset();		
		GMLIsoma::reset();
		TimeCounter = -1.0;
		rr = 0;
	
	}
	void calc(double x, double *y, double *f); 
};  

double MFIAF::V0 = -60;

double MFIAF::VON=40, MFIAF::VOFF=-60;
double MFIAF::spkdur=1, MFIAF::refact=2;


void MFIAF::calc(double x, double *y, double *f)
{
	q = (x - lastspk) - spkdur; //time since last spike ended
	TimeCounter -= DT;

	if (TimeCounter < -refact + TERR) {
		//if ( q>(refact-TERR) ) {			// refactory period over?
		if (y[0]>Vthr) {		// threshod reached?			
			y[0] = VON;
			lastspk = x;        // spike MLIunted
			TimeCounter = spkdur;
			spiketimes[Ca] = x;
			Ca++;
		}
		//}   else if ( q < TERR ) {		// spike still on
	}
	else if (TimeCounter > 0 - TERR) {	// spike still on
		y[0] = VON;
	}
	else {                  // spike off, refactory period on
		if (y[0] > 0) {	 	// turn spike off
			y[0] = VOFF;
			y[1] += gAHPbar;
		}
	}

	// IAF
	//Inoise = (GetRandomDoub(1.0) -1 ) * Imax;
	//if ( GetRandomDoub(1.0) < DT * BACKRATE ) {
	//	Inoise = -100;  
	//}   

#if TEST
#if MFSTIM
	if (Inoise < -1.) {
		y[0] = VON;
		lastspk = x;
	}
	else {
		y[0] = VOFF;
	}
	f[0] = 0.0;
	//	f[0] = (1./Cm) * ( -G_l*(y[0]-E_l) -y[1]*(y[0]-eAHP) - Inoise ); 
#else

#endif
#else
	///*
	//INPUT::calc(y[0], x);
	//f[0] = (1. / Cm) * (-G_l*(y[0] - E_l) - y[1] * (y[0] - eAHP) - Istim);
	//rr = 0.5 * Expondev::dev()*slow_invl + (1 - 0.5) * slow_invl;
	//printf("%4f\n", rr);
	//fprintf(f10, "%1f\n", rr);
#if POISSON
	//rr = rr - DT;
	//if (rr <= 0 || Inoise < -1.) {
	//	rr = Expondev::dev()*slow_invl;
	//	y[0] = VON;
	//	if(x>0) lastspk = x;
	//}
	//else {
	//	y[0] = VOFF;
	//}
	
	
		 if (Startburst) {
			 //y[0] = VOFF;
			INPUT::calc(y[0], x);
			if (Istim == 0) {
				y[0] = VOFF;
			}
			else {
				y[0] = VON;
				//(1. / Cm) * (-G_l*(y[0] - E_l) - y[1] * (y[0] - eAHP) - Istim);
			}
		 }else {
			 if (Inoise < -1. || rr == 1) {
				 INPUT::calc(y[0], x);
				 rr = 1;
				 if (Istim == 0) {
					 y[0] = VOFF;
				 }
				 else {
					 y[0] = VON;
					 //(1. / Cm) * (-G_l*(y[0] - E_l) - y[1] * (y[0] - eAHP) - Istim);
				 }
			 }		 		 
		 }
		f[0] = 0.0;

	
	//f[0] = 0.0;
#endif

#if MFSTIM
	if ( Inoise < -1. ) {
		y[0] = VON;
		lastspk = x;
	} else {
		y[0] = VOFF;
	}
	f[0] = 0.0;
#endif
#if BURST

		INPUT::calc(y[0], x);
		f[0] = (1. / Cm) * (-G_l*(y[0] - E_l) - y[1] * (y[0] - eAHP) - Istim);


#endif

	//*/
	//f[0] = (1./Cm) * ( -G_l*(y[0]-E_l) -y[1]*(y[0]-eAHP) - Inoise ); 
#endif

    f[1] = - y[1]/tauAHP;
 
	MGsoma::calc(lastspk, x);
	//MUsoma::calc(lastspk, x);
	GMLIsoma::calc(lastspk, x);
}



//-------------------Grangular cell model from ExpIAF by N. Brunel -----------------------------------------
class GCBRUNEL: public GCsoma,public GMLIsoma, public Normaldev {
	double V0, VOFF; 
	static double VON; 
	static double spkdur, refact; 
	

public:
	double Cm, G_l, E_l, taum, Ginh;
	int   Ca;  //Ca is spike number
	double spiketimes[SPKTIME], postB;
	double lastspk, TimeCounter;
    double Vthr, q;
	double tauAHP, gAHPbar, eAHP;

	double DeltaT, tauxAHP, xAHP, xAHPbar, G_ka, tauN, sigma;

	GCBRUNEL(): Normaldev ( 0.0, 1.0, rand()*rand() ){
		Ginh = 0.0;
		Cm = 4.9;       // pF
		taum = 7.0;   // ms			

		// I_Na
		DeltaT = 1.0; // mV
		Vthr = -50;
		E_l  = -90.0;       // mV
		V0 = -85.0;
		VOFF = -65.0;
        G_l  = 1.5;     // S  The maximum specific leakage MLInductance
		lastspk = -9999;
		TimeCounter = -1.0 * spkdur;
		Ca   = 0;
		postB = 0;

		tauAHP  = 3.0; //20 by niMLIlas
		gAHPbar = 1.0;   // nS; peak of AHP current
		eAHP = -90;
		tauxAHP = 1.0;
		xAHP = 0.0;
		xAHPbar = 1.0;
		G_ka = 0; //1.2; // nS

		//noise
		tauN = 1000.0;
		sigma = 0.12;


		q= 0;
		for (int i=0; i<SPKTIME; i++) spiketimes[i] = 0.;
	}

	void reset(){
		Ca = 0;
		for (int i=0; i<SPKTIME; i++) spiketimes[i] = 0.;
	}

	void init(double *y) {
        y[0] = V0;
        y[1] = 0.0; // xAHP;
		y[2] = 0.0; // zAHP;
		y[3] = 0.0; // noise
        //y[4] = 0; //EIsoma::R2;
        //y[5] = 0; //EIsoma::nmda;
		ampaGP = 0;
		ampaslGP = 0;
		nmdaGP = 0;
		ampaGMLI = 0;
		nmdaGMLI = 0;
		GCsoma::assignment();
		GCsoma::reset();
		
		GMLIsoma::reset();
		TimeCounter = -1.0 * spkdur;
	}

	void calc(double x, double *y, double *f); 
};  

double GCBRUNEL::VON=40;
double GCBRUNEL::spkdur=0.6, GCBRUNEL::refact=2;

void GCBRUNEL::calc(double x, double *y, double *f)
{
    q = (x- lastspk) - spkdur; //time since last spike ended
    TimeCounter -= DT;              
	
	if (TimeCounter < -refact+TERR) {
    //if ( q>(refact-TERR) ) {			// refactory period over?
		//if ( y[0]>=VON) {		// threshod reached?			
		if ( y[0]>=Vthr) {		// threshod reached?			
			y[0] = VON;
			lastspk = x;        // spike MLIunted
			TimeCounter = spkdur;
			spiketimes[Ca] = x;
			Ca ++;
			y[1] += xAHPbar;
		}
	//}   else if ( q < TERR ) {		// spike still on
    } else if (TimeCounter > 0-TERR) {	// spike still on
			y[0] = VON;
    } 	else {                  // spike off, refactory period on
		if ( y[0] > 0.0) {	 	// turn spike off
			y[0] = VOFF;
		}
	}

	//MLIut<<" vm is "<< y[0] <<" mv"<<endl;
	
	f[1] = - y[1]/tauxAHP;

	f[2] =  (1 - y[2]) * y[1] - y[2]/tauAHP;

	f[3] =  (- y[3] + sigma * sqrt(tauN) * Normaldev::dev()/sqrt(DT) )/tauN;

#if V_CLAMP
	y[0] = -65; // -75//mGluR2  -72//ex -70//inh
#else
	f[0] = ( -G_l*(y[0]-E_l)*exp(-(y[0]-E_l)/5.) - (G_ka+y[3])*y[0] -gAHPbar*y[2]*(y[0]-eAHP)  -TONICI/1000.0*(y[0]+75) ) / Cm;
	//f[0] = ( -G_l*(y[0]-E_l)*exp(-(y[0]-E_l)/5.) + G_l*DeltaT*exp((y[0]-Vthr)/DeltaT)  -gAHPbar*y[2]*(y[0]-eAHP) -TONICI/1000.0*(y[0]+75) ) / Cm;
	//f[0] = ( -G_l*(y[0]-E_l)*exp(-(y[0]-E_l)/5.) -gAHPbar*y[2]*(y[0]-eAHP) -TONICI/1000.0*(y[0]+75) ) / Cm;
	//f[0] = ( -G_l*(y[0]-E_l) - Inoise ) / Cm; 
#endif

	// from Stephane's data
    postB = (exp((y[0]+119.51)/38.427) + exp(-(y[0]+45.895)/28.357))
		/ (exp((y[0]+119.51)/38.427) + exp(-(y[0]+45.895)/28.357) + exp(-(y[0]-84.784)/38.427));

	// GC-MLI release
	GCsoma::calc(lastspk, x);
	GMLIsoma::calc(lastspk, x);

}


//-------------------MLIl model  -----------------------------------------
class MLIBRUNEL: public MLIsoma, public Normaldev {
	double V0, VOFF; 
	static double VON; 
	static double spkdur, refact; 
	

public:
	double Cm, G_l, E_l, taum;
	int   Ca;  //Ca is spike number
	double spiketimes[SPKTIME], postB;
	double lastspk, TimeCounter;
    double Vthr, q;
	double tauAHP, gAHPbar, eAHP;

	double DeltaT, tauxAHP, xAHP, xAHPbar, G_ka, tauN, sigma;

	MLIBRUNEL(): Normaldev ( 0.0, 1.0, rand()*rand() ){
		Cm = 20;// 60;       // pF
		taum = 20.0;   // ms			

		// I_Na
		DeltaT = 3.0; // mV
		Vthr = -45;
		E_l = -50;// -70.0;       // mV
		V0 = E_l;
		VOFF = E_l;
        G_l  = Cm / taum;     // S  The maximum specific leakage MLInductance
		lastspk = -9999;
        TimeCounter = -1.0 * spkdur;
		Ca   = 0;
		postB = 0;

		tauAHP  = 20.0; //20 by niMLIlas
		gAHPbar =  4.0;   // nS; peak of AHP current
		eAHP = -100;
		tauxAHP = 1.0;
		xAHP = 0.0;
		xAHPbar = 1.0;
		G_ka = 1.2; // nS

		//noise
		tauN = 1000.0;
		sigma = 0.12;


		q= 0;
		for (int i=0; i<SPKTIME; i++) spiketimes[i] = 0.;
	}

	void reset(){
		Ca = 0;
		for (int i=0; i<SPKTIME; i++) spiketimes[i] = 0.;
	}

	void init(double *y) {
        y[0] = V0;
        y[1] = 0.0; // xAHP;
		y[2] = 0.0; // zAHP;
		y[3] = 0.0; // noise
        //y[4] = 0; //EIsoma::R2;
        //y[5] = 0; //EIsoma::nmda;
		MLIsoma::reset();
		MLIsoma::assignment();

	}

	void calc(double x, double *y, double *f); 
};  

double MLIBRUNEL::VON=40;
double MLIBRUNEL::spkdur=1, MLIBRUNEL::refact=2;

void MLIBRUNEL::calc(double x, double *y, double *f)
{
    q = (x- lastspk) - spkdur; //time since last spike ended
    TimeCounter -= DT;              
	
	if (TimeCounter < -refact+TERR) {
    //if ( q>(refact-TERR) ) {			// refactory period over?
		if ( y[0]>=VON) {		// threshod reached?			
			y[0] = VON;
			lastspk = x;        // spike MLIunted
			TimeCounter = spkdur;
			spiketimes[Ca] = x;
			Ca ++;
			y[1] += xAHPbar;
		}
	//}   else if ( q < TERR ) {		// spike still on
    } else if (TimeCounter > 0-TERR) {	// spike still on
			y[0] = VON;
    } 	else {                  // spike off, refactory period on
		if ( y[0] > 0.0) {	 	// turn spike off
			y[0] = VOFF;
		}
	}

	//MLIut<<" vm is "<< y[0] <<" mv"<<endl;
	
	f[1] = - y[1]/tauxAHP;

	f[2] =  (1 - y[2]) * y[1] - y[2]/tauAHP;

	f[3] =  (- y[3] + sigma * sqrt(tauN) * Normaldev::dev()/sqrt(DT) )/tauN;

#if V_CLAMP
	y[0] = -65; // -75//mGluR2  -72//ex -70//inh
#else
	//f[0] = ( -G_l*(y[0]-E_l) + G_l*DeltaT*exp((y[0]-Vthr)/DeltaT) -gAHPbar*y[2]*(y[0]-eAHP) ) / Cm; 
	//f[0] = ( -G_l*(y[0]-E_l) + G_l*DeltaT*exp((y[0]-Vthr)/DeltaT) -gAHPbar*y[2]*(y[0]-eAHP) - (G_ka+y[3])*y[0] ) / Cm;
	f[0] = (-G_l*(y[0] - E_l) + G_l*DeltaT*exp((y[0] - Vthr) / DeltaT) - gAHPbar*y[2] * (y[0] - eAHP) - ( y[3])*y[0]) / Cm;
#endif

	postB = 1/(1+exp(-(y[0] - (-20))/13));
	// mGLR release
	MLIsoma::calc(lastspk, x);

}
//--------------------------One-MLImpartment-Purkinje cell------------
class Purkinjeone :  public Normaldev {
	double V0, VOFF;
	static double VON;
	static double spkdur, refact;


public:
	double Cm, G_l, E_l, taum;
	int   Ca;  //Ca is spike number
	double spiketimes[SPKTIME], postB;
	double lastspk, TimeCounter;
	double Vthr, q;
	double tauAHP, gAHPbar, eAHP;

	double DeltaT, tauxAHP, xAHP, xAHPbar, G_ka, tauN, sigma;

	Purkinjeone() : Normaldev(0.0, 1.0, rand()*rand()) {
		Cm = 250;// 60;       // pF
		taum = 20.0;   // ms			

					   // I_Na
		DeltaT = 3.0; // mV
		Vthr = -55;//-45
		E_l = -70;// -50.0;       // mV
		V0 = E_l;
		VOFF = E_l;
		G_l = Cm / taum;     // S  The maximum specific leakage MLInductance
		lastspk = -9999;
		TimeCounter = -1.0 * spkdur;
		Ca = 0;
		postB = 0;

		tauAHP = 20.0; //20 by niMLIlas
		gAHPbar = 4.0;   // nS; peak of AHP current
		eAHP = -100;
		tauxAHP = 1.0;
		xAHP = 0.0;
		xAHPbar = 1.0;
		G_ka = 1.2; // nS

					//noise
		tauN = 1000.0;
		sigma = 0.12;


		q = 0;
		for (int i = 0; i<SPKTIME; i++) spiketimes[i] = 0.;
	}

	void reset() {
		Ca = 0;
		for (int i = 0; i<SPKTIME; i++) spiketimes[i] = 0.;
	}

	void init(double *y) {
		y[0] = V0;
		y[1] = 0.0; // xAHP;
		y[2] = 0.0; // zAHP;
		y[3] = 0.0; // noise
		TimeCounter = -1.0 * spkdur;
					//y[4] = 0; //EIsoma::R2;
					//y[5] = 0; //EIsoma::nmda;
	}

	void calc(double x, double *y, double *f);
};

double Purkinjeone::VON = 40;
double Purkinjeone::spkdur = 1, Purkinjeone::refact = 2;

void Purkinjeone::calc(double x, double *y, double *f)
{
	q = (x - lastspk) - spkdur; //time since last spike ended
	TimeCounter -= DT;

	if (TimeCounter < -refact + TERR) {
		//if ( q>(refact-TERR) ) {			// refactory period over?
		if (y[0] >= VON) {		// threshod reached?			
			y[0] = VON;
			lastspk = x;        // spike MLIunted
			TimeCounter = spkdur;
			spiketimes[Ca] = x;
			Ca++;
			y[1] += xAHPbar;
		}
		//}   else if ( q < TERR ) {		// spike still on
	}
	else if (TimeCounter > 0 - TERR) {	// spike still on
		y[0] = VON;
	}
	else {                  // spike off, refactory period on
		if (y[0] > 0.0) {	 	// turn spike off
			y[0] = VOFF;
		}
	}

	//MLIut<<" vm is "<< y[0] <<" mv"<<endl;

	f[1] = -y[1] / tauxAHP;

	f[2] = (1 - y[2]) * y[1] - y[2] / tauAHP;

	f[3] = (-y[3] + sigma * sqrt(tauN) * Normaldev::dev() / sqrt(DT)) / tauN;

#if V_CLAMP
	y[0] = -65; // -75//mGluR2  -72//ex -70//inh
#else
	//f[0] = ( -G_l*(y[0]-E_l) + G_l*DeltaT*exp((y[0]-Vthr)/DeltaT) -gAHPbar*y[2]*(y[0]-eAHP) ) / Cm; 
	f[0] = (-G_l*(y[0] - E_l) + G_l*DeltaT*exp((y[0] - Vthr) / DeltaT) - gAHPbar*y[2] * (y[0] - eAHP) -  y[3]*y[0]) / Cm;

	//f[0] = (-G_l*(y[0] - E_l) - gAHPbar*y[2] * (y[0] - eAHP) - y[3] * y[0]) / Cm;
#endif

	postB = 1 / (1 + exp(-(y[0] - (-20)) / 13));
	// mGLR release
	//MLIsoma::calc(lastspk, x);

}
