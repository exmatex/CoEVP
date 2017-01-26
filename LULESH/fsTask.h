#ifndef __FSTASK_H__
#define __FSTASK_H__

#include "ElastoViscoPlasticity.h"

///TODO: Get some data to figure out what max state size shoudl be...
const size_t MAX_CM_STATE_SIZE = 256;

//Input for a NotTaylor call that is also the input of a VPSC call (for now)
struct TaylorInput
{
	///TODO: Figure out what these are
	int dummyVal;
	//Looks like args for advance
	//Big issue is "state size", but we can hackily set a max for now
};

//Output for a NotTaylor call that is also the output of a VPSC call (for now)
struct TaylorOutput
{
	///TODO: Figure out what these are
	int dummyVal;
};

//Initialize CM to perform FS Call. Optimally we don't do this every time
void initCM(Constitutive * cm, bool use_vpsc, double D[6], double W[3]);

TaylorOutput performFSCall(TaylorInput &input, Constitutive * cm);

TaylorOutput initAndPerformFSCall(TaylorInput &input, bool use_vpsc, double D[6], double W[3]);

void cleanUpCM(Constitutive *cm);

#endif
