#ifndef __FSTASK_H__
#define __FSTASK_H__

#include "ElastoViscoPlasticity.h"

///TODO: Get some data to figure out what max state size shoudl be...
const size_t MAX_CM_STATE_SIZE = 2048;

//Input for a NotTaylor call that is also the input of a VPSC call (for now)
struct TaylorInput
{
	double deltatime;
	//Tensor2Gen is a pointer to a class
	double tensor_a[9];
	double cm_vol_chng;
	//Big issue is "state size", but we can hackily set a max for now
	char cm_state[MAX_CM_STATE_SIZE];
};

//Output for a NotTaylor call that is also the output of a VPSC call (for now)
struct TaylorOutput
{
	double tensor_a[6];
	int num_models;
	int num_point_value_pairs;
	int num_Newton_iters;
	//Looks like Constitutive_Data and state
	char cm_state[MAX_CM_STATE_SIZE];
};

//Initialize CM to perform FS Call. Optimally we don't do this every time
Constitutive *  initCM(bool use_vpsc, double D[6], double W[3], void * cm_state);

TaylorOutput performFSCall(TaylorInput &input, Constitutive * cm);

TaylorOutput initAndPerformFSCall(TaylorInput &input, bool use_vpsc, double D[6], double W[3]);

void cleanUpCM(Constitutive *cm);

#endif
