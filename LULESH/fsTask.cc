#include "fsTask.h"

#include "MieGruneisen.h"
// Fine scale model options
#include "Taylor.h"        // the fine-scale plasticity model
#include "vpsc.h"

#include "domain.h"


TaylorOutput performFSCall(TaylorInput &input, Constitutive * cm)
{
	//Get inputs
	TaylorOutput output;
	Tensor2Gen tempTensor;
	for(int i = 0; i < 9; i++)
	{
		tempTensor.a[i] = input.tensor_a[i];
	}
	//Advance models
	ConstitutiveData dat = cm->advance(input.deltatime, tempTensor, input.cm_vol_chng, input.cm_state);
	//Copy outputs out
	output.num_models = dat.num_models;
	output.num_point_value_pairs = dat.num_point_value_pairs;
	output.num_Newton_iters = dat.num_Newton_iters;
	for(int i = 0; i < 6; i++)
	{
		output.tensor_a[i] = dat.sigma_prime.a[i];
	}
	//Copy the state out
	std::memcpy(&output.cm_state, &input.cm_state, MAX_CM_STATE_SIZE);

	return output;
}

void initCM(Constitutive * cm, bool use_vpsc, double D[6], double W[3])
{
	///TODO: Make this more robust and able to handle adaptive sampling
	// Current version doesn't support adaptive sampling as that is far more to consider than just an FS call
	
	//Assuming cm is a nullptr. Probably do something cleaner/less horrible

	///TODO: c_scaling should be an arg
	double c_scaling = 1.0;


	//Do we have to do anything with cm_global? It just seems to exist briefly
   ConstitutiveGlobal cm_global;

         Tensor2Gen L;

         L(1,1) = D[0];         // dxddx
         L(1,2) = D[5] - W[2];  // dyddx
         L(1,3) = D[4] + W[1];  // dzddx
         L(2,1) = D[5] + W[2];  // dxddy 
         L(2,2) = D[1];         // dyddy
         L(2,3) = D[3] - W[0];  // dzddy
         L(3,1) = D[4] - W[1];  // dxddz
         L(3,2) = D[3] + W[0];  // dyddz
         L(3,3) = D[2];         // dzddz


      double bulk_modulus = 1.94; // Tantallum (Mbar)
      double shear_modulus = 6.9e-1; // Tantallum (Mbar)


      // Construct the equation of state
      EOS* eos_model;
      {
         /* From Table 1 (converted from GPa to Mbar) in P. J. Maudlin et al.,
            "On the modeling of the Taylor cylinder impact test for orthotropic
            textured materials: experiments and simulations", Inter. J.
            Plasticity 15 (1999), pp. 139-166.
            */
         double k1 = 1.968;  // Mbar
         double k2 = 2.598;  // Mbar
         double k3 = 2.566;  // Mbar
         double Gamma = 1.60;  // dimensionless
         eos_model = (EOS*)(new MieGruneisen(k1, k2, k3, Gamma));
      }

      Plasticity* plasticity_model;

         // Construct the fine-scale plasticity model
         // These values are needed for both models now
         double D_0 = 1.e-2;
         double m = 1./20.;
         double g = 2.e-3; // (Mbar)
         //      double m = 1./2.;
         //      double g = 1.e-4; // (Mbar) Gives a reasonable looking result for m = 1./2.
         //      double m = 1.;
         //      double g = 2.e-6; // (Mbar) Gives a reasonable looking result for m = 1.
         //
      if (use_vpsc == true) {
         // New vpsc inititialization
         plasticity_model = (vpsc*) (new vpsc(D_0, m, g, c_scaling));
      } else {
         // Old Taylor initialization
         //

         plasticity_model = (Plasticity*)(new Taylor(D_0, m, g));
      }

/*
         size_t state_size;
         domain.cm(i) = (Constitutive*)(new ElastoViscoPlasticity(cm_global, ann, modelDB, L, bulk_modulus, shear_modulus, eos_model,
                  plasticity_model, sampling, state_size));
         domain.cm_state(i) = operator new(state_size);
         domain.cm(i)->getState(domain.cm_state(i));
*/
	size_t state_size;
	cm = (Constitutive*)new ElastoViscoPlasticity(cm_global, nullptr, nullptr, L, bulk_modulus, shear_modulus, eos_model, plasticity_model, false, state_size);
}

TaylorOutput initAndPerformFSCall(TaylorInput &input, bool use_vpsc, double D[6], double W[3])
{
	//Call initCM and use the resulting model to call performFSCall
	Constitutive * cm = nullptr;
	initCM(cm, use_vpsc, D, W);
	TaylorOutput output = performFSCall(input, cm);
	return output;
}

void cleanUpCM(Constitutive *cm)
{
	//ANything more than delete?
	delete cm;
}
