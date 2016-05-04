#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <stdlib.h>
#include "multinest.h"
#include "DMUtils.h"

double loglikelihood(double * params, int* num_hard, double *result);


/******************************************** loglikelihood routine ****************************************************/

// Now an example, sample an egg box likelihood

// Input arguments
// ndim 						= dimensionality (total number of free parameters) of the problem
// npars 						= total number of free plus derived parameters
// context						void pointer, any additional information
//
// Input/Output arguments
// Cube[npars] 						= on entry has the ndim parameters in unit-hypercube
//	 						on exit, the physical parameters plus copy any derived parameters you want to store with the free parameters
//	 
// Output arguments
// lnew 						= loglikelihood

void LogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{
	//Map from unit-cube to parameter cube
    double *param_min = new double[ndim];
	double *param_max = new double[ndim];

	param_min[0] = 1;
	param_max[0] = 3;
	
	param_min[1] = -40;
	param_max[1] = -38;
	
	/*
	param_min[2] = 50;
	param_max[2] = 350;
	
	param_min[3] = 50;
	param_max[3] = 350;
	*/
	
	for(int i = 2; i < ndim; i++)
	{
		param_min[i] = -20.0;
		param_max[i] = 20.0;
	}

	for(int i = 0; i < ndim; i++)
	{
		Cube[i] = param_min[i] + Cube[i]*(param_max[i] - param_min[i]);
		//std::cout << Cube[i] << std::endl;
	}
	int Nangbins = -1;
	
	if (Nangbins == 1)
	{
		Cube[2] = 7.0721303;
		Cube[3] = 2.20791021;
		Cube[4] = 0.11680584;
	}
	
	if (Nangbins == 2)
	{
		Cube[2] = 6.81751244;
		Cube[3] = 2.36025213;
		Cube[4] = 0.03740851;
		
		Cube[5] = 11.33756504;
		Cube[6] = 2.36025213;
		Cube[7] = 0.03740851;
		
	}

	if (Nangbins == 3)
	{
		Cube[2] = 6.507;
		Cube[3] = 2.4851;
		Cube[4] = 5.70509450e-03;
		
		Cube[5] = 9.07753874e+00;
		Cube[6] = 2.36025213e+00;
		Cube[7] = 3.74085115e-02;
		
		Cube[8] = 1.32870978e+01;
		Cube[9] = 2.48515295e+00;
		Cube[10] = 5.70509450e-03;
	}
	if (Nangbins == 4)
	{
		Cube[2] = 6.27549591e+00;
		Cube[3] = 2.53498230e+00;
		Cube[4] = 8.95737018e-04;
	
		Cube[5] = 7.98974799e+00;
		Cube[6] = 2.43056002e+00;
		Cube[7] = 1.59936997e-02;
	
		Cube[8] = 1.11859078e+01;
		Cube[9] = 2.43056002e+00;
		Cube[10] = 1.59936997e-02;
	
		Cube[11] = 1.39917083e+01;
		Cube[12] = 2.53498230e+00;
		Cube[13] = 8.95737018e-04;
	}
	
	loglikelihood(Cube, &npars, &lnew);
	//std::cout << std::endl;
	//std::cout <<" "<< npars << "\t" << lnew << std::endl;
	
	delete[] param_min;
	delete[] param_max;
	
}

/***********************************************************************************************************************/




/************************************************* dumper routine ******************************************************/

// The dumper routine will be called every updInt*10 iterations
// MultiNest doesn not need to the user to do anything. User can use the arguments in whichever way he/she wants
//
//
// Arguments:
//
// nSamples 						= total number of samples in posterior distribution
// nlive 						= total number of live points
// nPar 						= total number of parameters (free + derived)
// physLive[1][nlive * (nPar + 1)] 			= 2D array containing the last set of live points (physical parameters plus derived parameters) along with their loglikelihood values
// posterior[1][nSamples * (nPar + 2)] 			= posterior distribution containing nSamples points. Each sample has nPar parameters (physical + derived) along with the their loglike value & posterior probability
// paramConstr[1][4*nPar]:
// paramConstr[0][0] to paramConstr[0][nPar - 1] 	= mean values of the parameters
// paramConstr[0][nPar] to paramConstr[0][2*nPar - 1] 	= standard deviation of the parameters
// paramConstr[0][nPar*2] to paramConstr[0][3*nPar - 1] = best-fit (maxlike) parameters
// paramConstr[0][nPar*4] to paramConstr[0][4*nPar - 1] = MAP (maximum-a-posteriori) parameters
// maxLogLike						= maximum loglikelihood value
// logZ							= log evidence value from the default (non-INS) mode
// INSlogZ						= log evidence value from the INS mode
// logZerr						= error on log evidence value
// context						void pointer, any additional information

void dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr, double &maxLogLike, double &logZ, double &INSlogZ, double &logZerr, void *context)
{
	// convert the 2D Fortran arrays to C++ arrays
	
	
	// the posterior distribution
	// postdist will have nPar parameters in the first nPar columns & loglike value & the posterior probability in the last two columns
	
	/*
	int i, j;
	
	double postdist[nSamples][nPar + 2];
	for( i = 0; i < nPar + 2; i++ )
		for( j = 0; j < nSamples; j++ )
			postdist[j][i] = posterior[0][i * nSamples + j];
	
	
	
	// last set of live points
	// pLivePts will have nPar parameters in the first nPar columns & loglike value in the last column
	
	double pLivePts[nlive][nPar + 1];
	for( i = 0; i < nPar + 1; i++ )
		for( j = 0; j < nlive; j++ )
			pLivePts[j][i] = physLive[0][i * nlive + j];
	*/
}

/***********************************************************************************************************************/




/************************************************** Main program *******************************************************/



int main(int argc, char *argv[])
{
    std::ifstream pfile("MNparams.txt");
	//Read in the MultiNest sampling parameters

	int IS = read_param_int(&pfile, "IS");	// do Nested Importance Sampling?
	int mmodal = read_param_int(&pfile, "mmodal");	// do mode separation?
	int ceff = read_param_int(&pfile, "ceff");		// run in constant efficiency mode?
	
	int nlive = read_param_int(&pfile, "nlive");		// number of live points	
	double efr = read_param_double(&pfile, "efr");	// set the required efficiency
	double tol = read_param_double(&pfile, "tol");	// tol, defines the stopping criteria
	
	/*
	int ndims = read_param_int(&pfile, "ndims");		// dimensionality (no. of free parameters)
	int nPar = read_param_int(&pfile, "nPar");		// total no. of parameters including free & derived parameters
	int nClsPar = read_param_int(&pfile, "nClsPar");	// no. of parameters to do mode separation on
	*/
	
	
	int updInt = read_param_int(&pfile, "updInt");	// after how many iterations feedback is required & the output files should be updated
	// note: posterior files are updated & dumper routine is called after every updInt*10 iterations
	
	double Ztol = -1E90;				// all the modes with logZ < Ztol are ignored
	int maxModes = 100;				// expected max no. of modes (used only for memory allocation)
	
	char root[100] = "chains/data";	// root for output files
	
	int seed = -1;					// random no. generator seed, if < 0 then take the seed from system clock
	
	int fb = read_param_int(&pfile, "fb");		// need feedback on standard output?
	
	int resume = read_param_int(&pfile, "resume");    // resume from a previous job?
	
	int outfile = 1;				// write output files?
	
	int initMPI = 0;				// initialize MPI routines?, relevant only if compiling with MPI
							// set it to F if you want your main program to handle MPI initialization
	
	double logZero = -1E90;				// points with loglike < logZero will be ignored by MultiNest
	
	int maxiter = 0;				// max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it 
							// has done max no. of iterations or convergence criterion (defined through tol) has been satisfied
	
	void *context = 0;				// not required by MultiNest, any additional information user wants to pass

	
    std::ifstream fitfile("fitting.txt");
	//Read in the fitting parameters
	
	int ndims = 2;
	int nClsPar = ndims;
	int nPar = ndims;
	
	vmode = read_param_int(&fitfile, "vmode");
	if (vmode == 3)
	{
		N_ang = read_param_int(&fitfile, "N_ang");
		N_terms = read_param_int(&fitfile, "N_terms");
		ndims += N_ang*(N_terms-1); //Minus one for the fixed normalisation
		nPar = ndims;
		nClsPar = ndims;
		nPar += N_ang;
	}
	
	
	
	
	int pWrap[ndims];				// which parameters to have periodic boundary conditions?
	for(int i = 0; i < ndims; i++) pWrap[i] = 0;
	
	// calling MultiNest

	nested::run(IS, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI,
	logZero, maxiter, LogLike, dumper, context);
}

/***********************************************************************************************************************/
