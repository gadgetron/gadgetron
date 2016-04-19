#include "fatwater.h"

#include "hoMatrix.h"
#include "hoNDArray_linalg.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_reductions.h"


#define GAMMABAR 42.576 // MHz/T
#define PI 3.141592


namespace Gadgetron {
    hoNDArray< std::complex<float> > fatwater_separation(hoNDArray< std::complex<float> >& data, FatWaterParameters p, FatWaterAlgorithm a)
    {
		
	//Get some data parameters
	//7D, fixed order [X, Y, Z, CHA, N, S, LOC]
        uint16_t X = data.get_size(0);
        uint16_t Y = data.get_size(1);
        uint16_t Z = data.get_size(2);
        uint16_t CHA = data.get_size(3);
        uint16_t N = data.get_size(4);
        uint16_t S = data.get_size(5);
        uint16_t LOC = data.get_size(6);
	
	GDEBUG("Size of my array: %d, %d, %d .\n", X,Y,Z);

	hoNDArray< std::complex<float> > out(X,Y,Z,CHA,N,2,LOC); // S dimension gets replaced by water/fat stuff

	float fieldStrength = p.fieldStrengthT_;
        std::vector<float> echoTimes = p.echoTimes_;
	bool precessionIsClockwise = p.precessionIsClockwise_;
        for (auto& te: echoTimes) {
          te = te*0.001; // Echo times in seconds rather than milliseconds
        }

	GDEBUG("In toolbox - Field Strength: %f T \n", fieldStrength);
        for (auto& te: echoTimes) {
	  GDEBUG("In toolbox - Echo time: %f seconds \n", te);
        }
	GDEBUG("In toolbox - PrecessionIsClockwise: %d \n", precessionIsClockwise);


	//Get or set some algorithm parameters
	//Gadgetron::ChemicalSpecies w = a.species_[0];
	//Gadgetron::ChemicalSpecies f = a.species_[1];

	//	GDEBUG("In toolbox - Fat peaks: %f  \n", f.ampFreq_[0].first);
	//	GDEBUG("In toolbox - Fat peaks 2: %f  \n", f.ampFreq_[0].second);
	
	// Set some initial parameters so we can get going
	// These will have to be specified in the XML file eventually
	std::pair<float,float> range_r2star = std::make_pair(0.0,0.0);
	uint16_t num_r2star = 1;
	std::pair<float,float> range_fm = std::make_pair(-200.0,200.0);
	uint16_t num_fm = 101;
	uint16_t num_iterations = 40;
	uint16_t subsample = 1;
	float lmap_power = 2.0;
	float lambda = 0.02;
	float lambda_extra = 0.02;

	//Check that we have reasonable data for fat-water separation

	
	//Calculate residual
	//
	float relAmp, freq_hz;
	uint16_t npeaks;
	uint16_t nspecies = a.species_.size();
	uint16_t nte = echoTimes.size();
	GDEBUG("In toolbox - NTE: %d \n", nte);

	hoMatrix< std::complex<float> > phiMatrix(nte,nspecies);
	for( int k1=0;k1<nte;k1++) {
	  for( int k2=0;k2<nspecies;k2++) {
	    phiMatrix[k1,k2] = std::complex<float>(0.0,0.0);
	    npeaks = a.species_[k2].ampFreq_.size();
	    for( int k3=0;k3<npeaks;k3++) {
	      relAmp = a.species_[k2].ampFreq_[k3].first;
	      freq_hz = fieldStrength*GAMMABAR*a.species_[k2].ampFreq_[k3].second;
	      phiMatrix(k1,k2) += relAmp*std::complex<float>(cos(2*PI*echoTimes[k1]*freq_hz),sin(2*PI*echoTimes[k1]*freq_hz));
	    }
	    GDEBUG("Cur value phiMatrix = (%f,%f) \n", phiMatrix(k1,k2).real(), phiMatrix(k1,k2).imag());
	  }
	}
	
	hoMatrix< std::complex<float> > IdentMat(nte,nte);
	for( int k1=0;k1<nte;k1++) {
	  for( int k2=0;k2<nte;k2++) {
	    if( k1==k2 ) {
	      IdentMat(k1,k2) = std::complex<float>(1.0,0.0);
	    } else {
	      IdentMat(k1,k2) = std::complex<float>(0.0,0.0);
	    }
	  }
	}

	float fm;
	float fms[num_fm];
	fms[0] = range_fm.first;
	for(int k1=1;k1<num_fm;k1++) {
	  fms[k1] = range_fm.first + k1*(range_fm.second-range_fm.first)/(num_fm-1);
	}

	float r2star;
	float r2stars[num_r2star];
	r2stars[0] = range_r2star.first;
	for(int k2=1;k2<num_r2star;k2++) {
	  r2stars[k2] = range_r2star.first + k2*(range_r2star.second-range_r2star.first)/(num_r2star-1);
	}
	
	
	std::complex<float> curModulation;
	hoMatrix< std::complex<float> > tempM1(nspecies,nspecies);
	hoMatrix< std::complex<float> > tempM2(nspecies,nte);
	hoMatrix< std::complex<float> > psiMatrix(nte,nspecies);
	hoNDArray< std::complex<float> > Ps(nte,nte,num_fm,num_r2star);
	hoMatrix< std::complex<float> > P1(nte,nte);
	hoMatrix< std::complex<float> > P(nte,nte);
	
	for(int k3=0;k3<num_fm;k3++) {
	  fm = fms[k3];
	  for(int k4=0;k4<num_r2star;k4++) {
	    r2star = r2stars[k4];

	    for( int k1=0;k1<nte;k1++) {
	      curModulation = exp(-r2star*echoTimes[k1])*std::complex<float>(cos(2*PI*echoTimes[k1]*fm),sin(2*PI*echoTimes[k1]*fm));
	      for( int k2=0;k2<nspecies;k2++) {
		psiMatrix(k1,k2) = phiMatrix(k1,k2)*curModulation;
	      }
	    }
	    
	    gemm( tempM1, psiMatrix, true, psiMatrix, false );
	    
	    potri(tempM1);

	    gemm( tempM2, tempM1, false, psiMatrix, true );

	    gemm( P1, psiMatrix, false, tempM2, true );
	    
	    subtract(IdentMat,P1,P);
	    
	    
	    // Keep all projector matrices together
	    for( int k1=0;k1<nte;k1++) {
	      for( int k2=0;k2<nte;k2++) {
		Ps(k1,k2,k3,k4) = P(k1,k2);
	      }
	    }	    
	  }
	}
	
	
	// Need to check that S = nte
	// N should be the number of contrasts (eg: for PSIR)
	hoMatrix< std::complex<float> > tempResVector(S,N);
	hoMatrix< std::complex<float> > tempSignal(S,N);
	hoNDArray<float> residual(num_fm,X,Y);
	hoNDArray<uint16_t> r2starIndex(X,Y,num_fm);
	hoNDArray<uint16_t> fmIndex(X,Y);
	float curResidual, minResidual, minResidual2;
	for( int k1=0;k1<X;k1++) {
	  for( int k2=0;k2<Y;k2++) {
	    // Get current signal
	    for( int k4=0;k4<N;k4++) {
	      for( int k5=0;k5<S;k5++) {
		tempSignal(k5,k4) = data(k1,k2,0,0,k4,k5,0); 
	      }
	    }

	    minResidual2 = 1.0 + nrm2(&tempSignal);
	    
	    for(int k3=0;k3<num_fm;k3++) {

	      minResidual = 1.0 + nrm2(&tempSignal);

	      for(int k4=0;k4<num_r2star;k4++) {
		// Get current projector matrix
		for( int k5=0;k5<nte;k5++) {
		  for( int k6=0;k6<nte;k6++) {
		    P(k5,k6) = Ps(k5,k6,k3,k4);
		  }
		}	    
		
		// Apply projector
		gemm( tempResVector, P, false, tempSignal, false );

		curResidual = nrm2(&tempResVector);

		if (curResidual < minResidual) {
		  minResidual = curResidual;
		  r2starIndex(k1,k2,k3) = k4;
		}
	      }
	      residual(k3,k1,k2) = minResidual;
	      
	      if (minResidual < minResidual2) {
		minResidual2 = minResidual;
		fmIndex(k1,k2) = k3;
	      }
	      
	    }
	  }	
	}
	
	
	hoMatrix< std::complex<float> > curWaterFat(2,N);
	hoMatrix< std::complex<float> > AhA(2,2);
	// Do fat-water separation with current field map and R2* estimates
	for( int k1=0;k1<X;k1++) {
	  for( int k2=0;k2<Y;k2++) {
	    
	    // Get current signal
	    for( int k4=0;k4<N;k4++) {
	      for( int k5=0;k5<S;k5++) {
		tempSignal(k5,k4) = data(k1,k2,0,0,k4,k5,0); 
	      }
	    }
	    // Get current Psi matrix
	    fm = fms[fmIndex(k1,k2)];
	    r2star = r2stars[r2starIndex(k1,k2,fmIndex(k1,k2))];
	    for( int k3=0;k3<nte;k3++) {
	      curModulation = exp(-r2star*echoTimes[k3])*std::complex<float>(cos(2*PI*echoTimes[k3]*fm),sin(2*PI*echoTimes[k3]*fm));
	      for( int k4=0;k4<nspecies;k4++) {
		psiMatrix(k3,k4) = phiMatrix(k3,k4)*curModulation;
	      }
	    }
	    
	    // Solve for water and fat
	    gemm( curWaterFat, psiMatrix, true, tempSignal, false );
	    gemm( AhA, psiMatrix, true, psiMatrix, false );
	    hesv(AhA,curWaterFat);
	    for ( int k4=0;k4<N;k4++ ) {
	      for ( int k5=0;k5<2;k5++ ) { // 2 elements for water and fat currently
		out(k1,k2,0,0,k4,k5,0) = curWaterFat(k5,k4);
	      }
	    }
	    
	  }
	}

	
	//arma::Mat< std::complex<float> > arma_phiMatrix = as_arma_matrix( phiMatrix ); 

	
	//Do graph cut iterations

	
	//Do final calculations once the field map is done
	
	
	//Clean up as needed

	
        return out;
    }
}
