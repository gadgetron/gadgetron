#include "fatwater.h"

#include "hoMatrix.h"
#include "hoNDArray_linalg.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_reductions.h"
#include "hoArmadillo.h"

#include <boost/config.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/read_dimacs.hpp>
#include <boost/graph/graph_utility.hpp>


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
	std::pair<float,float> range_fm = std::make_pair(-80.0,80.0);
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
	    phiMatrix(k1,k2) = 0.0;
	    npeaks = a.species_[k2].ampFreq_.size();
	    for( int k3=0;k3<npeaks;k3++) {
	      relAmp = a.species_[k2].ampFreq_[k3].first;
	      freq_hz = fieldStrength*GAMMABAR*a.species_[k2].ampFreq_[k3].second;
	      phiMatrix(k1,k2) += relAmp*std::complex<float>(cos(2*PI*echoTimes[k1]*freq_hz),sin(2*PI*echoTimes[k1]*freq_hz));
	    }
	    GDEBUG("Cur value phiMatrix = (%f,%f) \n", phiMatrix(k1,k2).real(), phiMatrix(k1,k2).imag());
	  }
	}
	//auto a_phiMatrix = as_arma_matrix(&phiMatrix);
	//auto mymat2 = mymat.t()*mymat;

	for(int ka=0;ka<phiMatrix.get_size(0);ka++) {
	  for(int kb=0;kb<phiMatrix.get_size(1);kb++) {
	    GDEBUG("Check phiMatrix(%d,%d) = %f + i %f \n", ka,kb,phiMatrix(ka,kb).real(),phiMatrix(ka,kb).imag());
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
	//	auto a_phiMatrix = as_arma_matrix(&IdentMat);

	float fm;
	std::vector<float> fms(num_fm);
	fms[0] = range_fm.first;
	for(int k1=1;k1<num_fm;k1++) {
	  fms[k1] = range_fm.first + k1*(range_fm.second-range_fm.first)/(num_fm-1);
	}

	float r2star;
    std::vector<float> r2stars(num_r2star);
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

	    herk( tempM1, psiMatrix, 'L', true );
	    //	    tempM1.copyLowerTriToUpper();
	    for (int ka=0;ka<tempM1.get_size(0);ka++ ) {
	      for (int kb=ka+1;kb<tempM1.get_size(1);kb++ ) {
		tempM1(ka,kb) = conj(tempM1(kb,ka));
	      }
	    }

	    if(k3==50) {
	      for(int ka=0;ka<tempM1.get_size(0);ka++) {
		for(int kb=0;kb<tempM1.get_size(1);kb++) {
		  GDEBUG(" tempM1(%d,%d) = %f + i %f \n", ka,kb,tempM1(ka,kb).real(),tempM1(ka,kb).imag());
		}
	      }
	    }


	    potri(tempM1);
	    for (int ka=0;ka<tempM1.get_size(0);ka++ ) {
	      for (int kb=ka+1;kb<tempM1.get_size(1);kb++ ) {
		tempM1(ka,kb) = conj(tempM1(kb,ka));
	      }
	    }

	    if(k3==50) {
	      for(int ka=0;ka<tempM1.get_size(0);ka++) {
		for(int kb=0;kb<tempM1.get_size(1);kb++) {
		  GDEBUG(" inv tempM1(%d,%d) = %f + i %f \n", ka,kb,tempM1(ka,kb).real(),tempM1(ka,kb).imag());
		}
	      }
	    }


	    //GDEBUG(" (%d,%d) = (%d,%d) X (%d,%d) \n", tempM2.get_size(0),tempM2.get_size(1),tempM1.get_size(0),tempM1.get_size(1),psiMatrix.get_size(1),psiMatrix.get_size(0));
	    gemm( tempM2, tempM1, false, psiMatrix, true );

	    if(k3==50) {
	      for(int ka=0;ka<tempM2.get_size(0);ka++) {
		for(int kb=0;kb<tempM2.get_size(1);kb++) {
		  GDEBUG(" tempM2(%d,%d) = %f + i %f \n", ka,kb,tempM2(ka,kb).real(),tempM2(ka,kb).imag());
		}
	      }
	    }

	    //GDEBUG(" (%d,%d) = (%d,%d) X (%d,%d) \n", P1.get_size(0),P1.get_size(1),psiMatrix.get_size(0),psiMatrix.get_size(1),tempM2.get_size(0),tempM2.get_size(1));
	    gemm( P1, psiMatrix, false, tempM2, false );

	    if(k3==50) {
	      for(int ka=0;ka<P1.get_size(0);ka++) {
		for(int kb=0;kb<P1.get_size(1);kb++) {
		  GDEBUG(" P1(%d,%d) = %f + i %f \n", ka,kb,P1(ka,kb).real(),P1(ka,kb).imag());
		}
	      }
	    }


	    subtract(IdentMat,P1,P);

	    if(k3==50) {
	      for(int ka=0;ka<P.get_size(0);ka++) {
		for(int kb=0;kb<P.get_size(1);kb++) {
		  GDEBUG(" P(%d,%d) = %f + i %f \n", ka,kb,P(ka,kb).real(),P(ka,kb).imag());
		}
	      }

	    }





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
		if (k1==107 && k2==144) {
		  tempSignal(k5,k4) = std::complex<float>(1000.0,0.0);;
		  GDEBUG(" (%d,%d) -->  %f + i %f \n",k5,k4, tempSignal(k5,k4).real(),tempSignal(k5,k4).imag());
		}

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

	      if (k1==107 && k2==144) {
		GDEBUG(" %f -->  %f \n",fms[k3],minResidual);
	      }
	    }
	  }
	}



	//arma::Mat< std::complex<float> > arma_phiMatrix = as_arma_matrix( phiMatrix );


	//Do graph cut iterations
	using namespace boost;

	// create a typedef for the Graph type
	typedef adjacency_list<vecS, vecS, bidirectionalS> Graph;

	// Make convenient labels for the vertices
	enum { A, B, C, D, E };
	const int num_vertices = 5;
	const char* name = "ABCDE";

	// writing out the edges in the graph
	typedef std::pair<int, int> Edge;
	Edge edge_array[] =
	  { Edge(A,B), Edge(A,D), Edge(C,A), Edge(D,C),
	    Edge(C,E), Edge(B,D), Edge(D,E) };
	const int num_edges = sizeof(edge_array)/sizeof(edge_array[0]);

	// declare a graph object
	Graph g(edge_array, edge_array + sizeof(edge_array) / sizeof(Edge), num_vertices);


	/*
	property_map<Graph, edge_capacity_t>::type
	  capacity = get(edge_capacity, g);
	property_map<Graph, edge_reverse_t>::type
	  rev = get(edge_reverse, g);
	property_map<Graph, edge_residual_capacity_t>::type
	  residual_capacity = get(edge_residual_capacity, g);

	Traits::vertex_descriptor s, t;
	read_dimacs_max_flow(g, capacity, rev, s, t);

	long flow;
#if defined(BOOST_MSVC) && BOOST_MSVC <= 1300
	// Use non-named parameter version
	property_map<Graph, vertex_index_t>::type
	  indexmap = get(vertex_index, g);
	flow = push_relabel_max_flow(g, s, t, capacity, residual_capacity, rev, indexmap);
#else
	flow = push_relabel_max_flow(g, s, t);
#endif

	std::cout << "c  The total flow:" << std::endl;
	std::cout << "s " << flow << std::endl << std::endl;

	std::cout << "c flow values:" << std::endl;
	graph_traits<Graph>::vertex_iterator u_iter, u_end;
	graph_traits<Graph>::out_edge_iterator ei, e_end;
	for (boost::tie(u_iter, u_end) = vertices(g); u_iter != u_end; ++u_iter)
	  for (boost::tie(ei, e_end) = out_edges(*u_iter, g); ei != e_end; ++ei)
	    if (capacity[*ei] > 0)
	      std::cout << "f " << *u_iter << " " << target(*ei, g) << " "
			<< (capacity[*ei] - residual_capacity[*ei]) << std::endl;

	*/





	//Do final calculations once the field map is done
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
	    herk( AhA, psiMatrix, 'L', true );
	    //	    AhA.copyLowerTriToUpper();
	    for (int ka=0;ka<AhA.get_size(0);ka++ ) {
	      for (int kb=ka+1;kb<AhA.get_size(1);kb++ ) {
		AhA(ka,kb) = conj(AhA(kb,ka));
	      }
	    }

	    hesv(AhA,curWaterFat);
	    for ( int k4=0;k4<N;k4++ ) {
	      for ( int k5=0;k5<2;k5++ ) { // 2 elements for water and fat currently
		out(k1,k2,0,0,k4,k5,0) = curWaterFat(k5,k4);
	      }
	    }

	  }
	}



	//Clean up as needed


        return out;
    }
}
