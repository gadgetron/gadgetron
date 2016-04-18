#ifndef FATWATER_H
#define FATWATER_H

#include <vector>
#include <utility>

#include "fatwater_export.h"
#include "hoNDArray.h"

namespace Gadgetron
{

    struct FatWaterParameters
    {
        FatWaterParameters()
        : fieldStrengthT_(0)
        , precessionIsClockwise_(true)
        {}
        
        float fieldStrengthT_;
        std::vector<float> echoTimes_;
        bool precessionIsClockwise_;
    };


    struct ChemicalSpecies
    {
        std::vector< std::pair<float,float> > ampFreq_;
    };
    
    struct FatWaterAlgorithm
    {
        std::vector<ChemicalSpecies> species_;
    };
    
    EXPORTFATWATER hoNDArray< std::complex<float> > fatwater_separation(hoNDArray< std::complex<float> >& data, FatWaterParameters p, FatWaterAlgorithm a);

}

#endif //FATWATER_H
