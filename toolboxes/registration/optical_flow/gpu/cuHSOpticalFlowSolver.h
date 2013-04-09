/** \file cuHSOpticalFlowSolver.h
    \brief GPU-based Horn-Schunck optical flow registration solver.

    References to the solver implementation and some usage scenarios can be found in:

    An optimised multi-baseline approach for on-line MR-temperature monitoring on commodity graphics hardware
    BD de Senneville, KØ Noe, M Ries, M Pedersen, CTW Moonen, TS Sørensen.
    5th IEEE International Symposium on Biomedical Imaging: From Nano to Macro, 2008. ISBI 2008. pp. 1513-1516.

    Acceleration and validation of optical flow based deformable registration for image-guided radiotherapy.
    KØ Noe, BD de Senneville, UV Elstrøm, K Tanderup, TS Sørensen.
    Acta Oncologica 2008; 47(7): 1286-1293.

    Retrospective reconstruction of high temporal resolution cine images from real‐time MRI using iterative motion correction
    MS Hansen, TS Sørensen, AE Arai, P Kellman.
    Magnetic Resonance in Medicine 2012; 68(3): 741-750.
*/

#pragma once

#include "cuOpticalFlowSolver.h"

namespace Gadgetron{

  template<class REAL, unsigned int D> class EXPORTGPUREG cuHSOpticalFlowSolver 
    : public cuOpticalFlowSolver<REAL, D>
  {
  
  public:

    // Constructors / destructors
    //
  
    cuHSOpticalFlowSolver() : cuOpticalFlowSolver<REAL,D>(){ 
      alpha_ = REAL(0.1); 
    } 
  
    virtual ~cuHSOpticalFlowSolver() {}
  
    // Set the regularization weight
    //
  
    inline void set_alpha( REAL alpha ) { alpha_ = alpha; }
  
  protected:  
    virtual boost::shared_ptr< cuNDArray<REAL> > 
      core_solver( cuNDArray<REAL> *gradient_image, cuNDArray<REAL> *stencil_image );
    
  protected:
    REAL alpha_;
  };
}
