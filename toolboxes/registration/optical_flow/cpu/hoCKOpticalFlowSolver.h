/** \file hoCKOpticalFlowSolver.h
    \brief CPU-based Cornelius-Kanade optical flow registration solver.

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

#include "hoOpticalFlowSolver.h"

namespace Gadgetron{

  template<class T, unsigned int D> class hoCKOpticalFlowSolver
    : public hoOpticalFlowSolver<T, D>
  {

  public:

    // Constructors / destructors
    //

    hoCKOpticalFlowSolver() : hoOpticalFlowSolver<T,D>(){
      alpha_ = T(0.05);
      beta_ = T(1.0);
    }

    virtual ~hoCKOpticalFlowSolver() {}

    // Set the regularization weight
    //

    inline void set_alpha( T alpha ) { alpha_ = alpha; }
    inline void set_beta( T beta ) { beta_ = beta; }

  protected:
    virtual boost::shared_ptr< hoNDArray<T> >
      core_solver( hoNDArray<T> *gradient_image, hoNDArray<T> *stencil );

  protected:
    T alpha_;
    T beta_;
  };
}
