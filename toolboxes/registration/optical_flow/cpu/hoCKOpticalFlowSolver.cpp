#include "hoCKOpticalFlowSolver.h"
#include "vector_td_utilities.h"

#ifdef USE_OMP
#include <omp.h>
#endif

namespace Gadgetron {

// Helpers
//

template <unsigned int D>
inline bool is_border_pixel_for_stride(typename int64d<D>::Type stride, typename uint64d<D>::Type co,
                                       typename uint64d<D>::Type dims) {
    for (size_t d = 0; d < D; d++) {
        if (stride.vec[d] == -1) {
            if (co.vec[d] == 0) {
                return true;
            }
        } else if (stride.vec[d] == 1) {
            if (co.vec[d] == (dims.vec[d] - 1)) {
                return true;
            }
        }
    }
    return false;
}

template <size_t i, size_t j> struct Pow {
    enum { Value = i * Pow<i, j - 1>::Value };
};

template <size_t i> struct Pow<i, 1> {
    enum { Value = i };
};

//
// Implementation
//

template <class T, unsigned int D>
boost::shared_ptr<hoNDArray<T>> hoCKOpticalFlowSolver<T, D>::core_solver(hoNDArray<T>* _gradient_image,
                                                                         hoNDArray<T>* _stencil_image) {
    // Sanity checks
    //

    if (!_gradient_image) {
        throw std::runtime_error("hoCKOpticalFlowSolver::core_solver(): illegal input gradient image received.");
    }

    if (_gradient_image->get_number_of_dimensions() <= D) {
        throw std::runtime_error(
            "hoCKOpticalFlowSolver::core_solver(): number of gradient image dimensions is too small.");
    }

    // The dimensions of the displacement field should match the gradient field
    //

    std::vector<size_t> disp_dims = _gradient_image->get_dimensions();
    boost::shared_ptr<hoNDArray<T>> displacements_ping(new hoNDArray<T>(disp_dims));
    boost::shared_ptr<hoNDArray<T>> displacements_pong(new hoNDArray<T>(disp_dims));
    clear(displacements_ping.get());
    clear(displacements_pong.get());

    // We use "shared memory" to hold the averaged displacements
    boost::shared_ptr<hoNDArray<T>> _shared_mem(new hoNDArray<T>(disp_dims));
    T* shared_mem = _shared_mem->get_data_ptr();
    clear(_shared_mem.get());

    typename uint64d<D>::Type matrix_size = from_std_vector<size_t, D>(disp_dims);
    size_t number_of_elements = prod(matrix_size);
    size_t num_batches = 1;

    for (size_t d = D; d < _gradient_image->get_number_of_dimensions() - 1; d++) {
        num_batches *= _gradient_image->get_size(d);
    }

    // Get ready
    //

    size_t iteration_no = 0;
    hoNDArray<T>* ping = displacements_ping.get();
    hoNDArray<T>* pong = displacements_pong.get();

    if (this->output_mode_ >= hoOpticalFlowSolver<T, D>::OUTPUT_VERBOSE) {
        GDEBUG_STREAM(std::endl);
    }

    //
    // Main Jacobi loop
    //

    while (true) {

        if (this->output_mode_ >= hoOpticalFlowSolver<T, D>::OUTPUT_VERBOSE) {
            GDEBUG_STREAM("."; std::cerr.flush());
        }

        // Continuation flag used for early Jacobi termination
        size_t continue_flag = 0;

        // Number of elements per batch
        const size_t num_elements_per_batch = prod(matrix_size);

        // Number of elements per dim
        const size_t num_elements_per_dim = num_elements_per_batch * num_batches;

        T* in_disp = ping->get_data_ptr();
        T* out_disp = pong->get_data_ptr();
        T* gradient_image = _gradient_image->get_data_ptr();
        T* stencil_image = (_stencil_image) ? _stencil_image->get_data_ptr() : 0x0;

        //
        // Find the average velocities (shared memory)
        //

        for (size_t dim = 0; dim < D + 1; dim++) {
#ifdef USE_OMP
#pragma omp parallel for
#endif
            for (long long idx = 0; idx < (long long)num_elements_per_dim; idx++) {

                // Index to the shared memory
                const size_t shared_idx = dim * num_elements_per_dim + idx;

                // Batch idx (second slowest varying dimension)
                const size_t batch_idx = idx / num_elements_per_batch;

                // Local index to the image (or batch in our terminology)
                const size_t idx_in_batch = idx - batch_idx * num_elements_per_batch;

                if (stencil_image && stencil_image[idx_in_batch] > T(0))
                    continue;

                // Local co to the image
                const typename uint64d<D>::Type co = idx_to_co(idx_in_batch, matrix_size);
                const typename int64d<D>::Type zeros(0);
                const typename int64d<D>::Type ones(1);
                const typename int64d<D>::Type threes(3);

                const int num_neighbors = Pow<3, D>::Value;
                T num_contribs = T(0);

                shared_mem[shared_idx] = T(0);

                // Compute average of neighbors
                //

                for (long long i = 0; i < num_neighbors; i++) {

                    // Find the stride of the neighbor {-1, 0, 1}^D
                    const typename int64d<D>::Type stride = idx_to_co(i, threes) - ones;

                    size_t neighbor_idx;

                    const size_t base_offset = dim * num_elements_per_dim + batch_idx * num_elements_per_batch;

                    // Verify that the neighbor is not out of bounds (and not the thread itself)
                    if (!is_border_pixel_for_stride<D>(stride, co, matrix_size) && !(stride == zeros)) {
                        neighbor_idx = (size_t)co_to_idx(vector_td<long long, D>(co) + stride,
                                                         vector_td<long long, D>(matrix_size)) +
                                       base_offset;
                    } else {
                        neighbor_idx = idx_in_batch + base_offset;
                    }

                    shared_mem[shared_idx] += in_disp[neighbor_idx];
                    num_contribs += T(1);
                }

                // Normalize
                shared_mem[shared_idx] /= num_contribs;
            }
        }

        //
        // Update displacement field (Jacobi iteration)
        //

        const T disp_thresh_sqr = this->limit_ * this->limit_;

        for (size_t dim = 0; dim < D + 1; dim++) {
#ifdef USE_OMP
#pragma omp parallel for
#endif
            for (long long idx = 0; idx < num_elements_per_dim; idx++) {
                // Index to the shared memory
                const size_t shared_idx = dim * num_elements_per_dim + idx;

                // Batch idx (second slowest varying dimension)
                const size_t batch_idx = idx / num_elements_per_batch;

                // Local index to the image (or batch in our terminology)
                const size_t idx_in_batch = idx - batch_idx * num_elements_per_batch;

                if (stencil_image && stencil_image[idx_in_batch] > T(0))
                    continue;

                T phi = T(0);
                T norm = T(0);

                typename reald<T, D>::Type derivatives;

                // Contributions from the spatial dimensions
                //

                for (size_t d = 0; d < D; d++) {
                    derivatives.vec[d] = gradient_image[d * num_elements_per_dim + idx];
                    const size_t shared_idx_d = d * num_elements_per_dim + idx;
                    phi += (shared_mem[shared_idx_d] * derivatives.vec[d]);
                    norm += (derivatives.vec[d] * derivatives.vec[d]);
                }

                // Contributions from the temporal dimension
                //

                phi += gradient_image[D * num_elements_per_dim + idx];

                // Contribution from the intensity attentuation estimation
                //

                phi -= shared_mem[D * num_elements_per_dim + idx];

                // Normalize
                //

                phi /= ((alpha_ / beta_) * (alpha_ / beta_) + alpha_ * alpha_ + norm);

                // Form result displacement
                //

                T result;

                if (dim < D)
                    result = shared_mem[shared_idx] - derivatives.vec[dim] * phi;
                else
                    result = shared_mem[D * num_elements_per_dim + idx] + (alpha_ / beta_) * (alpha_ / beta_) * phi;

                // Clear the "termination" flag if the displacement field has changed above the threshold
                //

                T delta = result - in_disp[dim * num_elements_per_dim + idx];
                if (dim < D && delta * delta > disp_thresh_sqr)
                    continue_flag = 1;

                // Output result
                //

                out_disp[dim * num_elements_per_dim + idx] = result;
            }
        }

        // Swap in/out buffers
        //

        hoNDArray<T>* tmp = ping;
        ping = pong;
        pong = tmp;

        // Check termination criteria
        //

        if (continue_flag == 0) {
            if (this->output_mode_ >= hoOpticalFlowSolver<T, D>::OUTPUT_VERBOSE) {
                GDEBUG_STREAM(std::endl << "Break after " << iteration_no + 1 << " iterations" << std::endl);
            }
            break;
        }

        if (iteration_no > this->max_num_iterations_per_level_)
            break;

        iteration_no++;
    }

    if (ping == displacements_ping.get())
        return displacements_ping;
    else
        return displacements_pong;
}

//
// Template instantiation
//

template class hoCKOpticalFlowSolver<float, 1>;
template class hoCKOpticalFlowSolver<float, 2>;
template class hoCKOpticalFlowSolver<float, 3>;
template class hoCKOpticalFlowSolver<float, 4>;

template class hoCKOpticalFlowSolver<double, 1>;
template class hoCKOpticalFlowSolver<double, 2>;
template class hoCKOpticalFlowSolver<double, 3>;
template class hoCKOpticalFlowSolver<double, 4>;
} // namespace Gadgetron
