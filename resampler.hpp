#ifndef RESAMPLER_H_
#define RESAMPLER_H_

#include <algorithm>
#include <cassert>

#include "filter_design.hpp"

typedef float  complex_float[2];
typedef double complex_double[2];

template<typename T_complex>
inline void
complex_set_zero(T_complex out)
{
   out[0] = 0.0f;
   out[1] = 0.0f;
}

template<typename T_complex>
inline void
complex_copy(const T_complex src, T_complex dst)
{
   dst[0] = src[0];
   dst[1] = src[1];
}

template<typename T_complex, typename T>
inline void
complex_accumulate(T_complex accumulator,
                   const T_complex in,
                   T scalar)
{
   accumulator[0] += in[0] * scalar;
   accumulator[1] += in[1] * scalar;
}

// polyphase_resample
// Performs polyphase resampling on an input array of T_data samples
// using a full FIR filter. The function takes the input and desired output
// sampling rates (of type T), computes integer conversion factors, and uses
// the provided full FIR filter (designed via design_raised_cosine_filter) for
// filtering. Parameters:
//   input               : pointer to input complex signal (T_data*)
//   input_size          : number of input samples
//   input_sampling_rate : sampling rate of the input signal (Hz)
//   output_sampling_rate: desired output sampling rate (Hz)
//   full_filter         : full FIR filter coefficients
//   output_size         : (output) number of output samples computed
// Returns a dynamically allocated array of T_data samples (caller must
// delete[]).
template<typename T_data, typename T>
inline T_data*
polyphase_resample(const T_data*  input,
                   size_t                input_size,
                   T                     input_sampling_rate,
                   T                     output_sampling_rate,
                   const std::vector<T>& full_filter,
                   size_t&               output_size)
{
   // Compute integer conversion factors.
   // T local_epsilon = get_epsilon_for_type<T>();  // Only needed for
   // branch_sum
   up_down_factors conv = compute_conversion_factors<T>(input_sampling_rate,
                                                        output_sampling_rate);
   int             up_factor   = conv.up_factor;
   int             down_factor = conv.down_factor;

   int num_taps     = full_filter.size();
   int L            = up_factor;   // number of polyphase subfilters
   int phase_length = (num_taps + L - 1) / L;
   std::vector<std::vector<T>> poly_filters;
   poly_filters.resize(L);
   for(int p = 0; p < L; ++p)
   {
      poly_filters[p].resize(phase_length, T(0));
      for(int k = 0; k < phase_length; ++k)
      {
         int index = p + k * L;
         if(index < num_taps) { poly_filters[p][k] = full_filter[index]; }
         else { poly_filters[p][k] = T(0); }
      }
      /*
      // Normalize each phase so that the sum of coefficients is unity.
      T branch_sum = T(0);
      for (int k = 0; k < phase_length; ++k) {
         branch_sum += poly_filters[p][k];
      }
      if (std::abs(branch_sum) > local_epsilon) {
         for (int k = 0; k < phase_length; ++k) {
            poly_filters[p][k] /= branch_sum;
         }
      }
      */
   }

   // Group delay of the filter.
   int delay       = (num_taps - 1) / 2;
   int input_index = delay;
   int phase       = 0;

   size_t max_output_samples
      = (((input_size - delay) * up_factor) / down_factor) + 1;
   T_data* output_buffer = new T_data[max_output_samples];
   size_t         out_count     = 0;

   while(input_index < static_cast<int>(input_size))
   {
      T_data y;
      complex_set_zero<T_data>(y);
      for(int k = 0; k < phase_length; ++k)
      {
         int idx = input_index - k;
         if(idx < 0) { break; }
         complex_accumulate<T_data, T>(y,
                            input[idx],
                            static_cast<float>(poly_filters[phase][k]));
      }
      complex_copy<T_data>(y, output_buffer[out_count]);
      out_count++;

      phase += down_factor;
      int advance = phase / up_factor;
      phase       = phase % up_factor;
      input_index += advance;
   }
   output_size = out_count;
   return output_buffer;
}

#endif /* RESAMPLER_H_ */
