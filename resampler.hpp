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
complex_accumulate(T_complex accumulator, const T_complex in, T scalar)
{
   accumulator[0] += in[0] * scalar;
   accumulator[1] += in[1] * scalar;
}

// eases circular buffer emulation
inline const uint32_t
circular_index(size_t input_size, uint32_t index)
{
   uint32_t mod_index = index % static_cast<uint32_t>(input_size);
   if(mod_index < 0) { mod_index += static_cast<uint32_t>(input_size); }
   return mod_index;
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
polyphase_resample(const T_data*         input,
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
   uint32_t        up_factor   = conv.up_factor;
   uint32_t        down_factor = conv.down_factor;

   uint32_t num_taps     = full_filter.size();
   uint32_t L            = up_factor;   // number of polyphase subfilters
   uint32_t phase_length = (num_taps + L - 1) / L;
   std::vector<std::vector<T>> poly_filters;
   poly_filters.resize(L);
   for(uint32_t p = 0; p < L; ++p)
   {
      poly_filters[p].resize(phase_length, T(0));
      for(uint32_t k = 0; k < phase_length; ++k)
      {
         uint32_t index = p + k * L;
         if(index < num_taps) { poly_filters[p][k] = full_filter[index]; }
         else { poly_filters[p][k] = T(0); }
      }
      /*
      // Normalize each phase so that the sum of coefficients is unity.
      T branch_sum = T(0);
      for (uint32_t k = 0; k < phase_length; ++k) {
         branch_sum += poly_filters[p][k];
      }
      if (std::abs(branch_sum) > local_epsilon) {
         for (uint32_t k = 0; k < phase_length; ++k) {
            poly_filters[p][k] /= branch_sum;
         }
      }
      */
   }

   // Group delay of the filter.
   uint32_t delay       = (num_taps - 1) / 2;
   uint32_t input_index = delay;
   uint32_t phase       = 0;

   size_t max_output_samples
      = (((input_size - delay) * up_factor) / down_factor) + 1;
   T_data* output_buffer = new T_data[max_output_samples];
   size_t  out_count     = 0;

   while(input_index < static_cast<uint32_t>(input_size))
   {
      T_data y;
      complex_set_zero<T_data>(y);
      for(uint32_t k = 0; k < phase_length; ++k)
      {
         uint32_t idx = input_index - k;
         if(idx < 0) { break; }
         complex_accumulate<T_data, T>(y, input[idx], poly_filters[phase][k]);
         // static_cast<float>(poly_filters[phase][k]));
      }
      complex_copy<T_data>(y, output_buffer[out_count]);
      out_count++;

      phase += down_factor;
      uint32_t advance = phase / up_factor;
      phase            = phase % up_factor;
      input_index += advance;
   }
   output_size = out_count;
   return output_buffer;
}

// Pure decimation branch (up_factor == 1).
// In this branch, input_rate/output_rate is an integer
// Convolution is performed using circular indexing.
template<typename T_data, typename T>
inline T_data*
polyphase_decimate(const T_data*         input,
                   size_t                input_size,
                   T                     input_rate,
                   T                     output_rate,
                   const std::vector<T>& full_filter,
                   size_t&               output_size)
{
   uint32_t down_factor        = static_cast<uint32_t>(std::round(
      static_cast<double>(input_rate) / static_cast<double>(output_rate)));
   output_size                 = static_cast<size_t>(std::round(
      static_cast<double>(input_size) / static_cast<double>(down_factor)));
   uint32_t num_taps           = full_filter.size();
   uint32_t delay              = (num_taps - 1) / 2;
   double   double_down_factor = static_cast<double>(down_factor);

   T_data* output = new T_data[output_size];
   for(uint32_t i = 0; i < output_size; i++)
   {
      uint32_t index = delay
                       + static_cast<uint32_t>(std::round(
                          static_cast<double>(i) * double_down_factor));

      T_data y;
      complex_set_zero<T_data>(y);
      for(uint32_t k = 0; k < num_taps; k++)
      {
         uint32_t       idx          = index - k;
         const uint32_t sample_index = circular_index(input_size, idx);
         complex_accumulate<T_data, T>(y, input[sample_index], full_filter[k]);
      }
      complex_copy<T_data>(y, output[i]);
   }
   return output;
}

// Pure upsampling branch (down_factor == 1).
// In this branch, output_rate/input_rate is an integer
// Convolution is performed using circular indexing.
template<typename T_data, typename T>
inline T_data*
polyphase_upsample(const T_data*         input,
                   size_t                input_size,
                   T                     input_rate,
                   T                     output_rate,
                   const std::vector<T>& full_filter,
                   size_t&               output_size)
{
   uint32_t up_factor = static_cast<uint32_t>(std::round(
      static_cast<double>(output_rate) / static_cast<double>(input_rate)));
   output_size        = input_size * up_factor;
   uint32_t num_taps  = full_filter.size();
   uint32_t delay     = (num_taps - 1) / 2;

   uint32_t                    L            = up_factor;
   uint32_t                    phase_length = (num_taps + L - 1) / L;
   std::vector<std::vector<T>> poly_filters(
      L,
      std::vector<T>(phase_length, T(0)));
   for(uint32_t p = 0; p < L; p++)
   {
      for(uint32_t k = 0; k < phase_length; k++)
      {
         uint32_t index = p + k * L;
         if(index < num_taps) { poly_filters[p][k] = full_filter[index]; }
         else { poly_filters[p][k] = T(0); }
      }
   }

   T_data* output = new T_data[output_size];
   for(size_t i = 0; i < output_size; i++)
   {
      uint32_t base  = static_cast<uint32_t>(i / up_factor);
      uint32_t phase = static_cast<uint32_t>(i % up_factor);
      uint32_t index = delay + base;
      T_data   y;
      complex_set_zero<T_data>(y);
      for(uint32_t k = 0; k < phase_length; k++)
      {
         uint32_t       idx          = index - k;
         const uint32_t sample_index = circular_index(input_size, idx);
         complex_accumulate<T_data, T>(y,
                                       input[sample_index],
                                       poly_filters[phase][k]);
      }
      complex_copy<T_data>(y, output[i]);
   }
   return output;
}

// polyphase_resample_circular
// Performs polyphase resampling on a circular input array of T_data samples
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
circular_polyphase_resample(const T_data*         input,
                            size_t                input_size,
                            T                     input_sampling_rate,
                            T                     output_sampling_rate,
                            const std::vector<T>& full_filter,
                            size_t&               output_size)
{
   up_down_factors conv = compute_conversion_factors<T>(input_sampling_rate,
                                                        output_sampling_rate);
   uint32_t        up_factor   = conv.up_factor;
   uint32_t        down_factor = conv.down_factor;

   if(up_factor == 1 && down_factor > 1)
      return polyphase_decimate<T_data, T>(input,
                                           input_size,
                                           input_sampling_rate,
                                           output_sampling_rate,
                                           full_filter,
                                           output_size);
   if(down_factor == 1 && up_factor > 1)
      return polyphase_upsample<T_data, T>(input,
                                           input_size,
                                           input_sampling_rate,
                                           output_sampling_rate,
                                           full_filter,
                                           output_size);

   output_size = static_cast<size_t>(
      std::round(static_cast<double>(input_size)
                 * (static_cast<double>(output_sampling_rate)
                    / static_cast<double>(input_sampling_rate))));

   uint32_t num_taps     = full_filter.size();
   uint32_t delay        = (num_taps - 1) / 2;
   uint32_t L            = up_factor;   // number of polyphase subfilters
   uint32_t phase_length = (num_taps + L - 1) / L;
   std::vector<std::vector<T>> poly_filters(
      L,
      std::vector<T>(phase_length, T(0)));
   for(uint32_t p = 0; p < L; ++p)
   {
      for(uint32_t k = 0; k < phase_length; ++k)
      {
         uint32_t index = p + k * L;
         if(index < num_taps) { poly_filters[p][k] = full_filter[index]; }
         else { poly_filters[p][k] = T(0); }
      }
   }

   // Group delay of the filter.
   uint32_t input_index   = delay;
   uint32_t phase         = 0;
   T_data*  output_buffer = new T_data[output_size];

   for(size_t i = 0; i < output_size; i++)
   {
      T_data y;
      complex_set_zero<T_data>(y);
      for(uint32_t k = 0; k < phase_length; ++k)
      {
         uint32_t       idx          = input_index - k;
         const uint32_t circular_idx = circular_index(input_size, idx);
         complex_accumulate<T_data, T>(y,
                                       input[circular_idx],
                                       poly_filters[phase][k]);
      }
      complex_copy<T_data>(y, output_buffer[i]);

      phase += down_factor;
      uint32_t advance = static_cast<uint32_t>(std::floor(
         static_cast<double>(phase) / static_cast<double>(up_factor)));
      phase            = phase % up_factor;

      // Advance input_index circularly
      uint32_t relative_index = input_index - delay;
      relative_index
         = (relative_index + advance) % static_cast<uint32_t>(input_size);
      input_index = delay + relative_index;
   }
   return output_buffer;
}

#endif /* RESAMPLER_H_ */
