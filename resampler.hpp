#ifndef RESAMPLER_H_
#define RESAMPLER_H_

#include <algorithm>
#include <cassert>

typedef float  complex_float[2];
typedef double complex_double[2];

typedef struct
{
   int up_factor;
   int down_factor;
} up_down_factors;

inline void
complex_set_zero(complex_float out)
{
   out[0] = 0.0f;
   out[1] = 0.0f;
}

inline void
complex_copy(const complex_float src, complex_float dst)
{
   dst[0] = src[0];
   dst[1] = src[1];
}

inline void
complex_accumulate(complex_float       accumulator,
                   const complex_float in,
                   float               scalar)
{
   accumulator[0] += in[0] * scalar;
   accumulator[1] += in[1] * scalar;
}

inline int
gcd_int(int a, int b)
{
   while(b != 0)
   {
      int temp = a % b;
      a        = b;
      b        = temp;
   }
   return a;
}

template<typename T>
inline up_down_factors compute_conversion_factors(T input_rate, T target_rate)
{
   int input_rate_int = static_cast<int>(input_rate);
   int target_rate_int = static_cast<int>(target_rate);
   int gcd_val = gcd_int(target_rate_int, input_rate_int);
   up_down_factors conv;
   conv.up_factor = target_rate_int / gcd_val;
   conv.down_factor = input_rate_int / gcd_val;
   return conv;
}

// polyphase_resample
// Performs polyphase resampling on an input array of complex_float samples using a full FIR filter.
// The function takes the input and desired output sampling rates (of type T), computes integer conversion factors,
// and uses the provided full FIR filter (designed via design_raised_cosine_fir) for filtering.
// Parameters:
//   input               : pointer to input complex signal (complex_float*)
//   input_size          : number of input samples
//   input_sampling_rate : sampling rate of the input signal (Hz)
//   output_sampling_rate: desired output sampling rate (Hz)
//   full_filter         : full FIR filter coefficients
//   output_size         : (output) number of output samples computed
// Returns a dynamically allocated array of complex_float samples (caller must delete[]).
template<typename T>
inline complex_float* polyphase_resample(const complex_float* input, size_t input_size,
                                  T input_sampling_rate, T output_sampling_rate,
                                  const std::vector<T>& full_filter,
                                  size_t &output_size)
{
   // Compute integer conversion factors.
   up_down_factors conv = compute_conversion_factors<T>(input_sampling_rate, output_sampling_rate);
   int up_factor = conv.up_factor;
   int down_factor = conv.down_factor;

   int num_taps = full_filter.size();
   int L = up_factor;  // number of polyphase subfilters
   int phase_length = (num_taps + L - 1) / L;
   std::vector< std::vector<T> > poly_filters;
   poly_filters.resize(L);
   for (int p = 0; p < L; p++) {
      poly_filters[p].resize(phase_length, T(0));
      for (int k = 0; k < phase_length; k++) {
         int index = p + k * L;
         if (index < num_taps) {
            poly_filters[p][k] = full_filter[index];
         } else {
            poly_filters[p][k] = T(0);
         }
      }
   }

   // Group delay of the filter.
   int delay = (num_taps - 1) / 2;
   int input_index = delay;
   int phase = 0;

   size_t max_output_samples = (((input_size - delay) * up_factor) / down_factor) + 1;
   complex_float* output_buffer = new complex_float[max_output_samples];
   size_t out_count = 0;

   while (input_index < static_cast<int>(input_size))
   {
      complex_float y;
      complex_set_zero(y);
      for (int k = 0; k < phase_length; k++) {
         int idx = input_index - k;
         if (idx < 0) {
            break;
         }
         complex_accumulate(y, input[idx], static_cast<float>(poly_filters[phase][k]));
      }
      complex_copy(y, output_buffer[out_count]);
      out_count++;

      phase += down_factor;
      int advance = phase / up_factor;
      phase = phase % up_factor;
      input_index += advance;
   }
   output_size = out_count;
   return output_buffer;
}

#endif /* RESAMPLER_H_ */
