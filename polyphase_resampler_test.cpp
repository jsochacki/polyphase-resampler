#include "filter_design.hpp"
#include "resampler.hpp"

//g++ -std=c++17 -I/filter_design.hpp -O2 -o filter_design filter_design.cpp 

inline void polyphase_resampler_float_test(void)
{
   using T = float;
   // Set input and target sample rates.
   T input_rate = static_cast<T>(60e6);    // 60 MHz input
   T output_rate = static_cast<T>(1.92e6);   // 1.92 MHz desired output

   // Compute and print conversion factors.
   up_down_factors conv = compute_conversion_factors<T>(input_rate, output_rate);
   std::printf("Conversion ratio: up_factor = %d, down_factor = %d\n", conv.up_factor, conv.down_factor);

   // Filter design parameters.
   // For the filter design, use a sampling rate of: design_fs = input_rate * up_factor.
   T design_fs = input_rate * conv.up_factor;
   // Set cutoff equal to output_rate (i.e. design_cutoff = input_rate * up_factor / down_factor).
   T design_cutoff = input_rate * conv.up_factor / static_cast<T>(conv.down_factor);
   T rolloff = static_cast<T>(0.1);
   int filter_order = 32; // will be adjusted to odd if needed

   std::vector<T> full_filter = design_raised_cosine_filter<T>(filter_order, design_cutoff, rolloff, design_fs);

   // Generate an input signal: a pure tone at 1 MHz (within the passband).
   size_t num_samples = 10000;
   complex_float* input_signal = new complex_float[num_samples];
   T tone_freq = static_cast<T>(1e6); // 1 MHz tone
   for (size_t n = 0; n < num_samples; n++) {
      T t = static_cast<T>(n) / input_rate;
      input_signal[n][0] = std::cos(2.0 * M_PI * tone_freq * t);
      input_signal[n][1] = std::sin(2.0 * M_PI * tone_freq * t);
   }

   // Resample using the polyphase resampler.
   size_t output_size = 0;
   complex_float* output_signal = polyphase_resample<T>(input_signal, num_samples,
                                                          input_rate, output_rate,
                                                          full_filter, output_size);

   // Compute the average magnitude of the output signal.
   T sum_mag = 0;
   for (size_t i = 0; i < output_size; i++) {
      T mag = std::sqrt(output_signal[i][0] * output_signal[i][0] +
                        output_signal[i][1] * output_signal[i][1]);
      sum_mag += mag;
   }
   T avg_mag = sum_mag / static_cast<T>(output_size);
   std::printf("Average magnitude of output signal: %f\n", static_cast<float>(avg_mag));

   // Print the first 10 output samples.
   std::printf("First 10 output samples:\n");
   for (size_t i = 0; i < (output_size < 10 ? output_size : 10); i++) {
      std::printf("(%f, %f)\n", output_signal[i][0], output_signal[i][1]);
   }

   // Clean up.
   delete[] input_signal;
   delete[] output_signal;
}

int main(void)
{
   polyphase_resampler_float_test();
   return 0;
}
