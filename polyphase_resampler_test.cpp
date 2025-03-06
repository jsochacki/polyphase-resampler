#include "filter_design.hpp"
#include "resampler.hpp"

inline void
polyphase_resampler_float_test(void)
{
   using T          = float;
   using T_data     = complex_float;
   T   input_rate   = static_cast<T>(60e6);
   T   output_rate  = static_cast<T>(1.92e6);
   T   rolloff      = T(0.1);
   int filter_order = 1 << 12;   // 1 << 10 sufficient for RF sampling rates

   std::vector<T> full_filter = design_raised_cosine_filter(input_rate,
                                                            output_rate,
                                                            filter_order,
                                                            rolloff);

   FILE* file_coeff = std::fopen("filter_coefficients.txt", "w");
   if(file_coeff == nullptr)
   {
      std::perror("Error opening filter_coefficients.txt for writing");
      return;
   }
   for(size_t i = 0; i < full_filter.size(); ++i)
   {
      // Print with conditional formatting for the sign.
      if(full_filter[i] >= 0)
         std::fprintf(file_coeff, "+%f\n", full_filter[i]);
      else
         std::fprintf(file_coeff, "%f\n", full_filter[i]);
   }
   std::fclose(file_coeff);


   // Generate an input signal: a pure tone at 1.92 MHz / 2 (within the
   // passband).
   size_t  number_of_samples    = 10000;
   T_data* input_signal         = new T_data[number_of_samples];
   T       tone_freq_passband   = (output_rate / 2) - (output_rate / 316);
   T       tone_freq_rejectband = output_rate * 2;
   for(size_t n = 0; n < number_of_samples; ++n)
   {
      T t                = static_cast<T>(n) / input_rate;
      input_signal[n][0] = std::cos(2.0 * M_PI * tone_freq_passband * t);
      input_signal[n][1] = std::sin(2.0 * M_PI * tone_freq_passband * t);
      input_signal[n][0] += std::cos(2.0 * M_PI * tone_freq_rejectband * t);
      input_signal[n][1] += std::sin(2.0 * M_PI * tone_freq_rejectband * t);
   }

   // Resample using the polyphase resampler.
   size_t  output_size   = 0;
   T_data* output_signal = polyphase_resample<T_data, T>(input_signal,
                                                 number_of_samples,
                                                 input_rate,
                                                 output_rate,
                                                 full_filter,
                                                 output_size);


   FILE* file_in = std::fopen("input_signal.txt", "w");
   if(file_in == nullptr)
   {
      std::perror("Error opening input_signal.txt for writing");
      delete[] input_signal;
      delete[] output_signal;
      return;
   }
   for(size_t i = 0; i < number_of_samples; ++i)
   {
      // Print with conditional formatting for the sign.
      if(input_signal[i][1] >= 0)
         std::fprintf(file_in,
                      "%f+%fi\n",
                      input_signal[i][0],
                      input_signal[i][1]);
      else
         std::fprintf(file_in,
                      "%f%fi\n",
                      input_signal[i][0],
                      input_signal[i][1]);
   }
   std::fclose(file_in);

   // Write the output signal to a text file.
   FILE* file_out = std::fopen("output_signal.txt", "w");
   if(file_out == nullptr)
   {
      std::perror("Error opening output_signal.txt for writing");
      delete[] input_signal;
      delete[] output_signal;
      return;
   }
   for(size_t i = 0; i < output_size; ++i)
   {
      if(output_signal[i][1] >= 0)
         std::fprintf(file_out,
                      "%f+%fi\n",
                      output_signal[i][0],
                      output_signal[i][1]);
      else
         std::fprintf(file_out,
                      "%f%fi\n",
                      output_signal[i][0],
                      output_signal[i][1]);
   }
   std::fclose(file_out);

   // Compute the average magnitude of the signals
   T sum_mag_in  = 0;
   T sum_mag_out = 0;
   for(size_t i = 0; i < output_size; ++i)
   {
      T mag = std::sqrt(input_signal[i][0] * input_signal[i][0]
                        + input_signal[i][1] * input_signal[i][1]);
      sum_mag_in += mag;

      mag = std::sqrt(output_signal[i][0] * output_signal[i][0]
                      + output_signal[i][1] * output_signal[i][1]);
      sum_mag_out += mag;
   }
   /*
   T avg_mag = sum_mag_in / static_cast<T>(number_of_samples);
   std::printf("Average magnitude of input signal: %f\n",
               static_cast<float>(avg_mag));
   avg_mag = sum_mag_out / static_cast<T>(output_size);
   std::printf("Average magnitude of output signal: %f\n",
               static_cast<float>(avg_mag));


   // Print the first 10 output samples.
   std::printf("First 10 output samples:\n");
   for(size_t i = 0; i < (output_size < 10 ? output_size : 10); ++i)
   {
      std::printf("(%f, %f)\n", output_signal[i][0], output_signal[i][1]);
   }
   */

   // Clean up.
   delete[] input_signal;
   delete[] output_signal;
}


int
main(void)
{
   polyphase_resampler_float_test();
   return 0;
}
