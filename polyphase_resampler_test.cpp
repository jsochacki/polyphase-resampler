#include "filter_design.hpp"
#include "resampler.hpp"

// Just for timing
#include <chrono>
#include <string>

inline void
set_cf_t_value(cf_t& output, float real_value, float imag_value)
{
   __real__ output = real_value;
   __imag__ output = imag_value;
}

inline void
add_cf_t_value(cf_t& output, float real_value, float imag_value)
{
   __real__ output += real_value;
   __imag__ output += imag_value;
}

inline void
get_cf_t_parts(const cf_t& input, float& real_value, float& imag_value)
{
   real_value = __real__ input;
   imag_value = __imag__ input;
}

template<typename T, typename T_data>
inline void
test_polyphase_resampler(const std::string& prefix,
                         size_t             number_of_samples,
                         T                  input_rate,
                         T                  output_rate,
                         int                filter_order,
                         T                  rolloff)
{
   // Design the raised cosine filter.
   std::vector<T> full_filter = design_raised_cosine_filter(input_rate,
                                                            output_rate,
                                                            filter_order,
                                                            rolloff);

   auto start = std::chrono::high_resolution_clock::now();
   // Write filter coefficients to file.
   std::string coeff_file_name = prefix + "_filter_coefficients.txt";
   FILE*       file_coeff      = std::fopen(coeff_file_name.c_str(), "w");
   if(file_coeff == nullptr)
   {
      std::perror(
         ("Error opening " + coeff_file_name + " for writing").c_str());
      return;
   }
   for(size_t i = 0; i < full_filter.size(); ++i)
   {
      if(full_filter[i] >= 0)
         std::fprintf(file_coeff, "+%.24f\n", full_filter[i]);
      else
         std::fprintf(file_coeff, "%.24f\n", full_filter[i]);
   }
   std::fclose(file_coeff);
   auto end = std::chrono::high_resolution_clock::now();
   auto duration
      = std::chrono::duration_cast<std::chrono::microseconds>(end - start)
           .count();
   printf("%s_filter_coefficients_write duration: %lld microseconds\n",
          prefix.c_str(),
          (long long)duration);


   // Generate an input signal.
   // We generate a composite tone: one frequency in the passband and one in
   // the rejectband.
   T_data* input_signal         = new T_data[number_of_samples];
   T       tone_freq_passband   = (output_rate / 2) - 8 * (output_rate / 316);
   T       tone_freq_rejectband = (output_rate / 2) + 8 * (output_rate / 316);
   for(size_t n = 0; n < number_of_samples; ++n)
   {
      T t = static_cast<T>(n) / input_rate;
      if constexpr(std::is_same<T_data, cf_t>::value)
      {
         // Use helper functions for cf_t
         set_cf_t_value(input_signal[n],
                        std::cos(2.0 * M_PI * tone_freq_passband * t),
                        std::sin(2.0 * M_PI * tone_freq_passband * t));
         add_cf_t_value(input_signal[n],
                        std::cos(2.0 * M_PI * tone_freq_rejectband * t),
                        std::sin(2.0 * M_PI * tone_freq_rejectband * t));
      }
      else
      {
         input_signal[n][0] = std::cos(2.0 * M_PI * tone_freq_passband * t);
         input_signal[n][1] = std::sin(2.0 * M_PI * tone_freq_passband * t);
         input_signal[n][0] += std::cos(2.0 * M_PI * tone_freq_rejectband * t);
         input_signal[n][1] += std::sin(2.0 * M_PI * tone_freq_rejectband * t);
      }
   }

   // Resample using the polyphase_resample function and measure its duration.
   size_t output_size    = 0;
   start                 = std::chrono::high_resolution_clock::now();
   T_data* output_signal = polyphase_resample<T_data, T>(input_signal,
                                                         number_of_samples,
                                                         input_rate,
                                                         output_rate,
                                                         full_filter,
                                                         output_size);
   end                   = std::chrono::high_resolution_clock::now();
   duration
      = std::chrono::duration_cast<std::chrono::microseconds>(end - start)
           .count();
   printf("%s polyphase_resample Duration: %lld microseconds\n",
          prefix.c_str(),
          (long long)duration);
   printf("number of output samples: %ld \n", output_size);

   start = std::chrono::high_resolution_clock::now();
   // Write the input signal to file.
   std::string input_file_name = prefix + "_input_signal.txt";
   FILE*       file_in         = std::fopen(input_file_name.c_str(), "w");
   if(file_in == nullptr)
   {
      std::perror(
         ("Error opening " + input_file_name + " for writing").c_str());
      delete[] input_signal;
      delete[] output_signal;
      return;
   }
   for(size_t i = 0; i < number_of_samples; ++i)
   {
      if constexpr(std::is_same<T_data, cf_t>::value)
      {
         float real_part = __real__ input_signal[i];
         float imag_part = __imag__ input_signal[i];
         if(imag_part >= 0)
            std::fprintf(file_in, "%.24f+%.24fi\n", real_part, imag_part);
         else
            std::fprintf(file_in, "%.24f%.24fi\n", real_part, imag_part);
      }
      else
      {
         if(input_signal[i][1] >= 0)
            std::fprintf(file_in,
                         "%.24f+%.24fi\n",
                         input_signal[i][0],
                         input_signal[i][1]);
         else
            std::fprintf(file_in,
                         "%.24f%.24fi\n",
                         input_signal[i][0],
                         input_signal[i][1]);
      }
   }
   std::fclose(file_in);
   end = std::chrono::high_resolution_clock::now();
   duration
      = std::chrono::duration_cast<std::chrono::microseconds>(end - start)
           .count();
   printf("%s_input_signal write duration: %lld microseconds\n",
          prefix.c_str(),
          (long long)duration);

   start = std::chrono::high_resolution_clock::now();
   // Write the output signal to file.
   std::string output_file_name = prefix + "_output_signal.txt";
   FILE*       file_out         = std::fopen(output_file_name.c_str(), "w");
   if(file_out == nullptr)
   {
      std::perror(
         ("Error opening " + output_file_name + " for writing").c_str());
      delete[] input_signal;
      delete[] output_signal;
      return;
   }
   for(size_t i = 0; i < output_size; ++i)
   {
      if constexpr(std::is_same<T_data, cf_t>::value)
      {
         float real_part = __real__ output_signal[i];
         float imag_part = __imag__ output_signal[i];
         if(imag_part >= 0)
            std::fprintf(file_out, "%.24f+%.24fi\n", real_part, imag_part);
         else
            std::fprintf(file_out, "%.24f%.24fi\n", real_part, imag_part);
      }
      else
      {
         if(output_signal[i][1] >= 0)
            std::fprintf(file_out,
                         "%.24f+%.24fi\n",
                         output_signal[i][0],
                         output_signal[i][1]);
         else
            std::fprintf(file_out,
                         "%.24f%.24fi\n",
                         output_signal[i][0],
                         output_signal[i][1]);
      }
   }
   std::fclose(file_out);
   end = std::chrono::high_resolution_clock::now();
   duration
      = std::chrono::duration_cast<std::chrono::microseconds>(end - start)
           .count();
   printf("%s_output_signal write duration: %lld microseconds\n",
          prefix.c_str(),
          (long long)duration);


   // Clean up allocated memory.
   delete[] input_signal;
   delete[] output_signal;
}

// The following function tests the circular_polyphase_resample implementation.
// The structure is analogous to the previous function; however, we use
// parameters suited for circular buffering (for example, a different input
// sampling rate and a larger number of samples).
template<typename T, typename T_data>
inline void
test_circular_polyphase_resampler(const std::string& prefix,
                                  size_t             number_of_samples,
                                  T                  input_rate,
                                  T                  output_rate,
                                  int                filter_order,
                                  T                  rolloff)
{
   std::vector<T> full_filter = design_raised_cosine_filter(input_rate,
                                                            output_rate,
                                                            filter_order,
                                                            rolloff);

   auto        start           = std::chrono::high_resolution_clock::now();
   std::string coeff_file_name = prefix + "_filter_coefficients.txt";
   FILE*       file_coeff      = std::fopen(coeff_file_name.c_str(), "w");
   if(file_coeff == nullptr)
   {
      std::perror(
         ("Error opening " + coeff_file_name + " for writing").c_str());
      return;
   }
   for(size_t i = 0; i < full_filter.size(); ++i)
   {
      if(full_filter[i] >= 0)
         std::fprintf(file_coeff, "+%.24f\n", full_filter[i]);
      else
         std::fprintf(file_coeff, "%.24f\n", full_filter[i]);
   }
   std::fclose(file_coeff);
   auto end = std::chrono::high_resolution_clock::now();
   auto duration
      = std::chrono::duration_cast<std::chrono::microseconds>(end - start)
           .count();
   printf("%s_filter_coefficients_write duration: %lld microseconds\n",
          prefix.c_str(),
          (long long)duration);


   T_data* input_signal         = new T_data[number_of_samples];
   T       tone_freq_passband   = (output_rate / 2) - 8 * (output_rate / 316);
   T       tone_freq_rejectband = (output_rate / 2) + 8 * (output_rate / 316);
   for(size_t n = 0; n < number_of_samples; ++n)
   {
      T t = static_cast<T>(n) / input_rate;
      if constexpr(std::is_same<T_data, cf_t>::value)
      {
         set_cf_t_value(input_signal[n],
                        std::cos(2.0 * M_PI * tone_freq_passband * t),
                        std::sin(2.0 * M_PI * tone_freq_passband * t));
         add_cf_t_value(input_signal[n],
                        std::cos(2.0 * M_PI * tone_freq_rejectband * t),
                        std::sin(2.0 * M_PI * tone_freq_rejectband * t));
      }
      else
      {
         input_signal[n][0] = std::cos(2.0 * M_PI * tone_freq_passband * t);
         input_signal[n][1] = std::sin(2.0 * M_PI * tone_freq_passband * t);
         input_signal[n][0] += std::cos(2.0 * M_PI * tone_freq_rejectband * t);
         input_signal[n][1] += std::sin(2.0 * M_PI * tone_freq_rejectband * t);
      }
   }

   size_t output_size = 0;
   start              = std::chrono::high_resolution_clock::now();
   T_data* output_signal
      = circular_polyphase_resample<T_data, T>(input_signal,
                                               number_of_samples,
                                               input_rate,
                                               output_rate,
                                               full_filter,
                                               output_size);
   end = std::chrono::high_resolution_clock::now();
   duration
      = std::chrono::duration_cast<std::chrono::microseconds>(end - start)
           .count();
   printf("%s circular_polyphase_resample Duration: %lld microseconds\n",
          prefix.c_str(),
          (long long)duration);
   printf("number of output samples: %ld \n", output_size);

   start                       = std::chrono::high_resolution_clock::now();
   std::string input_file_name = prefix + "_input_signal.txt";
   FILE*       file_in         = std::fopen(input_file_name.c_str(), "w");
   if(file_in == nullptr)
   {
      std::perror(
         ("Error opening " + input_file_name + " for writing").c_str());
      delete[] input_signal;
      delete[] output_signal;
      return;
   }
   for(size_t i = 0; i < number_of_samples; ++i)
   {
      if constexpr(std::is_same<T_data, cf_t>::value)
      {
         float real_part = __real__ input_signal[i];
         float imag_part = __imag__ input_signal[i];
         if(imag_part >= 0)
            std::fprintf(file_in, "%.24f+%.24fi\n", real_part, imag_part);
         else
            std::fprintf(file_in, "%.24f%.24fi\n", real_part, imag_part);
      }
      else
      {
         if(input_signal[i][1] >= 0)
            std::fprintf(file_in,
                         "%.24f+%.24fi\n",
                         input_signal[i][0],
                         input_signal[i][1]);
         else
            std::fprintf(file_in,
                         "%.24f%.24fi\n",
                         input_signal[i][0],
                         input_signal[i][1]);
      }
   }
   std::fclose(file_in);
   end = std::chrono::high_resolution_clock::now();
   duration
      = std::chrono::duration_cast<std::chrono::microseconds>(end - start)
           .count();
   printf("%s_input_signal write duration: %lld microseconds\n",
          prefix.c_str(),
          (long long)duration);

   start                        = std::chrono::high_resolution_clock::now();
   std::string output_file_name = prefix + "_output_signal.txt";
   FILE*       file_out         = std::fopen(output_file_name.c_str(), "w");
   if(file_out == nullptr)
   {
      std::perror(
         ("Error opening " + output_file_name + " for writing").c_str());
      delete[] input_signal;
      delete[] output_signal;
      return;
   }
   for(size_t i = 0; i < output_size; ++i)
   {
      if constexpr(std::is_same<T_data, cf_t>::value)
      {
         float real_part = __real__ output_signal[i];
         float imag_part = __imag__ output_signal[i];
         if(imag_part >= 0)
            std::fprintf(file_out, "%.24f+%.24fi\n", real_part, imag_part);
         else
            std::fprintf(file_out, "%.24f%.24fi\n", real_part, imag_part);
      }
      else
      {
         if(output_signal[i][1] >= 0)
            std::fprintf(file_out,
                         "%.24f+%.24fi\n",
                         output_signal[i][0],
                         output_signal[i][1]);
         else
            std::fprintf(file_out,
                         "%.24f%.24fi\n",
                         output_signal[i][0],
                         output_signal[i][1]);
      }
   }
   std::fclose(file_out);
   end = std::chrono::high_resolution_clock::now();
   duration
      = std::chrono::duration_cast<std::chrono::microseconds>(end - start)
           .count();
   printf("%s_output_signal write duration: %lld microseconds\n",
          prefix.c_str(),
          (long long)duration);

   delete[] input_signal;
   delete[] output_signal;
}

int
main(void)
{
   uint16_t filter_size = 1 << 8;
   test_polyphase_resampler<float, complex_float>("polyphase_float",
                                                  10000,
                                                  static_cast<float>(60e6),
                                                  static_cast<float>(1.92e6),
                                                  filter_size,
                                                  0.01f);
   test_polyphase_resampler<double, complex_double>(
      "polyphase_double",
      10000,
      static_cast<double>(60e6),
      static_cast<double>(1.92e6),
      filter_size,
      0.01);

   test_circular_polyphase_resampler<float, complex_float>(
      "circular_float",
      3072000,
      static_cast<float>(30.72e6),
      static_cast<float>(1.92e6),
      filter_size,
      0.01f);
   test_circular_polyphase_resampler<double, complex_double>(
      "circular_double",
      3072000,
      static_cast<double>(30.72e6),
      static_cast<double>(1.92e6),
      filter_size,
      0.01);


   test_circular_polyphase_resampler<float, cf_t>("circular_float_cft",
                                                  3072000,
                                                  static_cast<float>(30.72e6),
                                                  static_cast<float>(1.92e6),
                                                  filter_size,
                                                  0.01f);
   test_circular_polyphase_resampler<double, cf_t>(
      "circular_double_cft",
      3072000,
      static_cast<double>(30.72e6),
      static_cast<double>(1.92e6),
      filter_size,
      0.01);

   test_circular_polyphase_resampler<double, complex_float>(
      "circular_float_data_double_filter",
      3072000,
      static_cast<float>(30.72e6),
      static_cast<float>(1.92e6),
      filter_size,
      0.01f);

   return 0;
}
