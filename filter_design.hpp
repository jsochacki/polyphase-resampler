#ifndef FILTER_DESIGN_H_
#define FILTER_DESIGN_H_

#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

template<typename T>
inline T
get_epsilon_for_type(void)
{
   if constexpr(std::is_same_v<T, float>) { return FLT_EPSILON; }
   else if constexpr(std::is_same_v<T, double>) { return DBL_EPSILON; }
   else
   {
      static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>,
                    "get_epsilon_for_type only supports float and double");
   }
}

typedef struct
{
   int up_factor;
   int down_factor;
} up_down_factors;

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
inline up_down_factors
compute_conversion_factors(T input_rate, T target_rate)
{
   int             input_rate_int  = static_cast<int>(input_rate);
   int             target_rate_int = static_cast<int>(target_rate);
   int             gcd_val         = gcd_int(target_rate_int, input_rate_int);
   up_down_factors conv;
   conv.up_factor   = target_rate_int / gcd_val;
   conv.down_factor = input_rate_int / gcd_val;
   return conv;
}

// Raised cosine FIR design function.
// Parameters:
//   filter_order : the number of filter taps (filter order)
//   cutoff  : cutoff frequency in Hz (should be symbol_rate/2)
//   rolloff : rolloff factor (alpha, typically between 0 and 1)
//   fs      : sampling rate in Hz (1 in normalized frequency corresponds to
//   fs)

template<typename T>
inline std::vector<T>
design_raised_cosine_filter(int filter_order, T cutoff, T rolloff, T fs)
{
   T local_epsilon = get_epsilon_for_type<T>();
   if((filter_order % 2) == 0) { filter_order += 1; }
   std::vector<T> h(filter_order, 0.0f);

   T center_tap    = T(filter_order - 1) / T(2);
   T symbol_period = T(1) / (T(2) * cutoff);

   // Compute the impulse response at each tap.
   for(uint16_t n = 0; n < filter_order; ++n)
   {
      // Time offset in seconds; note that sample spacing is 1/fs.
      T t     = (n - center_tap) / fs;
      T value = T(0);

      // Handle the t == 0 case.
      if(std::abs(t) < local_epsilon) { value = T(1); }
      // Handle the singularity at |t| == T/(2*rolloff)
      else if(std::abs(std::abs(t) - symbol_period / (T(2) * rolloff))
              < local_epsilon)
      {
         T temp = static_cast<T>(M_PI) / (T(2) * rolloff);
         value  = (rolloff / T(2)) * (std::sin(temp) / temp);
      }
      else
      {
         T t_over_T = t / symbol_period;
         // Sinc term: sin(π t/T)/(π t/T)
         T sinc_val = std::sin(static_cast<T>(M_PI) * t_over_T)
                      / (static_cast<T>(M_PI) * t_over_T);
         // Cosine term: cos(π·rolloff·t/T)
         T cos_val = std::cos(static_cast<T>(M_PI) * rolloff * t_over_T);
         // Denominator: 1 - 4·rolloff²·t²/T²
         T denom = T(1)
                   - T(4) * rolloff * rolloff * t * t
                        / (symbol_period * symbol_period);
         value = sinc_val * cos_val / denom;
      }
      h[n] = value;
   }

   // Normalize the filter coefficients so that their sum equals unity.
   T sum = T(0);
   for(int n = 0; n < filter_order; ++n) { sum += h[n]; }
   for(int n = 0; n < filter_order; ++n) { h[n] /= sum; }

   return h;
}

template<typename T>
inline std::vector<T>
design_raised_cosine_filter(T   input_rate,
                            T   output_rate,
                            int filter_order,
                            T   rolloff)
{
   up_down_factors conv
      = compute_conversion_factors<T>(input_rate, output_rate);
   // std::printf("Conversion ratio: up_factor = %d, down_factor = %d\n",
   // conv.up_factor, conv.down_factor);

   T upsampled_sample_rate = input_rate * conv.up_factor;
   T design_cutoff         = (input_rate * static_cast<T>(conv.up_factor))
                     / static_cast<T>(conv.down_factor);
   return design_raised_cosine_filter<T>(filter_order,
                                         design_cutoff,
                                         rolloff,
                                         upsampled_sample_rate);
}


#endif /* FILTER_DESIGN_H_ */
