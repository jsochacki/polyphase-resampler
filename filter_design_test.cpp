#include "filter_design.hpp"

int
main()
{
   // Example usage: these parameters mimic MATLAB's firrcos(2^12, 1920000/2,
   // 0.1, 60000000)
   int   filter_order = 1 << 12;
   float cutoff       = 1920000.0f / 2.0f;
   float rolloff      = 0.1f;
   float fs           = 60000000.0f;

   std::vector<float> h
      = design_raised_cosine_filter<float>(filter_order, cutoff, rolloff, fs);

   // for(float coef : h) { printf("%.24f\n", coef); }


   FILE* file_float = std::fopen("float_filter_coefficients.txt", "w");
   if(file_float == nullptr)
   {
      std::perror("Error opening float_filter_coefficients.txt for writing");
      return -1;
   }
   for(size_t i = 0; i < h.size(); ++i)
   {
      // Print with conditional formatting for the sign.
      if(h[i] >= 0)
         std::fprintf(file_float, "+%.24f\n", h[i]);
      else
         std::fprintf(file_float, "%.24f\n", h[i]);
   }
   std::fclose(file_float);


   //printf("That was float\n");

   double dcutoff  = 1920000.0d / 2.0d;
   double drolloff = 0.1d;
   double dfs      = 60000000.0d;

   std::vector<double> dh
      = design_raised_cosine_filter<double>(filter_order, cutoff, rolloff, fs);

   // for(double coef : dh) { printf("%.24f\n", coef); }

   FILE* file_double = std::fopen("double_filter_coefficients.txt", "w");
   if(file_double == nullptr)
   {
      std::perror("Error opening double_filter_coefficients.txt for writing");
      return -1;
   }
   for(size_t i = 0; i < dh.size(); ++i)
   {
      // Print with conditional formatting for the sign.
      if(dh[i] >= 0)
         std::fprintf(file_double, "+%.24f\n", dh[i]);
      else
         std::fprintf(file_double, "%.24f\n", dh[i]);
   }
   std::fclose(file_double);


   //printf("That was double\n");

   return 0;
}
