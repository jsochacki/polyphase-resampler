#include "filter_design.hpp"

int
main()
{
   // Example usage: these parameters mimic MATLAB's firrcos(32, 4, 0.1, 32)
   int   filter_order = 32;
   float cutoff       = 4.0f;
   float rolloff      = 0.1f;
   float fs           = 32.0f;

   std::vector<float> h
      = design_raised_cosine_filter<float>(filter_order, cutoff, rolloff, fs);

   for(float coef : h) { printf("%f\n", coef); }

   printf("That was float\n");

   double dcutoff       = 4.0d;
   double drolloff      = 0.1d;
   double dfs           = 32.0d;

   std::vector<double> dh
      = design_raised_cosine_filter<double>(filter_order, cutoff, rolloff, fs);

   for(double coef : dh) { printf("%f\n", coef); }

   printf("That was double\n");

   return 0;
}
