# Build Commands

## Build the polyphase resampler test executable

```bash
g++ -std=c++17 -I/filter_design.hpp -I/resampler.hpp -O2 -o polyphase_resampler_test polyphase_resampler_test.cpp
```

## Build the filter design test executable

```bash
g++ -std=c++17 -I/filter_design.hpp -O2 -o filter_design_test filter_design_test.cpp
```
