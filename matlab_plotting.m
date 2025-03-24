close all
clear all
clc

opts = delimitedTextImportOptions("NumVariables", 1);
opts.DataLines = [1, Inf];
opts.Delimiter = ",";
opts.VariableNames = "i1";
opts.VariableTypes = "double";
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
tbl = readtable("./polyphase_float_input_signal.txt", opts);
float_input_signal = tbl.i1;
clear tb1
tbl = readtable("./polyphase_float_output_signal.txt", opts);
float_output_signal = tbl.i1;
clear tb1
tbl = readtable("./polyphase_float_filter_coefficients.txt", opts);
test_h_float = tbl.i1;
clear tb1
tbl = readtable("./polyphase_double_input_signal.txt", opts);
double_input_signal = tbl.i1;
clear tb1
tbl = readtable("./polyphase_double_output_signal.txt", opts);
double_output_signal = tbl.i1;
clear tb1
tbl = readtable("./polyphase_double_filter_coefficients.txt", opts);
test_h_double = tbl.i1;

clear tb1
tbl = readtable("./circular_float_input_signal.txt", opts);
circular_float_input_signal = tbl.i1;
clear tb1
tbl = readtable("./circular_float_output_signal.txt", opts);
circular_float_output_signal = tbl.i1;
clear tb1
tbl = readtable("./circular_float_filter_coefficients.txt", opts);
circular_test_h_float = tbl.i1;
clear tb1
tbl = readtable("./circular_double_input_signal.txt", opts);
circular_double_input_signal = tbl.i1;
clear tb1
tbl = readtable("./circular_double_output_signal.txt", opts);
circular_double_output_signal = tbl.i1;
clear tb1
tbl = readtable("./circular_double_filter_coefficients.txt", opts);
circular_test_h_double = tbl.i1;

clear tb1
tbl = readtable("./circular_float_data_double_filter_input_signal.txt", opts);
circular_float_data_double_filter_input_signal = tbl.i1;
clear tb1
tbl = readtable("./circular_float_data_double_filter_output_signal.txt", opts);
circular_float_data_double_filter_output_signal = tbl.i1;
clear tb1
tbl = readtable("./circular_float_data_double_filter_filter_coefficients.txt", opts);
circular_test_h_float_data_double_filter = tbl.i1;

clear tb1
tbl = readtable("./circular_float_cft_input_signal.txt", opts);
circular_float_cft_input_signal = tbl.i1;
clear tb1
tbl = readtable("./circular_float_cft_output_signal.txt", opts);
circular_float_cft_output_signal = tbl.i1;
clear tb1
tbl = readtable("./circular_float_cft_filter_coefficients.txt", opts);
circular_test_h_float_cft = tbl.i1;
clear tb1
tbl = readtable("./circular_double_cft_input_signal.txt", opts);
circular_double_cft_input_signal = tbl.i1;
clear tb1
tbl = readtable("./circular_double_cft_output_signal.txt", opts);
circular_double_cft_output_signal = tbl.i1;
clear tb1
tbl = readtable("./circular_double_cft_filter_coefficients.txt", opts);
circular_test_h_double_cft = tbl.i1;
clear opts tbl

plot(linspace(0,60000000,length(float_input_signal)), 10*log10(abs(fft(float_input_signal, length(float_input_signal)))))
plot(linspace(0,1920000,length(float_output_signal)), 10*log10(abs(fft(float_output_signal, length(float_output_signal)))))
plot(10*log10(fftshift(abs(fft(float_output_signal, length(float_output_signal))))))
plot(real(float_output_signal))

plot(linspace(0,60000000,length(double_input_signal)), 10*log10(abs(fft(double_input_signal, length(double_input_signal)))))
plot(linspace(0,1920000,length(double_output_signal)), 10*log10(abs(fft(double_output_signal, length(double_output_signal)))))
plot(10*log10(fftshift(abs(fft(double_output_signal, length(double_output_signal))))))
plot(real(double_output_signal))


plot(linspace(0,30720000,length(circular_float_input_signal)), 10*log10(abs(fft(circular_float_input_signal, length(circular_float_input_signal)))))
plot(linspace(0,1920000,length(circular_float_output_signal)), 10*log10(abs(fft(circular_float_output_signal, length(circular_float_output_signal)))))
plot(10*log10(fftshift(abs(fft(circular_float_output_signal, length(circular_float_output_signal))))))
plot(real(circular_float_output_signal))

plot(linspace(0,30720000,length(circular_double_input_signal)), 10*log10(abs(fft(circular_double_input_signal, length(circular_double_input_signal)))))
plot(linspace(0,1920000,length(circular_double_output_signal)), 10*log10(abs(fft(circular_double_output_signal, length(circular_double_output_signal)))))
plot(10*log10(fftshift(abs(fft(circular_double_output_signal, length(circular_double_output_signal))))))
plot(real(circular_double_output_signal))


plot(linspace(0,30720000,length(circular_float_data_double_filter_input_signal)), 10*log10(abs(fft(circular_float_data_double_filter_input_signal, length(circular_float_data_double_filter_input_signal)))))
plot(linspace(0,1920000,length(circular_float_data_double_filter_output_signal)), 10*log10(abs(fft(circular_float_data_double_filter_output_signal, length(circular_float_data_double_filter_output_signal)))))
plot(10*log10(fftshift(abs(fft(circular_float_data_double_filter_output_signal, length(circular_float_data_double_filter_output_signal))))))
plot(real(circular_float_data_double_filter_output_signal))


plot(linspace(0,30720000,length(circular_float_cft_input_signal)), 10*log10(abs(fft(circular_float_cft_input_signal, length(circular_float_cft_input_signal)))))
plot(linspace(0,1920000,length(circular_float_cft_output_signal)), 10*log10(abs(fft(circular_float_cft_output_signal, length(circular_float_cft_output_signal)))))
plot(10*log10(fftshift(abs(fft(circular_float_cft_output_signal, length(circular_float_cft_output_signal))))))
plot(real(circular_float_cft_output_signal))

plot(linspace(0,30720000,length(circular_double_cft_input_signal)), 10*log10(abs(fft(circular_double_cft_input_signal, length(circular_double_cft_input_signal)))))
plot(linspace(0,1920000,length(circular_double_cft_output_signal)), 10*log10(abs(fft(circular_double_cft_output_signal, length(circular_double_cft_output_signal)))))
plot(10*log10(fftshift(abs(fft(circular_double_cft_output_signal, length(circular_double_cft_output_signal))))))
plot(real(circular_double_cft_output_signal))


float_mlresult = upfirdn(float_input_signal, test_h_float, 4, 125);
plot(linspace(0,1920000,length(float_mlresult)), 10*log10(abs(fft(float_mlresult, length(float_mlresult)))))
plot(10*log10(fftshift(abs(fft(float_mlresult, length(float_mlresult))))))
plot(real(float_mlresult))

double_mlresult = upfirdn(double_input_signal, test_h_double, 4, 125);
plot(linspace(0,1920000,length(double_mlresult)), 10*log10(abs(fft(double_mlresult, length(double_mlresult)))))
plot(10*log10(fftshift(abs(fft(double_mlresult, length(double_mlresult))))))
plot(real(double_mlresult))


hmatlab = firrcos(2^12, 60000000*4/125, 0.1, 60000000*4, 'rolloff');
fvtool(hmatlab)
fvtool(test_h_float)
fvtool(test_h_double)
fvtool(circular_test_h_float)
fvtool(circular_test_h_double)