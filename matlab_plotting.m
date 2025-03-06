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
tbl = readtable("./input_signal.txt", opts);
input_signal = tbl.i1;
clear tb1
tbl = readtable("./output_signal.txt", opts);
output_signal = tbl.i1;
clear tb1
tbl = readtable("./filter_coefficients.txt", opts);
h = tbl.i1;
clear opts tbl

plot(linspace(0,60000000,length(input_signal)), 10*log10(abs(fft(input_signal, length(input_signal)))))
plot(linspace(0,1920000,length(output_signal)), 10*log10(abs(fft(output_signal, length(output_signal)))))
plot(10*log10(fftshift(abs(fft(output_signal, length(output_signal))))))
plot(real(output_signal))

fvtool(h)
