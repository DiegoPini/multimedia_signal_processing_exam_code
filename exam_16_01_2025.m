close all
clc
clear all


%%In particular, f1 = 1KHz, f2 = 2KHz, f3 = 8KHz. Define the signal x related to this DFT, considering a number of temporal samples equal to 10 times its period. The sampling rate of the system is 80KHz.

f1 = 1*10^3; 
f2 = 2*10^3; 
f3 = 8*10^3;
Fs = 80 * 10^3;

P1_samples = Fs/f1;
P2_samples = Fs/f2;
P3_samples = Fs/f3;
samples = lcm(lcm(P1_samples, P2_samples), P3_samples);

time = 0:1/Fs:10*P_samples/Fs - 1/Fs;

x = cos(2*pi*f1*time) + cos(2*pi*f2*time) + cos(2*pi*f3*time);

X = fft(x);
N = length(x);
freq_axis = 0:Fs/N:Fs - Fs/N;
figure;
stem(freq_axis, abs(X), 'LineWidth', 2);

%% Design an FIR stop-band filter with the windowing method in order to remove the frequency f2 (but avoiding removing f1). Consider 89 filter samples and, as cut-off frequencies, the +20% and -20% of the frequency to remove.
N_filter = 89;

filter_order = N_filter - 1;

f2_norm = f2/Fs;

f_cutoff_1 = f2_norm - f2_norm*0.2;
f_cutoff_2 = f2_norm + f2_norm*0.2;

cutoffs = [f_cutoff_1, f_cutoff_2];

cutoffs_matlab = cutoffs * 2;

h = fir1(filter_order, cutoffs_matlab, 'stop');

y = filter(h, 1, x);

[H, omega] = freqz(h, 1, 1024, 'whole');
figure,
plot(omega./(2*pi)*Fs, abs(H));


%to check the signal
% Y = fft(y);
figure;
stem(freq_axis, abs(Y), 'LineWidth', 2);

%% We want to increase the sampling rate from 80KHz to 120KHz. Implement the sampling rate conversion of the signal y(n). You choose all the parameters needed. Define the new obtained signal as w(n).

y_upsampled = zeros(1, length(y) * 120);
y_upsampled(1:80:end) = y;

cutoff = min([1/(2*120), 1/(2*80)]);

cutoff_matlab = 2*cutoff;

h = 120*fir1(64, cutoff_matlab);

y_f = filter(h, 1, y_upsampled);

w = y_f(1:80:end);
N = length(w);
Fs = 120*10^3;
freq_axis = 0:Fs/N:Fs - Fs/N;

figure
plot(freq_axis, w);

%%

N = length(w);
W = fft(w);

freq_axis = 0:1/N:1 -1/N;
figure;
plot(freq_axis, abs(W));
