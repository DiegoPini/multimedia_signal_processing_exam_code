close all
clc
clear all

Ts = 62.5e-6;
Fs = 1/ Ts;
T = 50e-3;
n = 0: 1/Fs: T - 1/Fs;
x = cos(2 * pi* 0.125 * n) + cos(2 * pi* 0.25 * n);

%%
Fs_new = 12000;
[L, M] = rat(Fs_new / Fs);

x_up = zeros(1, length(x) * L);
x_up(1:L:end) = x;

cutoff = min([1/(2*L), 1/(2*M)]);
cutoff_filter = 2*cutoff;
h_multirate = L*fir1(128, cutoff_filter);
x_f = filter(h_multirate, 1, x_up);

y = x_f(1:M:end);
z = x_up(1:M:end);

N = length(y);
Y = fft(y);
Z = fft(z);
freq_axis = 0:1/N:1 -1/N;
figure;
plot(freq_axis, abs(Y));
grid;


figure;
plot(freq_axis, abs(Z));
grid;


%%

%fir
transition_band = 1/12;
cut_off_left = 1/6 - 1/12;
cut_off_right = 1/3 + 1/12;

h = fir1(64, [cut_off_left, cut_off_right]*2);

w = filter(h, 1, z);
W = fft(w);
figure;
plot(freq_axis, abs(W));
grid;

%notch
zeroes = [1, -1];
poles = [0.9, -0.9];
A = poly(poles);
B = poly(zeroes);

w = filter(B, A, z);
W = fft(w);

figure;
plot(freq_axis, abs(W));
grid;
