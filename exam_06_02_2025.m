close all
clc
clear all

%%The signal x(n) contains two sinusoidal contributions (with the same amplitude = 1) at the normalized frequencies 0.1 and 0.25. The period of x(n) is 1.25 [ms] and the duration is 100 [ms]. Define the signal x(n).
A = 1;
fn_1 = 0.1;
fn_2 = 0.25;
t = 1.25 * 10^-3;
d = 100 * 10^3;

samples = lcm(1/fn_2, 1/fn_1);
Fs = samples/t;
tSample = 0:1/Fs:d -1/Fs;

x = A * cos(2*pi*fn_1*Fs*tSample) + A *cos(2*pi*fn_2*Fs*tSample);

figure;
plot(tSample,x);

%% Define the filter H(z) 
B = conv([0.9025, 0, 1], [1, -1.8*cos(pi/5), 0.81]);
A = [1, 0, 0.9025];

[H, omega] = freqz(B,a, 1024, 'whole');
figure,
plot(omega./(2*pi), abs(H));
title('|DTFT| of the filter H(z)');
grid;
xlabel('f [norm]');

y = filter(B,A,x);

figure,
plot(tSample, y);
grid;


y_0 = x(1)*B(1);

%% Compute the all-pass minimum-phase decomposition of the filter H(z), defining H_ap(z) and H_min(z) as the two components.

B_ap = [0.9025, 0, 1];
A_ap = [1, 0, 0.9025];

B_min = [1, -1.8*cos(pi/5), 0.81];
A_min = 1;

figure;
zplane(B_ap, A_ap);
title('Z-plane of filter H_{ap}(z)');
grid;
[H, omega] = freqz(B_ap, A_ap, 1024, 'whole');

figure,
plot(omega./(2*pi), abs(H));
title('|DTFT| of the filter H_{ap}(z)');
grid;
xlabel('f [norm]');

figure;
zplane(B_min, A_min);
title('Z-plane of filter H_{min}(z)');
grid;
[H, omega] = freqz(B_min, A_min, 1024, 'whole');

figure,
plot(omega./(2*pi), abs(H));
title('|DTFT| of the filter H_{min}(z)');
grid;
xlabel('f [norm]');

y_ap = filter(B_ap, A_ap, x);
y_min = filter(B_min, A_min, x);

% Define the signal w as the arithmetic mean between x(n) and y_ap(n)
w = 0.5*x + 0.5*y_ap;

%Find the filter H_w such that W(z) = X(z) * H_w(z)
% W(z) = (X(z) + X(z)*H_ap(z)) / 2 = X(z) (H_ap(z) + 1)/2
B_w = (B_ap + A_ap)*0.5;
A_w = A_ap;

figure;
zplane(B_w, A_w);
title('Z-plane of filter H_w(z)');
grid;
[H, omega] = freqz(B_w, A_w, 1024, 'whole');
figure,
plot(omega./(2*pi), abs(H));
title('|DTFT| of the filter H_{w}(z)'); 
grid;
 xlabel('f [norm]');

 %% Compute the DFTs of the signals x(n), y(n), y_min(n), y_ap(n), w(n)

X = fft(x);
Y = fft(y);
Y_min = fft(y_min);
Y_ap = fft(y_ap);
W = fft(w);

N = length(x);
freq_axis = 0:1/N:1 - 1/N;

figure;
stem(freq_axis, abs(X));
title('Absolute value of the DFT of the signal x(n)');
grid;
xlabel('f [norm]');

figure;
stem(freq_axis, abs(Y));
title('Absolute value of the DFT of the signal x(n)');
grid;
xlabel('f [norm]');

figure;
stem(freq_axis, abs(W));
title('Absolute value of the DFT of the signal x(n)');
grid;
xlabel('f [norm]');