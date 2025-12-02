close all
clc
clear all

f1 = 500;
f2_n = 0.025;
n = 0:999;
period_n = 680;

P2_samples = 1/f2_n;
P1_samples = period_n/P2_samples;
f1_n= 1/P1_samples;
Fs = f1/f1_n;

x1 = cos(2*pi*f1_n*n);
x2 = cos(2*pi*f2_n*n);

delta = zeros(size(n));
delta(6)= 1;

x2 = conv(x2, delta);
n_conv = 0 : 2* n(end);
x2 = x2(n_conv <= n(end));
x = x1 + x2;
figure;
plot(n, x);
grid;

%%

Fs_new = 12750;
[L, M] = rat(Fs_new / Fs);

x_up = zeros(1, length(x) * L);
x_up(1:L:end) = x;

cutoff = min([1/(2*L), 1/(2*M)]);
cutoff_filter = 2*cutoff;
h = L*fir1(64, cutoff_filter);
x_f = filter(h, 1, x_up);

y = x_f(1:M:end);

N = 2048;
X = fft(x, N);
Y = fft(y, N);
freq_axis = 0:1/N:1 -1/N;

figure;
plot(freq_axis, abs(X), 'r');
hold on;
plot(freq_axis, abs(Y),'b');
grid on;

%% 
theta = 2*pi*f2_n/L*M;
rho = 0.95; 
a = 2;
b = -2*2*cos(theta);
c = 2;
d = 2*rho*cos(theta);
e = - rho^2;
B = [a, b, c];
A = [1, -d, -e];

[H, omega] = freqz(B, A, N, 'whole');

figure,
plot(omega./(2*pi), abs(H));
grid;
%% 

b = 1;
a = 1/2; 
B = [a, b];
A = [b, a];


[H, omega] = freqz(B, A, N, 'whole');
figure,
plot(omega./(2*pi), abs(H));
grid;

figure
plot(omega./(2*pi), angle(H));
grid;
