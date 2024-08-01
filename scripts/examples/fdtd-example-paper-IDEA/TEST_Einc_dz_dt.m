
%% TEST ë(ù) for HN medium -> paper IDEA...
clear; clc; close all;

% constants Havriliak - Negami
a = 0; b = 0.9; c = 0.8; tau = 7.234e-12;
e_inf = 1; e_s = 51; De = e_s-e_inf; 

eo = 8.854e-12; mu_o = 4*pi*1e-7;
co = 1 / sqrt(eo*mu_o);
fstart = 1e10; fend = 1e13; num_fr = 100000;

f = linspace(fstart,fend,num_fr);
wrad = 2*pi*f;

% Relative Permittivity å_r(ù)
eps_r = e_inf + De .* CP(i*wrad,tau,a,b,c) ;
rp = abs(sqrt(eps_r));

vel_medium = co ./ rp;   max_vel = max(vel_medium);
lamda = vel_medium ./ f;  
lmin = min(lamda);

% choose dz <= lmin/10
dz_max = lmin/10;
dz = 1e-6;
dt_max = dz/max_vel;

% choose dt <= dt_max
dt = 1.67e-15;

cfl = max_vel * dt / dz




%% TEST Einc - function
clear; clc; close all
dt = 1.67e-15;
fs = 1/dt;
To = 0.167e-12;
% To = dt*200
ao = 0.0634e-12;
N = 100;
n = 1:N;
t = -N*dt:dt:N*dt;
fc = 5000e9;
Einc = exp( - ( (t)/ao ) .^2 );

figure
plot(t*1e12,Einc);
xlabel('ps'); ylabel('E_i_n_c(t)');

% FFT of Einc(t)

L=length(Einc);
NFFT = 2^17;
FT_Einc = abs( fftshift(fft(Einc,NFFT)) );
f = (fs/1e12)*(-NFFT/2:NFFT/2-1)/NFFT; % Frequency Vector
e = ones(1,length(f)) * max(FT_Einc)/exp(1);

figure
semilogx(f,FT_Einc,'r'); hold on;
semilogx(f,e); hold off;
title('Magnitude of FFT');
xlabel('Frequency (THz)')
ylabel('Magnitude |Einc(f)|');
xlim([-30 30])
