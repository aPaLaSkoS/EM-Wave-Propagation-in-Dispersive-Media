
%% TEST ë(ù) for HN medium -> paper IDEA...
clear; clc; close all;

eo = 8.854e-12; mu_o = 4*pi*1e-7;
co = 1 / sqrt(eo*mu_o);
f = (0:0.0001:15)*1e9; wrad = 2*pi*f;

% dielectric
eps_r = 9;
vel_air = co;  max_vel = vel_air;
vel_d = co ./ eps_r;   
  
lmin = co / eps_r / f(end);

% choose dz <= lmin/10
dz_max = lmin/10;
dz = dz_max;
dt_max = dz/max_vel;

% choose dt <= dt_max
dt = dt_max/2;

cfl = max_vel * dt / dz




%% TEST Einc - function

fs = 1/dt;
spread = 1.26e10;
fe = 6e9;
N = 1000;
n = 1:N;
t = -N*dt:dt:N*dt;
Einc = exp( -spread^2 * t .^2 ) .* sin( 2*pi*fe*t  );

figure
plot(t,Einc);
title('Incident Wave Pulse');
xlabel('time (sec)')
ylabel('E');

%% FFT of Einc(t)

L=length(Einc);
NFFT = 2^13;
FT_Einc = fftshift(fft(Einc,NFFT));
f = (fs/1e9)*(-NFFT/2:NFFT/2-1)/NFFT; % Frequency Vector

figure
semilogx(f,abs(FT_Einc),'r');
plot(f,abs(FT_Einc).^2,'r');
title('Magnitude of FFT');
xlabel('Frequency (GHz)')
ylabel('Magnitude |Einc(f)|');
xlim([-20 20])
