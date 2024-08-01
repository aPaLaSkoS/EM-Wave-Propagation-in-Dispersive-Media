
%% TEST ë(ù) for RAICU medium
clear; clc; close all
a = 0.3; b = 0.65; c = 0.8;
tau = 153e-12; e_inf = 2; es = 50; De = es-e_inf;

eo = 8.854e-12; mu_o = 4*pi*1e-7; co = 1/sqrt(eo*mu_o);
f = (0:0.0001:18)*1e9; wrad = 2*pi*f;


% Relative Permittivity å_r(ù)
eps_r = e_inf + De .* CP(i*wrad,tau,a,b,c) ;
rp = abs(sqrt(eps_r));

vel_medium = co ./ rp;   max_vel = max(vel_medium);
lamda = vel_medium ./ f;  
lmin = min(lamda);

% choose dz <= lmin/10
dz = lmin/15;
dt_max = dz/max_vel;

% choose dt <= dt_max
dt = dt_max / 2;

cfl = vel_medium * dt / dz;
max(cfl)



%% TEST Einc - function
clear; clc; close all;
dt = 0.4e-12;
fs = 1/dt;
spread = 1.26e11;
fe = 1e9;
N = 1000;
n = 1:N;
t = -N*dt:dt:N*dt;
% Einc = exp( -spread^2 * t.^2 ) .* sin( 2*pi*fe*t  );
Einc = exp( -spread^2 * t.^2 );
% Einc = exp( -spread^2 * (t - 40/spread).^2 ) .* sin( 2*pi*fe*(t - 40/spread) );
% Einc = exp( - (t/spread) .^2 );

figure
plot(t,Einc);
title('Incident Wave Pulse');
xlabel('time (sec)')
ylabel('E');
% tm = 0.001;
% xlim([-tm tm])

% FFT of Einc(t)

L=length(Einc);
NFFT = 2^17;
FT_Einc = fftshift(abs(fft(Einc,NFFT)));
f = (fs/1e9)*(-NFFT/2:NFFT/2-1)/NFFT; % Frequency Vector

figure
plot(f,abs(FT_Einc),'r'); hold on;
plot(f, ( max(FT_Einc)/exp(1) ) * ones(1,length(f)) );
title('FFT Magnitude ');
xlabel('Frequency (GHz)')
ylabel('|Einc(f)|');
fm = 100;
xlim([-fm fm])

figure
semilogx(f,abs(FT_Einc),'r'); hold on;
semilogx(f, ( max(FT_Einc)/exp(1) ) * ones(1,length(f)) );
title('FFT Magnitude ');
xlabel('Frequency (GHz)')
ylabel('|Einc(f)|');

