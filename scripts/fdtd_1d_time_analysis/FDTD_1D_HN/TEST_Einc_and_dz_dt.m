
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
clear; clc; close all

Dt = 0.4e-12;  
fsD = 1/Dt;
spread = 1.26e10;
fe = 6e10;
ND = 700;
nD = 1:ND;
tD = -ND*Dt:Dt:ND*Dt;
Einc_D = exp( -spread^2 * tD .^2 ) .* sin( 2*pi * fe * tD  );
figure
plot(tD,Einc_D);
title('Incident Wave Pulse');
xlabel('time (sec)')
ylabel('E');


% dt = 15.3e-12; fsd = 1/dt;
% rt = dt/Dt;
% spread_pr = spread/rt;
% fe_pr = fe / rt ; 
% Nd = 700;
% nd = 1:Nd;
% td = -Nd*dt:dt:Nd*dt;
% Einc_d = exp( -(spread_pr * td) .^2 ) .* sin( 2*pi * fe_pr * td  );
% figure
% plot(td,Einc_d);
% title('Incident Wave Pulse');
% xlabel('time (sec)')
% ylabel('E');

% FFT of Einc(t)
LD=length(Einc_D);
NFFT_D = 2^17;
FT_Ein_D = fftshift(fft(Einc_D,NFFT_D));
fD = (fsD/1e9)*(-NFFT_D/2:NFFT_D/2-1)/NFFT_D; % Frequency Vector

figure
semilogx(fD,abs(FT_Ein_D),'r'); hold on;
semilogx(fD, ( max(abs(FT_Ein_D))/exp(1) ) * ones(1,length(fD)) );
title('FFT Magnitude ');
xlabel('Frequency (GHz)')
ylabel('|Einc(f)|');
% xlim([0.1 20])

% Ld=length(Einc_d);
% NFFT_d = 2^17;
% FT_Ein_d = fftshift(fft(Einc_d,NFFT_d));
% fd = (fsd/1e9)*(-NFFT_d/2:NFFT_d/2-1)/NFFT_d; % Frequency Vector
% 
% figure
% semilogx(fd,abs(FT_Ein_d),'r'); hold on;
% semilogx(fd, ( max(abs(FT_Ein_d))/exp(1) ) * ones(1,length(fd)) );
% title('FFT Magnitude ');
% xlabel('Frequency (GHz)')
% ylabel('|Einc(f)|');
% % xlim([0.1 20])
