%% 
clear; clc; close all;

eo = 8.854e-12; mu_o = 4*pi*1e-7;
co = 1 / sqrt(eo*mu_o);
fstart = 1e10; fend = 1e13; num_fr = 10000;
w = linspace(2*pi*fstart,2*pi*fend,num_fr);

% dielectric
er = 9 * ones(1,num_fr);
root_er = sqrt(er);
Tth = 2 ./ (1+root_er);
Rth = (1-root_er) ./ (1+root_er) ;

gamma = -i*w.*root_er/co;
dz = 0.1e-3; d = 60*dz;
Tdw = exp(d*gamma);

%% plots
figure
semilogx(w/(2*pi),real(er)); hold on;
semilogx(w/(2*pi),-imag(er));
xlabel('frequency (Hz)')
ylabel('å_r(ù)')
legend('real', 'imagine')
title('Havriliak - Negami')

figure
semilogx(w/(2*pi),abs(Rth)); hold on;
semilogx(w/(2*pi),abs(Tth)); hold on;
legend('R', 'T')
title('Havriliak - Negami')
xlabel('frequency (Hz)')

% figure
% semilogx(w/(2*pi),real(Tdw)); hold on;
% semilogx(w/(2*pi),imag(Tdw)); hold on;
% legend('real', 'imagine')
% title('Havriliak - Negami')
% ylabel('T(ù)')
% xlabel('frequency (GHz)')


