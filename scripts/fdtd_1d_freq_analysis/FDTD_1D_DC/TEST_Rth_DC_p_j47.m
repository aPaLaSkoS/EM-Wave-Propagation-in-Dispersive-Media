%% 
clear; clc; close all;

% constants Davidson - Cole -> p. j47
a = 0; b = 0.65; c = 1; tau = 153e-12;
e_inf = 2; e_s = 50; De = e_s-e_inf; 

eo = 8.854e-12; mu_o = 4*pi*1e-7;
co = 1 / sqrt(eo*mu_o);
f = (0:0.0001:15)*1e9; wrad = 2*pi*f;

er = zeros(1,length(wrad));
for k=1:length(f)
   er(k)= e_inf + De * CP(1i*wrad(k),tau,a,b,c); 
end
root_er = sqrt(er);
Tth = 2 ./ (1+root_er);
Rth = (1-root_er) ./ (1+root_er) ;

gamma = -i*wrad.*root_er/co;
dz = 0.1e-3; d = 60*dz;
Tdw = exp(d*gamma);

%% plots
% figure
% semilogx(w/(2*pi),real(er)); hold on;
% semilogx(w/(2*pi),-imag(er));
% xlabel('å_r(ù)')
% legend('real', 'imagine')
% title('Havriliak - Negami')

figure
semilogx(f/1e9,abs(Rth)); hold on;
semilogx(f/1e9,abs(Tth)); hold on;
legend('R', 'T')
title('Davidson - Cole')
xlabel('frequency (GHz)')

% figure
% semilogx(w/(2*pi),real(Tdw)); hold on;
% semilogx(w/(2*pi),imag(Tdw)); hold on;
% legend('real', 'imagine')
% title('Havriliak - Negami')
% ylabel('T(ù)')
% xlabel('frequency (GHz)')


