%% 
clear; clc; close all;

% constants Havriliak - Negami
a = 0; b = 0.9; c = 0.8; tau = 7.234e-12;
e_inf = 1; e_s = 51; De = e_s-e_inf; 

eo = 8.854e-12; mu_o = 4*pi*1e-7;
co = 1 / sqrt(eo*mu_o);
fstart = 0.01e12; fend = 10e12; num_fr = 100000;
f = linspace(fstart,fend,num_fr);
w = 2*pi*f;

er = zeros(1,length(w));
for k=1:num_fr
   er(k)= e_inf + De * CP(1i*w(k),tau,a,b,c); 
end
root_er = sqrt(er);
Tth = 2 ./ (1+root_er);
Rth = (1-root_er) ./ (1+root_er) ;

% gamma = -i*w.*root_er/co;
% dz = 0.1e-3; d = 60*dz;
% Tdw = exp(d*gamma);

%% plots
figure
semilogx(f/1e12,real(er)); hold on;
semilogx(f/1e12,-imag(er));
xlabel('frequency (THz)')
ylabel('å_r(ù)')
legend('real', 'imagine')
title('Havriliak - Negami')

figure
semilogx(f/1e12,abs(Rth)); hold on;
semilogx(f/1e12,abs(Tth)); hold on;
legend('R', 'T')
title('Havriliak - Negami')
xlabel('frequency (THz)')

% figure
% semilogx(w/(2*pi),real(Tdw)); hold on;
% semilogx(w/(2*pi),imag(Tdw)); hold on;
% legend('real', 'imagine')
% title('Havriliak - Negami')
% ylabel('T(ù)')
% xlabel('frequency (GHz)')


