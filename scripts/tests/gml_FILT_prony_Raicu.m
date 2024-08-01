%% implementation of FILT ,  PRONY , GML
clear; clc; 
close all;

% constants

% FILT
% Dt = 1.768e-12; Nt = 50000; % a = 0.1 (M=7) ,0.2 (M=6)
% Dt = 80e-12; Nt = 50000; % a = 0.3 (M=7)
% Dt = 3620e-12; Nt = 50000; % a = 0.4 (M=6)  Dt = 80e-12 
Dt = 1.768e-12; Nt = 50000; 

K = 21; alpha=5;
tn = (0.5:1:Nt+0.5);  % tn = t/Dt

% dispersive medium  - > LOOK PAPER j47
tau = 153e-12/Dt; tau_n = tau;
e_s = 1; e_inf = 0; 

% FILT

a = 0.1; b = 0.9; c = 0.8; p=8; 
tic
f_FILT = FILT(tn,K,p,alpha,tau,e_s,e_inf,a,b,c); 
toc

%% PRONY 

M = 7; % filter order
[fprony,w,poles] = myprony(f_FILT,tn,M);


%% Theoretical solution
alpha = c-a; 
beta = b*c;
gama = b; 
x = tn/tau_n; z = - x.^alpha;
tic
E = gml(z,alpha,beta,gama); 
g = (1/(tau_n))* x.^(beta-1); 
f_th = E.* g;
toc
%% PLOTS
nPlot = 1000; 
Dt = 1;
figure; 
plot(tn(1:nPlot),f_FILT(1:nPlot)/Dt); hold on ; 
plot(tn(1:nPlot),f_th(1:nPlot)/Dt); hold on;
% plot(tn(1:nPlot),fprony(1:nPlot)/Dt,'k'); hold off;
% legend('FILT','Theoretical',['FILT+Prony (M = ', num2str(M),')']);
legend('FILT','Theoretical');
title(['Raicu (á = ' ,num2str(a), ', â = ' ,num2str(b), ', ã = ' ,num2str(c), ')'] );
xlabel('Time Step Index'); ylabel('Impulse Response');

