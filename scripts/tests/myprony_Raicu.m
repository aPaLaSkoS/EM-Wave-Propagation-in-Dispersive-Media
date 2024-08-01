%% implementation of FILT + PRONY (Raicu)
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
tau = 153e-12/Dt;
e_s = 1; e_inf = 0; 

% FILT

a = 0.1; b = 0.9; c = 0.8; p=10; 
f = FILT(tn,K,p,alpha,tau,e_s,e_inf,a,b,c); 
plot(tn(1:3000),f(1:3000)); title('Havriliak-Negami (a = 0 , b = 0.9 , c = 0.8 , N_t = 10^5)');
xlabel('Time Step Index'); ylabel('Impulse Response'); hold on;

%% PRONY 
% filter order
M = 7; % μεχρι Μ = 10 λειτουργει η prony εδω
[fprony,w,x] = myprony(f,tn,M);


%%
figure
plot(tn(1:1000),f(1:1000)); hold on
plot(tn(1:1000),fprony(1:1000)); hold off
legend('FILT' , ['FILT+Prony (Μ = ' , num2str(M),')']);
title({'Raicu (α = 0.1 , β = 0.9 , γ = 0.8)' ; ['Number of points included in prony: ', num2str(length(fprony))]});
xlabel('Time Step Index'); ylabel('Impulse Response');
