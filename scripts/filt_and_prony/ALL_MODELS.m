clear; clc; close all;

%%constants
Dt = 1.768e-12; tau = 153e-12/Dt;
De = 1; Nt = 15000;
K = 21; alpha=5;
tn = (0.5:1:Nt+0.5);  % tn = t/Dt

%% implementation of FILT to various media

% Debye
a = 0; b = 1; c = 1; p=5;
fD = FILT(tn,K,p,alpha,tau,De,a,b,c);

% Cole - Cole
a = 0; b = 1; c = 0.8; p=6;
fCC = FILT(tn,K,p,alpha,tau,De,a,b,c);

% Davidson - Cole
a = 0; b = 0.9; c = 1; p=7;
fDC = FILT(tn,K,p,alpha,tau,De,a,b,c);

% Havriliak - Negami
a = 0; b = 0.9; c = 0.8; p=8;
fHN = FILT(tn,K,p,alpha,tau,De,a,b,c); 

% Raicu a=0.1
a = 0.1; b = 0.9; c = 0.8; p=10;
fR1 = FILT(tn,K,p,alpha,tau,De,a,b,c);

% % Raicu a=0.5
% a = 0.5; b = 0.9; c = 0.8; p=12;
% fR2 = FILT(tn,K,p,alpha,tau,De,a,b,c);

% Raicu a=0.9
a = 0.9; b = 0.9; c = 0.8; p=10;
fR3 = FILT(tn,K,p,alpha,tau,De,a,b,c);

%% plots
Nplot = 200;
figure
plot(tn(1:Nplot),fD(1:Nplot)); hold on
plot(tn(1:Nplot),fCC(1:Nplot)); hold on
plot(tn(1:Nplot),fDC(1:Nplot)); hold on
plot(tn(1:Nplot),fHN(1:Nplot)); hold on
% plot(tn(1:Nplot),fR1(1:Nplot)); hold on
% plot(tn(1:1000),fR2(1:1000)); hold on
% plot(tn(1:Nplot),fR3(1:Nplot)); hold on
% legend('Davidson-Cole' , 'Havriliak-Negami' , 'Raicu (á = 0.1)' , 'Raicu (á = 0.9)')
legend('Debye' , 'Cole-Cole' , 'Davidson-Cole' , 'Havriliak-Negami')
title('FILT - Various models')
xlabel('Time Step Index')
ylabel('Impulse Response')