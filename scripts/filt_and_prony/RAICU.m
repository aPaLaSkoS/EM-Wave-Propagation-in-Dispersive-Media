clear; clc;

%%constants
Dt = 1.768e-12; tau = 153e-12/Dt;
De = 1; Nt = 50000;
K = 21; alpha=5;
tn = (0.5:1:Nt+0.5);  % tn = t/Dt

%% implementation of FILT to RAICU media

% Raicu a=0.1
a = 0.1; b = 0.9; c = 0.8; p=10;
fR1 = FILT(tn,K,p,alpha,tau,De,a,b,c);

% Raicu a=0.35
a = 0.35; b = 0.9; c = 0.8; p=10;
fR2 = FILT(tn,K,p,alpha,tau,De,a,b,c);

% Raicu a=0.6
a = 0.6; b = 0.9; c = 0.8; p=10;
fR3 = FILT(tn,K,p,alpha,tau,De,a,b,c);

% Raicu a=0.9
a = 0.9; b = 0.9; c = 0.8; p=10;
fR4 = FILT(tn,K,p,alpha,tau,De,a,b,c);

%% plots
Nplot = 1000;
figure
plot(tn(1:Nplot),fR1(1:Nplot)); hold on
plot(tn(1:Nplot),fR2(1:Nplot)); hold on
plot(tn(1:Nplot),fR3(1:Nplot)); hold on
plot(tn(1:Nplot),fR4(1:Nplot));
legend('Raicu (á=0.1)', 'Raicu (á=0.35)', 'Raicu (á=0.6)', 'Raicu (á=0.9)')
title(['FILT - Raicu model (â = ',num2str(b),' , ã = ' , num2str(c),')'])
xlabel('Time Step Index')
ylabel('Impulse Response')
