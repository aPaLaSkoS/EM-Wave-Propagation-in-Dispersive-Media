clear; clc; close all;

% general constants
eo = 8.854e-12; mu_o = 4*pi*1e-7; co = 1/sqrt(eo*mu_o);

% STEP 1:  FILT

% dispersive medium constants
a = 0; b = 0.65; c = 1;
e_inf =2; e_s = 50; De = e_s-e_inf; 
tau = 153e-12;

% filt - constants
% Dt = 1.768e-12; 
Dt = 0.4e-12;
Nt = 50000;
tn = (0.5:1:Nt+0.5);  % tn = t/Dt
K = 21; alpha=5; p=8; 
f_FILT = FILT(tn,K,p,alpha,tau/Dt,e_s,e_inf,a,b,c); 

% STEP 2: PRONY 

M = 9; % filter order
[fprony,wDC,xDC] = myprony(f_FILT,tn,M);


% PLOTS
nPlot = 1000; 
figure; 
plot(tn(1:nPlot),f_FILT(1:nPlot)); hold on ; 
plot(tn(1:nPlot),fprony(1:nPlot),'k'); hold off;
legend('FILT',['FILT+Prony (M = ', num2str(M),')']);
title(['DC (á = ' ,num2str(a), ', â = ' ,num2str(b), ', ã = ' ,num2str(c), ')'] );
xlabel('Time Step Index'); ylabel('Impulse Response');

%%  STABILITY
clc;
lmin = 0.0064;
r = 15;

% Davidson - Cole
L = e_inf + De * sum(wDC);
B = - De * sum( xDC .* wDC ) / L 

J = 0:0.01:1;
R = sqrt( - L * (J-B) .* (J-1) / 4 ) ;

maxR = max(R)
minR = min(R);

dz = lmin/r;

dt_max = dz*maxR / (sind(180/r)*co)

% free space
Lfp = 1;
Bfp = 0; 

Jfp = 0:0.01:1;
Rfp = sqrt( - Lfp * (Jfp-Bfp) .* (Jfp-1) / 4 ) ;

maxRfp = max(Rfp)
minRfp = min(Rfp);

dz_fp = lmin/r;

dt_max_fp = dz_fp * maxRfp / (sind(180/r)*co)


