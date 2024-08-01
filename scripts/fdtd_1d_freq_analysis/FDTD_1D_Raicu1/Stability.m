clear; clc; close all;

% general constants
eo = 8.854e-12; mu_o = 4*pi*1e-7; co = 1/sqrt(eo*mu_o);

% STEP 1:  FILT

% dispersive medium constants (RAICU MODEL)
a = 0.3; b = 0.65; c = 0.8;
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

M = 10; % filter order
[fprony,wR,xR] = myprony(f_FILT,tn,M);


% PLOTS
nPlot = 1000; 
figure; 
plot(tn(1:nPlot),f_FILT(1:nPlot)); hold on ; 
plot(tn(1:nPlot),fprony(1:nPlot),'k'); hold off;
legend('FILT',['FILT+Prony (M = ', num2str(M),')']);
title(['DC (� = ' ,num2str(a), ', � = ' ,num2str(b), ', � = ' ,num2str(c), ')'] );
xlabel('Time Step Index'); ylabel('Impulse Response');

%%  STABILITY
clc;
r = 15;
J = 0.001:0.001:1;

% Raicu
lmin = 0.0049;

L = e_inf + De * sum(wR);
B = - De * sum( xR .* wR ) / L; 

R = sqrt( - L * (J-B) .* (J-1) / (4*J) ) ;

maxR = max(R);
minR = min(R);

dz = lmin/r

dt_max = dz*maxR / (sind(180/r)*co)

%% free space
f = (0:0.0001:18)*1e9;
lmin_fp = min(co./f);

Lfp = 1;
Bfp = 0; 

Rfp = sqrt( - Lfp * (J-Bfp) .* (J-1) / (4*J) ) ;

maxRfp = max(Rfp);
minRfp = min(Rfp);

dz_fp = lmin_fp/r

dt_max_fp = dz_fp * maxRfp / (sind(180/r)*co)