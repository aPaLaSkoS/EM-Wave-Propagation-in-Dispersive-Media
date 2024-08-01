%% Compare Time Execution of both FILT and GML
clear; clc; 
close all;

% constants
Num_of_ts = 10000:500:100000; 
TE_filt = zeros(1,length(Num_of_ts));
TE_gml = zeros(1,length(Num_of_ts));

% FILT
Dt = 1.768e-12; 

for m = 1:length(Num_of_ts)
Nt = Num_of_ts(m);
K = 21; alpha=5;
tn = (0.5:1:Nt+0.5);  % tn = t/Dt

% dispersive medium  - > LOOK PAPER j47
tau = 153e-12/Dt; tau_n = tau;
e_s = 1; e_inf = 0; 

% FILT

a = 0.1; b = 0.9; c = 0.8; p=8; 
tic
f_FILT = FILT(tn,K,p,alpha,tau,e_s,e_inf,a,b,c); 
TE_filt(m) = toc;

%% Theoretical solution
% Calculate 3-parameter Mittag-Leffler Function 
% using Roberto Girrappa's implementation of MLF

alpha = c-a; 
beta = b*c;
gama = b; 
x = tn/tau_n; z = - x.^alpha;
tic
E = gml(z,alpha,beta,gama); 
g = (1/(tau_n))* x.^(beta-1); 
f_th = E.* g;
TE_gml(m) = toc;

end


%% fit linear regression model

fitobject_TE_filt = fit(Num_of_ts',TE_filt','poly1');
p1_filt = fitobject_TE_filt.p1; p2_filt = fitobject_TE_filt.p2;
fitobject_TE_gml = fit(Num_of_ts',TE_gml','poly1');
p1_gml = fitobject_TE_gml.p1; p2_gml = fitobject_TE_gml.p2;

%% PLOTS

clc; close all;
figure 
plot(Num_of_ts,p1_filt*Num_of_ts+p2_filt,'b'); hold on;
plot(Num_of_ts,TE_filt,':r','LineWidth',1.5); hold on; 
plot(Num_of_ts,p1_gml*Num_of_ts+p2_gml,'k') ; hold on;
plot(Num_of_ts,TE_gml,':g','LineWidth',1.5); hold off;

legend(['T_e^F^I^L^T = ',num2str(p1_filt),'*N_t+',num2str(p2_filt)],'FILT' , ['T_e^G^M^L = ',num2str(p1_gml),'*N_t+',num2str(p2_gml)],'GML','Location','northwest');
title('Execution Time   vs   N_t');
xlabel('N_t (Number of Time Steps)'); ylabel('T_e (seconds)');
TEtot = sum(TE_filt+TE_gml)/60
