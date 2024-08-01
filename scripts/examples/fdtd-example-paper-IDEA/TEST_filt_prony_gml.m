%% STEP 1:  FILT
clear; clc; close all;

% dispersive medium constants
a = 0; c = 0.8; b = 0.9; 
e_inf = 1; e_s = 51; De = e_s-e_inf; 
tau = 7.234e-12;

% filt - constants
Dt = 1.67e-12; Nt = 4000; 
tn = (0.5:1:Nt+0.5);  % tn = t/Dt
K = 21; alpha=5; p=8; 
f_FILT = FILT(tn,K,p,alpha,tau/Dt,a,b,c); 

% STEP 2: PRONY 

M = 11; % filter order
[fprony,w,x] = myprony(f_FILT,tn,M);

% PLOTS
nPlot = 100; 
figure; 
plot(Dt*tn(1:nPlot),f_FILT(1:nPlot)); hold on ; 
plot(Dt*tn(1:nPlot),fprony(1:nPlot),'k'); 
legend('FILT',['FILT+Prony (M = ', num2str(M),')']);
title(['HN (á = ' ,num2str(a), ', â = ' ,num2str(b), ', ã = ' ,num2str(c), ')'] );
xlabel('Time Step Index'); ylabel('Impulse Response');

%%
% Nt_new = 400000; tnew = 0.5:1:Nt_new+0.5;
% dt = 1.67e-15;

% x_fdtd = x .^ (dt/Dt);
% fpr = zeros(1,length(tnew));
% for k=1:length(tnew)
%    for q=1:M
%        fpr(k) = fpr(k) + w(q)*(x_fdtd(q)^( tnew(k)-0.5 ) );
%    end
% end

dt = 1.67e-15;
fpr = zeros(1,length(tn));
for k=1:length(tn)
   mysum = 0;
   for q=1:M
       mysum = mysum + w(q)*( x(q)^(tn(k)-0.5) );
   end
    
   fpr(k) = mysum;
 
end

% PLOTS
nPlot = 100000; 
figure; 
% plot(tnew(1:nPlot),f_FILT(1:nPlot)); hold on ; 
plot(dt*tnew(1:nPlot),fpr(1:nPlot),'k'); 
% legend('FILT',['FILT+Prony (M = ', num2str(M),')']);
title(['HN (á = ' ,num2str(a), ', â = ' ,num2str(b), ', ã = ' ,num2str(c), ')'] );
xlabel('Time Step Index'); ylabel('Impulse Response');

%% Theoretical solution
alpha = c-a;
beta = b*c;
gama = b ;
tau_n = tau/Dt;
y = tn/tau_n; z = - y.^alpha;

E = gml(z,alpha,beta,gama); 
g = (1/(tau))* y.^(beta-1); 
f_th = E.* g;

% PLOTS
nPlot = 100; 
figure; 
plot(tn(1:nPlot),f_FILT(1:nPlot)); hold on ; 
plot(tn(1:nPlot),fprony(1:nPlot),'k'); 
plot(tn(1:nPlot),Dt*f_th(1:nPlot));
legend('FILT',['FILT+Prony (M = ', num2str(M),')'],'Analytical');
title(['HN (á = ' ,num2str(a), ', â = ' ,num2str(b), ', ã = ' ,num2str(c), ')'] );
xlabel('Time Step Index'); ylabel('Impulse Response');

figure; 
plot(tn(1:nPlot),f_th(1:nPlot));

