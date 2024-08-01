clear; clc; close all

% constants

% FILT
Dt = 0.01;  

K = 21; alpha=5;
tn = 0.02:Dt:1000;  % tn = t/tau

% dispersive medium  - > LOOK PAPER j47
tau = 1;
e_s = 1; e_inf = 0; 

% FILT
a = 0.4; b = 0.9; c = 0.65; p=8; 
tic
f_FILT = FILT(tn,K,p,alpha,tau,e_s,e_inf,a,b,c); 
toc

% GML TEST CODE
alpha = c-a; 
beta = b*c;
gama = b;
x = tn/tau; z = - x.^alpha;

tic
E = ml(z,alpha,beta,gama); 
g = (1/(tau))* x.^(beta-1); 
f_th = E.* g;
toc

figure; title('Havriliak-Negami');
xlabel('Time Step Index'); ylabel('Impulse Response');
plot(tn,f_FILT); hold on ; plot(tn,f_th); 
legend('FILT','Analytical');
f_th(end)