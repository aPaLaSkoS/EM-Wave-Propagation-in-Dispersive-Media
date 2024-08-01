clear; clc; close all

% FILT
% ...n = normalized
Dtn = 0.01;  
K = 21; alpha_filt=5;
tn = 0.02:Dtn:500;  % tn = t/tau

% dispersive medium  - > LOOK PAPER j47
tau = 153e-12; t = tau*tn; Dt = tau*Dtn;
tau_n = 1;
e_s = 1; e_inf = 0; 

% FILT
a = 0.3; b = 0.65; c = 0.8; p=8; 
f_FILT = FILT(tn,K,p,alpha_filt,tau_n,e_s,e_inf,a,b,c); 


%% PRONY 
N = length(f_FILT);  % f = [f(1) f(2) ... f(M)    f(M+1)... f(N)]
                %   = [f_0  f_1  ... f_(M-1) f_M   ... f_(N-1)]
% filter order
M = 4;
% prEnd -> number of points of f taken into account in "myprony"
prEnd = N;
f_FILT = f_FILT(1:prEnd);
tn = tn(1:prEnd); 

% STEP 1

% Solve for p coefficients -> Overdetermined Linear System
% Use of backslash operator <=> Least Squares Solution

% Create A - matrix  (Nt-M+1) x M
c = f_FILT(M+1:N);  
r = [f_FILT(N) f_FILT(M:-1:2)];
A = toeplitz(c,r);

% Create b - column Matrix (Nt-M+1) x 1
b = - f_FILT(1:N-M)' ;

% Solution
p = pinv(A)*b; % p = [p_M p_(M-1) ... p_1]
p = [p;1]; % p_0 = 1

% STEP 2

% Find the roots of Prony Polynomial
x = roots(p); % POLES of the IIR Filter

% STEP 3

X = zeros(N,M);
for k=1:N
   X(k,:) = x.^(k-1); 
end
w = pinv(X)*f_FILT'; % A coefficients of the IIR Filter

%% Re-create Impulse Response 
tn_prony = 0:1:length(tn);
fprony = zeros(1,length(tn));
for k=1:length(tn)
   mysum = 0;
   for q=1:M
       mysum = mysum + w(q)*(x(q)^tn_prony(k));%(tn(k)-tn(1)));
   end
    
   fprony(k) = mysum;
 
end
nPlot = 1000; 
% f_FILT = f_FILT/tau; fprony = fprony/tau;

figure; title('Havriliak-Negami');
xlabel('Time Step Index'); ylabel('Impulse Response');
plot(tn(1:nPlot),f_FILT(1:nPlot)); hold on ; 
plot(tn(1:nPlot),fprony(1:nPlot)); 
legend('FILT','Prony');

%% Theoretical solution
alpha = c-a; 
beta = b*c;
gama = b; 
x = tn/tau_n; z = - x.^alpha;
E = gml(z,alpha,beta,gama); 
g = (1/(tau_n))* x.^(beta-1); 
f_th = E.* g;

%% figures
nPlot = 1000; f_FILT = f_FILT/tau; f_th = f_th/tau;

figure; title('Havriliak-Negami');
xlabel('Time Step Index'); ylabel('Impulse Response');
plot(tn(1:nPlot),f_FILT(1:nPlot)); hold on ; 
plot(tn(1:nPlot),f_th(1:nPlot)); 
legend('FILT','Analytical');

%%
% % constants
% Dt = 1.768e-12; Nt = 15000; 
% tau = 153e-12;
% tau_n = 153e-12/Dt;

%% some comparisons
% t1 = t(1)
% te = t(end)
% 
% Dtnew = 1.768e-12; Nt = 100000; 
% tnew = Dtnew*(0.5:1:Nt+0.5);  
% tnew1 = tnew(1)
% tnew_e = tnew(end)