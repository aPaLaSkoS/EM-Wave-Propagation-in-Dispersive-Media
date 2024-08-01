%% implementation of FILT + PRONY (Havriliak - Negami)
clear; clc; 
% close all;

% constants

% FILT
Dt = 1.768e-12; Nt = 50000; 

K = 21; alpha=5;
tn = (0.5:1:Nt+0.5);  % tn = t/Dt

% dispersive medium  - > LOOK PAPER j47
tau = 153e-12/Dt;
e_s = 1; e_inf = 0; 

% FILT

a = 0; b = 0.9; c = 0.8; p=8; 
f = FILT(tn,K,p,alpha,tau,e_s,e_inf,a,b,c); 
% figure; plot(tn,f); title('Havriliak-Negami (a = 0 , b = 0.9 , c = 0.8 , N_t = 10^5)');
% xlabel('Time Step Index'); ylabel('Impulse Response');

%% PRONY 
N = length(f);  % f = [f(1) f(2) ... f(M)    f(M+1)... f(N)]
                %   = [f_0  f_1  ... f_(M-1) f_M   ... f_(N-1)]
% filter order
M = 8;
% prEnd -> number of points of f taken into account in "myprony"
prEnd = N;
f = f(1:prEnd);
tn = tn(1:prEnd); 

%% STEP 1

% Solve for p coefficients -> Overdetermined Linear System
% Use of backslash operator <=> Least Squares Solution

% Create A - matrix  (Nt-M+1) x M
c = f(M+1:N);  
r = [f(N) f(M:-1:2)];
A = toeplitz(c,r);

% Create b - column Matrix (Nt-M+1) x 1
b = - f(1:N-M)' ;

% Solution
p = pinv(A)*b; % p = [p_M p_(M-1) ... p_1]
p = [p;1]; % p_0 = 1

% STEP 2

% Find the roots of Prony Polynomial
x = roots(p); % POLES of the IIR Filter

%% STEP 3

X = zeros(N,M);
for k=1:N
   X(k,:) = x.^(k-1); 
end
w = pinv(X)*f'; % A coefficients of the IIR Filter

%% Re-create Impulse Response 
fprony = zeros(1,length(tn));
for k=1:length(tn)
   mysum = 0;
   for q=1:M
       mysum = mysum + w(q)*(x(q)^(tn(k)-0.5));
   end
    
   fprony(k) = mysum;
 
end

%%
figure
plot(tn(1:1000),f(1:1000)); hold on
plot(tn(1:1000),fprony(1:1000)); hold off
title({'Havriliak-Negami' ; ['Number of points included in prony: ', num2str(length(fprony))]});
xlabel('Time Step Index'); ylabel('Impulse Response');
legend('FILT' , 'FILT+Prony');

