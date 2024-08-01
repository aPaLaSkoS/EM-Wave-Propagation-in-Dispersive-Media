%% implementation of FILT + PRONY (Havriliak - Negami)
clear; clc; 
close all;

% constants

% FILT
% Dt = 1.768e-12; Nt = 50000;
Dt = 1.768e-12; Nt = 50000;
K = 21; alpha=5;
tn = (0.5:1:Nt+0.5);  % tn = t/Dt

% dispersive medium  - > LOOK PAPER j47
tau = 153e-12/Dt;
e_s = 1; e_inf = 0; 

% FILT

clc;
a = 0.3; b = 0.65; c = 0.8; p=8; 
x = FILT(tn,K,p,alpha,tau,e_s,e_inf,a,b,c); 
% figure; plot(tn,x); title('Havriliak-Negami (a = 0 , b = 0.9 , c = 0.8 , N_t = 10^5)');
% xlabel('Time Step Index'); ylabel('Impulse Response');

%% PRONY 
Nt = length(x);  % f = [f(1) f(2) ... f(M)    f(M+1)... f(Nt)]
                 %   = [f_0  f_1  ... f_(M-1) f_M   ... f_(Nt-1)]
% filter order
M = 6;
% prEnd -> number of points of f taken into account in "myprony"
prEnd = Nt;
x = x(1:prEnd);
tn = tn(1:prEnd);

% STEP 1

% Solve for p coefficients -> Overdetermined Linear System
% Use of backslash operator <=> Least Squares Solution

% Create A - matrix  (prEnd-M) x M
c = x(M:prEnd-1);  
r = [x(prEnd-1) x(M-1:-1:1)];
X = toeplitz(c,r);

% Anew = zeros(prEnd-M,M);
% for k=0:prEnd-M-1
%    Anew(k+1,:) = flip(f(k+1:k+M));  
% end

% Create b - column Matrix (prEnd-M) x 1
B = - x(M+1:prEnd)' ;

% Solution
gamma = pinv(X)*B; 
% gamma2 = lsqminnorm(X,B);

% STEP 2
gamma = [1;gamma];
% Find the roots of Prony Polynomial
mu = roots(gamma); % POLES of the IIR Filter
% mu = flip(mu);

%% STEP 3

Mar = zeros(prEnd,M);
for k=1:prEnd
   Mar(k,:) = mu.^(k-1); 
end
avec = pinv(Mar)*x(1:prEnd)'; % A coefficients of the IIR Filter
% avec2 = lsqminnorm(Mar,x');


%% Re-create Impulse Response 
fprony = zeros(1,length(tn));
for k=1:length(tn)
   mysum = 0;
   for q=1:M
       mysum = mysum + avec(q)*(mu(q)^(tn(k)-0.5));
   end
    
   fprony(k) = mysum;
 
end

%
figure
plot(tn(1:1000),x(1:1000)); hold on
plot(tn(1:1000),fprony(1:1000)); hold off
title({'Havriliak-Negami' ; ['Number of points included in prony: ', num2str(prEnd)]});
xlabel('Time Step Index'); ylabel('Impulse Response');
legend('Original' , 'Prony');

