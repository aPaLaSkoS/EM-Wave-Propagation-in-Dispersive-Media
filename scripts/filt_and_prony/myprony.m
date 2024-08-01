function [fprony,w,x] = myprony(f,tn,M)

%========================= INPUTS =========================================
% M  -> Filter Order
% f  -> values of the function f, on which Prony is applied (see p.38 -
%       eq.(3.3))
% tn -> Normalized Time (t/Dt), where Dt: the time step used for sampling f       
%==========================================================================

%========================= OUTPUTS ========================================
% fprony  -> sum of decaying exponentials approximating the function f
% x       -> exponentials of the Prony method (see p.39 eq.(3.7) calculated
%            in STEP 2
% w       -> weights of the Prony method (see p.38 eq.(3.3) and p.43 eq.
%           (3.20)) calculated in STEP 3
%=========================================================================

% NOTE: here N corresponds to Nt in pdf
N = length(f);  % f = [f(1) f(2) ... f(M)    f(M+1)... f(N)]
                %   = [f_0  f_1  ... f_(M-1) f_M   ... f_(N-1)]

% prEnd -> number of points of f taken into account in "myprony"
prEnd = N;

f = f(1:prEnd);
tn = tn(1:prEnd);

%% STEP 1

% Solve for p coefficients -> Overdetermined Linear System
% Use of backslash operator <=> Least Squares Solution

% Create A - matrix  (N-M) x M
c = f(M+1:N);  
r = [f(N) f(M:-1:2)];
A = toeplitz(c,r);

% Create b - column Matrix (N-M) x 1
b = - f(1:N-M)' ;

% Solution
p = pinv(A)*b; % p = [p_M p_(M-1) ... p_1]
p = [p;1]; % p_0 = 1

%% STEP 2

% Find the roots of Prony Polynomial
x = roots(p); % POLES of the IIR Filter

%% STEP 3

X = zeros(N,M);
for k=1:N
   X(k,:) = x.^(k-1); 
end
w = pinv(X)*f'; % A coefficients of the IIR Filter

%% Create Impulse Response 

fprony = zeros(1,length(tn));
for k=1:length(tn)
   mysum = 0;
   for q=1:M
       mysum = mysum + w(q)*(x(q)^(tn(k)-0.5));
   end
    
   fprony(k) = mysum;
 
end


end

