clear; clc; close all

%%constants
Dt = 1.67*10^(-15); tau = 7.234*10^(-12)/Dt;
De = 50; Nt = 4000;
K = 21; alpha=5;
tn = (0.5:1:Nt+0.5);  % tn = t/Dt

% Raicu a=0.5
pf=9; pl=9;
pVec = pf:pl; Tr = zeros(1,length(pVec)); Tbf = zeros(1,length(pVec));
a = 0.5; b = 0.9; c = 0.8; 



for p=pVec(1):pVec(end)
    
% 1st option -> RECURSION

tic;
Apq_r= [zeros(1,p-1) 1 ];

for q = 1 : p-1
    if q==1 
        Cpq = 1; Bpq = p;
    else
        Cpq = (p-q+1)/q;
        Bpq = Cpq * Bpq;
    end
   Apq_r(p-q) =  Apq_r(p-q+1) +Bpq;
end

% calculate time-domain response
f_r = zeros(1,length(tn));
for k=1:length(tn)
    Fn = 0;
    for n=1:K
       Fn = Fn + (-1)^n * imag( CP( (alpha+i*pi*(n-0.5))/tn(k) , tau, De, a, b, c ) ) ;  
    end
    A_r=0;
    for q=1:p
        A_r = A_r + Apq_r(q) * (-1)^(K+q) * imag( CP( (alpha+i*pi*(K+q-0.5))/tn(k) , tau, De, a, b, c ) ) ;
    end
    
f_r(k) = (exp(alpha)/tn(k)) * (Fn + 2^(-p)*A_r);
    
end
Tr(p-pf+1) = toc;


% 2nd option -> BRUTE FORCE

tic;
fbf = zeros(1,length(tn));
for k=1:length(tn)
    
    Fn = 0;
    for n=1:K
       Fn = Fn + (-1)^n * imag( CP( (alpha+i*pi*(n-0.5))/tn(k) , tau, De, a, b, c ) ) ;  
    end
    
    sum = 0;
    for n = 0:p-1
        F = 0;
        for w = 0:n
           x = n-w+K+1;
           F = F + nchoosek(n,w) *  (-1)^x * imag( CP( (alpha+i*pi*(x-0.5))/tn(k) , tau, De, a, b, c ) );
        end
        sum = sum + 2^(-n-1) * F ; 
    end
    
fbf(k) = (exp(alpha)/tn(k)) * (Fn + sum);
    
end 
Tbf(p-pf+1) = toc;


end


%% PLOTS
if length(pVec)==1  % μια τιμη του p για ελεγχο ομοιοτητας των 2 μεθοδων επιλυσης
    figure 
    plot(tn(1:1000),f_r(1:1000)); hold on
    plot(tn(1:30:1000),fbf(1:30:1000),'o')
    legend('Recursion' , 'Brute-Force');
    xlabel('Time Step Index'); ylabel('Impulse Response'); title('FILT - Raicu model(α = 0.5 p = 9)');

else  % για διαφορες τιμες του p, συγκρινω το Execution Time των 2 αλγοριθμων
    figure
    plot(pVec,Tr); hold on
    plot(pVec,Tbf);
    legend('Recursion' , 'Brute-Force', 'Location','northwest');
    xlabel('p'); ylabel('Execution Time (secs)'); title('FILT Implementation');
end