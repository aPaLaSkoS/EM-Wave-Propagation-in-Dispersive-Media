function f = FILT(t,K,p,alpha,tau,e_s,e_inf,a,b,c)

Apq= [zeros(1,p-1) 1 ];

% calculate Apq coefficients recursively
for q = 1 : p-1
    if q==1 
        Cpq = 1; Bpq = p;
    else
        Cpq = (p-q+1)/q;
        Bpq = Cpq * Bpq;
    end
   Apq(p-q) =  Apq(p-q+1) +Bpq;
end

% calculate time-domain response
f = zeros(1,length(t));
for k=1:length(t)
    Fn = 0;
    for n=1:K
       Fn = Fn + (-1)^n * imag( CP( (alpha+i*pi*(n-0.5))/t(k) , tau, a, b, c ) ) ;  
    end
    A=0;
    for q=1:p
        A = A + Apq(q) * (-1)^(K+q) * imag( CP( (alpha+i*pi*(K+q-0.5))/t(k) , tau, a, b, c ) ) ;
    end
    
f(k) = (exp(alpha)/t(k)) * (Fn + 2^(-p)*A);
    
end

end

