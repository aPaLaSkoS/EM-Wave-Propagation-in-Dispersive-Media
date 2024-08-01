function F = CP(s,tau,a,b,c)

F = 1 ./ ( (tau*s).^a + (tau*s).^c ) .^ b ;

end

