% 1. FREE SPACE simulation -> get Ex_i(z1)
% 2. FREE SPACE + DIELECTRIC simulation -> get Ex_t(z1)

clear; clc; close all;

% general constants
eo = 8.854e-12; mu_o = 4*pi*1e-7; co = 1/sqrt(eo*mu_o);

% STEP 1:  FILT

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
[fprony,wHN,xHN] = myprony(f_FILT,tn,M);

% PLOTS
nPlot = 100; 
figure; 
plot(Dt*tn(1:nPlot),f_FILT(1:nPlot)); hold on ; 
plot(Dt*tn(1:nPlot),fprony(1:nPlot),'k'); 
legend('FILT',['FILT+Prony (M = ', num2str(M),')']);
title(['HN (á = ' ,num2str(a), ', â = ' ,num2str(b), ', ã = ' ,num2str(c), ')'] );
xlabel('Time Step Index'); ylabel('Impulse Response');

%% STEP 3: FDTD - 1D - Simulation of a pulse hitting a 
%                      Havriliak - Negami MEDIUM

% Initialization
ke = 7000;
Ex = zeros(1,ke); 
Hy = zeros(1,ke);
Dx = zeros(1,ke);
A = zeros(M,ke);
L = ones(1,ke);

% FDTD constants
dt = 1.67e-15;
dt_r = dt/Dt;

r = 0.5; 
dz = co*dt/r ; 

% calculate new wHN && xHN
wHN = wHN .* dt_r;
xHN = xHN .^ dt_r;

% input_wave constants
To = 0.167e-12;
ao = 0.0634e-12;
pos_inc = 3000; % position of incident wave

% MEDIUM starts at cell with number: k_start
k_start = 4000;
L(k_start:end) = e_inf + De * sum(wHN);

% MAIN FDTD Loop

% n -> TIME
% k -> SPACE
N = 6000;
z1 = k_start - 1;
Ex_tot = zeros(1,N);
for n = 1:N
    
    % update Dx
    for k=2:ke
       Dx(k) = Dx(k) + (co*dt/dz) * ( Hy(k-1)-Hy(k) ); 
    end
    
    % create a gaussian pulse 
    pulse = exp( - ( (n*dt-To)/ao ) .^2 );
    Dx(pos_inc) = pulse + Dx(pos_inc);
    
    for k = 2:ke 
       
       % update Ex
       Ex(k) = Dx(k) ;
       if k >= k_start
           for j=1:M
               Ex(k) = Ex(k) - De*wHN(j)*xHN(j)*A(j,k);
           end
       end
       Ex(k) = Ex(k) / L(k);
       
       % update A
       if k >= k_start
           for j = 1:M
               A(j,k) = xHN(j)*A(j,k) + Ex(k); 
           end
       end
       
    end
    
    % magnetic field - Hy
    for k=1:ke-1
       Hy(k) = Hy(k) + (co*dt/dz) * ( Ex(k)-Ex(k+1) ); 
    end 
    
    % calculate Total field 1-space-step BEFORE Raicu medium
    Ex_tot(n) = Ex(z1);
end

% PLOT wave propagation in FREE SPACE
fplot  = [zeros(1,k_start-1) max(Ex)*ones(1,ke-k_start+1)];
figure;  plot(Ex); 
hold on; plot(1:ke,fplot,'--');

%% STEP 4: FDTD - 1D - Simulation of a pulse propagating in FREE SPACE

% Initialization
ke = 7000;
Ex = zeros(1,ke); 
Hy = zeros(1,ke);

% FDTD constants
dt = 1.67e-15;
r = 0.5;
dz = co*dt/r ;

% input_wave constants
To = 0.167e-12;
ao = 0.0634e-12;
pos_inc = 3000; % position of incident wave

% MEDIUM starts at cell with number: k_start
k_start = 4000;
z1 = k_start - 1;

% MAIN FDTD Loop

% n -> TIME
% k -> SPACE
N = 6000;
Ex_i = zeros(1,N);
for n = 1:N
    
    % update Dx
    for k=2:ke
       Ex(k) = Ex(k) + (co*dt/dz) * ( Hy(k-1)-Hy(k) ); 
    end
    
    % create a gaussian pulse 
    pulse = exp( - ( (n*dt-To)/ao )^2 );
    Ex(pos_inc) = pulse + Ex(pos_inc);
    
    % update Hy
    for k=1:ke-1
       Hy(k) = Hy(k) + (co*dt/dz) * ( Ex(k)-Ex(k+1) ); 
    end 
    
    % calculate Incident field 1-space-step BEFORE Raicu medium
    Ex_i(n) = Ex(z1);
end

% PLOT wave propagation in FREE SPACE
fplot  = [zeros(1,k_start-1) max(Ex)*ones(1,ke-k_start+1)];
figure;  plot(Ex); 
hold on; plot(1:ke,fplot,'--');


%% calculate Reflection pulse wave
Ex_r = Ex_tot - Ex_i;

%% FFT
nfft = 2^17;
Ex_f_i = fftshift( abs(fft(Ex_i,nfft)) );
Ex_f_r = fftshift( abs(fft(Ex_r,nfft)) );
f = (1/dt/1e12) * (-nfft/2:nfft/2-1)/nfft; % f -> THz

% keep frequencies 0.01 - 10 THz
posf = find(f>=0.01 & f<=10);
Ex_f_i = Ex_f_i(posf);
Ex_f_r = Ex_f_r(posf);
f = f(posf);
% f = 1e9*f;

% PLOT Incident and Reflection waves
% in FREQUENCY domain
figure; semilogx(f,Ex_f_i); hold on; 
semilogx(f,Ex_f_r); hold off;
legend('Incident','Reflection')
xlabel('GHz')

% Reflection Coefficient
R_fdtd = zeros(1,length(posf));
for k=1:length(R_fdtd)
   R_fdtd(k) = Ex_f_r(k) / Ex_f_i(k);
end

% Theoretical Reflection Coefficient - HN medium 
clc;  
f = f*1e12;
wrad = 2*pi*f; num_fr = length(wrad);

er = zeros(1,length(wrad));
for k=1:num_fr
   er(k)= e_inf + De * CP(1i*wrad(k),tau,a,b,c); 
end
root_er = sqrt(er);
Rth = abs( (1-root_er) ./ (1+root_er) );

% PLOT Reflection Coefficients
figure
f = f/1e12;
semilogx(f, R_fdtd); hold on;
semilogx(f, Rth);
legend('R_f_d_t_d','R_t_h')
title('Havriliak - Negami')
xlabel('frequency (THz)')
% xlim([0.7 1.7])

