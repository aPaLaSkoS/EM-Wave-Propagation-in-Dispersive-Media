% 1. FREE SPACE simulation -> get Ex_i(z1)
% 2. FREE SPACE + DIELECTRIC simulation -> get Ex_t(z1)

clear; clc; close all;

% general constants
eo = 8.854e-12; mu_o = 4*pi*1e-7; co = 1/sqrt(eo*mu_o);

% STEP 1:  FILT

% dispersive medium constants
a = 0; b = 0.9; c = 0.8;
e_inf =2; e_s = 50; De = e_s-e_inf; 
tau = 153e-12;

% filt - constants
% Dt = 1.768e-12; 
Dt = 0.4e-12;
Nt = 50000;
tn = (0.5:1:Nt+0.5);  % tn = t/Dt
K = 21; alpha=5; p=8; 
f_FILT = FILT(tn,K,p,alpha,tau/Dt,e_s,e_inf,a,b,c); 

% STEP 2: PRONY 

M = 10; % filter order
[fprony,wHN,xHN] = myprony(f_FILT,tn,M);


% PLOTS
nPlot = 1000; 
figure; 
plot(tn(1:nPlot),f_FILT(1:nPlot)); hold on ; 
plot(tn(1:nPlot),fprony(1:nPlot),'k'); hold off;
legend('FILT',['FILT+Prony (M = ', num2str(M),')']);
title(['HN (á = ' ,num2str(a), ', â = ' ,num2str(b), ', ã = ' ,num2str(c), ')'] );
xlabel('Time Step Index'); ylabel('Impulse Response');

%% STEP 3: FDTD - 1D - Simulation of a pulse propagating in FREE SPACE

e_inf =1; e_s = 1; De = e_s-e_inf;
M = 4; wfp = zeros(M,1);  xfp = zeros(M,1);

% Initialization
ke = 7000;
Ex = zeros(1,ke); 
Hy = zeros(1,ke);
Dx = zeros(1,ke);
P = zeros(M,ke);
L = ones(1,ke);

% FDTD-constants Ät and Äz
dt = Dt;  dz = 0.41e-3;

% input_wave constants
spread = 1.26e10;
fe = 6e9;
pos_inc = 3000; % position of incident wave

% MEDIUM starts at cell with number k_start
k_start = 4000;
L(k_start:end) = e_inf + De * sum(wfp);

% MAIN FDTD Loop

% n -> TIME
% k -> SPACE
N = 6000;
z1 = k_start - 1;
Ex_i = zeros(1,N);
for n = 1:N
    
    % update Dx
    for k=2:ke
       Dx(k) = Dx(k) + (dt/dz) * ( Hy(k-1)-Hy(k) ); 
    end
    
    % create a gaussian pulse 
    pulse = exp( -spread^2 * (n*dt - 4/spread)^2 ) * sin( 2*pi*fe*(n*dt - 4/spread) );
    Dx(pos_inc) = pulse + Dx(pos_inc);
    
    for k = 2:ke 
       
       % update Ex
       Ex(k) = Dx(k)/eo ;
       if k >= k_start
           % As n changes , P(j,k) changes, so xP changes
           % because xP = x(1,...,M) .* P[(1,...,M) k]
           xP = zeros(M,1);

           for j=1:M
               xP(j) = xfp(j) * P(j,k);
           end
           Ex(k) = Ex(k) - De * sum(xP);
       end
       Ex(k) = Ex(k) / L(k);
       
       % update P
       if k >= k_start
           for j = 1:M
              P(j,k) = xP(j) + wfp(j)*Ex(k); 
           end          
       end

    end
    
    % magnetic field - Hy
    for k=1:ke-1
       Hy(k) = Hy(k) + ( dt / (dz*mu_o) ) * ( Ex(k)-Ex(k+1) ); 
    end 
    Ex_i(n) = Ex(z1);
end

fplot  = [zeros(1,k_start-1) max(Ex)*ones(1,ke-k_start+1)];
figure;  plot(Ex); 
hold on; plot(1:ke,fplot,'--');


%% STEP 4: FDTD - 1D - Simulation of a pulse hitting a 
%                      Davidson - Cole MEDIUM

e_inf =2; e_s = 50; De = e_s-e_inf;
M = length(xHN);

% Initialization
ke = 7000;
Ex = zeros(1,ke); 
Hy = zeros(1,ke);
Dx = zeros(1,ke);
P = zeros(M,ke);
L = ones(1,ke);

% FDTD-constants Ät and Äz
dt = Dt;  dz = 0.41e-3;

% input_wave constants
spread = 1.26e10;
fe = 6e9;
pos_inc = 3000; % position of incident wave

% MEDIUM starts at cell with number k_start
k_start = 4000;
L(k_start:end) = e_inf + De * sum(wHN);

% MAIN FDTD Loop

% n -> TIME
% k -> SPACE
N = 6000;
Ex_tot = zeros(1,N);
z1 = k_start - 1;
for n = 1:N
    
    % update Dx
    for k=2:ke
       Dx(k) = Dx(k) + (dt/dz) * ( Hy(k-1)-Hy(k) ); 
    end
    
    % create a gaussian pulse 
    pulse = exp( -spread^2 * (n*dt - 4/spread)^2 ) * sin( 2*pi*fe*(n*dt - 4/spread) );
    Dx(pos_inc) = pulse + Dx(pos_inc);
    
    for k = 2:ke 
       
       % update Ex
       Ex(k) = Dx(k)/eo ;
       if k >= k_start
           % As n changes , P(j,k) changes, so xP changes
           % because xP = x(1,...,M) .* P[(1,...,M) k]
           xP = zeros(M,1);
           
           for j=1:M
               xP(j) = xHN(j) * P(j,k);
           end
           Ex(k) = Ex(k) - De * sum(xP);
       end
       Ex(k) = Ex(k) / L(k);
       
       % update P
       if k >= k_start
           for j = 1:M
              P(j,k) = xP(j) + wHN(j)*Ex(k); 
           end          
       end

    end
    
    % magnetic field - Hy
    for k=1:ke-1
       Hy(k) = Hy(k) + ( dt / (dz*mu_o) ) * ( Ex(k)-Ex(k+1) ); 
    end 
    Ex_tot(n) = Ex(z1);
end
%
fplot  = [zeros(1,k_start-1) max(Ex)*ones(1,ke-k_start+1)];
figure;  plot(Ex); 
hold on; plot(1:ke,fplot,'--');

% calculate Reflection pulse wave
Ex_r = Ex_tot - Ex_i;

%% FFT
nfft = 2^17;
Ex_f_i = fftshift( abs(fft(Ex_i,nfft)) );
Ex_f_r = fftshift( abs(fft(Ex_r,nfft)) );
f = (1/dt/1e9) * (-nfft/2:nfft/2-1)/nfft; % f -> GHz

% figure; plot(f,Ex_f_i); hold on; 
% plot(f,Ex_f_r); 
% legend('Incident','Reflection')
% xlabel('GHz')

% keep frequencies 0 - 18 GHz
f0 = find(f==0);
flag = f0;
while f(flag) < 18
    flag = flag + 1;
end

posf = f0:flag;
Ex_f_i = Ex_f_i(posf);
Ex_f_r = Ex_f_r(posf);
f = f(posf);
f = 1e9*f;

% PLOT Incident and Reflection waves
% in FREQUENCY domain
figure; semilogx(f/1e9,Ex_f_i); hold on; 
semilogx(f/1e9,Ex_f_r); hold off;
legend('Incident','Reflection')
xlabel('GHz')

% Reflection Coefficient
R_fdtd = zeros(1,length(posf));
for k=1:length(R_fdtd)
   R_fdtd(k) = Ex_f_r(k) / Ex_f_i(k);
end



% Theoretical Reflection Coefficient
clc; 

% Davidson - Cole medium - > paper j47
a = 0; b = 0.9; c = 0.8; tau = 153e-12;
e_inf = 2; e_s = 50; De = e_s-e_inf; 

w = 2*pi*f; num_fr = length(w);

er = zeros(1,length(w));
for k=1:num_fr
   er(k)= e_inf + De * CP(1i*w(k),tau,a,b,c); 
end
root_er = sqrt(er);
Rth = abs( (1-root_er) ./ (1+root_er) );

% gamma = -i*w.*root_er/co;
% dz = 0.1e-3; d = 60*dz;
% Tdw = exp(d*gamma);

% PLOT Reflection Coefficients
figure
% fvec = linspace(f(1),f(end),80);
% fvec = 1:3:length(f);
% semilogx(f(fvec)/1e9, R_fdtd(fvec),'o'); hold on;
semilogx(f/1e9, R_fdtd); hold on;
semilogx(f/1e9, Rth);
legend('R_f_d_t_d','R_t_h')
title('Havriliak - Negami')
xlabel('frequency (GHz)')
% xlim([0.7 1.7])

DR = abs( R_fdtd-Rth );
max(DR)
