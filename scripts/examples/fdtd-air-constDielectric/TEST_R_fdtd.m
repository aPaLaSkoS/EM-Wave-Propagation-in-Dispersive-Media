% 1. FREE SPACE simulation -> get Ex_i(z1)
% 2. FREE SPACE + DIELECTRIC simulation -> get Ex_t(z1)

%% 
clear; clc; close all;
% STEP 1: FDTD - 1D - Simulation of a pulse propagating in free space

e_inf =1; e_s = 1; De = e_s-e_inf;
M = 5; wfp = zeros(M,1);  xfp = zeros(M,1);

% Initialization
ke = 7000;
Ex = zeros(1,ke); 
Hy = zeros(1,ke);
Dx = zeros(1,ke);
P = zeros(M,ke);
L1 = ones(1,ke);


% general constants
eo = 8.854e-12; mu_o = 4*pi*1e-7; co = 1/sqrt(eo*mu_o);
f = (0.1:0.01:10)*1e9; wrad = f/(2*pi);

% find FDTD-constants Ät and Äz
% dt = 1.768e-12; dz = 1.1e-3;
dt = 3.7e-13; dz = 2.22e-4;

% input_wave constants
spread = 1.26e10;
fe = 6e9;
pos_inc = 3000; % position of incident wave

% MEDIUM starts at cell with number k_start
k_start = 4000;
L1(k_start:end) = e_inf + De * sum(wfp);

% MAIN FDTD Loop

% n -> TIME
% k -> SPACE
N = 6000;

Ex_i = zeros(1,N);
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
               xP(j) = xfp(j) * P(j,k);
           end
           Ex(k) = Ex(k) - De * sum(xP);
       end
       Ex(k) = Ex(k) / L1(k);
       
       
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


%% STEP 2: FDTD - 1D - Simulation of a pulse hitting a CONSTANT dielectric

e_inf =9; e_s = 9; De = e_s-e_inf; 
M=5; wd = zeros(M,1);  xd = zeros(M,1);

% Initialization
ke=7000;
Ex = zeros(1,ke); 
Hy = zeros(1,ke);
Dx = zeros(1,ke);
P = zeros(M,ke);
L1 = ones(1,ke);


% general constants
eo = 8.854e-12; mu_o = 4*pi*1e-7; co = 1/sqrt(eo*mu_o);
% dt = 1.768e-12; dz = 1.1e-3;
dt = 3.7e-13; dz = 2.22e-4;

% input_wave constants
spread = 1.26e10;
fe = 6e9;
pos_inc = 3000; % position of incident wave

% MEDIUM starts at cell with number k_start
k_start = 4000;
L1(k_start:end) = e_inf + De * sum(wd);

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
               xP(j) = xd(j) * P(j,k);
           end
           Ex(k) = Ex(k) - De * sum(xP);
       end
       Ex(k) = Ex(k) / L1(k);
       
       
       % update P
       if k >= k_start
           for j = 1:M
              P(j,k) = xP(j) + wd(j)*Ex(k); 
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


%%
Ex_r = Ex_tot - Ex_i;

%% FFT
nfft = 2^15;
Ex_f_i = fftshift( abs(fft(Ex_i,nfft)) );
Ex_f_r = fftshift( abs(fft(Ex_r,nfft)) );
f = (1/dt/1e9) * (-nfft/2:nfft/2-1)/nfft;

% figure; plot(f,Ex_f_i); hold on; 
% plot(f,Ex_f_r); 
% legend('Incident','Reflection')
% xlabel('GHz')

f0 = find(f==0);
flag = f0;
while f(flag) < 16
    flag = flag + 1;
end

posf = f0:flag;
Ex_f_i = Ex_f_i(posf);
Ex_f_r = Ex_f_r(posf);
f = f(posf);
f = 1e9*f;
figure; semilogx(f/1e9,Ex_f_i); hold on; 
semilogx(f/1e9,Ex_f_r); hold off;
legend('Incident','Reflection')
xlabel('GHz')

% Reflection Coefficient
Rw = zeros(1,length(posf));
for k=1:length(Rw)
   Rw(k) = Ex_f_r(k) / Ex_f_i(k);
end



%% Theoretical Reflection Coefficient
clc; 

% dielectric
er = 9 * ones(1,length(f));
root_er = sqrt(er);
Rth = (1-root_er) ./ (1+root_er) ;

% plots
f = f/1e9;
figure
semilogx(f, Rw); hold on;
semilogx(f,abs(Rth)); hold on; 
legend('R_f_d_t_d','R_t_h')
title('Havriliak - Negami')
xlabel('frequency (GHz)')
ylim([0.47 0.57])

