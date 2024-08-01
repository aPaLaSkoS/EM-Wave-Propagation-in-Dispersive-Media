clear; clc; close all;

% FDTD - 1D - Simulation of a pulse hitting a DEBYE (b=1) or DC (b~=0) MEDIUM

% Initialization
ke=3000;
Ex = zeros(1,ke); 
Hy = zeros(1,ke);
Dx = zeros(1,ke);
Ix = zeros(1,ke);
Sx = zeros(1,ke);
gax = ones(1,ke);
gbx = zeros(1,ke);
gcx = zeros(1,ke);

% general constants
co = 3e8; eo = 8.854e-12;
f = (0.1:0.01:10)*1e9; wrad = f/(2*pi);

% create dielectric profile
e_inf =2; es = 50; De = es-e_inf; sigma = 0;
b = 0.65;   tau = 153e-12;

% find FDTD-constants Ät and Äz
rp = real( e_inf + De ./ (1 + 1i*tau*wrad).^b );  % Relative Permittivity å_r(ù)
vel_medium = co ./ rp;   max_vel = max(vel_medium);
lamda = vel_medium ./ f;  lmin = min(lamda);
dz = lmin/10; % Äz = Äz_min
dt_max = dz/max_vel;
dt = 1.5e-12; % < dt_max

% input_wave constants
spread = 1.26e10;
fe = 6e10;
pos_inc = 1800; % position of incident wave

% MEDIUM starts at cell with number k_start
k_start = 2300;
gax(k_start:end) = 1 / (e_inf + (sigma*dt/eo) + De*dt/tau);
gbx(k_start:end) = sigma*dt/eo ;
gcx(k_start:end) = De*dt/tau ;
del_exp = exp(-dt/tau);

%% MAIN FDTD Loop

% n -> TIME
% k -> SPACE
N = 2500;
for n = 1:N+1
   
    % electric field - Ex
    for k=2:ke
       Dx(k) = Dx(k) + 0.5 * ( Hy(k-1)-Hy(k) ); 
    end
    
    % create a modulated gaussian pulse 
    pulse = exp( -spread^2 * (n*dt - 4/spread).^2 ) .* sin( 2*pi*fe*(n*dt - 4/spread) );
    Dx(pos_inc) = pulse + Dx(pos_inc);
    
    for k = 2:ke
       Ex(k) = gax(k) * ( Dx(k)-Ix(k)-del_exp*Sx(k) );
       Ix(k) = Ix(k) + gbx(k) * Ex(k);
       Sx(k) = del_exp * Sx(k) + gcx(k) * Ex(k);
    end
    
    % magnetic field - Hy
    for k=1:ke-1
       Hy(k) = Hy(k) + 0.5* ( Ex(k)-Ex(k+1) ); 
    end
    
    if n == 1100
        ex1 = Ex; N1 = n;
    end
    
    if n == 1408
        ex2 = Ex; N2 = n;
    end
    
    if n==2000
        ex3 = Ex; N3 = n;
    end
    
    if n == 2500
        ex4 = Ex; N4 = n;
    end
    

end
fplot  = [zeros(1,k_start-1) ones(1,ke-k_start+1)];

%% PLOTS
figure

subplot(2,2,1)
plot(ex1); ylabel('E_x'); xlabel('FDTD cells'); 
ylim([-1.2 1.2]); title(['FDTD - 1D - dielectric medium (T = ', num2str(N1) , ')' ]);
hold on; plot(1:3000,fplot,'--');

subplot(2,2,2)
plot(ex2); ylabel('E_x'); xlabel('FDTD cells'); 
ylim([-1.2 1.2]); title(['FDTD - 1D - dielectric medium (T = ', num2str(N2) , ')' ]);
hold on; plot(1:3000,fplot,'--'); 

subplot(2,2,3)
plot(ex3); ylabel('E_x'); xlabel('FDTD cells'); 
ylim([-1.2 1.2]); title(['FDTD - 1D - dielectric medium (T = ', num2str(N3) , ')' ]);
hold on; plot(1:3000,fplot,'--');

subplot(2,2,4)
plot(ex4); ylabel('E_x'); xlabel('FDTD cells'); 
ylim([-1.2 1.2]); title(['FDTD - 1D - dielectric medium (T = ', num2str(N4) , ')' ]);
hold on; plot(1:3000,fplot,'--');



% figure
% plot(ex1100); ylabel('E_x'); xlabel('FDTD cells'); 
% ylim([-1.2 1.2]); title(['FDTD - 1D - dielectric medium (T = ', num2str(N1100) , ')' ]);
% hold on; plot(1:3000,fplot,'--');
% 
% figure
% plot(ex); ylabel('E_x'); xlabel('FDTD cells'); 
% ylim([-1.2 1.2]); title(['FDTD - 1D - dielectric medium (T = ', num2str(N) , ')' ]);
% hold on; plot(1:3000,fplot,'--');