clear; clc; close all;

% FDTD - 1D - Simulation of a pulse hitting a DEBYE MEDIUM

ke=3000;
ex = zeros(1,ke); 
hy = zeros(1,ke);
dx = zeros(1,ke);
ix = zeros(1,ke);
sx = zeros(1,ke);

ddx = 0.01;
dt = ddx / 6e8;

t0 = 40;
spread = 10;

% create dielectric profile
epsz = 8.854e-12;
epsr =2;
sigma = 0.01;
tau = 0.001 * 1e-6;
chi = 2;
k_start = 2300;

gax = ones(1,ke);
gbx = zeros(1,ke);
gcx = zeros(1,ke);
gax(k_start:end) = 1 / (epsr + (sigma*dt/epsz) + chi*dt/tau);
gbx(k_start:end) = sigma*dt/epsz ;
gcx(k_start:end) = chi*dt/tau ;
del_exp = exp(-dt/tau);

% MAIN FDTD Loop
% n -> TIME
% k -> SPACE
N = 1500;
for n = 1:N+1
   
    % electric field - Ex
    for k=2:ke
       dx(k) = dx(k) + 0.5 * ( hy(k-1)-hy(k) ); 
    end
    
    % create a gaussian pulse 
    pulse = exp( -0.5 * ((t0-n)/spread)^2 );
    dx(1800) = pulse + dx(1800);
    
    for k = 2:ke
       ex(k) = gax(k) * ( dx(k)-ix(k)-del_exp*sx(k) );
       ix(k) = ix(k) + gbx(k) * ex(k);
       sx(k) = del_exp * sx(k) + gcx(k) * ex(k);
    end
    
    % magnetic field - Hy
    for k=1:ke-1
       hy(k) = hy(k) + 0.5* ( ex(k)-ex(k+1) ); 
    end
    
    if n == 600
        ex600 = ex; N600 = n;
    end
    
    if n==800
        ex800 = ex; N800 = n;
    end
    
    if n == 1000
        ex1000 = ex; N1000 = n;
    end
    
    if n == 1100
        ex1100 = ex; N1100 = n;
    end
end
fplot  = [zeros(1,k_start-1) ones(1,ke-k_start+1)];
%%
figure

subplot(2,2,1)
plot(ex600); ylabel('E_x'); xlabel('FDTD cells'); 
ylim([-1.2 1.2]); title(['FDTD - 1D - dielectric medium (T = ', num2str(N600) , ')' ]);
hold on; plot(1:3000,fplot,'--'); 

subplot(2,2,2)
plot(ex800); ylabel('E_x'); xlabel('FDTD cells'); 
ylim([-1.2 1.2]); title(['FDTD - 1D - dielectric medium (T = ', num2str(N800) , ')' ]);
hold on; plot(1:3000,fplot,'--');

subplot(2,2,3)
plot(ex1000); ylabel('E_x'); xlabel('FDTD cells'); 
ylim([-1.2 1.2]); title(['FDTD - 1D - dielectric medium (T = ', num2str(N1000) , ')' ]);
hold on; plot(1:3000,fplot,'--');

subplot(2,2,4)
plot(ex1100); ylabel('E_x'); xlabel('FDTD cells'); 
ylim([-1.2 1.2]); title(['FDTD - 1D - dielectric medium (T = ', num2str(N1100) , ')' ]);
hold on; plot(1:3000,fplot,'--');

figure
plot(ex1100); ylabel('E_x'); xlabel('FDTD cells'); 
ylim([-1.2 1.2]); title(['FDTD - 1D - dielectric medium (T = ', num2str(N1100) , ')' ]);
hold on; plot(1:3000,fplot,'--');

figure
plot(ex); ylabel('E_x'); xlabel('FDTD cells'); 
ylim([-1.2 1.2]); title(['FDTD - 1D - dielectric medium (T = ', num2str(N) , ')' ]);
hold on; plot(1:3000,fplot,'--');