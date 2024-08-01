clear; clc; close all;

% FDTD - 1D - Simulation of a pulse hitting a lossy medium 

ke=3000;
ex = zeros(1,ke); hy = zeros(1,ke);

ddx = 0.01;
dt = ddx / 6e8;
% freq_in = 200e6;

% pulse parameters
kc = 2*ke/3;
t0 = 40;
spread = 12;

% create dielectric profile
epsz = 8.854e-12;
epsilon = 4;
sigma = 0.04;

ca = ones(1,ke);
cb = 0.5 * ones(1,ke);
cb_start = 2300;

eaf = dt * sigma / (2*epsz*epsilon);
ca(cb_start:end) = (1-eaf) / (1+eaf);
cb(cb_start:end) = 0.5 / (epsilon * (1+eaf) );

% MAIN FDTD Loop
% n -> TIME
% k -> SPACE
N = 1500;
for n = 1:N+1
   
    % electric field - Ex
    for k=2:ke
       ex(k) = ca(k) * ex(k) + cb(k)* ( hy(k-1)-hy(k) ); 
    end
    
    % put a gaussian pulse in the middle
    pulse = exp( -0.5 * ((t0-n)/spread)^2 );
    ex(1800) = pulse + ex(1800);
%     pulse = sin( 2*pi*freq_in*dt*n);
%     ex(1800) = pulse + ex(1800);

    
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
    
    if n == 1200
        ex1200 = ex; N1200 = n;
    end
end
fplot  = [zeros(1,cb_start-1) ones(1,ke-cb_start+1)];
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
plot(ex1200); ylabel('E_x'); xlabel('FDTD cells'); 
ylim([-1.2 1.2]); title(['FDTD - 1D - dielectric medium (T = ', num2str(N1200) , ')' ]);
hold on; plot(1:3000,fplot,'--');

% figure
% plot(hy); ylabel('H_y'); xlabel('FDTD cells'); 
% ylim([-1.2 1.2]); title(['FDTD - 1D - free space (T = ', num2str(N) , ')' ]);



