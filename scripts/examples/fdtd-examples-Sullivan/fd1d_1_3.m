clear; clc; close all;

% FDTD - 1D - Simulation of a pulse hitting a DIELECTRIC MEDIUM

ke=3000;
ex = zeros(1,ke); hy = zeros(1,ke);

% pulse parameters
kc = 2*ke/3;
t0 = 40;
spread = 12;

% create dielectric profile
cb = ones(1,ke);
cb = 0.5*cb;
cb_start = 2300;
epsilon = 4;
cb(cb_start:end) = 0.5/epsilon;

% MAIN FDTD Loop
% n -> TIME
% k -> SPACE
N = 800;
for n = 1:N+1
   
    for k=2:ke
       ex(k) = ex(k) + cb(k)* ( hy(k-1)-hy(k) ); 
    end
    
    % put a gaussian pulse in the middle
    pulse = exp( -0.5 * ((t0-n)/spread)^2 );
    ex(1800) = pulse + ex(1800);
    
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
end
for k=1:1800
    ex800(k) = 0;
end
fplot  = [zeros(1,cb_start-1) ones(1,ke-cb_start+1)];

%%
n = 1:N;
pulse = exp( -0.5 * ((t0-n)/spread).^2 );
pos_max = find(pulse == max(pulse))
figure; plot(n,pulse);

figure
plot(ex800); ylabel('E_x'); xlabel('FDTD cells'); 
ylim([-1.2 1.2]); title(['FDTD - 1D - dielectric medium (T = ', num2str(N800) , ')' ]);
hold on; plot(1:3000,fplot,'--');
k_max = find(ex800 == max(ex800))
d = k_max - 1800 + 20
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
plot(ex); ylabel('E_x'); xlabel('FDTD cells'); 
ylim([-1.2 1.2]); title(['FDTD - 1D - dielectric medium (T = ', num2str(N) , ')' ]);
hold on; plot(1:3000,fplot,'--');





