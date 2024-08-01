clear; clc; close all;

% FDTD - 1D - FREE SPACE

ke=200;
ex = zeros(1,ke); hy = zeros(1,ke);

% pulse parameters
kc = ke/2;
t0 = 40;
spread = 12;
N = 100;

% MAIN FDTD Loop
% n -> TIME
% k -> SPACE
for n = 1:N+1
   
    for k=2:ke
       ex(k) = ex(k) + 0.5* ( hy(k-1)-hy(k) ); 
    end
    
    % put a gaussian pulse in the middle
    pulse = exp( -0.5 * ((t0-n)/spread)^2 );
    ex(kc) = pulse;
    
    for k=1:ke-1
       hy(k) = hy(k) + 0.5* ( ex(k)-ex(k+1) ); 
    end
    
end

%%
figure
plot(ex); ylabel('E_x'); xlabel('FDTD cells'); 
ylim([-1.2 1.2]); title(['FDTD - 1D - free space (T = ', num2str(N) , ')' ]);
figure
plot(hy); ylabel('H_y'); xlabel('FDTD cells'); 
ylim([-1.2 1.2]); title(['FDTD - 1D - free space (T = ', num2str(N) , ')' ]);



