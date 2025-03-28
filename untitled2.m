% Finite Difference code for 1-D advection-diffusion problem.

% Parameters
n = 200; % Number of grid points
dx = 0.500; % Spatial step size
dt = 0.2500; % Time step size
u = 0.10; % Advection velocity
alpha = 0.25; % Diffusion coefficient
mstep = 1600; % Number of time steps

% Initialize arrays
fo = zeros(1, n+1); % Old concentration values
f = zeros(1, n+1); % New concentration values

% Initial and boundary conditions
fo(1) = 1.0; % Left boundary
f(1) = 1.0; % Left boundary
fo(n+1) = 0.0; % Right boundary
f(n+1) = 0.0; % Right boundary

% Time evolution
for kk = 1:mstep
    for i = 2:n % Loop over interior points
        adv = dt * u * (fo(i) - fo(i-1)) / dx; % Advection term
        f(i) = fo(i) + dt * alpha * (fo(i+1) - 2*fo(i) + fo(i-1)) / (dx^2) - adv; % Update concentration
    end
    fo = f; % Update old concentration values
end

% Output results
x = 0:dx:n*dx; % Spatial coordinates
figure(1)
plot(x, f)
xlabel('X')
ylabel('Concentration')
title('1-D Advection-Diffusion Solution')
