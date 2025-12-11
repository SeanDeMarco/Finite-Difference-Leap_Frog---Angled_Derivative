% Initialising given parameters

x_min = -1;         % Min. x value
x_max = 1;          % Max. x value
sigma = 0.6;        % Courant number
a = 0.75;           % Propagation velocity
T = 0.5;            % Final time
dx = 0.005;         % Spatial step size
dt = sigma * dx / a;  % Time step based on Courant number

N = round((2 / dx) + 1);  % Number of grid points
x = linspace(x_min, x_max, N);        % Creating the spatial domain

% N.B: _lfs and _ads refer to the leap_frog_scheme and
% angled_derivative_scheme respectively in the code

% Initial condition for square wave at t = 0 - this will be the n-1 th
% value for u
u = zeros(1, N);                     
u(x >= -0.2 & x <= 0.2) = 1;          

% Initialise u for both schemes
u_lfs = u;
u_ads = u;

u_dt = zeros(1, N);
x_shifted_dt = x-a*dt;
u_dt(x_shifted_dt >= -0.2 & x_shifted_dt <= 0.2) = 1;

% Initialise u for t = dt using the analytical solution - this will
% progress as the n th value for u
u_n_lfs = u_dt;                          
u_n_ads = u_dt;                          

num_steps = round(T / dt);

% Time stepping section
for n = 2:num_steps
    % Initialising an array of u for the next time step - the n + 1 th value of u 
    u_next_lfs = zeros(1, N);
    u_next_ads = zeros(1, N);

    % Solving the schemes in the internal region
    for i = 2:N-1
        % Implementing the leap frog scheme
        u_next_lfs(i) = u_lfs(i) - sigma * (u_n_lfs(i+1) - u_n_lfs(i-1));
        
        % Implementing the angled derivative scheme
        u_next_ads(i) = u_ads(i-1) + (1 - 2 * sigma) * (u_n_ads(i) - u_n_ads(i-1));
    end

    % Solving the periodic boundary conditions
    u_next_lfs(1) = u_lfs(1) - sigma * (u_n_lfs(2) - u_n_lfs(N));
    u_next_lfs(N) = u_lfs(N) - sigma * (u_n_lfs(1) - u_n_lfs(N-1));
    
    u_next_ads(1) = u_ads(N) + (1 - 2 * sigma) * (u_n_ads(1) - u_n_ads(N));
    u_next_ads(N) = u_ads(N-1) + (1 - 2 * sigma) * (u_n_ads(N) - u_n_ads(N-1));

    % Plotting the solutions at each time step
    plot(x, u_next_lfs, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Leapfrog Scheme');
    hold on;
    plot(x, u_next_ads, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Angled Derivative Scheme');
    
    xlabel('t');
    ylabel('u(x, t)');
    
    title(['Time step: ', num2str(n), ', Time: ', num2str(n * dt, '%.3f')]);
    legend('Leapfrog Scheme', 'Angled Derivative Scheme');
    axis([x_min x_max -0.5 1.5]);
    hold off;
    drawnow;

    % Updating the u and the u_n arrays for the next time step
     u_lfs = u_n_lfs;             %n-1 --> n
     u_n_lfs = u_next_lfs;            %n ---> n+1
     u_ads = u_n_ads;             %n-1 --> n
     u_n_ads = u_next_ads;            %n ---> n+1    
end

% Calculating the exact solution at time T
u_exact = zeros(1, N);
x_shifted = x-a*T;
u_exact(x_shifted >= -0.2 & x_shifted <= 0.2) = 1;

% Plotting the exact solution against the numerical solutions at T
figure;
plot(x, u_next_lfs, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Leapfrog Scheme');
hold on;
plot(x, u_next_ads, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Angled Derivative Scheme');
plot(x, u_exact, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Exact Solution');

xlabel('x');
ylabel('u(x, T)');
title(['Solutions at Final Time T = ', num2str(T)]);
legend('Leapfrog Scheme', 'Angled Derivative Scheme', 'Exact Solution');
axis([x_min x_max -0.5 1.5]);
hold off;
