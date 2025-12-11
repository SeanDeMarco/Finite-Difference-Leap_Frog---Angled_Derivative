% Initialising given parameters

x_min = -1;         % Min. x value
x_max = 1;          % Max. x value
sigma = 0.6;        % Courant number
a = 0.75;           % Proagation velocity
T = 1;              % Final time

% N.B: _lfs and _ads refer to the leap_frog_scheme and
% angled_derivative_scheme respectively in the code

% Initialising an arrary for the dx values to be able to loop over
x_list = [0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001];          

% Initialising a zeros matrix to store the max. error of each dx
errors_lfs = zeros(size(dx_list));          
errors_ads = zeros(size(dx_list)); 

% Selecting dx value
for i_dx = 1:length(dx_list)
    dx = dx_list(i_dx);
    dt = (sigma * dx) / a;          % Calculating time step based on courant number
    N = round((2/dx) + 1);         % Number of grid points using given equation
    x = linspace(x_min, x_max, N);          % Creating the spatial domain 
    
    % Solving the initial condition and t = 0 - this will be the n-1 th value for u   
    u = sin(2 * pi * x);       
    %Solving the analytical solution at t = dt - this will progress as the n th value for u
    u_n = sin(2 * pi * (x - a * dt));            
    
    % Assigning the initial conditions for the leap frog scheme
    u_lfs = u;
    u_n_lfs = u_n;
    
    % Assigning the initial conditions for the angled derivative scheme
    u_ads = u;
    u_n_ads = u_n;

    % Time stepping section
    num_steps = round(T/dt);            
    for n = 2:num_steps         % Starts from n = 2 as n = 0 and 1 done previously
        % Initialising an array of u for the next time step - the n + 1 th value of u
        u_next_lfs = zeros(1,N);            
        u_next_ads = zeros(1,N);
        
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

        %  % lotting at each time step so the numerical solution vs
        %  % analytical solution can be seen propagating in time visually
        %  plot(x, u_next_lfs, 'b', 'LineWidth', 2, 'DisplayName', 'Leapfrog Scheme');
        %    hold on;
        % plot(x, u_next_ads, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Angled Derivative Scheme');
        % plot(x, sin(2 * pi * (x - a * n * dt)), 'r--', 'DisplayName', 'Exact Solution'); 
        % hold off;        
        % xlabel('t');
        % ylabel('u(x, t)');
        % title(['Time step: ', num2str(n), ', Time: ', num2str(n * dt, '%.3f'), ', dx: ', num2str(dx)]);
        % legend('show'); 
        % axis([x_min x_max -1 1]);
        % drawnow;

        % Updating the u and u_n arrays for the next time step
        u_lfs = u_n_lfs;             %n-1 --> n
        u_n_lfs = u_next_lfs;            %n ---> n+1
        u_ads = u_n_ads;             %n-1 --> n
        u_n_ads = u_next_ads;            %n ---> n+1        
    end

    % Solving analytically at t = T to get the discrepancy between numerical and
    % exact solution    
    exact_solution = sin(2 * pi * (x - a * T));         % Where t = T = 1

    % Computing maximum error at final time T
    errors_lfs(i_dx) = max(abs(u_n_lfs - exact_solution));
    errors_ads(i_dx) = max(abs(u_n_ads - exact_solution));
end


%% 
% Plotting the error on a log-log scale for both schemes
figure;
loglog(dx_list, errors_lfs, '-o', 'DisplayName', 'Leapfrog Scheme');  
hold on;
loglog(dx_list, errors_ads, '-s', 'DisplayName', 'Angled Derivative Scheme'); 
hold off;
xlabel('\Delta x');
ylabel('Maximum Error');
title('Error Analysis for Leapfrog and Angled Derivative Schemes');
legend('show'); 
grid on;

