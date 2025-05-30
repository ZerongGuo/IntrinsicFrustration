clc;
clear;

% Parameter settings
N = 1000;             % Number of oscillators in the network
dt = 0.001;           % Time step
t_total = 10000;      % Total simulation time
K = 130;              % Coupling strength
omega0 = 50;          % Frequency shift

% Model parameters
n = 3;                % Order for the coupling function

% Initialize theta and omega
D = 1;  % Standard deviation for the Lorentzian distribution

% Generate Lorentzian-distributed natural frequencies centered at 0
lorentz = generate_sorted_lorentz_array(N/2, 0, D);

% Construct the full omega array with symmetric frequency shift
omega = [lorentz + omega0, -lorentz - omega0];

% Initialize phase distribution:
% 10% of the oscillators are initialized at 2π/3, the rest at 0
unite1 = [0, 0, 0, 0, 0, 0, 0, 2*pi/n, 0, 0];  % Fixed phase pattern with strict 9:1:0 ratio
theta1 = repmat(unite1, 1, N/20);  % Initial phases for the first group
%     theta1 = (rand(N/2, 1) < 0.1) * 2/3*pi;  % (Alternative) Random initial condition
theta2 = flip(theta1);  % Initial phases for the second group (mirror of the first)
                         % The two groups have identical initial phases but asymmetric dynamics
theta = [theta1, theta2];  % Full phase vector

% Preallocate variables
dtheta = zeros(N, 1);
timesteps = 0:dt:t_total;

% Create figure for real-time animation
figure;
h1 = scatter(0, 0, 'bo');  % Scatter plot for the first half of the oscillators
hold on;
h2 = scatter(0, 0, 'bo');  % Scatter plot for the second half
h3 = scatter(0, 0, 'o');   % Marker for the global order parameter

title('Phase Evolution Animation');
xlabel('X');
ylabel('Y');
axis equal;         % Ensure equal scaling for X and Y
axis([-1 1 -1 1]);

grid on;

% Main simulation loop
for t_idx = 1:length(timesteps)
    % Compute order parameters
    z = (1/N) * sum(exp(1i * theta));
    z2 = (1/N) * sum(exp(1i * 2 * theta));
    R = abs(z);
    phase = angle(z);

    % Phase update rule for higher-order Kuramoto model
    dtheta = omega - K * R^n * sin(n * (theta - phase));
    
    % Update phases
    theta = theta + dt * dtheta;

    % Recompute order parameter after phase update
    z = (1/N) * sum(exp(1i * theta));
    R = abs(z);
    phase = angle(z);

    % Convert phases to Cartesian coordinates
    x = cos(theta);
    y = sin(theta);

    % Update scatter plots
    set(h1, 'XData', x(1:500), 'YData', y(1:500));
    set(h2, 'XData', x(501:end), 'YData', y(501:end));
    set(h1, 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'red');    
    set(h3, 'XData', R * cos(phase), 'YData', R * sin(phase));
    set(h3, 'MarkerFaceColor', 'green', 'MarkerEdgeColor', 'green');    

    % Refresh figure for animation
    drawnow;
end
