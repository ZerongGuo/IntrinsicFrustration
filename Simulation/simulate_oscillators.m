function [z1, z2] = simulate_oscillators(N, omega, t_max, dt, K)
    % Parameter settings
    n = 3;  % Exponent of the nonlinearity
    t_steps = t_max / dt;  % Number of time steps
    unite1 = [0, 0, 0, 0, 0, 0, 0, 2*pi/n, 0, 0];  % Fixed phase pattern with strict 9:1:0 ratio
    theta1 = repmat(unite1, 1, N/20)';  % Initial phases for the first group
%     theta1 = (rand(N/2, 1) < 0.1) * 2/3*pi;  % (Alternative) Random initial condition
    theta2 = flip(theta1);  % Initial phases for the second group (mirror of the first)
                             % The two groups have identical initial phases but asymmetric dynamics
    theta = [theta1; theta2];  % Full phase vector

    % Define the 4th-order Runge-Kutta method
    function dtheta = theta_dot(theta, omega, K, N, n)
        z = (1/N) * sum(exp(1i * theta));  % Compute the global order parameter z
        R = abs(z);  % Magnitude of z
        phase = angle(z);  % Global phase of z
        dtheta = omega - (K * R^n) * sin(n * (theta - phase));  % Phase velocity
    end

    % Time evolution
    for t = 1:t_steps
        % 4th-order Runge-Kutta update
        k1 = dt * theta_dot(theta, omega, K, N, n);
        k2 = dt * theta_dot(theta + 0.5 * k1, omega, K, N, n);
        k3 = dt * theta_dot(theta + 0.5 * k2, omega, K, N, n);
        k4 = dt * theta_dot(theta + k3, omega, K, N, n);
        theta = theta + (k1 + 2*k2 + 2*k3 + k4) / 6;  % Update phases
    end

    % Compute z1 and z2
    z1 = (2/N) * sum(exp(1i * theta(1:N/2)));       % Partial order parameter for forward-rotating group
    z2 = (2/N) * sum(exp(1i * theta(N/2+1:N)));     % Partial order parameter for backward-rotating group
end
