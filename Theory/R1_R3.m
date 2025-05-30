% Define parameter ranges
K_values = linspace(0.01, 150, 1000);  % Range of K values
Delta = 1;  % Does not affect the functional relationship, only the value range. Should not be too large or too small.
gamma0 = 0.9;
gamma1 = 0.1;
gamma2 = 0;  % Proportions of oscillators in the three synchronized groups

zf = gamma0 * sum(exp(1i * 0)) + gamma1 * sum(exp(1i * 2*pi/3)) + gamma2 * sum(exp(1i * 4*pi/3));
pf = 3 * angle(zf);  % Intrinsic frustration phase

% Initialize result matrices
solutions = zeros(length(K_values), 1);     % Store R1 values
solutions3 = zeros(length(K_values), 1);    % Store R3 values

% Compute for each K
for i = 1:length(K_values)
    KR13 = K_values(i);
    omega0 = 0;
    
    % Define the Lorentzian distribution
    g = @(omega) (Delta/pi)./((omega - omega0).^2 + Delta^2);   

    % Iterative computation
    % Approximate Omega_1 as KR13*sin(pf), i.e., the angular velocity at θ = Θ_3/3 for oscillators with ω_i = 0
    solutions(i) = 3 * KR13 * integral(@(theta) cos(theta).*cos(3*theta - pf) .* ...
        g(KR13 * sin(3*theta - pf) + KR13 * sin(pf)), -pi/6 + pf/3, pi/6 + pf/3); 
    solutions3(i) = 3 * KR13 * integral(@(theta) cos(3*theta).*cos(3*theta - pf) .* ...
        g(KR13 * sin(3*theta - pf) + KR13 * sin(pf)), -pi/6 + pf/3, pi/6 + pf/3);
end

% Create a new figure window
figure;

% Plot the relationship between R_{1,p} and R_{3,p} in one synchronized region under quasi-static conditions
plot(solutions, solutions3, 'b-', 'LineWidth', 2, 'DisplayName', 'R1 vs R3');

% Add axis labels
% title('Plot of R1 vs R3');
xlabel('R_{1,p}', 'FontSize', 15, 'FontName', 'Arial');
ylabel('R_{3,p}', 'FontSize', 15, 'FontName', 'Arial');
ax = gca;  % Get current axes
ax.FontSize = 15;  % Set axis font size

% save('R1_910.mat', 'solutions');  % Save R1 data to file
% save('R3_910.mat', 'solutions3');  % Save R3 data to file
