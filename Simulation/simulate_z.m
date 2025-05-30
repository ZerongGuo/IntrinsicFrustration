clear;
clc;

N = 1000;  % Number of oscillators
t_max = 20;
dt = 0.0001;

if isempty(gcp('nocreate'))
    parpool(24);  % Create a parallel pool with 24 workers if none exists
end

% Initialize parameter ranges
K_values = linspace(0, 200, 100);       % Coupling strength K: 100 values from 0 to 200
omega0_values = linspace(0, 150, 100);  % Distance between two frequency peaks: 100 values from 0 to 150
z_diff_matrix = zeros(length(K_values), length(omega0_values));  % (optional) Stores |z+ - z-|
z_1 = zeros(length(K_values), length(omega0_values));  % Stores z+ (order parameter for one group)
z_2 = zeros(length(K_values), length(omega0_values));  % Stores z- (order parameter for the other group)
D = 1;  % Standard deviation of the Lorentzian distribution
lorentz = generate_sorted_lorentz_array(N/2, 0, D)';  % Lorentzian distribution centered at 0 with std D

SL = length(omega0_values);
% Loop through each combination of K and omega0
parfor i = 1:length(K_values)
    for j = 1:SL
        K = K_values(i);           % Current coupling strength
        omega0 = omega0_values(j);  % Current peak separation

        % Generate natural frequency array omega
        omega = [lorentz + omega0; -lorentz - omega0];

        % Simulate the oscillator system
        [z1, z2] = simulate_oscillators(N, omega, t_max, dt, K);  % Compute boundary between chimera and standing wave states

        % Store results
%         z_diff_matrix(i, j) = abs(abs(z2) - abs(z1));  % Optional: store difference between order parameters
        z_1(i, j) = z1;
        z_2(i, j) = z2;
    end
end

z_p = z_1;
z_n = z_2;

save('z1ChimeraSim910D1.mat', 'z_p');  % Save data to file
save('z2ChimeraSim910D1.mat', 'z_n');  % Save data to file

delete(gcp('nocreate'))  % Shut down parallel pool if it exists
