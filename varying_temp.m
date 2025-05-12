% File for ploting the bird population with temperature

% get the time gap first

% these are in daily^-1

% Parameters for 1950
% rI is incubation rate, and it becomes death rate after the peak
params_1950 = struct('M0', 200, 'R0', 200, 'I0', 1000, ...
                            'KI', 27000, 'rR', 0.00027, 'rM', 0.00027, 'rI', 0.08, ...  
                            'gammaR', 0.02, 'gammaM', 0.02, ...
                            'thetaR', 0.6, 'thetaM', 0.4, ...
                            'alpha', 1, 'beta', 1, 't', 12, 'T', 10, ...
                            'Ea', 60, 'R', 0.008314, 'B', 0.000025, 'a', 0.65); 
% 60 days
time_span = [0 60];

% Updated population dynamics function with insects
function dYdt = population_dynamics(Y, params)
    M = Y(1); % Pied Flycatcher population
    R = Y(2); % Great Tit population
    I = Y(3); % Insect population


    % Differential equations
    dMdt = params.rM * M * (1 - (M + params.alpha * R)/ (params.KI + I * params.thetaM * (M / (M + R))));
    dRdt = params.rR * R * (1 - (R + params.beta * M) / (params.KI + I * params.thetaR * (R / (M + R))));
   
    % Arrhenius term for temperature dependence of growth rate
    temp_effect = exp(-params.Ea * params.B^(params.a) / (params.R * (params.T + 273.15))); 
    dIdt = params.rI * I * temp_effect;

    if I >= 5000
        dIdt = params.rI * I - (params.gammaM * M * I + params.gammaR * R * I);
    end
    % Return derivatives
    dYdt = [dMdt; dRdt; dIdt];
end

% Initial conditions now include insect population
initial_conditions_1950 = [params_1950.M0, params_1950.R0, params_1950.I0]; % Populations (M, R, I)

temperatures = [10, 17, 23, 29];

max_insect_populations = zeros(size(temperatures));
time_of_max_populations = zeros(size(temperatures));

% Loop over each temperature
for i = 1:length(temperatures)
    params_1950.T = temperatures(i); % Update temperature

    % Solve ODE for current temperature
    [t1950, Y1950] = ode23s(@(t, Y) population_dynamics(Y, params_1950), time_span, initial_conditions_1950);

    % Find the maximum insect population and corresponding time
    [maxI, idx] = max(Y1950(:, 3));
    max_insect_populations(i) = maxI;
    time_of_max_populations(i) = round(t1950(idx));
end

% Display results for each temperature
for i = 1:length(temperatures)
    fprintf('At %d°C, the maximum insect population is %.2f, occurring at day %.2f.\n', ...
            temperatures(i), max_insect_populations(i), time_of_max_populations(i));
end


% plot bird with time gap

% Parameters for 1950
params_1950 = struct('M0', 200, 'R0', 200, 'I0', 1000, ...
                            'KI', 27000, 'rR', 0.3, 'rM', 0.3, 'rI', 3, ...  
                            'gammaR', 200, 'gammaM', 200, ...
                            'thetaR', 0.6, 'thetaM', 0.4, ...
                            'alpha', 1, 'beta', 1, 'T', 10, ...
                            'Ea', 60, 'R', 0.008314, 'peak', 28, ...
                            'B', 0.000025, 'a', 0.65); 

time_span = [0 70];

function dYdt = population_dynamics_actual(Y, params)
    M = Y(1); % Pied Flycatcher population
    R = Y(2); % Great Tit population
    I = Y(3); % Insect population

    arrival_time = 41;
    t_T = arrival_time - params.peak;

    time_factor_M = exp(-0.001 * t_T); % Exponentially decaying time factor for Pied Flycatchers
    time_factor_R = exp(0.001 * t_T);  % Exponentially increasing time factor for Great Tits


    % Differential equations
    dMdt = params.rM * M * (1 - (M + params.alpha * R) / (params.KI + I * params.thetaM * (M / (M + R)))) * time_factor_M;
    dRdt = params.rR * R * (1 - (R + params.beta * M) / (params.KI + I * params.thetaR * (R / (M + R)))) * time_factor_R;

    % Arrhenius term for temperature dependence of growth rate
    temp_effect = exp(-params.Ea * params.B^(params.a) / (params.R * (params.T + 273.15))); 
    dIdt = params.rI * I * temp_effect;

    if I >= 5000
        dIdt = params.rI * I - (params.gammaM * M * I + params.gammaR * R * I);
    end

    % Return derivatives
    dYdt = [dMdt; dRdt; dIdt];
end


% Initial conditions
initial_conditions_1950 = [params_1950.M0, params_1950.R0, params_1950.I0]; % Initial populations (M, R, I)

% Set up colors for different time gaps
colors = lines(length(time_of_max_populations)); % Using distinct colors for each time gap

% Plot populations for both bird species with different time gaps
figure('Position', [100, 100, 600, 400]); % Adjusted figure size
set(gcf, 'PaperPosition', [0.25 0.25 6 4]); % Adjust the paper position for printing

% Subplot for Pied Flycatcher population
% subplot(1, 2, 1);
hold on;
for i = 1:length(time_of_max_populations)
    params_1950.peak = time_of_max_populations(i); % Update the arrival time

    % Solve ODE
    [t, Y] = ode23s(@(t, Y) population_dynamics_actual(Y, params_1950), time_span, initial_conditions_1950);

    % Extract Pied Flycatcher population
    M_population = Y(:, 1);

    % Plot with different colors for each time gap
    plot(t, M_population, 'Color', colors(i, :), 'DisplayName', sprintf('Temperature = %d °C', temperatures(i)), 'LineWidth', 1.5);
end
xlabel('Time (years)');
ylabel('Pied Flycatcher Population (million)');
title('Pied Flycatcher Population Dynamics with varying temperatures');
legend;
grid on;


figure('Position', [100, 100, 600, 400]); % Adjusted figure size
set(gcf, 'PaperPosition', [0.25 0.25 6 4]); % Adjust the paper position for printing
% Subplot for Great Tit population
% subplot(1, 2, 2);
hold on;
for i = 1:length(time_of_max_populations)
    params_1950.peak = time_of_max_populations(i); % Update the arrival time

    % Solve ODE
    [t, Y] = ode23s(@(t, Y) population_dynamics_actual(Y, params_1950), time_span, initial_conditions_1950);

    % Extract Great Tit population
    R_population = Y(:, 2);

    % Plot with different colors for each time gap
    plot(t, R_population, 'Color', colors(i, :), 'DisplayName', sprintf('Temperature = %d °C', temperatures(i)), 'LineWidth', 1.5);
end
xlabel('Time (years)');
ylabel('Great Tit Population (million)');
title('Great Tit Population Dynamics with varying temperatures');
legend;
grid on;

