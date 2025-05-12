% file for time series of bird species


% get the time gap first

% these are in daily^-1

% Parameters for 1950
params_1950 = struct('M0', 200, 'R0', 200, 'I0', 1000, ...
                            'KI', 27000, 'rR', 0.00027, 'rM', 0.00027, 'rI', 0.08, ...  
                            'gammaR', 0.02, 'gammaM', 0.02, ...
                            'thetaR', 0.6, 'thetaM', 0.4, ...
                            'alpha', 1, 'beta', 1, 'T', 10, ...
                            'Ea', 60, 'R', 0.008314, 'B', 0.000025, 'a', 0.65); 

params_2018 = struct('M0', 200, 'R0', 200, 'I0', 1000, ...
                            'KI', 27000, 'rR', 0.00027, 'rM', 0.00027, 'rI', 0.08, ...  
                            'gammaR', 0.02, 'gammaM', 0.02, ...
                            'thetaR', 0.6, 'thetaM', 0.4, ...
                            'alpha', 1.5, 'beta', 0.5, 'T', 16, ...
                            'Ea', 60, 'R', 0.008314, 'B', 0.000025, 'a', 0.65); 

% yearly^-1
% Parameters for 1950 
params_1950_real = struct('M0', 200, 'R0', 200, 'I0', 1000, ...
                            'KI', 27000, 'rR', 0.3, 'rM', 0.3, 'rI', 3, ...  
                            'gammaR', 200, 'gammaM', 200, ...
                            'thetaR', 0.6, 'thetaM', 0.4, ...
                            'alpha', 1, 'beta', 1, 'T', 10, ...
                            'Ea', 60, 'R', 0.008314, 'peak', 37, ...
                            'B', 0.000025, 'a', 0.65); 

% Parameters for 2018 
params_2018_real = struct('M0', 200, 'R0', 200, 'I0', 1000, ...
                            'KI', 27000, 'rR', 0.3, 'rM', 0.3, 'rI', 3, ...
                            'gammaR', 200, 'gammaM', 200, ...
                            'thetaR', 0.6, 'thetaM', 0.4, ...
                            'alpha', 1.5, 'beta', 0.5, 'T', 16, ...
                            'Ea', 60, 'R', 0.008314, 'peak', 26, ...
                            'B', 0.000025, 'a', 0.65);  


time_span = [0 70];

% population dynamics function with insects
function dYdt = population_dynamics_test(Y, params)
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

temperatures = [13, 16];

max_insect_populations = zeros(size(temperatures));
time_of_max_populations = zeros(size(temperatures));

% Loop over each temperature
for i = 1:length(temperatures)
    params_1950.T = temperatures(i); % Update temperature

    % Solve ODE for current temperature
    [t1950, Y1950] = ode23s(@(t, Y) population_dynamics_test(Y, params_1950), time_span, initial_conditions_1950);

    % Find the maximum insect population and corresponding time
    [maxI, idx] = max(Y1950(:, 3));
    max_insect_populations(i) = maxI;
    time_of_max_populations(i) = round(t1950(idx));
end


% population dynamics function with insects
function dYdt = population_dynamics(Y, params)
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

    temp_effect = exp(-params.Ea * params.B^(params.a) / (params.R * (params.T + 273.15)));  
    % dIdt = params.rI * I * temp_effect - (params.gammaM * M * I + params.gammaR * R * I); % Insect dynamics
    dIdt = params.rI * I * temp_effect;

    if I >= 5000
        dIdt = params.rI * I - (params.gammaM * M * I + params.gammaR * R * I);
    end
    % Return derivatives
    dYdt = [dMdt; dRdt; dIdt];
end

params_1950_real.peak = time_of_max_populations(1);
params_2018_real.peak = time_of_max_populations(2);


% Initial conditions now include insect population
initial_conditions_1950 = [params_1950_real.M0, params_1950_real.R0, params_1950_real.I0]; % Populations (M, R, I)
initial_conditions_2018 = [params_2018_real.M0, params_2018_real.R0, params_2018_real.I0]; % Populations (M, R, I)


% Solve ODE for 1950
[t1950, Y1950] = ode23s(@(t, Y) population_dynamics(Y, params_1950_real), time_span, initial_conditions_1950);

% Solve ODE for 2018
[t2018, Y2018] = ode23s(@(t, Y) population_dynamics(Y, params_2018_real), time_span, initial_conditions_2018);

% Plot results for Great Tits
figure('Position', [100, 100, 600, 400]); % Adjusted figure size
set(gcf, 'PaperPosition', [0.25 0.25 6 4]); % Adjust the paper position for printing
% subplot(2, 1, 1);
plot(t1950, Y1950(:,2), 'b-', 'LineWidth', 1.5); hold on;
plot(t2018, Y2018(:,2), 'r--', 'LineWidth', 1.5);
xlabel('Years');
ylabel('Population (millions)');
title('Population of Great Tits (1950 vs. 2018)');
legend('Starting from 1950', 'Starting from 2018');
grid on;

% Plot results for Pied Flycatchers
figure('Position', [100, 100, 600, 400]); % Adjusted figure size
set(gcf, 'PaperPosition', [0.25 0.25 6 4]); % Adjust the paper position for printing
% subplot(2, 1, 2);
plot(t1950, Y1950(:,1), 'b-', 'LineWidth', 1.5); hold on;
plot(t2018, Y2018(:,1), 'r--', 'LineWidth', 1.5);
xlabel('Years');
ylabel('Population (millions)');
title('Population of Pied Flycatchers (1950 vs. 2018)');
legend('Starting from 1950', 'Starting from 2018');
grid on;

% Plot results for 1950
figure('Position', [100, 100, 600, 400]); % Adjusted figure size
set(gcf, 'PaperPosition', [0.25 0.25 6 4]); % Adjust the paper position for printing
% subplot(2, 1, 2);
plot(t1950, Y1950(:,2), 'b-', 'LineWidth', 1.5); hold on;
plot(t1950, Y1950(:,1), 'r--', 'LineWidth', 1.5);
xlabel('Years');
ylabel('Population (million)');
title('Population of both species from 1950 to 70 years after');
legend('Great Tits', 'Pied Flycatcher');
grid on;

% Plot results for 2018
figure('Position', [100, 100, 600, 400]); % Adjusted figure size
set(gcf, 'PaperPosition', [0.25 0.25 6 4]); % Adjust the paper position for printing
% subplot(2, 1, 2);
plot(t2018, Y2018(:,2), 'b-', 'LineWidth', 1.5); hold on;
plot(t2018, Y2018(:,1), 'r--', 'LineWidth', 1.5);
xlabel('Years');
ylabel('Population (million)');
title('Population of both species from 2018 to 70 years after');
legend('Great Tits', 'Pied Flycatcher');
grid on;