% extensions: insect peak & competition coef

% get the time gap first

% these are in daily^-1

% Parameters 
% rI is incubation rate, and it becomes death rate after the peak
params_1950 = struct('M0', 200, 'R0', 200, 'I0', 1000, ...
                            'KI', 27000, 'rR', 0.00027, 'rM', 0.00027, 'rI', 0.08, ...  
                            'gammaR', 0.02, 'gammaM', 0.02, ...
                            'thetaR', 0.6, 'thetaM', 0.4, ...
                            'alpha', 1, 'beta', 1, 't', 12, 'T', 13, ...
                            'Ea', 60, 'R', 0.008314, 'IP', 5000, ...
                            'B', 0.000025, 'a', 0.65); 


% Parameters for 1950 for years
params_1950_actual = struct('M0', 200, 'R0', 200, 'I0', 1000, ...
                            'KI', 27000, 'rR', 0.3, 'rM', 0.3, 'rI', 3, ...  
                            'gammaR', 200, 'gammaM', 200, ...
                            'thetaR', 0.6, 'thetaM', 0.4, ...
                            'alpha', 1, 'beta', 1, 'T', 13, ...
                            'Ea', 60, 'R', 0.008314, 'peak', 38, ...
                            'B', 0.000025, 'a', 0.65); 


% days
time_span = [0 70];

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

    if I >= params.IP
        dIdt = params.rI * I - (params.gammaM * M * I + params.gammaR * R * I);
    end
    % Return derivatives
    dYdt = [dMdt; dRdt; dIdt];
end


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

IP_values = [0.25 * params_1950.IP, params_1950.IP, 1.75 * params_1950.IP];
alpha_values = [0.25 * params_1950_actual.alpha, params_1950_actual.alpha, 1.75 * params_1950_actual.alpha];

figure;

% Loop to generate subplots for each IP and alpha value
for i = 1:3
    % Set the IP value
    params_1950.IP = IP_values(i);
    initial_conditions_1950 = [params_1950.M0, params_1950.R0, params_1950.I0]; % Initial populations (M, R, I)

    % Solve ODE for the population dynamics
    [t1950, Y1950] = ode23s(@(t, Y) population_dynamics(Y, params_1950), time_span, initial_conditions_1950);

    % Find the maximum insect population and corresponding time
    [maxI, idx] = max(Y1950(:, 3)); % Maximum insect population and the index where it occurs
    max_day = round(t1950(idx));  % Corresponding day for max insect population
    params_1950_actual.peak = max_day;

    initial_conditions_1950 = [params_1950_actual.M0, params_1950_actual.R0, params_1950_actual.I0]; 


    % Solve for actual population dynamics
    [t_actual, Y_actual] = ode23s(@(t, Y) population_dynamics_actual(Y, params_1950_actual), time_span, initial_conditions_1950);

    % Extract populations from the solution
    M_actual = Y_actual(:, 1); % Pied Flycatcher population
    R_actual = Y_actual(:, 2); % Great Tit population

    % Plot each subplot
    subplot(2, 3, i);
    plot(t_actual, M_actual, 'r--', 'DisplayName', 'Pied Flycatcher', 'LineWidth', 1.5); 
    hold on;
    plot(t_actual, R_actual, 'b-', 'DisplayName', 'Great Tit', 'LineWidth', 1.5);
    xlabel('Time (year)');
    ylabel('Population (million)');
    ylim([0, 30000]);
    title(sprintf('Insect Population Peak = %.2f billion', IP_values(i) / 1000));
    grid on;
end

for i = 1:3
    % Set the alpha value
    params_1950_actual.alpha = alpha_values(i);
    params_1950_actual.beta = 2 - alpha_values(i);

    params_1950.IP = 5000;
    initial_conditions_1950 = [params_1950.M0, params_1950.R0, params_1950.I0]; % Initial populations (M, R, I)

    % Solve ODE for the population dynamics
    [t1950, Y1950] = ode23s(@(t, Y) population_dynamics(Y, params_1950), time_span, initial_conditions_1950);

    % Find the maximum insect population and corresponding time
    [maxI, idx] = max(Y1950(:, 3)); % Maximum insect population and the index where it occurs
    max_day = round(t1950(idx));  % Corresponding day for max insect population
    params_1950_actual.peak = max_day;

    % Solve for actual population dynamics with modified alpha and beta
    [t_actual, Y_actual] = ode23s(@(t, Y) population_dynamics_actual(Y, params_1950_actual), time_span, initial_conditions_1950);

    % Extract populations from the solution
    M_actual = Y_actual(:, 1); % Pied Flycatcher population
    R_actual = Y_actual(:, 2); % Great Tit population

    % Plot each subplot
    subplot(2, 3, i + 3); % Plot on the bottom row
    plot(t_actual, M_actual, 'r--', 'DisplayName', 'Pied Flycatcher', 'LineWidth', 1.5); 
    hold on;
    plot(t_actual, R_actual, 'b-', 'DisplayName', 'Great Tit', 'LineWidth', 1.5);
    xlabel('Time (year)');
    ylabel('Population (million)');
    title(['Alpha = ' num2str(alpha_values(i)), ', Beta = ' num2str(2-alpha_values(i))])
    grid on;

    if i == 2
        ylim([0, 30000]); % Set y-axis limit for the 5th graph
    end

end

legend('Pied Flycatcher', 'Great Tit', 'Location', 'northeast', 'Orientation', 'vertical', 'FontSize', 10);
