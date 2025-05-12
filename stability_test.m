% file for stability test around eq point

% Parameters

params = struct('M0', 200, 'R0', 200, 'I0', 1000, ...
                'KI', 27000, 'rR', 0.3, 'rM', 0.3, 'rI', 3, ...
                'gammaR', 200, 'gammaM', 200, ...
                'thetaR', 0.6, 'thetaM', 0.4, ...
                'alpha', 1.5, 'beta', 0.5, 'T', 16, ...
                'Ea', 60, 'R', 0.008314, 'peak', 26, ...
                'B', 0.000025, 'a', 0.65);

arrival_time = 41;
t_T = arrival_time - params.peak;

time_factor_M = exp(-0.001 * t_T);
time_factor_R = exp(0.001 * t_T);

all_stable = true;
all_unstable = true;

stable_equilibria = [];

% Dynamical equations with current I
syms M R
dM_dt = params.rM * M * (1 - (M + params.alpha * R) / (params.KI + I * params.thetaM * (M / (M + R)))) * time_factor_M;
dR_dt = params.rR * R * (1 - (R + params.beta * M) / (params.KI + I * params.thetaR * (R / (M + R)))) * time_factor_R;

% Solve the system for equilibrium points (dM/dt = 0, dR/dt = 0)
equations = [dM_dt == 0, dR_dt == 0];
[sol_M, sol_R] = solve(equations, [M, R]);


% Loop through equilibrium points
for i = 1:length(sol_M)
    M_star = sol_M(i);
    R_star = sol_R(i);

    % Choose a Lyapunov function
    V = 0.5 * (M - M_star)^2 + 0.5 * (R - R_star)^2; % quadratic form

    % sub in all possible I values
    for I = 1000:100:5000
        % Compute time derivative of the Lyapunov function
        dV_dt = diff(V, M) * dM_dt + diff(V, R) * dR_dt;
        
        % Simplify for interpretation
        dV_dt_simplified = simplify(dV_dt);
        
        % Convert symbolic derivative to a function
        dV_dt_func = matlabFunction(dV_dt_simplified);
        
        % Check if the derivative is negative at the equilibrium
        dV_dt_value = dV_dt_func(double(M_star), double(R_star)); % Numerically evaluate the derivative at the equilibrium
        
        % Check stability based on sign of dV/dt
        if real(dV_dt_value) >= 0
            all_stable = false; % If any equilibrium is unstable, set to false
        elseif real(dV_dt_value) < 0
            all_unstable = false; % If any equilibrium is stable, set to false
        end
    end
end

% Display final result
if all_stable
    disp('All equilibria are Lyapunov stable.');
elseif all_unstable
    disp('All equilibria are unstable.');
else
    disp('Some equilibria are unstable.');
end
