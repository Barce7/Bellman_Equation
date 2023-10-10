%%%%%%%%%%%%%%%%% Parameter Setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha       = 0.3;      % Capital share
beta        = 0.95;     % Discounted rate
A           = 2;
delta       = 0.1;     % Depreciation rate of capital stock
tol         = 1e-5;     % Convergence tolerance
error       = 1;        % Initial value of error
iter        = 1;        % Number of loops
% Set the iterations for which you want to save the value function
iterations_to_save = [0, 1, 2, 10, 50, 100, 226];
%%%%%%%%%%%%%%%%% Grid And Initial Guess %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the steady state of capital
kss     = (((1/beta)+delta-1)/(A*alpha))^(1/(alpha-1));
kmin    = 0.01;
kmax    = 15;

% Generate grids
nk      = 1500;  % Choose No. of grids to cover kss
kgrid   = linspace(kmin, kmax, nk)';

% Initial guess of value function
v       = zeros(nk,1);
tv      = zeros(nk,1);
dr      = zeros(nk,1);
%%%%%%%%%%%%%%%%% Main Loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
consumption = 2*kgrid.^alpha + (1- delta).*kgrid - kgrid'; 
consumption(consumption < 0) = NaN;

utility = log(consumption);
utility(isnan(utility)) = -inf;
% Initialize variables to store value functions at specific iterations
v0 = zeros(nk,1);
v1 = zeros(nk,1);
v2 = zeros(nk,1);
v10 = zeros(nk,1);
v50 = zeros(nk,1);
v100 = zeros(nk,1);
v226 = zeros(nk,1);
kp0 = zeros(nk,1);
kp1 = zeros(nk,1);
kp2 = zeros(nk,1);
kp10 = zeros(nk,1);
kp50 = zeros(nk,1);
kp100 = zeros(nk,1);
kp226 = zeros(nk,1);
while error > tol
    [tv, dr] = max(utility + beta*v', [], 2);
    error    = max(abs(tv - v));
    v        = tv;
    
    % Check if the current iteration is in the list to save
    if ismember(iter, iterations_to_save)
        % Store the value function in the corresponding variable
        switch iter
            case 0
                kp0 = kgrid(dr);
                v0 = v;
            case 1
                kp1 = kgrid(dr);
                v1 = v;
            case 2
                kp2 = kgrid(dr);
                v2 = v;
            case 10
                kp10 = kgrid(dr);
                v10 = v;
            case 50
                kp50 = kgrid(dr);
                v50 = v;
            case 100
                kp100 = kgrid(dr);
                v100 = v;
            case 226
                kp226 = kgrid(dr);
                v226 = v;
        end
    end

    iter     = iter + 1;
end

% Policy function
kp = kgrid(dr);
cp = 2*kgrid.^alpha + (1 - delta).*kgrid - kp;
% Create a figure with a white background
figure('Color', 'w');

% Plot v0 with a thicker line
plot(kgrid, v0, 'LineWidth', 2);
hold on;

% Plot v1 with a thicker line
plot(kgrid, v1, 'LineWidth', 2);

% Plot v2 with a thicker line
plot(kgrid, v2, 'LineWidth', 2);

% Plot v10 with a thicker line
plot(kgrid, v10, 'LineWidth', 2);

% Plot v50 with a thicker line
plot(kgrid, v50, 'LineWidth', 2);

% Plot v100 with a thicker line
plot(kgrid, v100, 'LineWidth', 2);

% Plot v226 with a thicker line
plot(kgrid, v226, 'LineWidth', 2);

% Plot converged v with a thicker line
plot(kgrid, v, 'LineWidth', 2);

% Label the axes with thicker text
xlabel('Capital Stock k today', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Value Function', 'FontSize', 12, 'FontWeight', 'bold');
title('True and Approximated Value Function', 'FontSize', 14, 'FontWeight', 'bold');
legend('v0', 'v1', 'v2', 'v10', 'v50', 'v100', 'v226', 'Converged v', 'FontSize', 10, 'FontWeight', 'bold');

% Add gridlines
grid on;

% Hold off to stop overlaying additional plots
hold off;

% Create a figure with a white background for the policy function
figure('Color', 'w');

% Plot kp0 with a thicker line
plot(kgrid, kp0, 'LineWidth', 2);
hold on;

% Plot kp1 with a thicker line
plot(kgrid, kp1, 'LineWidth', 2);

% Plot kp2 with a thicker line
plot(kgrid, kp2, 'LineWidth', 2);

% Plot kp10 with a thicker line
plot(kgrid, kp10, 'LineWidth', 2);

% Plot kp50 with a thicker line
plot(kgrid, kp50, 'LineWidth', 2);

% Plot kp100 with a thicker line
plot(kgrid, kp100, 'LineWidth', 2);

% Plot kp226 with a thicker line
plot(kgrid, kp226, 'LineWidth', 2);

% Plot converged kp with a thicker line
plot(kgrid, kp, 'LineWidth', 2);

% Label the axes with thicker text
xlabel('Capital Stock k today', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Policy Function', 'FontSize', 12, 'FontWeight', 'bold');
title('True and Approximated Policy Function', 'FontSize', 14, 'FontWeight', 'bold');
legend('v0', 'v1', 'v2', 'v10', 'v50', 'v100', 'v226', 'Converged v', 'FontSize', 10, 'FontWeight', 'bold');

% Add gridlines
grid on;

% Hold off to stop overlaying additional plots
hold off;

%%%%%%%%%%%%%%%%% steady state %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a figure with a white background
figure('Color', 'w');

% Plot the policy function with a thicker line
plot(kgrid, kp, 'LineWidth', 2);
hold on;

% Plot the 45-degree line with a thicker red line
plot(kgrid, kgrid, '--r', 'LineWidth', 2);

% Label the axes with a larger font size and thicker font
xlabel('Capital stock today (k)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Policy Function', 'FontSize', 16, 'FontWeight', 'bold');
title('Policy Function and Steady State', 'FontSize', 18, 'FontWeight', 'bold');

% Add gridlines
grid on;

% Find the steady-state level of capital
steady_state = interp1(kgrid, kp, kgrid, 'linear', 'extrap');

% Highlight the steady state with a smaller green circle
plot(steady_state, steady_state, 'og', 'MarkerSize', 4);

% Find the intersection with the 45-degree line
intersection = fminbnd(@(k) abs(k - interp1(kgrid, kp, k)), min(kgrid), max(kgrid));

% Highlight the intersection with a red "X"
plot(intersection, intersection, 'rx', 'MarkerSize', 8);

% Connect the intersection point to the axes with thicker dotted lines
plot([intersection, intersection], [0, intersection], ':k', 'LineWidth', 1.5);
plot([0, intersection], [intersection, intersection], ':k', 'LineWidth', 1.5);

% Display the steady-state value and intersection coordinates with a larger font size and thicker font
fprintf('Steady-State Level of Capital: %.4f\n', steady_state);
fprintf('Intersection coordinates: (%.4f, %.4f)\n', intersection, intersection);

% Set axis labels to be in bold and larger font
set(gca, 'FontWeight', 'bold', 'FontSize', 12);

% Legend with a larger font size and thicker font
legend('Policy Function', '45-degree line', 'Steady State', 'Intersection', 'FontSize', 10, 'FontWeight', 'bold');

% Hold off to stop overlaying additional plots
hold off;



%%%%%%%%%%%%%%%%% question(2) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = 20; % Number of periods
k0_values = [0.1, 1, 10]; % Initial values of k0 to test

% Initialize a matrix to store the time series of capital for each k0
k_time_series = zeros(T+1, length(k0_values));

for i = 1:length(k0_values)
    % Set the initial value of k0
    k0 = k0_values(i);
    
    % Initialize an array to store the time series of capital
    k_t = zeros(T+1, 1);
    
    % Set the initial capital value
    k_t(1) = k0;
    
    % Simulate the time series of capital using the policy function
    for t = 1:T
        % Use the policy function to determine the next period's capital
        k_t(t+1) = interp1(kgrid, kp, k_t(t), 'spline');
    end
    
    % Store the time series of capital in the matrix
    k_time_series(:, i) = k_t;
end
% Create a figure with a white background
figure('Color', 'w');

for i = 1:length(k0_values)
    subplot(length(k0_values), 1, i);
    plot(0:T, k_time_series(:, i), 'LineWidth', 2); % Adjust LineWidth as desired
    xlabel('Time Period', 'FontWeight', 'bold'); % Set FontWeight to 'bold'
    ylabel('Capital', 'FontWeight', 'bold'); % Set FontWeight to 'bold'
    title(['Time Series of Capital (k0 = ' num2str(k0_values(i)) ')'], 'FontWeight', 'bold'); % Set FontWeight to 'bold'
    grid on;
end

% Add a legend
legend('Location', 'Best', 'FontWeight', 'bold'); % Set FontWeight to 'bold'

% Adjust the subplot layout
sgtitle('Time Series of Capital for Different Initial Values of k0', 'FontWeight', 'bold'); % Set FontWeight to 'bold'
% Your previous code for plotting the figures

% Create a figure with a white background
figure('Color', 'w');

% Initialize a vector to store the steady-state line
steady_state_line = kss * ones(T+1, 1);

% Create a single subplot for all lines
subplot(1, 1, 1);

for i = 1:length(k0_values)
    plot(0:T, k_time_series(:, i), 'LineWidth', 2); % Adjust LineWidth as desired
    hold on; % Hold the plot for adding the steady-state line
end

% Plot the steady-state line
plot(0:T, steady_state_line, 'r--', 'LineWidth', 1); % Dotted line for steady state

xlabel('Time Period', 'FontWeight', 'bold'); % Set FontWeight to 'bold'
ylabel('Capital', 'FontWeight', 'bold'); % Set FontWeight to 'bold'
title('Time Series of Capital and Steady State', 'FontWeight', 'bold'); % Set FontWeight to 'bold'
grid on;

% Add a legend
legend('k0 = 0.1', 'k0 = 1', 'k0 = 10', 'Steady State', 'Location', 'Best', 'FontWeight', 'bold'); % Set FontWeight to 'bold'



% Save the first figure as a PNG file
figure(1);
saveas(gcf, 'D:\Econ\CUHKSZ\Macro I\code\Bellman_Equation\graph\figure1.png');

% Save the second figure as a PNG file
figure(2);
saveas(gcf, 'D:\Econ\CUHKSZ\Macro I\code\Bellman_Equation\graph\figure2.png');

% Save the third figure as a PNG file
figure(3);
saveas(gcf, 'D:\Econ\CUHKSZ\Macro I\code\Bellman_Equation\graph\figure3.png');

% Save the fourth figure as a PNG file
figure(4);
saveas(gcf, 'D:\Econ\CUHKSZ\Macro I\code\Bellman_Equation\graph\figure4.png');

% Save the fifth figure as a PNG file
figure(5);
saveas(gcf, 'D:\Econ\CUHKSZ\Macro I\code\Bellman_Equation\graph\figure5.png');

