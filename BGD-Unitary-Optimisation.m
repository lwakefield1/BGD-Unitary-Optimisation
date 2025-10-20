% Function to compute cost
computeCost = @(Uref, Utest) norm(Uref - Utest, 'fro')^2;

% Set a fixed seed for the random number generator
rng_seed = 1111;
rng(rng_seed);

% Prompt the user to enter the value of d (dimensions of matrix)
d = input("Enter the value of d (dimensions of matrix): ");

% Calculate the matrix size
n = d;

% Generate a random complex matrix
A = randn(n) + 1i * randn(n);
% Compute the QR decomposition of A
[Q, R] = qr(A);
% Construct the reference unitary matrix Uref
Uref = Q * diag(diag(R) ./ abs(diag(R)));

% Set the number of iterations
num_iter = 50;

% Prompt the user to enter the number of simulations
num_simulations = input("Enter the number of simulations: ");

learning_rate = 0.35; % Learning rate for gradient descent

% Initialise arrays to store results
avg_iterations = zeros(1, num_simulations);
best_costs = zeros(1, num_simulations);
best_iterations = zeros(1, num_simulations);
iteration_counts = zeros(1, num_simulations);
best_U_all = cell(1, num_simulations); % Store best unitary matrices
all_cost_history = zeros(num_iter, num_simulations);

% -------------------------------------------------------------------------
% Main simulation loop
for sim = 1:num_simulations

    % Generate a random real matrix P
    P = randn(d);

    % Initialise variables for this simulation
    best_cost = inf;
    best_U = [];
    best_iter = 0;
    cost_history = zeros(1, num_iter);
    cost_function_evaluations = 0; % Initialise the cost function evaluation counter

    % Loop over iterations to find the best unitary matrix
    for i = 1:num_iter
        % Compute Utest from P
        Utest = UC(P);

        % Compute cost
        cost = computeCost(Uref, Utest);

        % Update best cost and matrix
        if cost < best_cost
            best_cost = cost;
            best_U = Utest;
            best_iter = i;
        end

        % Store cost in history
        cost_history(i) = cost;

        % Compute the gradient of the cost function using finite differences
        epsilon = 1e-6; % Small epsilon for finite differences
        gradient = zeros(size(P));
        for row = 1:d
            for col = 1:d
                % Perturb the current element of P (real part)
                P_plus = P;
                P_plus(row, col) = P_plus(row, col) + epsilon;
                P_minus = P;
                P_minus(row, col) = P_minus(row, col) - epsilon;

                % Compute the corresponding Utest values (real part)
                Utest_plus = UC(P_plus);
                Utest_minus = UC(P_minus);

                % Compute the partial derivative approximation (real part)
                gradient(row, col) = (computeCost(Uref, Utest_plus) - computeCost(Uref, Utest_minus)) / (2 * epsilon);
            end
        end

        % Update P using gradient descent with learning rate
        P = P - learning_rate * gradient;

        % Update the cost function evaluation counter
        cost_function_evaluations = cost_function_evaluations + 1;

        % Check if cost goes below 0.01
        if cost < 0.01
            iteration_counts(sim) = cost_function_evaluations; % Store the number of cost function evaluations
            break;
        end
    end

    % Store results for this simulation
    avg_iterations(sim) = mean(find(cost_history > 0));
    best_costs(sim) = best_cost;
    best_iterations(sim) = best_iter;
    best_U_all{sim} = best_U; % Save best U for this simulation
    all_cost_history(:, sim) = cost_history;
end

% -------------------------------------------------------------------------
% Plot cost history for each simulation and the average
figure;
set(gca, 'FontSize', 14);
hold on;

for sim = 1:num_simulations
    if sim == 1
        plot(1:num_iter, all_cost_history(:, sim), '-o', 'Color', [0.7 0.7 0.7]);
    else
        plot(1:num_iter, all_cost_history(:, sim), '-o', 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');
    end
end
avg_cost_history = mean(all_cost_history, 2);
plot(1:num_iter, avg_cost_history, 'k-', 'LineWidth', 2);
hold off;
xlabel('Iteration', 'FontSize', 16);
ylabel('Cost', 'FontSize', 16);
title('Cost History for Each Simulation (BGD)', 'FontSize', 18);
legend('Seed 1111', sprintf('Average (%d Simulations)', num_simulations), 'Location', 'Best', 'FontSize', 14);
grid on;

% -------------------------------------------------------------------------
% Compute and plot percentiles
figure;
set(gca, 'FontSize', 14);
hold on;

percentile_10_history = prctile(all_cost_history, 10, 2);
percentile_50_history = prctile(all_cost_history, 50, 2);
percentile_90_history = prctile(all_cost_history, 90, 2);
best_percentile_history = min(all_cost_history, [], 2);
worst_percentile_history = max(all_cost_history, [], 2);

plot(1:num_iter, avg_cost_history, 'k-', 'LineWidth', 2);
for sim = 1:num_simulations
    plot(1:num_iter, all_cost_history(:, sim), '-o', 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');
end
plot(1:num_iter, percentile_10_history, 'r--', 'LineWidth', 1.5);
plot(1:num_iter, percentile_90_history, 'b--', 'LineWidth', 1.5);
plot(1:num_iter, best_percentile_history, 'g-.', 'LineWidth', 1.5);
plot(1:num_iter, worst_percentile_history, 'm-.', 'LineWidth', 1.5);
plot(1:num_iter, percentile_50_history, 'c-', 'LineWidth', 1.5);
hold off;

xlabel('Iteration', 'FontSize', 16);
ylabel('Value', 'FontSize', 16);
title('Average Cost and Percentile Range (BGD)', 'FontSize', 18);
legend('Seed 1111', 'Average Cost', '10th Percentile', '90th Percentile', ...
       'Best Percentile', 'Worst Percentile', 'Median', 'Location', 'Best', 'FontSize', 14);
grid on;

% -------------------------------------------------------------------------
% Plot iteration counts for each simulation
figure;
set(gca, 'FontSize', 14);
bar(iteration_counts);

avg_iteration = mean(iteration_counts);
avg_text = sprintf('Average: %.2f\nα=%.3f\nSeed=%d', avg_iteration, learning_rate, rng_seed);
annotation('textbox', [0.85 0.6 0.1 0.3], 'String', avg_text, ...
           'FitBoxToText', 'on', 'BackgroundColor', 'w', 'EdgeColor', 'k', 'FontWeight', 'bold');

xlabel('Simulation', 'FontSize', 16);
ylabel('Cost Function Evaluations', 'FontSize', 16);
title('Number of Cost Function Evaluations until Cost < 0.01 (BGD)', 'FontSize', 18);
grid on;

% -------------------------------------------------------------------------
% Display the best unitary matrix from all simulations
[~, best_sim_idx] = min(best_costs);
fprintf('\nBest simulation: #%d\n', best_sim_idx);
fprintf('Best cost achieved: %.10f\n', best_costs(best_sim_idx));

best_U_overall = best_U_all{best_sim_idx};
disp('Best Unitary Matrix (U_best):');
disp(best_U_overall);

% Optionally save the best unitary matrix
save('best_unitary_matrix.mat', 'best_U_overall');
fprintf('Saved best unitary matrix to best_unitary_matrix.mat\n');

% -------------------------------------------------------------------------
% Save cost history as CSV
data_to_save = [1:num_iter; avg_cost_history'];
csv_filename = 'cost_history_dataBGD5.csv';
csvwrite(csv_filename, data_to_save);
fprintf('\nData has been saved to %s.\n', csv_filename);

fprintf('Cost after %d iterations (average simulation): %.5e\n', num_iter, avg_cost_history(end));


