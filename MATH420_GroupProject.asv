%% AMSC 420 Group Homework 2
% Group: Robert "Eddie" Bull, Alexander Klein
clearvars
%% Question 1 Initializers
T = readtable("project5_data.xlsx");

Y = table2array(T(3, 13:1103));
V = table2array(T(2, 13:1103));


% Initializing parameters
N = 909327;
Tau_0 = 7;
V_min = 5;
p = [1, 2, inf];
c = [0, 1; 1, 1];

% Establishing Omega
R_0 = 0.8:0.05:2.2;
alpha = 0.05:0.01:0.4;
s_R = size(R_0, 2);
s_a = size(alpha, 2);

offset = 0;
% Actually creating Omega
Om = zeros([s_a * s_R, 2]);
for i = 1:s_a
    for j = 1:s_R
        Om(i + (offset * (s_R - 1)) + j - 1, 1) = alpha(i);
        Om(i + (offset * (s_R - 1)) + j - 1, 2) = alpha(i) * R_0(j);
    end
    offset = offset + 1;
end


%% SIR Model
% Establish important params first
format shortG
h = 0.01;
initials = [N, I_t(1), 0];
C_s = [0 1; 1 1];
T_max = 120;


% Setting infections and deaths
Y_t = Y(52:(T_max + 51));
V_t = V(52:(T_max + 51));

% Setting infection rates
I_t = V((52:(T_max + 51)) + Tau_0) - V((52:(T_max + 51)) - Tau_0);

% Make lists for the J values so we can get the minimum when it's important
J_1s = zeros([size(Om, 1), size(C_s, 1)]);
J_2s = zeros([size(Om, 1), size(C_s, 1)]);
J_infs = zeros([size(Om, 1), size(C_s, 1)]);

% Initialize gamma and rho vectors
gammas = zeros(size(Om,1),3);
rhos = zeros(size(Om,1),3);

for i=1:size(Om, 1)
    % syms f(gamma)
    % First, get the results for the simulation.
    set = Om(i, :);
    r = euler_SIR(set(1), set(2), initials, T_max, h, N);

    % Then, for each pair, we calculate rho and gamma for each p value and
    % in the same loop, we calculate three J values for each p value

    % We'll do gamma first
    gammas(i, :) = optimize_hat_vals(Y_t, r(:, 3), T_max);

    % Now we'll optimize rho
    rhos(i, :) = optimize_hat_vals(I_t, r(:, 2), T_max);

    % Precalculating the rho and gamma summations for later
    rho_sum = [sum(abs(I_t - rhos(i, 1) * r(:, 2))),...
        sum((I_t - rhos(i, 2) * r(:, 2)).^2), ...
        max(abs(I_t - rhos(i, 3) * r(:, 2)))];

    gamma_sum = [sum(abs(Y_t - gammas(i, 1) * r(:, 3))), ...
        sum((Y_t - gammas(i, 2) * r(:, 3)).^2), ...
        max(abs(Y_t - gammas(i, 3) * r(:, 3)))];

    % Then we calculate the two different cI and cY
    % Cool method that lets us redefine C_s with no consequence!
    % Multiplies rho and gamma sums by the cI and cY values defined
    % above in C_s and turns it into a row vector.
    J_1 = sum(C_s .* [rho_sum(1), gamma_sum(1)], 2)';
    J_2 = sum(C_s .* [rho_sum(2), gamma_sum(2)], 2)';
    J_inf = sum(C_s .* [rho_sum(3), gamma_sum(3)], 2)';
    
    J_1s(i, :) = J_1;
    J_2s(i, :) = J_2;
    J_inf(i, :) = J_inf;
end

% After the loop, we have all that we need to start getting data.
% So now let's find the minimum value and index of the J values.

% Returns two row vectors denoting the minimum values and indices of J_p
[val_1, ind_1] = min(J_1s);
[val_2, ind_2] = min(J_2s);
[val_inf, ind_inf] = min(J_infs);

% Now we have all the data we need for T_max = 120

rho_hat_1 = rhos(ind_1(1), 1);
res_1 = euler_SIR(Om(ind_1(1)), Om(ind_1(2)), initials, T_max, h, N);

%plot_graph(res_1(2), I_t(:), res_1(3), Y_t, rhos, gammas, ind_1, "1")

%% Functions

% We'll use this function to calculate the optimal values for rho hat and
% gamma hat to ensure our code is prettier.
function optimals = optimize_hat_vals(observed, simulated, T_max)
options = optimoptions('linprog','Display','none');
% Row vector of optimized values given observed and simulated arrays
optimals = zeros(1, 3);
r_k = zeros(T_max,1);
f_k = zeros(T_max,1);
for k = 1:(T_max)
    r_k(k) = observed(k)/simulated(k);
    % Don't consider values of r that aren't in [0,1]
    if r_k(k) >= 0 && r_k(k) <= 1
        f_k(k,1) = sum(abs(transpose(observed(1,1:T_max)) - r_k(k)*simulated(1:T_max)));
    else
        f_k(k,1) = inf;
    end
end
[~,index] = min(f_k);
optimals(1) = r_k(index);

% Trying closed form solution with I_t and I_sim
optimals(2) = min(sum(observed(1:T_max)'.*simulated(1:T_max))./sum(simulated(1:T_max).^2), 1);

% And then the linear program with I_t. This doesn't seem to work.
A = [-1 * ones(2*(T_max), 1), zeros(2*(T_max), 1)];
b = zeros(2*(T_max), 1);
f = [1;0];
for j = 1:(size(A, 1))
    A(j,2) = (-1)^mod(j,2)...
        *simulated(ceil(j/2));
    b(j,1) = (-1)^mod(j,2)*observed(ceil(j/2));
end
x = linprog(f,A,b, [], [], [0; 0], [inf; 1], options);
optimals(3) = x(2);

end

% This is the euler function we'll use for the different alphas and betas
% Should take inits as S I R and functions as S I R
function results = euler_SIR(alpha, beta, inits, T_max, step, N)
% Initialize diffeqs to pretty up code
dS = @(a, b, S, I) -1 * b * S * (I/N);
dI = @(a, b, S, I) b * S * (I/N) - a * I;
dR = @(a, b, S, I) a * I;

% S_sim, I_sim, R_sim
results = zeros([T_max/step, 3]);
results(1, :) = inits;
index = 1;


for t=step:step:T_max-1
    index = index + 1;
    % Assign to a temporary variable for easy reading...
    S = results(index - 1, 1);
    I = results(index - 1, 2);
    R = results(index - 1, 3);

    % Calculate using the defined function above, and store to array
    results(index, 1) = S + step * dS(alpha, beta, S, I);
    results(index, 2) = I + step * dI(alpha, beta, S, I);
    results(index, 3) = R + step * dR(alpha, beta, S, I);
end

results = downsample(results, 1/step);
end

function plot_graph(sim_field_I, field_I, R_sim, Y, rhos, gammas, indices, p_val)

rho = rhos(indices(1));
gamma = gammas(indices(1));

figure
subplot(2, 1, 1)
hold on
% This plots the case where it's (0, 1)
plot(rho * sim_field_I(indices(1), :), 'r-')
plot(field_I(:), 'b-')
% Now print the Y(t) and R_sim functions
plot(Y, 'g-')
plot(gamma * R_sim(indices(1), :), 'm-')
axis tight

hold off
title("I_{sim} vs I for p = " + p_val + ", (0, 1)")

rho = rhos(indices(2));
gamma = gammas(indices(2));

subplot(2, 1, 2)
hold on
% This plots the case where it's (1, 1)
plot(rho * sim_field_I(indices(2), :), 'r-')
plot(field_I(:), 'b-')
% Now print the Y(t) and R_sim functions
plot(Y, 'g-')
plot(gamma * R_sim(indices(2), :), 'm-')
axis tight

hold off
title("I_{sim} vs I for p = " + p_val + ", (1, 1)")
end
