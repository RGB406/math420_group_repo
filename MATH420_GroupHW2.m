%% AMSC 420 Group Homework 2
% Group: Robert "Eddie" Bull, Alexander Klein
clearvars
%% Question 1 Initializers
T = readtable("project5_data.xlsx");
T_max = 120;

Y = table2array(T(3, 13:1103));
V = table2array(T(2, 13:1103));
Y_t = Y(52:(T_max + 51));
V_t = V(52:(T_max + 51));

% Initializing parameters
N = 909327;
Tau_0 = 7;
V_min = 5;
p = [1, 2, inf];
c = [0, 1; 1, 1];

% Setting I(t)
I_t = V((52:(T_max + 51)) + Tau_0) - V((52:(T_max + 51)) - Tau_0);

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

%% Question 1
% First, establish parameters.
format short
h = 0.01;
initials = [N, I_t(1), 0];

S_sim = zeros(size(Om, 1), T_max);
I_sim = zeros(size(Om, 1), T_max);
R_sim = zeros(size(Om, 1), T_max);

% Get SIR sim
for i=1:size(Om, 1)
    set = Om(i, :);
    results = euler_SIR(set(1), set(2), initials, T_max, h, N) * h;
    r = downsample(results, 1/h);

    S_sim(i, :) = r(:, 1);
    I_sim(i, :) = r(:, 2);
    R_sim(i, :) = r(:, 3);
end


% First, estimate gamma and rho for each p
for p_i=p
    min_func_gamma = @(g) norm(Y_t - g * R_sim(:, :), p_i);
    fminsearch(min_func_gamma, 0)

    min_func_rho = @(rho) norm(I_t - rho * I_sim(:, :), p_i);
    fminsearch(min_func_rho, 0)

end



% This is the euler function we'll use for the different alphas and betas
% Should take inits as S I R and functions as S I R
function results = euler_SIR(alpha, beta, inits, T_max, step, N) 
dS = @(a, b, S, I) -1 * b * S * (I/N);
dI = @(a, b, S, I) b * S * (I/N) - a * I;
dR = @(a, b, S, I) a * I;

% S_sim, I_sim, R_sim
results = zeros([T_max/step, 3]);
results(1, :) = inits;

for t=1.02:step:T_max+1
    index = round((t - 1)/step);
    results(index, 1) = results(index - 1, 1) ...
        + step * dS(alpha, beta, results(index - 1, 1), results(index - 1, 2));

    results(index, 2) = results(index - 1, 2) ...
        + step * dI(alpha, beta, results(index - 1, 1), results(index - 1, 2));

    results(index, 3) = results(index - 1, 3) ...
        + step * dR(alpha, beta, results(index - 1, 1), results(index - 1, 2));
end
end
