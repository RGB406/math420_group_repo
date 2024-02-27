%% AMSC 420 Group Homework 2
% Group: Robert "Eddie" Bull, Alexander Klein

%% Question 1 Initializers
T = readtable("project5_data.xlsx");
T_max = 120;

Y = table2array(T(3, 13:1103));
V = table2array(T(2, 13:1103));
Y_t = Y(52:(T_max + 52));
V_t = V(52:(T_max + 52));

% Initializing parameters
N = 909327;
Tau_0 = 7;
V_min = 5;
p = [1, 2, inf];
c = [0, 1; 1, 1];

% Setting I(t)
I_t = V((52:(T_max + 52)) + Tau_0) - V((52:(T_max + 52)) - Tau_0);

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
        Om(i + (offset * (s_R - 1)) + j, 1) = alpha(i);
        Om(i + (offset * (s_R - 1)) + j, 2) = alpha(i) * R_0(j);
    end
    offset = offset + 1;
end

%% Question 1
% First, establish the diffeqs.
h = 0.01;
initials = [N, I_t(1), 0];

results = euler_SIR(0.2, 0.4, initials, T_max, h, N) * h;
r = downsample(results, 1/h);


subplot(3, 1, 1)
plot(r(:, 1), 'r-')
subplot(3, 1, 2)
hold on
plot(r(:, 2), 'b-')
plot(I_t(:), 'g-')
hold off
subplot(3, 1, 3)
plot(r(:, 3), 'g-')

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
    R = results(index - 1, 1) ...
        + step * dS(alpha, beta, results(index - 1, 1), results(index - 1, 2));
    results(index, 1) = R;

    results(index, 2) = results(index - 1, 2) ...
        + step * dI(alpha, beta, results(index - 1, 1), results(index - 1, 2));

    results(index, 3) = results(index - 1, 3) ...
        + step * dR(alpha, beta, results(index - 1, 1), results(index - 1, 2));
end
end

