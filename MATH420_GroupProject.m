%% AMSC 420 Group Homework 2
% Group: Robert "Eddie" Bull, Alexander Klein
clearvars
%% Question 1 Initializers
T = readtable("project5_data.xlsx");
% Change this to change the overall amount of data we use!
T_max = 40;

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


%% SIR Model
% Establish important params first
format shortG
h = 0.01;
initials = [N, I_t(1), 0];

% Make lists for the J values so we can get the minimum when it's important
J_1s = zeros([size(Om, 1), 2]);
J_2s = zeros([size(Om, 1), 2]);
J_infs = zeros([size(Om, 1), 2]);

% Initialize gamma and rho vectors
gammas = zeros(size(Om,1),3);
rhos = zeros(size(Om,1),3);

% First, optimize the alpha and beta
for i=1:size(Om, 1)
    syms f(gamma)
    % First, get the results for the simulation.
    set = Om(i, :);
    results = euler_SIR(set(1), set(2), initials, T_max, h, N);
    r = downsample(results, 1/h);
    
    % Then, for each pair, we calculate rho and gamma for each p value and
    % in the same loop, we calculate three J values for each p value

    % We'll do gamma first
    % Using the I^1 estimator for gamma_1:
    r_k = zeros(T_max-Tau_0,1);
    f_k = zeros(T_max-Tau_0,1);
    for k = 1:(T_max-Tau_0)
        r_k(k) = Y_t(k + Tau_0)/r(k,3);
        % Don't consider values of r that aren't in [0,1]
        if r_k(k) >= 0 && r_k(k) <= 1
            f_k(k,1) = sum(abs(transpose(Y_t(1,Tau_0:T_max)) - r_k(k)*r(1:T_max-Tau_0 + 1,3)));
        else
            f_k(k,1) = inf;
        end
    end
    [m,index] = min(f_k);
    gammas(i,1) = r_k(index);

    % Closed form solution for gamma when p = 2
    gammas(i,2) = sum(Y_t(Tau_0:T_max)'.*r(Tau_0:T_max, 3))./sum(r(Tau_0:T_max, 3).^2);

    % Slides have this using V_t, but I get the impression we should be
    % using R_sim.
    % Solving for gamma when p = inf. We will use a linear programming
    % approach.
    A = zeros(2*(T_max-Tau_0+1), 2);
    b = zeros(2*(T_max-Tau_0+1), 1);
    f = [1;0];
    for j = 1:size(A)
        A(j,1) = -1;
        A(j,2) = (-1)^mod(j,2)*r(ceil(j/2),3);
        b(j,1) = (-1)^mod(j,2)*Y_t(Tau_0 + ceil(j/2) - 1);
    end
    x = linprog(f,A,b);
    gammas(i,3) = x(2);


    % Now we'll optimize rho
    % rho_1 =
    % rho_2 =
    % rho_inf =

    % Precalculating the rho and gamma summations for later
    %rho_sum = [sum((I_t - rho_1 * I_sim(i, :))),...
    %sum((I_t - rho_2 * I_sim(i, :)).^2), ...
    %max((I_t - rho_inf * I_sim(i, :)))];

    %gamma_sum = [sum(Y_t - gamma_1 * R_sim(i, :)), ...
    %sum((Y_t - gamma_2 * R_sim(i, :)).^2), ...
    %max(Y_t - gamma_3 * R_sim(i, :))];

    % Then we calculate the two different cI and cY
    % J_1(1) is (0, 1) and J_1(2) is (1, 1)
    %J_1 = [gamma_sum(1), rho_sum(1) + gamma_sum(1)];
    %J_2 = [gamma_sum(2), rho_sum(2) + gamma_sum(2)];
    %J_inf = [gamma_sum(3), rho_sum(3) + gamma_sum(3)];
end

%% Functions

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

for t=(1 + step + step):step:T_max+1
    index = round((t - 1)/step);
    results(index, 1) = results(index - 1, 1) ...
        + step * dS(alpha, beta, results(index - 1, 1), results(index - 1, 2));

    results(index, 2) = results(index - 1, 2) ...
        + step * dI(alpha, beta, results(index - 1, 1), results(index - 1, 2));

    results(index, 3) = results(index - 1, 3) ...
        + step * dR(alpha, beta, results(index - 1, 1), results(index - 1, 2));
end
end
