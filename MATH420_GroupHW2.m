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
format shortG
h = 0.01;
initials = [N, I_t(1), 0];

S_sim = zeros(size(Om, 1), T_max);
I_sim = zeros(size(Om, 1), T_max);
R_sim = zeros(size(Om, 1), T_max);

% Get SIR simulated
for i=1:size(Om, 1)
    set = Om(i, :);
    results = euler_SIR(set(1), set(2), initials, T_max, h, N);
    r = downsample(results, 1/h);

    % SIR saved in individual variables for ease of access later.
    S_sim(i, :) = r(:, 1);
    I_sim(i, :) = r(:, 2);
    R_sim(i, :) = r(:, 3);
end

rhos = zeros(3, 1);
gammas = zeros(3, 1);
% First, estimate gamma and rho for each p
for j=1:3
    % We use fminsearch here to get the variables.
    % I'm sorry. I really tried to use the numeric methods for this.
    % Symbolic variables absolutely ATE my processor and took eons to load.
    % I don't regret this.
    min_func_gamma = @(g) norm(Y_t - g * R_sim(:, :), p(j));
    gammas(j) = fminsearch(min_func_gamma, 0);

    min_func_rho = @(rho) norm(I_t - rho * I_sim(:, :), p(j));
    % Restrain the coefficient to be maxed at 1.
    % Since our N is so absolutely massive for our dataset, this kinda
    % tracks.
    rhos(j) = min(fminsearch(min_func_rho, 0), 1);
end
[rhos, gammas];

% Now, let's get the J function using 0, 1 and 1, 1
% We'll use multiple lists to find the minimum.
% This will be ugly.
J_1s = zeros([size(Om, 1), 2]);
J_2s = zeros([size(Om, 1), 2]);
J_infs = zeros([size(Om, 1), 2]);


for j=1:size(Om, 1)
% First, calculate the sums so we can reference them later for when we need
% to use all three p-vals.
rho_sum = [sum((I_t - rhos(1) * I_sim(j, :))),...
    sum((I_t - rhos(2) * I_sim(j, :)).^2), ...
    max((I_t - rhos(3) * I_sim(j, :)))];
% We'll just precompute these for ease of reading.
gamma_sum = [sum(Y_t - gammas(1) * R_sim(j, :)), ...
    sum((Y_t - gammas(2) * R_sim(j, :)).^2), ...
    max(Y_t - gammas(3) * R_sim(j, :))];

% Then we calculate the two different cI and cY
% J_1(1) is (0, 1) and J_1(2) is (1, 1)
J_1 = [gamma_sum(1), rho_sum(1) + gamma_sum(1)];
J_2 = [gamma_sum(2), rho_sum(2) + gamma_sum(2)];
J_inf = [gamma_sum(3), rho_sum(3) + gamma_sum(3)];

% Putting the J functions for each pair into the matrix. This way, we
% should be able to store 
J_1s(j, :) = abs(J_1);
J_2s(j, :) = abs(J_2);
J_infs(j, :) = abs(J_inf);
end


% Now plot some of the surfaces!   
figure
subplot(2, 1, 1)
surf(reshape(J_1s(:, 1), 29, 36))
title('Surface for J_1 using p = 1 and (0, 1)')

subplot(2, 1, 2)
surf(reshape(J_1s(:, 2), 29, 36))
title('Surface for J_1 using p = 1 and (1, 1)')

figure
subplot(2, 1, 1)
surf(reshape(J_2s(:, 1), 29, 36))
title('Surface for J_2 using p = 2 and (0, 1)')

subplot(2, 1, 2)
surf(reshape(J_2s(:, 2), 29, 36))
title('Surface for J_2 using p = 2 and (1, 1)')

figure
subplot(2, 1, 1)
surf(reshape(J_infs(:, 1), 29, 36))
title('Surface for J_inf using p = inf and (0, 1)')

subplot(2, 1, 2)
surf(reshape(J_infs(:, 2), 29, 36))
title('Surface for J_inf using p = inf and (1, 1)')

% Now lets find alphas and betas for the min values
% This probably returns more than one value for the minimum. If it does,
% we'll just use the minimum later.
[~, min_1] = min(J_1s);
[~, min_2] = min(J_2s);
[~, min_inf] = min(J_infs);
% Print the data for 1.1
[m,i] = min(J_1s(:,1));
disp("[min_J,a_hat,b_hat,r0_hat,g_hat,rho_hat] for (c_I,c_Y) = (0,1), p = 1 given by: ")
[m, Om(i,1), Om(i,2), Om(i,2)/Om(i,1), gammas(1), rhos(1)]

[m,i] = min(J_2s(:,1));
disp("[min_J,a_hat,b_hat,r0_hat,g_hat,rho_hat] for (c_I,c_Y) = (0,1), p = 2 given by: ")
[m, Om(i,1), Om(i,2), Om(i,2)/Om(i,1), gammas(2), rhos(2)]

[m,i] = min(J_infs(:,1));
disp("[min_J,a_hat,b_hat,r0_hat,g_hat,rho_hat] for (c_I,c_Y) = (0,1), p = inf given by: ")
[m, Om(i,1), Om(i,2), Om(i,2)/Om(i,1), gammas(3), rhos(3)]

[m,i] = min(J_1s(:,2));
disp("[min_J,a_hat,b_hat,r0_hat,g_hat,rho_hat] for (c_I,c_Y) = (1,1), p = 1 given by: ")
[m, Om(i,1), Om(i,2), Om(i,2)/Om(i,1), gammas(1), rhos(1)]

[m,i] = min(J_2s(:,2));
disp("[min_J,a_hat,b_hat,r0_hat,g_hat,rho_hat] for (c_I,c_Y) = (1,1), p = 2 given by: ")
[m, Om(i,1), Om(i,2), Om(i,2)/Om(i,1), gammas(2), rhos(2)]

[m,i] = min(J_infs(:,2));
disp("[min_J,a_hat,b_hat,r0_hat,g_hat,rho_hat] for (c_I,c_Y) = (1,1), p = inf given by: ")
[m, Om(i,1), Om(i,2), Om(i,2)/Om(i,1), gammas(3), rhos(3)]


plot_graph(I_sim, I_t, R_sim, Y_t, rhos(1), gammas(1), min_1, "1")
plot_graph(I_sim, I_t, R_sim, Y_t, rhos(2), gammas(2), min_2, "2")
plot_graph(I_sim, I_t, R_sim, Y_t, rhos(3), gammas(3), min_inf, "inf")


%% Problem 2 Initialization
clearvars
format shortG
T = readtable("project5_data.xlsx");
Tau_0 = 7;
T_max = 120;
V_min = 5;
p = [1, 2, inf];
c = [0, 1; 1, 1];
Y = table2array(T(3, 13:1103));
V = table2array(T(2, 13:1103));
Y_t = Y(52:(T_max + 51));
V_t = V(52:(T_max + 51));
I_t = V((52:(T_max + 51)) + Tau_0) - V((52:(T_max + 51)) - Tau_0);
N = 909327;

%Set initial conditions
S_0 = N;
R_0 = 0;
I_0 = I_t(1);
E_0 = I_0;

%Define variables that make up the set Omega.
r_0 = 0.8:0.05:2.2;
alpha = 0.05:0.01:0.4;
delta = 0.05:0.01:0.4;
gammas = zeros(3,1);
rhos = zeros(3,1);
%% Problem 2 part 1 

% Compute numerical approximations for S,E,I,R
Om = zeros(size(alpha,2)*size(delta,2)*size(r_0,2),3);
S_sim = zeros(size(Om,1),120);
E_sim = zeros(size(Om,1),120);
I_sim = zeros(size(Om,1),120);
R_sim = zeros(size(Om,1),120);
%Populate the Om matrix and SEIR sim matrices
vec = zeros(4,T_max);
index = 1;
for a = 1:36
    for r = 1:29
        beta = r_0(r) * alpha(a);
        for d = 1:36
            Om(index,:) = [alpha(a),beta,delta(d)];
            % Store the results of the euler scheme into our sim matrices
            % Each row represents a different alpha,beta,delta pair from
            % the set
            vec = euler_SIER(alpha(a),beta,delta(d),T_max,0.05,N,[S_0,E_0,I_0,R_0]);
            S_sim(index,:) = vec(:,1);
            E_sim(index,:) = vec(:,2);
            I_sim(index,:) = vec(:,3);
            R_sim(index,:) = vec(:,4);
            index = index + 1;
        end
    end
end


% Find optimal gamma and rho for each p-value:for i = 1:3
for i = 1:3
    min_func_gamma = @(g) norm(Y_t - g * R_sim(:, :), p(i));
    gammas(i) = fminsearch(min_func_gamma, 0);

    min_func_rho = @(rho) norm(I_t - rho * I_sim(:, :), p(i));
    rhos(i) = fminsearch(min_func_rho, 0);
end
[rhos, gammas];

% Compute the values for Objective function J:
% c_I,c_Y = (0,1)
J_1 = zeros(size(Om,1),3);
% c_I,c_Y = (1,1)
J_2 = zeros(size(Om,1),3);

for i = 1:size(Om)
    J_1(i,1) = norm(Y_t-gammas(1)*R_sim(i,:),p(1));
    J_1(i,2) = norm(Y_t-gammas(2)*R_sim(i,:),p(2));
    J_1(i,3) = norm(Y_t-gammas(3)*R_sim(i,:),p(3));

    J_2(i,1) = norm(I_t-rhos(1)*I_sim(i,:),p(1)) + norm(Y_t-gammas(1)*R_sim(i,:),p(1));
    J_2(i,2) = norm(I_t-rhos(2)*I_sim(i,:),p(2)) + norm(Y_t-gammas(2)*R_sim(i,:),p(2));
    J_2(i,3) = norm(I_t-rhos(3)*I_sim(i,:),p(3)) + norm(Y_t-gammas(3)*R_sim(i,:),p(3));
end

%store the optimal values for later
hat_1 = zeros(3,3);
hat_2 = zeros(3,3);
%Printing information for problem 2.1
for j = 1:3
    disp("[min_J,a_hat,b_hat,d_hat,r0_hat,g_hat,rho_hat] for (c_I,c_Y) = (0,1), p = " + j + " given by: ")
    [m i] = min(J_1(:,j));
    hat_1(j,:) = [Om(i,1), Om(i,2), Om(i,3)];
    [m, Om(i,1), Om(i,2), Om(i,3), Om(i,2)/Om(i,1), gammas(j), rhos(j)]
end

for j = 1:3
    disp("[min_J,a_hat,b_hat,d_hat,r0_hat,g_hat,rho_hat] for (c_I,c_Y) = (1,1), p = " + j + " given by: ")
    [m i] = min(J_2(:,j));
    hat_2(j,:) = [Om(i,1), Om(i,2), Om(i,3)];
    [m, Om(i,1), Om(i,2), Om(i,3), Om(i,2)/Om(i,1), gammas(j), rhos(j)]
end


%% Problem 2, part 2
% Plotting the surfaces
% First we need to compute sub-matrices of Om that store the alpha-beta
% pairs we want:
J_1_ab = zeros(size(alpha,2),size(r_0,2));
J_2_ab = zeros(size(alpha,2),size(r_0,2));
tolerance = 0.0001;
% Surfing over the function (alpha,beta) --> J
for j = 1:3
    for i = 1:size(Om,1)
        %build J_ab for each p-value (1 to 3)
        if(Om(i,3) == hat_1(j,3))
            J_1_ab(find(alpha == Om(i,1)), find(abs(r_0 - (Om(i,2)/Om(i,1))) < tolerance)) = J_1(i,j);
        end
        if(Om(i,3) == hat_2(j,3))
            J_2_ab(find(alpha == Om(i,1)), find(abs(r_0 - (Om(i,2)/Om(i,1))) < tolerance)) = J_2(i,j);
        end
    end
    %surf over each J_ab
    figure();
    subplot(2, 1, 1)
    surf(J_1_ab)
    title("(alpha,beta) --> J for p = " + p(j) + " (c_I,c_Y) = (0,1)")
    subplot(2, 1, 2)
    surf(J_2_ab)
    title("(alpha,beta) --> J for p = " + p(j) + " (c_I,c_Y) = (1,1)")
end

% Surfing over the function (alpha,delta) --> J
J_1_ad = zeros(size(alpha,2),size(delta,2));
J_2_ad = zeros(size(alpha,2),size(delta,2));
for j = 1:3
    for i = 1:size(Om,1)
        %build J_ad for each p-value (1 to 3)
        if(abs(Om(i,2)/Om(i,1) - hat_1(j,2)/hat_1(j,1)) < tolerance)
            J_1_ad(find(alpha == Om(i,1)), find(delta == Om(i,3))) = J_1(i,j);
        end
        if(abs(Om(i,2)/Om(i,1) - hat_2(j,2)/hat_2(j,1)) < tolerance)
            J_2_ad(find(alpha == Om(i,1)), find(delta == Om(i,3))) = J_2(i,j);
        end
    end
    %surf over each J_ad
    figure();
    subplot(2, 1, 1)
    surf(J_1_ad)
    title("(alpha,delta) --> J for p = " + p(j) + " (c_I,c_Y) = (0,1)")
    subplot(2, 1, 2)
    surf(J_2_ad)
    title("(alpha,delta) --> J for p = " + p(j) + " (c_I,c_Y) = (1,1)")
end

% Surfing over the function (beta,delta) --> J
J_1_bd = zeros(size(r_0,2),size(delta,2));
J_2_bd = zeros(size(r_0,2),size(delta,2));
for j = 1:3
    for i = 1:size(Om,1)
        %build J_bd for each p-value (1 to 3)
        if(Om(i,1) == hat_1(j,1))
            J_1_bd(find(abs(r_0 - (Om(i,2)/Om(i,1))) < tolerance), find(delta == Om(i,3))) = J_1(i,j);
        end
        if(Om(i,1) == hat_2(j,1))
            J_2_bd(find(abs(r_0 - (Om(i,2)/Om(i,1))) < tolerance), find(delta == Om(i,3))) = J_2(i,j);
        end
    end
    %surf over each J_bd
    figure();
    subplot(2, 1, 1)
    surf(J_1_bd)
    title("(beta,delta) --> J for p = " + p(j) + " (c_I,c_Y) = (0,1)")
    subplot(2, 1, 2)
    surf(J_2_bd)
    title("(beta,delta) --> J for p = " + p(j) + " (c_I,c_Y) = (1,1)")
end

%% Problem 2, part 3
% Iterate through all 3 p-values:
for j = 1:3
    % Plot the I(t), Y(t) data for (0,1) case
    [m i] = min(J_1(:,j));
    figure()
    subplot(2, 1, 1);
    plot(I_t)
    hold on;
    plot(I_sim(i,:));
    title("(0,1) plot of I(t) vs I_sim for p = " + p(j))
    subplot(2, 1, 2);
    plot(Y_t);
    hold on;
    plot(gammas(j)*R_sim(i,:));
    title("(0,1) plot of Y(t) vs Y_sim for p = " + p(j))

    % (1,1) case
    [m i] = min(J_2(:,j));
    figure()
    subplot(2, 1, 1);
    plot(I_t)
    hold on;
    plot(I_sim(i,:));
    title("(1,1) plot of I(t) vs I_sim for p = " + p(j))
    subplot(2, 1, 2);
    plot(Y_t);
    hold on;
    plot(gammas(j)*R_sim(i,:));
    title("(1,1) plot of Y(t) vs Y_sim for p = " + p(j))
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

% Euler function for problem 2 (SIER modelling)
function results = euler_SIER(alpha,beta,delta,T_max,step,N,inits)
% Initialize the diffeqs and results matrix
S = inits(1);
E = inits(2);
I = inits(3);
R = inits(4);
results = zeros(T_max,4);

% Run the euler scheme
for t = 1:step:T_max
    S = S + step*((-1 * beta * S * I / N));
    E = E + step*((beta * S * I / N) - (delta * E));
    I = I + step*((delta * E) - (alpha * I));
    R = R + step*(alpha*I);
    results(round(t),1) = S;
    results(round(t),2) = E;
    results(round(t),3) = I;
    results(round(t),4) = R;
end
end

% This function should be used for the last part of question 1.
% Mostly just to cut down on space.
function plot_graph(sim_field_I, field_I, R_sim, Y, rho, gamma, indices, p_val)
    
figure
    subplot(2, 1, 1)
    hold on
    % This plots the case where it's (0, 1)
    plot(rho * sim_field_I(indices(1), :), 'r-')
    plot(field_I(:), 'b-')
    % Now print the Y(t) and R_sim functions
    plot(Y, 'g-')
    plot(gamma * R_sim(indices(1), :), 'm-')

    hold off
    title("I_{sim} vs I for p = " + p_val + ", (0, 1)")

    subplot(2, 1, 2)
    hold on
    % This plots the case where it's (1, 1)
    plot(rho * sim_field_I(indices(2), :), 'r-')
    plot(field_I(:), 'b-')
    % Now print the Y(t) and R_sim functions
    plot(Y, 'g-')
    plot(gamma * R_sim(indices(1), :), 'm-')

    hold off
    title("I_{sim} vs I for p = " + p_val + ", (1, 1)")
end
