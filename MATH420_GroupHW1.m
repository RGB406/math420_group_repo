%% AMSC 420 Group Homework 1
% Group: Robert "Eddie" Bull, Alexander Klein

%% Question 1 Initializers
T = readtable("project5_data.xlsx");
V = table2array(T(2, 13:1103));
V_t = V(52:1091);

%% Question 1
N_max = 909327;
N_min = 1 + V_t(1040);
disp("t0 is 52")

% Below is Algorithm 1
B_hat_max = (6/(1040 * (1040 + 1) * (2*1040 + 1)));
B_hat_min = B_hat_max;

% Sum of I(t)
sum_max = 0;
sum_min = 0;
c_max = (N_max - V_t(1));
c_min = (N_min - V_t(1));
for t = 1:1040
    % Calculating the summation in both B_hat approximations.
    sum_max = sum_max + t * (log((V_t(t) * c_max) / (V_t(1) * (N_max - V_t(t)))));
    sum_min = sum_min + t * (log((V_t(t) * c_min) / (V_t(1) * (N_min - V_t(t)))));
end

% This gives us the final B_hat value for both N_max and N_min.
B_hat_max = B_hat_max * sum_max
B_hat_min = B_hat_min * sum_min

% Precalculate some of the expression to reduce runtime.
b_max = log(V_t(1) / (N_max - V_t(1)));
b_min = log(V_t(1) / (N_min - V_t(1)));
J_max = 0;
J_min = 0;

for t = 1:1040
    % Calculating the objective function summation.
    J_max = J_max + (B_hat_max * t - log(V_t(t) / (N_max - V_t(t))) - b_max);
    J_min = J_min + (B_hat_min * t - log(V_t(t) / (N_min - V_t(t))) - b_min);
end

% Display the calculated values.
J_max
J_min

NI_max = N_max * V_t(1);
NI_min = N_min * V_t(1);
hold on

plot(0:1039, V_t(1:1040), 'r-')
plot(0:1039, NI_max * (1./(V_t(1) + (N_max - V_t(1)).*exp(-1 * B_hat_max * (0:1039)))), 'g-')
plot(0:1039, NI_min * (1./(V_t(1) + (N_min - V_t(1)).*exp(-1 * B_hat_min * (0:1039)))), 'b-')
legend('Real', 'Estimated N_{max}', 'Estimated N_{min}')
title('Graph of I(t)')
hold off


%% Question 2