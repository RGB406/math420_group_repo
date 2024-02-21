%% AMSC 420 Group Homework 1
% Group: Robert "Eddie" Bull, Alexander Klein

%% Question 1 Initializers
T = readtable("project5_data.xlsx");
V = table2array(T(2, 13:1103));
V_t = V(52:(119 + 51));

%% Question 1, part a
N_max = 909327;
N_min = 1 + V_t(119);
a = (6/(119 * (119 + 1) * (2*119 + 1)));
disp("t0 is 52")

% Below is Algorithm 1
B_hat_max = a;
B_hat_min = B_hat_max;

% Sum of I(t)
sum_max = 0;
sum_min = 0;
c_max = (N_max - V_t(1));
c_min = (N_min - V_t(1));
for t = 1:119
    % Calculating the summation in both B_hat approximations.
    sum_max = sum_max + t * (log((V_t(t) * c_max) / (V_t(1) * (N_max - V_t(t)))));
    sum_min = sum_min + t * (log((V_t(t) * c_min) / (V_t(1) * (N_min - V_t(t)))))
end

% This gives us the final B_hat value for both N_max and N_min.
B_hat_max = B_hat_max * sum_max
B_hat_min = B_hat_min * sum_min

% Precalculate some of the expression to reduce runtime.
b_max = log(V_t(1) / (N_max - V_t(1)));
b_min = log(V_t(1) / (N_min - V_t(1)));
J_max = 0;
J_min = 0;

for t = 1:119
    % Calculating the objective function summation.
    J_max = J_max + (B_hat_max * t - log(V_t(t) / (N_max - V_t(t))) - b_max)^2;
    J_min = J_min + (B_hat_min * t - log(V_t(t) / (N_min - V_t(t))) - b_min)^2;
end

% Display the calculated values.
J_max
J_min

% Precalculate NI values for N_max and N_min
NI_max = N_max * V_t(1);
NI_min = N_min * V_t(1);

hold on
% First, plot I(t)
plot(0:118, V_t(1:119), 'r-')
% Then, plot the estimated value using the N_max pop
plot(0:118, NI_max * (1./(V_t(1) + (N_max - V_t(1)).*exp(-1 * B_hat_max * (0:118)))), 'g-')
% Finally, plot with the N_min estimation
plot(0:118, NI_min * (1./(V_t(1) + (N_min - V_t(1)).*exp(-1 * B_hat_min * (0:118)))), 'b-')
axis tight
legend({'Real', 'Estimated N_{max}', 'Estimated N_{min}'}, 'Location', 'northwest')
title('Graph of I(t)')
hold off

%% Part b
N = N_min;
J_old = inf;
J_N = 0;
% These two variables are just used to keep data from this part.
J_N_prev = 0;
N_prev = 0;

init = V_t(1);
a = (6/(119 * (119 + 1) * (2*119 + 1)));
% Assigning T to be a vector.
t = 1:119;
J_N_vals = 0;

while true
    % This is algorithm 2, step 2.1.
    first = sum( abs( log( V_t(t) * (N - init) ./ ( init * (N - V_t(t)) ) )).^2);
    second = sum( t .* log( V_t(t) * (N - init) ./ ( init * (N - V_t(t)) ) ))^2;
    J_N = first - a* second;

    if (J_N >= J_old)
        break;
    end

    % Keep track of J_N values for later.
    J_N_vals = [J_N_vals J_N];
    J_old = J_N;
    N = N+1;
end

% Plot all of the intermediary J_N values
plot(0:(size(J_N_vals, 2)-2), J_N_vals(2:(size(J_N_vals, 2))))
title('Intermediate J(N) values')

% Calculate beta hat
B_hat = a * sum(t * (log((V_t(t) * c_max) / (V_t(1) * (N - V_t(t))))));

% Print N hat (N estimate) and B hat
disp(append("B_hat estimate is: ", string(B_hat)))
disp(append("N estimate is: ", string(N)))

% Get the objective function J(beta, N) below
b = log(V_t(1) / (N - V_t(1)));
J_obj = sum((B_hat * t - log(V_t(t) / (N - V_t(t))) - b).^2);
disp(append("Objective function has value: ", string(J_obj)))

% Plot I(t) vs our estimate
figure
hold on
NI = N * init;
% Plot I(t) first...
plot(0:118, V_t(1:119), 'r-')
% ...and our estimate second!
plot(0:118, NI * (1./(V_t(1) + (N - V_t(1)).*exp(-1 * B_hat * (0:118)))), 'g-')
legend({'Real', 'Estimated'}, 'Location', 'southeast')
axis tight
title('Graph of I(t)')
hold off

% Saving these for later parts.
J_N_prev = J_N;
N_prev = N;

%% Part b, IV
N = N_min;
J_old = inf;
J_N = 0;
length = N_max - N;
J_N_vals = zeros(1, length);


while N < N_max
    % This is algorithm 2, step 2.1.
    first = sum( abs( log( V_t(t) * (N - init) ./ ( init * (N - V_t(t)) ) )).^2);
    second = sum( t .* log( V_t(t) * (N - init) ./ ( init * (N - V_t(t)) ) ))^2;
    J_N = first - a* second;

    if (J_N < J_old)
        J_old = J_N;
    end

    % Keep track of J_N values for later.
    J_N_vals(N - N_max + length + 1) = J_N;

    N = N+1;
end

% Plot all of the intermediary J_N values.
plot(0:(size(J_N_vals, 2)-2), J_N_vals(2:(size(J_N_vals, 2))))
title('Intermediate J(N) values')

% Printing and checking the minimum J_N value.
disp(append("Global minimum for J_N values: ", string(min(J_N_vals))))
disp(append("The previous estimated min was: ", string(J_N_prev)))
% The results seem to check out with the last part, so there doesn't seem
% to be much of a change, no.

%% Question One, part 2, sub c
N_vals = N_min:N_prev;
B_hat_vals = zeros(1, size(N_vals, 2));
for b=1:size(N_vals, 2)
    B_hat_vals(b) =  a * sum(t * (log((V_t(t) * c_max) / (V_t(1) * (N_vals(b) - V_t(t))))));
end

I_vals = zeros(1, size(N_vals, 2));
for b=1:size(N_vals, 2)
    I_vals(b) = min(sum(abs((N_vals(b) * init) * (1./(V_t(1) + (N_vals(b) - V_t(1)).*exp(-1 * B_hat_vals(b) * (0:119)))))));
end

hold on
plot(N_vals, I_vals, 'r-')
axis tight
title("Function N to I(N, B_N)")

%% Question 2 Initializers
clearvars
T = readtable("project5_data.xlsx");
V = table2array(T(2, 13:1103));
Y = table2array(T(3, 13:1103));
V_t = V(52:171);
Y_t = Y(52:171);
% Setting up I like the problem suggests
I = V((52:171) + 7) - V((52:171) - 7);

%% Question 2, 1
% Plot the rates of detected infections and cumulative deaths
figure
plot(I, 'r-')
title('Calculated Rates of Detected Infections')

figure
plot(Y, 'g-')
title('Cumulative Deaths')

%% Question 2.2 + 2.3
clearvars
%Euler scheme for solving the SIR problem:
a = (0.05:0.01:0.2);
r_ratio = (0.8:0.05:2.2);
size(a)
size(r_ratio)
index = 1;
%Calculating the euler approximation for all values of alpha and beta
%described in the set.
for i = 1:16
    for r = 1:29
        alpha = a(i);
        beta = r_ratio(r)*alpha;
        [I(index,1),R(index,1)] = euler(909327,6,0,.01,119,alpha,beta);
        index = index + 1;
    end
end

[I,R]

%Function that implements the euler scheme (Q 2.2)
function [I,R] = euler(S_0,I_0,R_0, h,t_max,alpha,beta)
    S = S_0;
    I = I_0;
    R = R_0;
    S_old = S;
    I_old = I;
    for t = 0:h:t_max
        S = S + ((-1)*beta*S_old*(I_old/909327))*h;
        I = I + (beta*S_old*(I_old/909327)-alpha*(I_old))*h;
        R = R + (alpha*I_old)*h;
        S_old = S;
        I_old = I;
    end
end
