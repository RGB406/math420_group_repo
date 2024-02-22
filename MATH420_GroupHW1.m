%% AMSC 420 Group Homework 1
% Group: Robert "Eddie" Bull, Alexander Klein

%% Question 1 Initializers
T = readtable("project5_data.xlsx");
V = table2array(T(2, 13:1103));
V_t = V(52:(119 + 52));

%% Question 1, part a
N_max = 909327;
N_min = 1 + V_t(120);
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
plot(0:119, V_t(1:120), 'r-')
% Then, plot the estimated value using the N_max pop
plot(0:119, NI_max * (1./(V_t(1) + (N_max - V_t(1)).*exp(-1 * B_hat_max * (0:119)))), 'g-')
% Finally, plot with the N_min estimation
plot(0:119, NI_min * (1./(V_t(1) + (N_min - V_t(1)).*exp(-1 * B_hat_min * (0:119)))), 'b-')
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
t = 1:120;
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
plot(0:119, V_t(1:120), 'r-')
% ...and our estimate second!
plot(0:119, NI * (1./(V_t(1) + (N - V_t(1)).*exp(-1 * B_hat * (0:119)))), 'g-')
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
%2.3 part (a):
alph = (0.05:0.01:0.2);
r_ratio = (0.8:0.05:2.2);
index = 1;
%Calculating the euler approximation for all values of alpha and beta
%described in the set.
I_sim = zeros(464,120);
R_sim = zeros(464,120);
for i = 1:16
    for r = 1:29
        alpha = alph(i);
        beta = r_ratio(r)*alpha;
        for t = 1:120
            [I_sim(index,t),R_sim(index,t)] = euler(909327,6,0,.01,t,alpha,beta);
        end
        index = index + 1;
    end
end
%Each row in I_sim stores a time series of euler approximations of I(t) for 
% a different alpha beta pair. Time series ranges from t = t0 to t = t_max.
%R_sim is structured similarly. (464 rows, 464 alpha-beta pairs).


%2.3 part (b):
%Optimize 'p' for each alpha,beta pair:
syms f(p)
p_vec = zeros(464,1);
for i = 1:464
    f(p) = (I-p*I_sim(i,:))*transpose(I-p*I_sim(i,:));
    eqn = diff(f) == 0;
    p_vec(i,1) = solve(eqn,p);
end
%Each row of p_vec stores the optimizing p-value for the alpha-beta pair
%encoded in the same row of I_sim and R_sim (464 rows, 464 alpha-beta
%pairs).

%Optimize gamma for each alpha,beta pair:
syms g(gamma)
gamma_vec = zeros(464,1);
for i = 1:464
    g(gamma) = (Y_t-gamma*R_sim(i,:))*transpose(Y_t-gamma*R_sim(i,:));
    eqn = diff(g) == 0;
    gamma_vec(i,1) = solve(eqn,gamma);
end
%Each row of gamma_vec stores the optimizing gamma-value for the 
% alpha-beta pair encoded in the same row of I_sim and R_sim 
% (464 rows, 464 alpha-beta pairs).


%% Question 2.4
%Part 3.c for (c_I,c_Y) = (0,1)
index = 1;
J_1 = zeros(16,29);
for a = 1:16
    for r = 1:29
        J_1(a,r) = (Y_t-gamma_vec(index)*R_sim(index,:))*transpose(Y_t-gamma_vec(index)*R_sim(index,:));
        index = index + 1;
    end
end
figure();
surf(J_1)
title('J(alpha,beta) for (c_I,c_Y) = (0,1)')
%axes are based on J indices, not raw alpha,beta values.
disp("[alpha, beta] that minimize J given by:")
[row, col] = find(J_1 == min(J_1(:)));
[alph(row), r_ratio(col)*alph(row)]
disp("Minimum value obtained at J = " + J_1(row,col))
%J_1 is a 16x29 matrix, with the (i,j) entry corresponding to the value of 
% J (J parameters given by (c_I,c_Y) = (0,1)), for alpha = a(i) and 
% beta = b(j). Where 'a' and 'b' are vectors that store the range of alpha 
% and beta.




%Part 3.c for (c_I,c_Y) = (1,1)
index = 1;
J_2 = zeros(16,29);
for a = 1:16
    for r = 1:29
        J_2(a,r) = (I-p_vec(index)*I_sim(index,:))*transpose(I-p_vec(index)*I_sim(index,:)) + (Y_t-gamma_vec(index)*R_sim(index,:))*transpose(Y_t-gamma_vec(index)*R_sim(index,:));
        index = index + 1;
    end
end
figure();
surf(J_2)
%axes are based on J indices, not raw alpha,beta values.
title('J(alpha,beta) for (c_I,c_Y) = (1,1)')
disp("[alpha, beta] that minimize J given by:")
[row, col] = find(J_2 == min(J_2(:)));
[alph(row), r_ratio(col)*alph(row)]
disp("Minimum value obtained at J = " + J_2(row,col))
%J_2 is similar to J_1 as described above, with J being given by parameters
% (c_I,c_Y) = (1,1) instead this time.


% Adjusting the I(0) could lead to a better fit. I(0) = 6 was used for all
% euler approximations, since that was the first day were cumulative
% infections was greater than 5. It is important that I(0) > 0 for an
% epidemic to develop, and a greater value will result in a more rapid
% increase of the approximated time series I_sim(t). Choosing an I(0)
% greater than 6 could lead to a model that better predicted the rapid rate
% of infection. For instance, choosing I(0) = 26, which was 8 days later 
% than the chosen time for I(0), would have resulted in a J_min value 
% 15% smaller than we observed with parameters (c_I,c_Y) = (1,1).


%Function that implements the euler scheme (Q 2.2)
function [I,R] = euler(S_0,I_0,R_0, h,t_max,alpha,beta)
    S = S_0;
    I = I_0;
    R = R_0;
    S_old = S;
    I_old = I;
    for t = 0:h:t_max
        S = S + (-1)*beta*S_old*(I_old/909327)*h;
        I = I + (beta*S_old*(I_old/909327)-alpha*(I_old))*h;
        R = R + (alpha*I_old)*h;
        S_old = S;
        I_old = I;
    end
end
