% Assignment 2 - Parameter identification

load('output.mat')
load('input.mat')

% Now we do some identification!
% First of all lets split the data in estimation and validation sets
% (half and half)

N = length(y);  % number of data
uest = u(1:N/2);
yest = y(1:N/2);
uval = u(N/2+1:end);
yval = y(N/2+1:end);

% ARX model 32a
PHI_a = zeros(3,N/2);
PHI_a(:, 1) = [uest(1);0;0];
PHI_a(:, 2) = [uest(2);yest(1);0];

for i=3:N/2
    PHI_a(:,i) = [uest(i); yest(i-1); yest(i-2)];
end

theta_hat_a = inv(PHI_a*PHI_a')*PHI_a*yest;


% ARX model 32b
PHI_b = zeros(4,N/2);
PHI_b(:, 1) = [uest(1);0;0;0];
PHI_b(:, 2) = [uest(2);uest(1);yest(1);0];

for i=3:N/2
    PHI_b(:,i) = [uest(i);uest(i-1); yest(i-1); yest(i-2)];
end

theta_hat_b = inv(PHI_b*PHI_b')*PHI_b*yest;

% ARX model 32c
PHI_c = zeros(4,N/2);
PHI_c(:, 1) = [0;0;0;0];
PHI_c(:, 2) = [uest(1);yest(1);0;0];
PHI_c(:,3) = [uest(2);yest(2);yest(1);0];

for i=4:N/2
    PHI_c(:,i) = [uest(i-1);yest(i-1); yest(i-2); yest(i-3)];
end

theta_hat_c = inv(PHI_c*PHI_c')*PHI_c*yest;

%% --------- Predictions and Simulations ---------

y_sim_a = zeros(N/2,1);
y_sim_b = zeros(N/2,1);
y_sim_c = zeros(N/2,1);

% ARX model 32c
y_sim_a(1) = [uval(1) 0 0] * theta_hat_a;
y_sim_a(2) = [uval(2) uval(1) 0] * theta_hat_a;

y_pred_a(1) = [uval(1) 0 0] * theta_hat_a;
y_pred_a(2) = [uval(2) uval(1) 0] * theta_hat_a;

for i=3:N/2
    y_sim_a(i) = [uval(i) y_sim_a(i-1) y_sim_a(i-2)] * theta_hat_a;
    y_pred_a(i) = [uval(i) yval(i-1) yval(i-2)]  * theta_hat_a;
end

% ARX model 32b
y_sim_b(1) = [uval(1) 0 0 0] * theta_hat_b;
y_sim_b(2) = [uval(2) uval(1) y_sim_b(1) 0] * theta_hat_b;

y_pred_b(1) = [uval(1) 0 0 0] * theta_hat_b;
y_pred_b(2) = [uval(2) uval(1) yval(1) 0] * theta_hat_b;

for i=3:N/2
    y_sim_b(i) = [uval(i) uval(i-1) y_sim_b(i-1) y_sim_b(i-2)] * theta_hat_b;
    y_sim_b(i) = [uval(i) uval(i-1) yval(i-1) yval(i-2)] * theta_hat_b;
end

% ARX model 32c
y_sim_c(1) = [0 0 0 0] * theta_hat_c;
y_sim_c(2) = [uval(1) y_sim_c(1) 0 0] * theta_hat_c;
y_sim_c(3) = [uval(2) y_sim_c(2) y_sim_c(1) 0] * theta_hat_c;

y_pred_c(1) = [0 0 0 0] * theta_hat_c;
y_pred_c(2) = [uval(1) yval(1) 0 0] * theta_hat_c;
y_pred_c(3) = [uval(2) yval(2) yval(1) 0] * theta_hat_c;

for i=4:N/2
    y_sim_c(i) = [uval(i-1) y_sim_c(i-1) y_sim_c(i-2) y_sim_c(i-3)] * theta_hat_c;
    y_sim_c(i) = [uval(i-1) yval(i-1) yval(i-2) yval(i-3)] * theta_hat_c;  
end

% Prediction errors
pred_error_a = y_pred_a - yval;
pred_RMSE_a = rms(pred_error_a);

pred_error_b = y_pred_b - yval;
pred_RMSE_b = rms(pred_error_b);

pred_error_c = y_pred_c - yval;
pred_RMSE_c = rms(pred_error_c);

% Simulation errors
sim_error_a = y_sim_a - yval;
sim_RMSE_a = rms(sim_error_a);

sim_error_b = y_sim_b - yval;
sim_RMSE_b = rms(sim_error_b);

sim_error_c = y_sim_c - yval;
sim_RMSE_c = rms(sim_error_c);

cov_a = 0.01*inv(PHI_a*PHI_a');
cov_b = 0.01*inv(PHI_b*PHI_b');
cov_c = 0.01*inv(PHI_c*PHI_c');
