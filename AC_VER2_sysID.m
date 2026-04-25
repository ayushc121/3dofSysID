function [y_out, p_out, y_true, p_true, u_cmd, u_act, gust] = AC_VER2_sysID(filename1, filename2)

% 3 DoF LONGITUDINAL SYSTEM IDENTIFICATION

%% importing data

data = readmatrix(filename1);
paramdata = readtable(filename2);

% state variable time series
time     = data(:, 1);
V        = data(:, 2);
alpha    = data(:, 3);
gamma    = data(:, 4);
Q        = data(:, 5);

y_true = [V, alpha, gamma, Q, time];

u_cmd = data(:, 6:7);
u_act = data(:, 8:9);

u_thrust_cmd = u_cmd(:, 1);
u_elev_cmd = u_cmd(:, 2);

gust = data(:, 10);

% timestep
dt = time(2) - time(1);

N = length(time);

% constants
g    = paramdata.value(strcmp(paramdata.parameter, 'g'));
rho  = paramdata.value(strcmp(paramdata.parameter, 'rho'));
S    = paramdata.value(strcmp(paramdata.parameter, 'S'));
cbar = paramdata.value(strcmp(paramdata.parameter, 'cbar'));
m    = paramdata.value(strcmp(paramdata.parameter, 'm'));
Jy   = paramdata.value(strcmp(paramdata.parameter, 'Jy'));

% true parameters
gt_CL0      = paramdata.value(strcmp(paramdata.parameter, 'C_L0'));
gt_CL_alpha = paramdata.value(strcmp(paramdata.parameter, 'C_La'));
gt_CD0      = paramdata.value(strcmp(paramdata.parameter, 'C_D0'));
gt_k_CD     = paramdata.value(strcmp(paramdata.parameter, 'k_CLCD'));
gt_CM0      = paramdata.value(strcmp(paramdata.parameter, 'C_M0'));
gt_CM_alpha = paramdata.value(strcmp(paramdata.parameter, 'C_Ma'));
gt_CM_Q     = paramdata.value(strcmp(paramdata.parameter, 'C_MQ'));
gt_CM_e     = paramdata.value(strcmp(paramdata.parameter, 'C_Me'));

p_true = [gt_CL0, gt_CL_alpha, gt_CD0, gt_k_CD, ...
          gt_CM0, gt_CM_alpha, gt_CM_Q, gt_CM_e];

%% data smoothing

% not needed for this simulated data

%% differentiation
Vdot     = deriv(V,     dt);
gammadot = deriv(gamma, dt);
Qdot     = deriv(Q,     dt);

%% aerodynamic terms
% dynamic pressure [N x 1]
% NOTE: S is applied here so all Phi matrices are consistent
q_dyn_times_S = 0.5 * rho * V.^2 * S;


%% linear regression

% =========================================================================
% PITCHING MOMENT, identifies [CM0, CM_alpha, CM_Q, CM_e]
% =========================================================================

z_M   = Jy * Qdot;
Phi_M = (q_dyn_times_S) .* [ones(N,1), alpha, Q, u_elev_cmd];

% Define bounds for [CM0, CM_alpha, CM_Q, CM_e]
% CM_alpha should be negative (static stability)
% CM_Q should be negative (pitch damping)
% CM_e should be strictly negative (positive elevator = pitch down)

lb_M = [-Inf; -Inf; -Inf; -Inf]; % Lower bounds
ub_M = [ Inf;   -0.01;   -0.01; -0.01]; % Upper bounds forces negatives

options = optimoptions('lsqlin', 'Display', 'off');
theta_M_hat = lsqlin(Phi_M, z_M, [], [], [], [], lb_M, ub_M, [], options);

CM0_hat      = theta_M_hat(1);
CM_alpha_hat = theta_M_hat(2);
CM_Q_hat     = theta_M_hat(3);
CM_e_hat     = theta_M_hat(4);


% =========================================================================
% LIFT, identifies [CL0, CL_alpha]
%   Note: q_dyn already includes S (defined above)
% =========================================================================

z_L   = m * V .* gammadot - u_thrust_cmd .* sin(alpha) + m * g * cos(gamma);

Phi_L = q_dyn_times_S .* [ones(N,1), alpha];

[~, theta_L_hat, Cbb_L, S2_L] = lesq(Phi_L, z_L);

CL0_hat      = theta_L_hat(1);
CL_alpha_hat = theta_L_hat(2);


% =========================================================================
% DRAG, identifies [CD0, k_CD]
% =========================================================================

z_D   = -(m * Vdot + m * g * sin(gamma) - u_thrust_cmd .* cos(alpha));

% Using the CL_hat from your Lift regression
CL_hat = CL0_hat + CL_alpha_hat .* alpha;
Phi_D = (q_dyn_times_S) .* [ones(N,1), CL_hat.^2];

% Define Lower Bounds (lb) and Upper Bounds (ub) for [CD0, k_CD]
% We force CD0 to be at least 0.005, and k_CD to be at least 0.01
lb_D = [0.005; 0.01]; 
ub_D = [Inf; Inf];    % No upper limit

% Solve using Constrained Least Squares
% Syntax: lsqlin(C, d, A, b, Aeq, beq, lb, ub)
options = optimoptions('lsqlin', 'Display', 'off');
theta_D_hat = lsqlin(Phi_D, z_D, [], [], [], [], lb_D, ub_D, [], options);

CD0_hat  = theta_D_hat(1);
k_CD_hat = theta_D_hat(2);


%% output error optimization

p0 = [CL0_hat; CL_alpha_hat; CD0_hat; k_CD_hat; CM0_hat; CM_alpha_hat; CM_Q_hat; CM_e_hat];
% parameter initial guesses

c  = [g; rho; S; cbar; m; Jy];  % constants

x0 = [V(1); alpha(1); gamma(1); Q(1)];  % initial states

u_oe = [u_thrust_cmd, u_elev_cmd];   % [Nx2] - inputs

z_oe = [V, alpha, gamma, Q];  % [Nx4] - all four measured states

[y_out, p_out, ~, ~, ~] = oe('AC_VER2_statespace', p0, u_oe, time, x0, c, z_oe, 1);

% reading estimates
CL0_oe      = p_out(1);
CL_alpha_oe = p_out(2);
CD0_oe      = p_out(3);
k_CD_oe     = p_out(4);
CM0_oe      = p_out(5);
CM_alpha_oe = p_out(6);
CM_Q_oe     = p_out(7);
CM_e_oe     = p_out(8);

%% results display

% LINEAR REGRESSION
err_CL0      = 100 * (CL0_hat      - gt_CL0)      / gt_CL0;
err_CL_alpha = 100 * (CL_alpha_hat - gt_CL_alpha)  / gt_CL_alpha;
err_CD0      = 100 * (CD0_hat      - gt_CD0)       / gt_CD0;
err_k_CD     = 100 * (k_CD_hat     - gt_k_CD)      / gt_k_CD;
err_CM0      = 100 * (CM0_hat      - gt_CM0)       / gt_CM0;
err_CM_alpha = 100 * (CM_alpha_hat - gt_CM_alpha)   / gt_CM_alpha;
err_CM_Q     = 100 * (CM_Q_hat     - gt_CM_Q)      / gt_CM_Q;
err_CM_e     = 100 * (CM_e_hat     - gt_CM_e)      / gt_CM_e;

fprintf('\n========== INITIAL PARAMETER GUESSES==========\n');
fprintf('  %-12s = %10.5f  (err = %8.3f%%)\n', 'CL0',      CL0_hat,      err_CL0);
fprintf('  %-12s = %10.5f  (err = %8.3f%%)\n', 'CL_alpha', CL_alpha_hat, err_CL_alpha);
fprintf('  %-12s = %10.5f  (err = %8.3f%%)\n', 'CD0',      CD0_hat,      err_CD0);
fprintf('  %-12s = %10.5f  (err = %8.3f%%)\n', 'k_CD',     k_CD_hat,     err_k_CD);
fprintf('  %-12s = %10.5f  (err = %8.3f%%)\n', 'CM0',      CM0_hat,      err_CM0);
fprintf('  %-12s = %10.5f  (err = %8.3f%%)\n', 'CM_alpha', CM_alpha_hat, err_CM_alpha);
fprintf('  %-12s = %10.5f  (err = %8.3f%%)\n', 'CM_Q',     CM_Q_hat,     err_CM_Q);
fprintf('  %-12s = %10.5f  (err = %8.3f%%)\n', 'CM_e',     CM_e_hat,     err_CM_e);

% OUTPUT ERROR
err_CL0_oe      = 100 * (CL0_oe      - gt_CL0)      / gt_CL0;
err_CL_alpha_oe = 100 * (CL_alpha_oe - gt_CL_alpha)  / gt_CL_alpha;
err_CD0_oe      = 100 * (CD0_oe      - gt_CD0)       / gt_CD0;
err_k_CD_oe     = 100 * (k_CD_oe     - gt_k_CD)      / gt_k_CD;
err_CM0_oe      = 100 * (CM0_oe      - gt_CM0)       / gt_CM0;
err_CM_alpha_oe = 100 * (CM_alpha_oe - gt_CM_alpha)   / gt_CM_alpha;
err_CM_Q_oe     = 100 * (CM_Q_oe     - gt_CM_Q)      / gt_CM_Q;
err_CM_e_oe     = 100 * (CM_e_oe     - gt_CM_e)      / gt_CM_e;

fprintf('\n========== OUTPUT ERROR OPTIMIZED PARAMETERS ==========\n');
fprintf('  %-12s = %10.5f  (err = %8.3f%%)\n', 'CL0',      CL0_oe,      err_CL0_oe);
fprintf('  %-12s = %10.5f  (err = %8.3f%%)\n', 'CL_alpha', CL_alpha_oe, err_CL_alpha_oe);
fprintf('  %-12s = %10.5f  (err = %8.3f%%)\n', 'CD0',      CD0_oe,      err_CD0_oe);
fprintf('  %-12s = %10.5f  (err = %8.3f%%)\n', 'k_CD',     k_CD_oe,     err_k_CD_oe);
fprintf('  %-12s = %10.5f  (err = %8.3f%%)\n', 'CM0',      CM0_oe,      err_CM0_oe);
fprintf('  %-12s = %10.5f  (err = %8.3f%%)\n', 'CM_alpha', CM_alpha_oe, err_CM_alpha_oe);
fprintf('  %-12s = %10.5f  (err = %8.3f%%)\n', 'CM_Q',     CM_Q_oe,     err_CM_Q_oe);
fprintf('  %-12s = %10.5f  (err = %8.3f%%)\n', 'CM_e',     CM_e_oe,     err_CM_e_oe);