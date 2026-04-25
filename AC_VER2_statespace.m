function y = AC_VER2_statespace(p, u, t, x0, c)
% AC_3DoF_statespace  3DoF longitudinal aircraft model for SIDPAC oe.m
%
% Integrates the nonlinear longitudinal equations of motion using RK4
% and returns all four states as the model output.
%
% Usage (called internally by oe.m):
%   y = AC_3DoF_statespace(p, u, t, x0, c)
%
% Inputs:
%   p   = initial parameter guess vector [8x1]:
%   u   = input matrix [Nx2]:  col 1 = delta_T [N], col 2 = delta_e [rad]
%   t   = time vector [Nx1]
%   x0  = initial state [4x1]: [V0, alpha0, gamma0, Q0]
%   c   = constants vector [6x1]:
%
% Output:
%   y   = model output matrix [Nx4]: [V, alpha, gamma, Q]
%         each column is one state time history

% constants
g    = c(1);
rho  = c(2);
S    = c(3);
cbar = c(4);
mass = c(5);
Jy   = c(6);

% initial parameter guesses
CL0      = p(1);
CL_alpha = p(2);
CD0      = p(3);
k_CD     = p(4);
CM0      = p(5);
CM_alpha = p(6);
CM_Q     = p(7);
CM_e     = p(8);

N  = length(t);

% allocate state matrix
% rows = [V; alpha; gamma; Q], one column per time step
x = zeros(4, N);
x(:, 1) = x0(:);

% RK4 integration
for i = 1:N-1

    dt = t(i+1) - t(i);
    u_i = u(i, :);

    k1 = dynamics(x(:,i), u_i, g, rho, S, cbar, mass, Jy, ...
                  CL0, CL_alpha, CD0, k_CD, CM0, CM_alpha, CM_Q, CM_e);

    k2 = dynamics(x(:,i)+0.5*dt*k1, u_i, g, rho, S, cbar, mass, Jy, ...
                  CL0, CL_alpha, CD0, k_CD, CM0, CM_alpha, CM_Q, CM_e);

    k3 = dynamics(x(:,i)+0.5*dt*k2, u_i, g, rho, S, cbar, mass, Jy, ...
                  CL0, CL_alpha, CD0, k_CD, CM0, CM_alpha, CM_Q, CM_e);

    k4 = dynamics(x(:,i)+dt*k3, u_i, g, rho, S, cbar, mass, Jy, ...
                  CL0, CL_alpha, CD0, k_CD, CM0, CM_alpha, CM_Q, CM_e);

    x(:,i+1) = x(:,i) + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);

    % NEW SMART NAN GUARD
    if any(isnan(x(:,i+1))) || any(isinf(x(:,i+1)))
        % The optimizer guessed an unstable parameter. 
        % Fill the rest of the simulation with an artificially massive number.
        % This creates a massive cost function penalty without crashing oe.m.
        x(:, i+1:end) = 1e8; 
        break;
    end

end
y = x';   % transpose to [Nx4]: each row is one time step

end


% =========================================================================
% Local function — nonlinear equations of motion
%
%   State:  x = [V, alpha, gamma, Q]'
%   Input:  u = [delta_T, delta_e]
%
%   CL = CL0 + CL_alpha * alpha
%   CD = CD0 + k_CD * CL^2
%   CM = CM0 + CM_alpha*alpha + CM_Q*qhat + CM_e*delta_e
%        where qhat = cbar*Q / (2*V)
%
%   q_dyn = 0.5 * rho * V^2
%   L     = q_dyn * S * CL
%   D     = q_dyn * S * CD
%   M     = q_dyn * S * cbar * CM
%
%   Vdot     = (-D + T*cos(alpha) - m*g*sin(gamma)) / m
%   gammadot = ( L + T*sin(alpha) - m*g*cos(gamma)) / (m*V)
%   Qdot     = M / Jy
%   alphadot = Q - gammadot
% =========================================================================

function xdot = dynamics(x, u, g, rho, S, cbar, mass, Jy, ...
                          CL0, CL_alpha, CD0, k_CD, ...
                          CM0, CM_alpha, CM_Q, CM_e)

% =========================
% STATE EXTRACTION
% =========================
V     = x(1);
alpha = x(2);
gamma = x(3);
Q     = x(4);

delta_T = u(1);
delta_e = u(2);

% =========================
% SAFE VELOCITY
% =========================
V_safe = max(V, 0.5);

% =========================
% AERODYNAMIC COEFFS
% =========================
CL = CL0 + CL_alpha * alpha;

CD = CD0 + k_CD * CL^2;

CM = CM0 + CM_alpha * alpha + CM_Q * Q + CM_e * delta_e;

% =========================
% FORCES / MOMENTS
% =========================
q_dyn = 0.5 * rho * V^2;

L = q_dyn * S * CL;
D = q_dyn * S * CD;
M = q_dyn * S * CM;  

T = delta_T;

% =========================
% EQUATIONS OF MOTION
% =========================
Vdot = (-D + T*cos(alpha) - mass*g*sin(gamma)) / mass;

gammadot = (L + T*sin(alpha) - mass*g*cos(gamma)) / (mass * V_safe);

Qdot = M / Jy;

alphadot = Q ...
    + (-T*sin(alpha) - q_dyn*S*CL)/(mass*V_safe) ...
    + g*cos(gamma)/V_safe;

xdot = [Vdot; alphadot; gammadot; Qdot];

end