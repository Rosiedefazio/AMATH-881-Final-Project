function out = CancerModel
% Cancer-Immune ODE system for MatCont continuation
%
% States      : x = [C, T, M]
% Parameters  : p = [r, b, gamma, alpha, lam, beta, delta, rho, eta, omega]
% Index:             1  2    3      4     5     6      7     8    9     10

out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobianp;
out{5} = [];   % hessians   (unused)
out{6} = [];   % hessiansp  (unused)
out{7} = [];   % 3rd deriv  (unused)
out{8} = [];
out{9} = [];

% ── ODE right-hand side ────────────────────────────────────────────────────
function dydt = fun_eval(~, x, r, b, gamma, alpha, lam, ...
                         beta, delta, rho, eta, omega)
    C = x(1);  T = x(2);  M = x(3);
    dydt = [
        r*C*(1 - b*C)*(C - gamma) - alpha*C*T;
        lam*C + beta*M*T - delta*T;
        rho - eta*M - omega*C*M
    ];
end

% ── Default init (overridden by runner) ────────────────────────────────────
function [tspan, y0, options] = init
    options = odeset();
    tspan   = [0, 10];
    y0      = [0.5; 0.5; 0.5];
end

% ── Jacobian  ∂f/∂x ────────────────────────────────────────────────────────
function J = jacobian(~, x, r, b, gamma, alpha, lam, ...
                      beta, delta, rho, eta, omega)
    C = x(1);  T = x(2);  M = x(3);
    J = [
        r*(2*C - gamma - 3*b*C^2 + 2*b*gamma*C) - alpha*T, ...
            -alpha*C,            0;
        lam,                     beta*M - delta,   beta*T;
        -omega*M,                0,               -eta - omega*C
    ];
end

% ── Parameter Jacobian  ∂f/∂p ──────────────────────────────────────────────
% Columns ordered: [r, b, gamma, alpha, lam, beta, delta, rho, eta, omega]
function Jp = jacobianp(~, x, r, b, gamma, alpha, lam, ...
                        beta, delta, rho, eta, omega)
    C = x(1);  T = x(2);  M = x(3);
    Jp = [
        C*(1-b*C)*(C-gamma), -r*C^2*(C-gamma), -r*C*(1-b*C), -C*T, 0,   0,   0,  0,  0,     0;
        0,                    0,                 0,             0,   C,  M*T, -T,  0,  0,     0;
        0,                    0,                 0,             0,   0,   0,   0,  1, -M,  -C*M
    ];
end

end % CancerModel
