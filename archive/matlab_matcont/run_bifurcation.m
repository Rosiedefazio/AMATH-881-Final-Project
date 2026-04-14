function run_bifurcation(p_default, vary_idx, x0, p1_range, p2_range, options)
% Run MatCont continuation and two-parameter bifurcation diagram.
%
% INPUTS
%   p_default  : 1×10 default parameter vector
%                [r, b, gamma, alpha, lam, beta, delta, rho, eta, omega]
%   vary_idx   : 1×2 vector of parameter indices to vary, e.g. [3, 4]
%                (index 3 = gamma, index 4 = alpha)
%   x0         : 3×1 initial state [C0; T0; M0]
%   p1_range   : [min, max] for first parameter (1D sweep)
%   p2_range   : [min, max] for second parameter (in 2D continuation)
%   options    : (optional) struct with fields:
%                  .MaxNumPoints  (default 300)
%                  .MaxStepsize   (default 0.05)
%                  .MinStepsize   (default 1e-6)
%                  .Singularities (default 1)
%
% EXAMPLE — vary alpha (idx 4) and gamma (idx 3):
%   p = [0.5, 0.1, 0.2, 0.3, 0.4, 0.1, 0.5, 0.2, 0.3, 0.1];
%   run_bifurcation(p, [4, 3], [0.5;0.1;0.3], [0.05, 1.5], [0.05, 0.8]);

% ── Default options ─────────────────────────────────────────────────────────
if nargin < 6 || isempty(options)
    options = struct();
end
MaxPts  = getfield_default(options, 'MaxNumPoints', 300);
MaxStep = getfield_default(options, 'MaxStepsize',  0.05);
MinStep = getfield_default(options, 'MinStepsize',  1e-6);
Sing    = getfield_default(options, 'Singularities', 1);

param_names = {'r','b','gamma','alpha','lam','beta','delta','rho','eta','omega'};
p1_idx = vary_idx(1);
p2_idx = vary_idx(2);
fprintf('Varying: %s (idx %d) and %s (idx %d)\n', ...
    param_names{p1_idx}, p1_idx, param_names{p2_idx}, p2_idx);

% ── Step 1: find equilibrium near x0 ────────────────────────────────────────
p_col = p_default(:);   % ensure column vector

rhs = @(x) model_rhs(x, p_col);
eq_opts = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-10);
[x_eq, ~, flag] = fsolve(rhs, x0, eq_opts);
if flag <= 0
    warning('fsolve did not converge. Using x0 as initial guess.');
    x_eq = x0;
end
fprintf('Equilibrium found: C=%.4f, T=%.4f, M=%.4f\n', x_eq(1), x_eq(2), x_eq(3));

% ── Step 2: MatCont options ──────────────────────────────────────────────────
opt = contset;
opt = contset(opt, 'MaxNumPoints',  MaxPts);
opt = contset(opt, 'MaxStepsize',   MaxStep);
opt = contset(opt, 'MinStepsize',   MinStep);
opt = contset(opt, 'Singularities', Sing);
opt = contset(opt, 'Backward',      0);

% ── Step 3: 1D equilibrium continuation along p1 ────────────────────────────
fprintf('\n--- 1D continuation along %s ---\n', param_names{p1_idx});

[x0c, v0] = init_EP_EP(@CancerModel, x_eq, p_col, p1_idx);
[x1d, v1d, s1d, ~, ~] = cont(@equilibrium, x0c, v0, opt);

plot_1d(x1d, s1d, p1_idx, p1_range, param_names);

% ── Step 4: 2D continuation of bifurcation curves ───────────────────────────
fprintf('\n--- 2D continuation in [%s, %s] ---\n', ...
    param_names{p1_idx}, param_names{p2_idx});

figure; hold on; grid on;
xlabel(param_names{p1_idx}); ylabel(param_names{p2_idx});
title(sprintf('Bifurcation diagram: %s vs %s', ...
    param_names{p1_idx}, param_names{p2_idx}));

found_any = false;
for k = 1:length(s1d)
    lbl = s1d(k).label;

    if strcmp(lbl, 'LP')   % Fold / Limit Point
        try
            [x0_lp, v0_lp] = init_LP_LP(@CancerModel, x1d, s1d(k), ...
                                         p_col, [p1_idx, p2_idx]);
            [xlp, ~, slp] = cont(@limitpoint, x0_lp, v0_lp, opt);
            plot_2d_curve(xlp, slp, p1_idx, p2_idx, 'b', 'LP');
            found_any = true;
        catch ME
            warning('LP continuation failed at point %d: %s', k, ME.message);
        end

    elseif strcmp(lbl, 'H')  % Hopf
        try
            [x0_h, v0_h] = init_H_H(@CancerModel, x1d, s1d(k), ...
                                     p_col, [p1_idx, p2_idx]);
            [xh, ~, sh] = cont(@hopf, x0_h, v0_h, opt);
            plot_2d_curve(xh, sh, p1_idx, p2_idx, 'r', 'Hopf');
            found_any = true;
        catch ME
            warning('Hopf continuation failed at point %d: %s', k, ME.message);
        end
    end
end

if ~found_any
    warning('No LP or Hopf points detected in 1D sweep. Try adjusting p1_range.');
end
legend show;
end

% ════════════════════════════════════════════════════════════════════════════
%  HELPER FUNCTIONS
% ════════════════════════════════════════════════════════════════════════════

function dydt = model_rhs(x, p)
% Thin wrapper for fsolve (sets t=0)
    dydt = feval(@CancerModel);
    dydt = dydt{2}(0, x, p(1), p(2), p(3), p(4), p(5), ...
                          p(6), p(7), p(8), p(9), p(10));
end

function plot_1d(x, s, p_idx, p_range, pnames)
    figure; hold on; grid on;
    n_states = 3;
    state_names = {'C','T','M'};
    colors = lines(n_states);
    p_vals = x(end, :);   % parameter is last row in x

    for st = 1:n_states
        plot(p_vals, x(st,:), 'Color', colors(st,:), ...
             'DisplayName', state_names{st});
    end

    % Mark bifurcation points
    for k = 1:length(s)
        lbl = s(k).label;
        if ismember(lbl, {'LP','H','BP','PD'})
            idx = s(k).index;
            plot(x(end, idx), x(1, idx), 'ko', 'MarkerSize', 8, ...
                 'DisplayName', lbl);
            text(x(end, idx), x(1, idx), sprintf('  %s', lbl));
        end
    end

    xlabel(pnames{p_idx}); ylabel('State at equilibrium');
    title(sprintf('1D Continuation along %s', pnames{p_idx}));
    xlim(p_range); legend show;
end

function plot_2d_curve(x, s, p1_idx, p2_idx, clr, lbl)
    % In 2D continuation x has states + 2 parameters
    % The last two rows are the two continuation parameters
    n = size(x, 1);
    p1_row = n - 1;
    p2_row = n;
    plot(x(p1_row,:), x(p2_row,:), 'Color', clr, ...
         'LineWidth', 1.5, 'DisplayName', lbl);
end

function val = getfield_default(s, field, default)
    if isfield(s, field)
        val = s.(field);
    else
        val = default;
    end
end
