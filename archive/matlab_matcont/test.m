function run_bifurcation(p_default, vary_idx, x0, p1_range, p2_range, options)
%RUN_BIFURCATION  Run MatCont continuation and produce a two‑parameter
%                 bifurcation diagram for the Cancer model.
%
% INPUTS
%   p_default : 1×10 default parameter vector
%               [r, b, gamma, alpha, lam, beta, delta, rho, eta, omega]
%   vary_idx  : 1×2 vector of parameter indices to vary, e.g. [3,4]
%               (index 3 = gamma, index 4 = alpha)
%   x0        : 3×1 initial state [C0; T0; M0]
%   p1_range  : [min, max] for the first continuation parameter
%   p2_range  : [min, max] for the second continuation parameter
%   options   : (optional) struct with fields
%                .MaxNumPoints   (default 300)
%                .MaxStepsize    (default 0.05)
%                .MinStepsize    (default 1e‑6)
%                .Singularities  (default 1)
%
% EXAMPLE
%   p = [0.5 0.1 0.2 0.3 0.4 0.1 0.5 0.2 0.3 0.1];
%   run_bifurcation(p, [4 3], [0.5;0.1;0.3], [0.05 1.5], [0.05 0.8]);

% -------------------------------------------------------------------------
% 1️⃣  Set default options
% -------------------------------------------------------------------------
if nargin < 6 || isempty(options)
    options = struct();
end
MaxPts  = getfield_default(options, 'MaxNumPoints',   300);
MaxStep = getfield_default(options, 'MaxStepsize',    0.05);
MinStep = getfield_default(options, 'MinStepsize',    1e-6);
Sing    = getfield_default(options, 'Singularities',  1);

param_names = {'r','b','gamma','alpha','lam','beta','delta','rho','eta','omega'};
p1_idx = vary_idx(1);
p2_idx = vary_idx(2);
fprintf('Varying: %s (idx %d) and %s (idx %d)\n', ...
    param_names{p1_idx}, p1_idx, param_names{p2_idx}, p2_idx);

% -------------------------------------------------------------------------
% 2️⃣  Find an equilibrium close to the supplied initial state
% -------------------------------------------------------------------------
p_col = p_default(:);               % ensure column vector
rhs   = @(x) model_rhs(x, p_col);   % wrapper for fsolve
eq_opt = optimoptions('fsolve','Display','off','TolFun',1e-10);
[x_eq,~,flag] = fsolve(rhs, x0, eq_opt);
if flag <= 0
    warning('fsolve did not converge – using supplied x0 as equilibrium.');
    x_eq = x0;
end
fprintf('Equilibrium found: C=%.4f, T=%.4f, M=%.4f\n', x_eq(1), x_eq(2), x_eq(3));

% -------------------------------------------------------------------------
% 3️⃣  MatCont continuation options
% -------------------------------------------------------------------------
opt = contset;
opt = contset(opt,'MaxNumPoints',   MaxPts);
opt = contset(opt,'MaxStepsize',    MaxStep);
opt = contset(opt,'MinStepsize',    MinStep);
opt = contset(opt,'Singularities',  Sing);
opt = contset(opt,'Backward',       0);

% -------------------------------------------------------------------------
% 4️⃣  1‑D equilibrium continuation along the first parameter (p1)
% -------------------------------------------------------------------------
fprintf('\n--- 1‑D continuation along %s ---\n', param_names{p1_idx});

[x0c, v0] = init_EP_EP(@CancerModel, x_eq, p_col, p1_idx);
[x1d, v1d, s1d, ~, ~] = cont(@equilibrium, x0c, v0, opt);

plot_1d(x1d, s1d, p1_idx, p1_range, param_names);

% -------------------------------------------------------------------------
% 5️⃣  2‑D continuation of codimension‑2 bifurcation curves
% -------------------------------------------------------------------------
fprintf('\n--- 2‑D continuation in [%s, %s] ---\n', ...
    param_names{p1_idx}, param_names{p2_idx});

figure; hold on; grid on;
xlabel(param_names{p1_idx}); ylabel(param_names{p2_idx});
title(sprintf('Two‑parameter bifurcation diagram: %s vs %s', ...
    param_names{p1_idx}, param_names{p2_idx}));

found_any = false;   % flag to know whether any codim‑2 curve was traced

for k = 1:numel(s1d)
    lbl = s1d(k).label;   % e.g. 'LP', 'H', 'GH', 'BT', 'SN', ...

    % --------------------------------------------------------------
    % 5.1  Fold / Saddle‑node (LP or SN)
    % --------------------------------------------------------------
    if strcmp(lbl,'LP') || strcmp(lbl,'SN')
        try
            [x0_lp, v0_lp] = init_LP_LP(@CancerModel, x1d, s1d(k), ...
                                         p_col, [p1_idx, p2_idx]);
            [xlp, ~, slp] = cont(@limitpoint, x0_lp, v0_lp, opt);
            plot_2d_curve(xlp, slp, p1_idx, p2_idx, 'b', 'LP/SN');
            % mark the start point of the curve
            mark_bif_point(xlp(end-1,1), xlp(end,1), 'LP/SN', 'b');
            found_any = true;
        catch ME
            warning('LP/SN continuation failed at point %d: %s',k,ME.message);
        end

    % --------------------------------------------------------------
    % 5.2  Hopf (H)
    % --------------------------------------------------------------
    elseif strcmp(lbl,'H')
        try
            [x0_h, v0_h] = init_H_H(@CancerModel, x1d, s1d(k), ...
                                      p_col, [p1_idx, p2_idx]);
            [xh, ~, sh] = cont(@hopf, x0_h, v0_h, opt);
            plot_2d_curve(xh, sh, p1_idx, p2_idx, 'r', 'Hopf');
            mark_bif_point(xh(end-1,1), xh(end,1), 'Hopf', 'r');
            found_any = true;
        catch ME
            warning('Hopf continuation failed at point %d: %s',k,ME.message);
        end

    % --------------------------------------------------------------
    % 5.3  Generalized Hopf (GH) – Bautin point
    % --------------------------------------------------------------
    elseif strcmp(lbl,'GH')
        try
            [x0_gh, v0_gh] = init_GH_GH(@CancerModel, x1d, s1d(k), ...
                                         p_col, [p1_idx, p2_idx]);
            [xgh, ~, sgh] = cont(@generalizedHopf, x0_gh, v0_gh, opt);
            plot_2d_curve(xgh, sgh, p1_idx, p2_idx, 'm', 'Gen‑Hopf');
            mark_bif_point(xgh(end-1,1), xgh(end,1), 'Gen‑Hopf', 'm');
            found_any = true;
        catch ME
            warning('Generalized Hopf continuation failed at point %d: %s',k,ME.message);
        end

    % --------------------------------------------------------------
    % 5.4  Bogdanov–Takens (BT)
    % --------------------------------------------------------------
    elseif strcmp(lbl,'BT')
        try
            [x0_bt, v0_bt] = init_BT_BT(@CancerModel, x1d, s1d(k), ...
                                         p_col, [p1_idx, p2_idx]);
            [xbt, ~, sbt] = cont(@bogtakens, x0_bt, v0_bt, opt);
            plot_2d_curve(xbt, sbt, p1_idx, p2_idx, 'c', 'Bog‑Tak');
            mark_bif_point(xbt(end-1,1), xbt(end,1), 'Bog‑Tak', 'c');
            found_any = true;
        catch ME
            warning('Bogdanov–Takens continuation failed at point %d: %s',k,ME.message);
        end
    end
end

if ~found_any
    warning('No codimension‑2 points (LP/SN, Hopf, GH, BT) detected in the 1‑D sweep.');
end

legend('show','Location','best');
end

% -------------------------------------------------------------------------
%  HELPER FUNCTIONS
% -------------------------------------------------------------------------

function dydt = model_rhs(x, p)
% Wrapper that evaluates the RHS of the Cancer model at time t = 0.
% The CancerModel function is expected to return a cell array where the
% second element is a function handle @(t,x,...) that computes the ODE RHS.
    tmp = feval(@CancerModel);
    rhs = tmp{2};
    dydt = rhs(0, x, p(1), p(2), p(3), p(4), p(5), ...
                     p(6), p(7), p(8), p(9), p(10));
end

function plot_1d(x, s, p_idx, p_range, pnames)
% Plot the result of a 1‑D equilibrium continuation.
%   x      – matrix of continuation data (states + parameter row)
%   s      – structure array of detected singularities
%   p_idx  – index of the continuation parameter
%   p_range – limits for the x‑axis
%   pnames – cell array with parameter names
    figure; hold on; grid on;
    n_states = 3;
    state_names = {'C','T','M'};
    colors = lines(n_states);
    p_vals = x(end, :);               % last row holds the continuation parameter

    for st = 1:n_states
        plot(p_vals, x(st, :), 'Color', colors(st, :), ...
             'DisplayName', state_names{st});
    end

    % Mark detected bifurcation points on the 1‑D plot
    for k = 1:numel(s)
        lbl = s(k).label;
        if ismember(lbl, {'LP','SN','H','GH','BT'})
            idx = s(k).index;
            plot(p_vals(idx), x(1, idx), 'ko', 'MarkerSize', 8, ...
                 'DisplayName', lbl);
            text(p_vals(idx), x(1, idx), sprintf(' %s', lbl));
        end
    end

    xlabel(pnames{p_idx}); ylabel('Equilibrium state');
    title(sprintf('1‑D continuation along %s', pnames{p_idx}));
    xlim(p_range);
    legend('show','Location','best');
end

function plot_2d_curve(x, s, p1_idx, p2_idx, clr, lbl)
% Plot a 2‑D continuation curve (codimension‑2 bifurcation).
%   x      – matrix of continuation data (states + 2 parameter rows)
%   s      – singularity structure (unused here but kept for consistency)
%   p1_idx – index of the first continuation parameter
%   p2_idx – index of the second continuation parameter
%   clr    – colour for the curve
%   lbl    – label for the legend
    n = size(x, 1);
    p1_row = n - 1;   % second‑last row holds the first parameter
    p2_row = n;       % last row holds the second parameter
    plot(x(p1_row, :), x(p2_row, :), 'Color', clr, ...
         'LineWidth', 1.5, 'DisplayName', lbl);
end

function mark_bif_point(p1, p2, lbl, clr)
% Add a coloured marker and a short text label for a codimension‑2 point.
    plot(p1, p2, 'o', 'MarkerEdgeColor', clr, 'MarkerFaceColor', clr, ...
         'MarkerSize', 8, 'HandleVisibility','off');
    text(p1, p2, sprintf(' %s', lbl), 'Color', clr, 'FontSize', 9, ...
         'VerticalAlignment','bottom','HorizontalAlignment','left');
end

function val = getfield_default(s, field, default)
% Return s.(field) if it exists, otherwise return DEFAULT.
    if isfield(s, field)
        val = s.(field);
    else
        val = default;
    end
end