%% main_bifurcation.m
% Driver script for cancer-immune bifurcation analysis using MatCont.
%
% Requires:
%   - CancerModel.m       (system definition)
%   - run_bifurcation.m   (continuation runner)
%   - MatCont on MATLAB path





clear; clc; close all;

%% ── MatCont path (edit this to your installation) ────────────────────────
addpath(genpath('/Users/rosiedefazio/Desktop/AMATH881/Final Project/Code/matlab_matcont'));   % <-- change to your MatCont folder
opts.MaxNumPoints  = 300;
opts.MaxStepsize   = 0.05;
opts.MinStepsize   = 1e-6;
opts.Singularities = 1;

p = [0.5, 0.1, 0.2, 0.3, 0.4, 0.1, 0.5, 0.2, 0.3, 0.1];
run_bifurcation(p, [4, 3], [0.5;0.1;0.3], [0.05, 1.5], [0.05, 0.8], opts);
%% ── Parameter vector ─────────────────────────────────────────────────────
%    [r,   b,    gamma, alpha, lam,  beta, delta, rho,  eta,  omega]
%p = [0.08, 0.01,  0.20,  0.10,  0.20, 0.01590026216975, 1.0,  0.0, 0.0, 1.0];

%% ── Initial state guess ──────────────────────────────────────────────────
%x0 = [5; 5; 5];   % [C0; T0; M0]

%% ── Continuation options (optional, remove to use defaults) ──────────────
opts.MaxNumPoints  = 300;
opts.MaxStepsize   = 0.05;
opts.MinStepsize   = 1e-6;
opts.Singularities = 1;

%% ── Run 1: vary alpha (idx 4) and gamma (idx 3) ──────────────────────────
%fprintf('=== Run 1: alpha vs gamma ===\n');
%run_bifurcation(p, [4, 3], x0, [0.01, 0.5], [0.0, 30], opts);

%% ── Run 2: vary rho (idx 8) and delta (idx 7) ───────────────────────────
%fprintf('=== Run 2: rho vs delta ===\n');
%run_bifurcation(p, [8, 7], x0, [0.01, 1.0], [0.10, 1.5], opts);

%% ── Add more runs below as needed ────────────────────────────────────────
% run_bifurcation(p, [p1_idx, p2_idx], x0, [p1_min, p1_max], [p2_min, p2_max]);
