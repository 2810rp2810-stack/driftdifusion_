%EXAMPLE_CHARACTERIZATION End-to-end perovskite bilayer characterization example.

clear; close all;

% Load parameters and equilibrate
par = pinParams;
[~, ~, ~, sol_i_eq_SR] = equilibrate_minimal(par);

% Run illuminated J-V sweep
JVsol = doJV(sol_i_eq_SR, par, 1e-2, 100, 1, 0, 1.5);

% Optional dark sweep for overlay in Figure 1
JVsol.dark = doJV(sol_i_eq_SR, par, 1e-2, 100, 0, 0, 1.5);

% Comprehensive characterization
results = characterize_perovskite_bilayer(JVsol);

disp(results.metrics);
