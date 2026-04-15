%CHARACTERIZATION_UNIT_TEST Basic validation for characterization helpers.

% synthetic JV solution structure
x = linspace(0, 1e-4, 40);
t = linspace(0, 1, 60);
Vapp = linspace(0, 1.2, numel(t))';
Jn = 20 - 30*Vapp;

sol = zeros(numel(t), numel(x), 4);
for i = 1:numel(t)
    sol(i,:,1) = 1e10 * exp(-x / max(x)) + i;
    sol(i,:,2) = 2e10 * exp(-(max(x)-x) / max(x));
    sol(i,:,3) = 1e15;
    sol(i,:,4) = linspace(0, Vapp(i), numel(x));
end

p = pinParams;
p.Ana = 1;
p.OC = 0;
p.JV = 1;

JVsol = struct('sol', sol, 'x', x, 't', t, 'p', p, 'Vapp', Vapp, 'Jn', Jn);

results = characterize_perovskite_bilayer(JVsol);
assert(isstruct(results));
assert(isfield(results, 'metrics'));
assert(isfield(results.metrics, 'Voc'));
assert(numel(results.CF.f) == numel(results.CF.Cf));
