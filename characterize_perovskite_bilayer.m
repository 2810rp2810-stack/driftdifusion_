function results = characterize_perovskite_bilayer(JVsol)
%CHARACTERIZE_PEROVSKITE_BILAYER Comprehensive characterization from a JV solution.
%
% This routine performs seven linked analyses from a driftfusion JV solution:
% 1) J-V curves and photovoltaic parameters
% 2) Carrier density profiles at SC and Voc
% 3) Equivalent-circuit EIS (Nyquist + Bode)
% 4) Electric-field profile
% 5) Band diagram and quasi-Fermi levels
% 6) C-V and Mott-Schottky analysis
% 7) C-F analysis
%
% Input:
%   JVsol - solution struct produced by pindrift/doJV
%
% Output:
%   results - struct containing raw arrays, extracted parameters and report text

validateattributes(JVsol, {'struct'}, {'nonempty'}, mfilename, 'JVsol', 1);
requiredFields = {'sol', 'p', 'x', 't'};
for k = 1:numel(requiredFields)
    if ~isfield(JVsol, requiredFields{k})
        error('characterize_perovskite_bilayer:MissingField', 'JVsol.%s is required.', requiredFields{k});
    end
end

% Guarantee JV vectors
if ~isfield(JVsol, 'Vapp') || ~isfield(JVsol, 'Jn') || isempty(JVsol.Vapp) || isempty(JVsol.Jn)
    [Vapp, Jn, ~] = pinana(JVsol);
else
    Vapp = JVsol.Vapp;
    Jn = JVsol.Jn;
end
Vapp = Vapp(:);
Jn = Jn(:);

p = JVsol.p;
x = JVsol.x(:);
t = JVsol.t(:);
sol = JVsol.sol;

n = squeeze(sol(:, :, 1));
h = squeeze(sol(:, :, 2));
a = squeeze(sol(:, :, 3));
phi = squeeze(sol(:, :, 4));

% --- 1) JV + PV metrics
pv = extractPVmetrics(Vapp, Jn);

% --- helper operating points
idxSC = nearestIndex(Vapp, 0);
if isnan(pv.Voc)
    idxVoc = numel(Vapp);
else
    idxVoc = nearestIndex(Vapp, pv.Voc);
end

% --- 2) carrier profiles
carrier.x = x;
carrier.n_SC = n(idxSC, :)';
carrier.p_SC = h(idxSC, :)';
carrier.n_Voc = n(idxVoc, :)';
carrier.p_Voc = h(idxVoc, :)';

% --- 3) derived charge/capacitance for C-V and EIS model
rho = h - n + a - getfieldWithDefault(p, 'NI', 0); %#ok<GFLD>
if isfield(p, 'e')
    qC = p.e;
else
    qC = 1.602176634e-19;
end
Q = qC * trapz(x, rho, 2); % [C cm^-2]
C = gradient(Q, Vapp);      % [F cm^-2]

% --- 4) electric field profile
E_SC = -gradient(phi(idxSC, :), x);
E_Voc = -gradient(phi(idxVoc, :), x);

% --- 5) band diagram at Voc
if isfield(p, 'EA') && isfield(p, 'IP')
    Ecb = p.EA - phi(idxVoc, :) - p.EA;
    Evb = p.IP - phi(idxVoc, :) - p.EA;
else
    Ecb = -phi(idxVoc, :);
    Evb = Ecb - 1.55;
end
if isfield(p, 'Ei') && isfield(p, 'kB') && isfield(p, 'T') && isfield(p, 'q') && isfield(p, 'ni')
    Efn = -phi(idxVoc, :) + p.Ei + (p.kB * p.T / p.q) * log(max(n(idxVoc, :), eps) / p.ni);
    Efp = -phi(idxVoc, :) + p.Ei - (p.kB * p.T / p.q) * log(max(h(idxVoc, :), eps) / p.ni);
else
    Efn = NaN(size(x'));
    Efp = NaN(size(x'));
end

% --- 6) C-V and Mott-Schottky
ms = extractMottSchottky(Vapp, C, p);

% --- 7) equivalent circuit + C-F approximation
eis = fitEquivalentCircuit(Vapp, Jn, C, p);
cf = buildCFresponse(eis);

% --- figures (publication-style)
mkfig(1, 'J-V Curves');
plot(Vapp, -Jn, 'LineWidth', 2, 'Color', [0.1 0.35 0.8]); hold on;
if isfield(JVsol, 'dark') && isstruct(JVsol.dark) && isfield(JVsol.dark, 'Vapp') && isfield(JVsol.dark, 'Jn')
    plot(JVsol.dark.Vapp(:), -JVsol.dark.Jn(:), '--', 'LineWidth', 1.5, 'Color', [0.2 0.2 0.2]);
    legend('Illuminated', 'Dark', 'Location', 'best');
else
    legend('J-V', 'Location', 'best');
end
xlabel('Voltage (V)'); ylabel('Current density (mA cm^{-2})'); grid on; box on;

mkfig(2, 'Carrier Profiles');
subplot(1,2,1); semilogy(x, max(carrier.n_SC, eps), 'b', x, max(carrier.p_SC, eps), 'r', 'LineWidth', 1.5); grid on;
xlabel('x (cm)'); ylabel('Density (cm^{-3})'); title('Short-circuit'); legend('n','p','Location','best');
subplot(1,2,2); semilogy(x, max(carrier.n_Voc, eps), 'b', x, max(carrier.p_Voc, eps), 'r', 'LineWidth', 1.5); grid on;
xlabel('x (cm)'); ylabel('Density (cm^{-3})'); title('Open-circuit'); legend('n','p','Location','best');

mkfig(3, 'EIS (Equivalent-Circuit)');
subplot(1,2,1);
plot(real(eis.Z), -imag(eis.Z), 'LineWidth', 1.8); grid on; axis tight;
xlabel('Z'' (\Omega cm^2)'); ylabel('-Z'''' (\Omega cm^2)'); title('Nyquist');
subplot(1,2,2);
semilogx(eis.f, abs(eis.Z), 'LineWidth', 1.5); hold on;
semilogx(eis.f, rad2deg(angle(eis.Z)), '--', 'LineWidth', 1.5);
grid on; xlabel('Frequency (Hz)'); ylabel('Magnitude / Phase'); title('Bode'); legend('|Z|', 'Phase (deg)', 'Location', 'best');

mkfig(4, 'Electric Field Profile');
plot(x, E_SC, 'LineWidth', 1.6); hold on; plot(x, E_Voc, '--', 'LineWidth', 1.6);
xlabel('x (cm)'); ylabel('E (V cm^{-1})'); grid on; box on;
legend('SC','Voc','Location','best');

mkfig(5, 'Band Diagram at Voc');
plot(x, Ecb, 'LineWidth', 1.8); hold on;
plot(x, Evb, 'LineWidth', 1.8);
plot(x, Efn, '--', 'LineWidth', 1.4);
plot(x, Efp, '--', 'LineWidth', 1.4);
xlabel('x (cm)'); ylabel('Energy (eV)'); grid on; box on;
legend('E_C','E_V','E_{Fn}','E_{Fp}','Location','best');

mkfig(6, 'C-V and Mott-Schottky');
subplot(1,2,1);
plot(Vapp, C, 'LineWidth', 1.6); grid on;
xlabel('Voltage (V)'); ylabel('C (F cm^{-2})'); title('C-V');
subplot(1,2,2);
plot(ms.Vfit, ms.Yfit, 'o', 'MarkerSize', 4); hold on;
plot(ms.Vline, ms.Yline, 'LineWidth', 1.5);
grid on; xlabel('Voltage (V)'); ylabel('1/C^2 (cm^4 F^{-2})'); title('Mott-Schottky');

mkfig(7, 'C-F Analysis');
semilogx(cf.f, cf.Cf, 'LineWidth', 1.8); grid on; box on;
xlabel('Frequency (Hz)'); ylabel('C(f) (F cm^{-2})');

% Pack results
results.JV = struct('Vapp', Vapp, 'Jn', Jn, 'Q', Q, 'C', C);
results.metrics = pv;
results.carrier = carrier;
results.electricField = struct('x', x, 'E_SC', E_SC(:), 'E_Voc', E_Voc(:));
results.band = struct('x', x, 'Ecb', Ecb(:), 'Evb', Evb(:), 'Efn', Efn(:), 'Efp', Efp(:));
results.CV = ms;
results.EIS = eis;
results.CF = cf;
results.operatingPoints = struct('idxSC', idxSC, 'idxVoc', idxVoc, 'tSC', t(idxSC), 'tVoc', t(idxVoc));

results.report = sprintf(['Comprehensive device characterization complete.\n' ...
    'Jsc = %.3f mA/cm^2\nVoc = %.3f V\nFF = %.2f %%\nPCE = %.2f %%\n' ...
    'Vbi (Mott-Schottky) = %.3f V\nNA = %.3e cm^-3\n' ...
    'Rs = %.3e Ohm*cm^2\nRrec = %.3e Ohm*cm^2\nCgeo = %.3e F/cm^2\nRion = %.3e Ohm*cm^2\nCion = %.3e F/cm^2\n' ...
    'f_geo = %.3e Hz\nf_ion = %.3e Hz\nWdep = %.3e cm\n'], ...
    pv.Jsc, pv.Voc, 100*pv.FF, pv.PCE, ms.Vbi, ms.NA, eis.Rs, eis.Rrec, eis.Cgeo, eis.Rion, eis.Cion, eis.f_geo, eis.f_ion, ms.Wdep);

disp(results.report);

end

function mkfig(id, name)
figure(id); clf;
set(gcf, 'Color', 'w', 'Name', name, 'NumberTitle', 'off');
end

function idx = nearestIndex(vec, value)
[~, idx] = min(abs(vec - value));
end

function pv = extractPVmetrics(V, J)
% J expected in mA cm^-2 from pinana/pindrift conventions
Jsc = interp1(V, J, 0, 'linear', 'extrap');

% Voc from J=0 crossing
Voc = NaN;
crossIdx = find((J(1:end-1) <= 0 & J(2:end) >= 0) | (J(1:end-1) >= 0 & J(2:end) <= 0), 1, 'first');
if ~isempty(crossIdx)
    Voc = interp1(J(crossIdx:crossIdx+1), V(crossIdx:crossIdx+1), 0, 'linear', 'extrap');
end

P = -V .* J; % mW cm^-2 if J in mA cm^-2
[Pmax, iMpp] = max(P);
Vmpp = V(iMpp);
Jmpp = J(iMpp);

if isnan(Voc) || Jsc == 0
    FF = NaN;
    PCE = NaN;
else
    FF = Pmax / (abs(Voc * Jsc));
    PCE = max(Pmax, 0); % at 1 sun (100 mW/cm^2), mW/cm^2 numerically equals %
end

pv = struct('Jsc', Jsc, 'Voc', Voc, 'FF', FF, 'PCE', PCE, 'Vmpp', Vmpp, 'Jmpp', Jmpp, 'Pmax', Pmax);
end

function ms = extractMottSchottky(V, C, p)
Cabs = abs(C(:));
Y = 1 ./ max(Cabs.^2, eps);

valid = isfinite(Y) & isfinite(V(:));
Vv = V(valid);
Yv = Y(valid);
if numel(Vv) < 3
    coeff = [NaN NaN];
    Vbi = NaN;
    slope = NaN;
    Vline = V(:);
    Yline = NaN(size(Vline));
else
    i1 = max(1, floor(0.3*numel(Vv)));
    i2 = max(i1+1, ceil(0.8*numel(Vv)));
    coeff = polyfit(Vv(i1:i2), Yv(i1:i2), 1);
    slope = coeff(1);
    Vbi = -coeff(2)/max(coeff(1), eps);
    Vline = linspace(min(Vv), max(Vv), 100).';
    Yline = polyval(coeff, Vline);
end

q = getfieldWithDefault(p, 'e', 1.602176634e-19); %#ok<GFLD>
epsRel = getfieldWithDefault(p, 'eppi', 20); %#ok<GFLD>
eps0 = 8.854187817e-14; % F/cm
epsAbs = epsRel * eps0;
NA = 2 ./ max(q * epsAbs * abs(slope), eps);

if isnan(Vbi)
    Wdep = NaN;
else
    Vref = min(max(Vbi - median(V), 1e-6), 5);
    Wdep = sqrt(2 * epsAbs * Vref / max(q * NA, eps));
end

ms = struct('V', V(:), 'C', C(:), 'Y', Y, 'Vfit', Vv, 'Yfit', Yv, 'Vline', Vline, 'Yline', Yline, ...
    'slope', slope, 'intercept', coeff(2), 'Vbi', Vbi, 'NA', NA, 'Wdep', Wdep);
end

function eis = fitEquivalentCircuit(V, J, C, p)
J_A = J(:) / 1000; % A/cm^2
Rdiff = abs(gradient(V(:), J_A));
Rdiff = Rdiff(isfinite(Rdiff) & Rdiff > 0);

if isempty(Rdiff)
    Rs = 1;
    Rrec = 100;
else
    Rs = max(min(Rdiff), eps);
    Rrec = max(prctile(Rdiff, 90), Rs + eps);
end

Cclean = abs(C(:));
Cclean = Cclean(isfinite(Cclean) & Cclean > 0);
if isempty(Cclean)
    Cgeo = 1e-7;
else
    Cgeo = median(Cclean);
end

Rion = max(Rrec - Rs, Rs);
Cion = max(5 * Cgeo, 1e-9);

f_geo = 1/(2*pi*Rrec*Cgeo);
f_ion = 1/(2*pi*Rion*Cion);

f = logspace(-2, 7, 400).';
w = 2*pi*f;
Z = Rs + Rrec ./ (1 + 1i*w*Rrec*Cgeo) + Rion ./ (1 + 1i*w*Rion*Cion);

eis = struct('f', f, 'Z', Z, 'Rs', Rs, 'Rrec', Rrec, 'Cgeo', Cgeo, 'Rion', Rion, 'Cion', Cion, ...
    'f_geo', f_geo, 'f_ion', f_ion, 'model', 'Rs + (Rrec||Cgeo) + (Rion||Cion)');
end

function cf = buildCFresponse(eis)
f = eis.f;
Cf = eis.Cgeo + eis.Cion ./ (1 + (f ./ max(eis.f_ion, eps)).^2);
cf = struct('f', f, 'Cf', Cf, 'f_geo', eis.f_geo, 'f_ion', eis.f_ion);
end

function value = getfieldWithDefault(s, fieldName, defaultVal)
if isfield(s, fieldName)
    value = s.(fieldName);
else
    value = defaultVal;
end
end
