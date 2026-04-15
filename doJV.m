function JVsol = doJV(sol, varargin)
%DOJV Configure and run a linear J-V sweep with driftfusion.
%
% Syntax:
%   JVsol = doJV(sol, par, scan_rate, scan_points, Int, Vstart, Vend)
%   JVsol = doJV(..., calcJ)
%   JVsol = doJV(..., calcJ, BC)
%   JVsol = doJV(sol, scan_rate, scan_points, Int, calcJ, Vstart, Vend)
%   JVsol = doJV(sol, scan_rate, scan_points, Int, calcJ, Vstart, Vend, BC)
%
% Inputs:
%   sol         - initial solution structure (e.g. from equilibrate)
%   par         - parameter struct (e.g. pinParams); use [] to reuse sol.p
%   scan_rate   - scan rate [V s^-1]
%   scan_points - number of points in voltage sweep
%   Int         - illumination intensity (0 dark, 1 one-sun)
%   Vstart      - start voltage [V]
%   Vend        - end voltage [V]
%   calcJ       - optional current method (default: par.calcJ or 4)
%   BC          - optional boundary condition (default: par.BC)
%
% Output:
%   JVsol       - pindrift output with guaranteed JVsol.Vapp and JVsol.Jn

if ~isstruct(sol) || ~isfield(sol, 'sol')
    error('doJV:InvalidInput', 'sol must be a driftfusion solution struct containing field .sol.');
end

% Parse signatures:
%   New: doJV(sol, par, scan_rate, scan_points, Int, Vstart, Vend, [calcJ], [BC])
%   Legacy: doJV(sol, scan_rate, scan_points, Int, calcJ, Vstart, Vend, [BC])
if numel(varargin) < 6
    error('doJV:InvalidInput', 'Not enough input arguments for doJV.');
end

if isstruct(varargin{1}) || isempty(varargin{1})
    par = varargin{1};
    scan_rate = varargin{2};
    scan_points = varargin{3};
    Int = varargin{4};
    Vstart = varargin{5};
    Vend = varargin{6};
    opt = varargin(7:end);
else
    par = [];
    scan_rate = varargin{1};
    scan_points = varargin{2};
    Int = varargin{3};
    Vstart = varargin{5};
    Vend = varargin{6};
    opt = varargin(7:end);
    opt = [{varargin{4}}, opt]; % legacy 5th input is calcJ
end

if isempty(par)
    if isfield(sol, 'p')
        p = sol.p;
    else
        p = pinParams;
    end
else
    p = par;
end

if ~isstruct(p)
    error('doJV:InvalidParameters', 'This driftfusion version expects par to be a parameter struct (e.g. from pinParams).');
end

if ~exist('pindrift', 'file')
    error('doJV:MissingDependency', 'pindrift.m was not found in the MATLAB path.');
end

if scan_rate <= 0
    error('doJV:InvalidInput', 'scan_rate must be > 0.');
end
if scan_points < 2
    error('doJV:InvalidInput', 'scan_points must be >= 2.');
end

calcJ = pickOptional(opt, 1, getfieldWithDefault(p, 'calcJ', 4)); %#ok<GFLD>
BC = pickOptional(opt, 2, getfieldWithDefault(p, 'BC', 3)); %#ok<GFLD>

% Configure JV scan parameters
p.figson = 0;
p.Ana = 1;
p.pulseon = 0;
p.OC = 0;
p.JV = 1;
p.Int = Int;
p.calcJ = calcJ;
p.BC = BC;
p.Vstart = Vstart;
p.Vend = Vend;
p.JVscan_rate = scan_rate;
p.JVscan_pnts = scan_points;
p.tmax = abs(Vend - Vstart) / scan_rate;
p.t0 = 0;
p.tmesh_type = 1;
p.tpoints = scan_points;

% Run simulation
JVsol = pindrift(sol, p);

% Ensure J-V arrays are available
if ~isfield(JVsol, 'Vapp') || ~isfield(JVsol, 'Jn') || isempty(JVsol.Vapp) || isempty(JVsol.Jn)
    [Vapp, Jn, ~] = pinana(JVsol);
    JVsol.Vapp = Vapp;
    JVsol.Jn = Jn;
end

JVsol.JVmeta.scan_rate = scan_rate;
JVsol.JVmeta.scan_points = scan_points;
JVsol.JVmeta.Int = Int;
JVsol.JVmeta.Vstart = Vstart;
JVsol.JVmeta.Vend = Vend;
JVsol.JVmeta.calcJ = calcJ;
JVsol.JVmeta.BC = BC;

end

function val = pickOptional(args, idx, defaultVal)
if numel(args) >= idx && ~isempty(args{idx})
    val = args{idx};
else
    val = defaultVal;
end
end

function value = getfieldWithDefault(s, fieldName, defaultVal)
if isfield(s, fieldName)
    value = s.(fieldName);
else
    value = defaultVal;
end
end
