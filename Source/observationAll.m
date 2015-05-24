function obs = observationAll(tF)

obs.tF = tF;
obs.Complex = true;
obs.Simulation = @simulation;
obs.Sensitivity = @sensitivity;
obs.Curvature = @curvature;

obs = pastestruct(observationZero(), obs);

    function sim = simulation(int)
        sim.Type = 'Simulation.System.All';
        sim.Name = int.Name;
        sim.t    = int.t;
        sim.x    = @(t, varargin)evaluate_states(int, t, varargin{:});
        sim.u    = @(t, varargin)evaluate_inputs(int, t, varargin{:});
        sim.y    = @(t, varargin)evaluate_outputs(int, t, varargin{:});
        sim.int  = int;
    end

    function sim = sensitivity(int)
        sim = simulation(int);
        
        sim.Type = 'Simulation.Sensitivity.All';
        
        sim.dxdT = @(t, varargin)evaluate_state_sensitivity(int, t, varargin{:});
        sim.dudT = @(t, varargin)evaluate_input_sensitivity(int, t, varargin{:});
        sim.dydT = @(t, varargin)evaluate_output_sensitivity(int, t, varargin{:});
    end

    function sim = curvature(int)
        sim = sensitivity(int);
        
        sim.Type = 'Simulation.Curvature.All';
        
        sim.d2xdT2 = @(t, varargin)evaluate_state_curvature(int, t, varargin{:});
        sim.d2udT2 = @(t, varargin)evaluate_input_curvature(int, t, varargin{:});
        sim.d2ydT2 = @(t, varargin)evaluate_output_curvature(int, t, varargin{:});
    end
end

function val = evaluate_states(int, t, ind)
val = int.x(t);
if nargin >= 3
    val = val(ind,:);
end
end

function val = evaluate_inputs(int, t, ind)
val = int.u(t);
if nargin >= 3
    val = val(ind,:);
end
end

function val = evaluate_outputs(int, t, ind)
val = int.y(t);
if nargin >= 3
    val = val(ind,:);
end
end

function val = evaluate_state_sensitivity(int, t, ind)
nx = int.nx;
nT = int.nT;
nt = numel(t);

val = int.dxdT(t);
if nargin >= 3
    val = reshape(val, nx,nT*nt); % x_Tt
    val = val(ind,:); % x_Tt chopped out rows
    val = reshape(val, nnz(ind)*nT,nt); % xT_t
end
end

function val = evaluate_input_sensitivity(int, t, ind)
nu = int.nu;
nT = int.nT;
nt = numel(t);

val = int.dudT(t);
if nargin >= 3
    val = reshape(val, nu,nT*nt); % u_Tt
    val = val(ind,:); % u_Tt chopped out rows
    val = reshape(val, nnz(ind)*nT,nt); % uT_t
end
end

function val = evaluate_output_sensitivity(int, t, ind)
ny = int.ny;
nT = int.nT;
nt = numel(t);

val = int.dydT(t);
if nargin >= 3
    val = reshape(val, ny,nT*nt); % y_Tt
    val = val(ind,:); % y_Tt chopped out rows
    val = reshape(val, nnz(ind)*nT,nt); % yT_t
end
end

function val = evaluate_state_curvature(int, t, ind)
nx = int.nx;
nT = int.nT;
nt = numel(t);

val = int.d2xdT2(t);
if nargin >= 3
    val = reshape(val, nx,nT*nT*nt); % x_TTt
    val = val(ind,:); % x_TTt chopped out rows
    val = reshape(val, nnz(ind)*nT*nT,nt); % xTT_t
end
end

function val = evaluate_input_curvature(int, t, ind)
nu = int.nu;
nT = int.nT;
nt = numel(t);

val = int.d2udT2(t);
if nargin >= 3
    val = reshape(val, nu,nT*nT*nt); % u_TTt
    val = val(ind,:); % u_TTt chopped out rows
    val = reshape(val, nnz(ind)*nT*nT,nt); % uTT_t
end
end

function val = evaluate_output_curvature(int, t, ind)
ny = int.ny;
nT = int.nT;
nt = numel(t);

val = int.d2ydT2(t);
if nargin >= 3
    val = reshape(val, ny,nT*nT*nt); % y_TTt
    val = val(ind,:); % y_TTt chopped out rows
    val = reshape(val, nnz(ind)*nT*nT,nt); % yTT_t
end
end
