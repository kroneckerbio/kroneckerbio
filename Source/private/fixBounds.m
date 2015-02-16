function bounds = fixBounds(bounds, UseParams, UseSeeds, UseInputControls, UseDoseControls)

% Constants
ns = size(UseSeeds, 1);
nk = numel(UseParams);
nq = numel(cat(1,UseInputControls{:}));
nh = numel(cat(1,UseDoseControls{:}));
nTk = nnz(UseParams);
nTs = nnz(UseSeeds);
nTq = nnz(cat(1,UseInputControls{:}));
nTh = nnz(cat(1,UseDoseControls{:}));
nT = nTk + nTs + nTq + nTh;
nCon = size(UseSeeds, 2);


% Bounds can be nk+ns*nCon+nq+nh, nT, nk, nk+ns*nCon, or 1
l = numel(bounds);

if l == 1
    bounds = zeros(nT,1) + bounds;
elseif l == nk+ns*nCon+nq+nh
    bounds = bounds([UseParams; vec(UseSeeds); cat(1,UseInputControls{:})]);
elseif l == nT
    %bounds = bounds;
elseif l == nk && nTs == 0 && nTq == 0 && nTh == 0
    bounds = bounds(UseParams);
elseif l == nk+ns*nCon && nTq == 0 && nTh == 0
    bounds = bounds([UseParams; vec(UseSeeds)]);
else
    error('KroneckerBio:BoundSize', ...
        'LowerBound and UpperBound must be vectors the length of m.nk+m.ns*numel(con)+m.nq, number of variable parameters, m.nk if there are no variable ICs or controls, m.nk+m.ns*numel(con) if there are no variable controls, or scalar')
end
