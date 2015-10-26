function ic = StateInitialCondition(ic)
% Allowed Analytic model state initial condition inputs:
%   - empty -> set to default initial amount = 0, no seed
%   - scalar number -> value set to initial amount, no seed
%   - string expression -> taken as is, may contain seeds

if isempty(ic)
    ic = 0;
end

if ischar(ic)
    % pass
elseif isnumeric(ic) && isscalar(ic)
    ic = num2str(ic);
else
    error('KroneckerBio:FieldValidator:StateInitialCondition', 'Initial condition type not recognized')
end