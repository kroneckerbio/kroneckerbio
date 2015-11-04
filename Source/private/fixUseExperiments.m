function UseExperiments = fixUseExperiments(UseExperiments, n_obs, n_con)
%useExperiments = fixUseExperiments(useExperiments)

if isnumeric(UseExperiments)
    temp = false(n_obs, n_con);
    temp(UseExperiments) = true;
    UseExperiments = temp;
elseif islogical(UseExperiments)
    assert(ismatrix(UseExperiments) && all(size(UseExperiments) == [n_obs, n_con]), 'KroneckerBio:fixUseExperiments:SizeOfUseExperiments', 'If UseExperiments is provided as a logical array, it must be the size of n_obs by n_con')
else
    error('KroneckerBio:fixUseExperiments:InvalidClass', 'UseExperiments must be a logical array or vector of indexes')
end
