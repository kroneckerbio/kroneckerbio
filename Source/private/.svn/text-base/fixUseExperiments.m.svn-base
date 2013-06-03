function useExperiments = fixUseExperiments(useExperiments, nObj, nCon)
%useExperiments = fixUseExperiments(useExperiments)

if isnumeric(useExperiments)
    temp = false(nObj, nCon);
    temp(useExperiments) = true;
    useExperiments = temp;
elseif islogical(useExperiments)
    assert(ndims(useExperiments) == 2 && all(size(useExperiments) == [nObj, nCon]), 'KroneckerBio:fixUseExperiments:SizeOfUseExperiments', 'If UseExperiments is provided as a logical array, it must be the size of nObj by nCon')
else
    error('KroneckerBio:fixUseExperiments:InvalidClass', 'UseExperiments must be a logical array or vector of indexes')
end
