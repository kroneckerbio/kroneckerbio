function scurr = rng(seed)
% Function written for old versions of MATLAB that don't have the rng
% function. This function only implements changing the seed of the random
% number generator either by providing a seed integer or a struct returned
% by this function.

if nargin < 1
    seed = [];
end

% Get current RandStream gentype
sdefault = RandStream.getDefaultStream;
scurr.Type = sdefault.Type;
scurr.Seed = sdefault.Seed;
scurr.State = sdefault.State;

if ~isempty(seed)
    % Generate a new stream with the provided seed
    if isnumeric(seed)
        snew = RandStream(sdefault.Type, 'Seed', seed);
    elseif isstruct(seed)
        snew = RandStream(seed.Type, 'Seed', seed.Seed);
        snew.State = seed.State;
    end
    % Set that stream as the default
    RandStream.setDefaultStream(snew);
end

end