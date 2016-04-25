function NumWorkers = setParallelPoolSize(NumWorkers)
% setParallelPoolSize(NumWorkers)
% Sets the number of parallel pool workers to the provided number.
%
% Input arguments:
%   NumWorkers [ positive integer scalar {[]} ]
%       Number of parallel workers to use. If [], the value is set to
%       the maximum number of workers specified in the 'local' cluster
%       profile, which is typically the number of cores on the machine.
%       If a parallel pool is not initialized, a pool will be
%       initialized on the local cluster with the provided number of
%       workers. If a pool is initialized, but has the wrong number of
%       workers, the existing pool will be closed, and a new pool will
%       be opened with the provided number of workers. If a pool is
%       already open with the provided number of workers, no new pool
%       will be opened.

if isempty(NumWorkers)
    % Get default number of max workers from local cluster profile
    localCluster = parcluster('local');
    NumWorkers = localCluster.NumWorkers;
end

try
    % Creating new temporary cluster profile allowing
    % specified number of workers
    myCluster = parcluster('local');
    myCluster.NumWorkers = NumWorkers;
    parpool(myCluster, NumWorkers);
catch ME
    if strcmp(ME.identifier, 'parallel:convenience:ConnectionOpen')
        % If a pool is already open, check that it has the
        % correct number of workers
        p = gcp;
        if p.NumWorkers ~= NumWorkers
            % Restart pool with correct number of workers
            p.delete;
            parpool(myCluster, NumWorkers);
        end
    else
        rethrow(ME);
    end
end

end