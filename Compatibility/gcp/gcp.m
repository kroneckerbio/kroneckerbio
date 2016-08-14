function poolStruct = gcp( input_args )
%GCP Reproduces some of the behaviors of the gcp function for versions of
%MATLAB older than R2013b.

NumWorkers = matlabpool('size');

if NumWorkers > 0
    poolStruct.NumWorkers = NumWorkers;
else
    matlabpool open
    NumWorkers = matlabpool('size');
    poolStruct.NumWorkers = NumWorkers;
end

end

