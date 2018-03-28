function eve = eventRisesAbove(m, output, threshold)
%eventRisesAbove Detect when an output drops below a threshold
%
%   eve = eventRisesAbove(m, output, threshold)
%
%   Given an output index and a threshold, detect any time that output rises
%   above the threshold.
%
%   Inputs
%   m: [ model struct scalar ]
%       The KroneckerBio model for which the simulation will be run.
%   output: [ positive integer ]
%       An index into the model outputs.
%   threshold: [ numeric scalar ]
%       The number which will trigger the event.
%
%   Outputs
%   eve: [ Event scalar ]
%       A KroneckerBio event struct.
%
%   For the meanings of the fields of eve see "help Event"

% (c) 2018 David R Hagen and Bruce Tidor
% This work is released under the MIT license.


% real() call is needed for imaginary finite differences to work with this
% event
eve = Event(1, @(t,y)real(y(output))-threshold);
