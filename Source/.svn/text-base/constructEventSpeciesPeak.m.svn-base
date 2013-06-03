function eventFcn = constructEventSpeciesPeak(m, speciesIndex, maxOrMin)
% 
% 
% function eventFcn = constructEventSpeciesPeak(m, speciesIndex, maxOrMin)
% 
% This function constructs an events function that looks for local maxima
% of the model output.  It stops integration when it finds one.
% 
% Inputs
%       m               -   The model from which output is collected
%       speciesIndex    -   The index of the output for which the ISE will
%                           be computed
%       maxOrMin        -   Specifies whether to collect maxima or minima.
%                           This argument takes 3 values:   1 = maxima
%                                                           0 = all extrema
%                                                          -1 = minima
% 
% See MATLAB's ODE solver documentation (for instance, type 'doc ode15s')
% for more information on events functions.
% 
% See also:
%       constructEventOutputPeak
%       ode15s

% This work is licensed under the Creative Commons
% Attribution-Noncommercial-No Derivative 
% Works 3.0 United States License. To view a copy of this license, visit 
% http://creativecommons.org/licenses/by-nc-nd/3.0/us/ or send a letter to Creative Commons,
% 171 Second Street, Suite 300, San Francisco, California, 94105, USA.

% argument processing and error checking
if nargin<2 || isempty(speciesIndex)
    speciesIndex = 1;
end

speciesIndex = speciesIndex(:);

if nargin<3 || isempty(maxOrMin)
    maxOrMin = ones(size(speciesIndex));   % by default, collect only maxima.
end

assert(isnumeric(speciesIndex), '', 'The second argument to constructEventSpeciesPeak must be a state index.');
assert(all(speciesIndex <= m.nX), '', 'The state index %d is greater than the total number of states %d.', speciesIndex, m.nX);

% return event function
eventFcn.events = @events;
eventFcn.update = @update;

    function eventFcn = update(mNew)
        eventFcn = constructEventSpeciesPeak(mNew, speciesIndex, maxOrMin);
    end

%% The events function that finds peaks in the state trajectory
    function [val isTerm direction] = events(t, x, u)
        val       =  real(m.f(t, x, u));
        val       =  val(speciesIndex);
        isTerm    =  ones(size(speciesIndex));
        direction = -1*maxOrMin;
    end
end