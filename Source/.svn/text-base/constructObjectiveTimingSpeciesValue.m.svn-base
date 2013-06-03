function obj = constructObjectiveTimingSpeciesValue(m, speciesIndex, eventNumber)
% 
% 
% function obj = constructObjectiveOutputTimingSpeciesValue(m, speciesIndex, eventNumber)
% 
% This function constructs an objective function that returns the value of
% the timing of an event, as well as its derivatives.
% 
% Inputs
%       m               -   the model from which output is collected.
%       speciesIndex    -   the index of the species whose peak was
%                           detected by the event described in eventIndex.
%       eventNumber     -   the number of the event of interest.  For
%                           example, if the event of interest is a peak in
%                           concentration, setting eventNumber = 3 will
%                           compute the sensitivity of the third pulse timing.
% 
% Outputs
%       obj             -   an objective function structure containing
%                           three functions:
%                           obj.G is the objective function value
%                           obj.dGdx is the derivative with respect to species
%                           obj.dGdp is the derivative with respect to parameters
% 
% See also:
%   constructObjectiveOutputValue
%   constructObjectiveOutputISE
%   constructObjectiveOutputSSE
%   constructObjectiveSpeciesValue
%   constructObjectiveSpeciesISE
%   constructObjectiveSpeciesSSE
%   constructObjectiveTimingSpeciesSSE
%   constructObjectiveTimingOutputValue
%   constructObjectiveTimingOutputSSE

% This work is licensed under the Creative Commons Attribution-Noncommercial-No Derivative 
% Works 3.0 United States License. To view a copy of this license, visit 
% http://creativecommons.org/licenses/by-nc-nd/3.0/us/ or send a letter to Creative Commons,
% 171 Second Street, Suite 300, San Francisco, California, 94105, USA.


assert(length(speciesIndex) == 1, '', 'speciesIndex must be of length 1.');

% Compute the matrix that selects out the interesting species.
C = eye(m.nX);
C = C(speciesIndex,:);

obj.update      = @update;
obj.G           = @G;
obj.dGdx        = @dGdx;
obj.dGdp        = @dGdp;
obj.d2Gdx2      = @d2Gdx2;
obj.d2Gdp2      = @d2Gdp2;
obj.d2Gdxdp     = @d2Gdxdp;
obj.d2Gdpdx     = @d2Gdpdx;

    function objNew = update(mNew)
        objNew = constructObjectiveTimingSpeciesValue(mNew, speciesIndex, eventNumber);
    end

    function [val stopTimes] = G(t, xSol, u, nSim)
        % (1) Event detection
        val = inf;
        stopTimes = [];
        if isempty(xSol(nSim).ie), return, end      % No events detected
        
        ei = find(xSol(nSim).ie); % == eventIndex
        if length(ei) < eventNumber, return, end    % Not enough events detected
        
        % (2) Computing objective function value
        val       =  xSol(nSim).xe(ei(eventNumber));
        stopTimes =  xSol(nSim).xe(ei(eventNumber));
    end

    function val = dGdx(t, xSol, u, nSim)
        val = zeros(1, m.nX);
        if isempty(xSol(nSim).ie), return, end      % No events detected
        
        ei = find(xSol(nSim).ie); % == eventIndex
        if length(ei) < eventNumber, return, end    % Not enough events detected
        
        if t == xSol(nSim).xe(ei(eventNumber))
            ye      = deval(xSol(nSim), t, 1:m.nX);
            f       = m.f(t, ye, u);
            dfdx    = m.dfdx(t, ye, u);
            dfidx   = C*dfdx;
            val     = -1/(dfidx*f)*dfidx;
        end
    end

    function val = dGdp(t, xSol, u, nSim)
        val = zeros(1, m.nP);
        
        if isempty(xSol(nSim).ie), return, end      % No events detected
        
        ei = find(xSol(nSim).ie); % == eventIndex
        if length(ei) < eventNumber, return, end    % Not enough events detected
        
        if t == xSol(nSim).xe(ei(eventNumber))
            ye      = deval(xSol(nSim), t, 1:m.nX);
            f       = m.f(t, ye, u);
            dfdx    = m.dfdx(t, ye, u);
            dfidx   = C*dfdx;
            dfdp    = m.dfdp(t, ye, u);
            dfidp   = C*dfdp;
            val     = -1/(dfidx*f)*dfidp;
        end
    end

    function val = d2Gdx2(t, xSol, u, nSim)
        val = zeros(m.nX, m.nX);
        if isempty(xSol(nSim).ie), return, end      % No events detected
        
        ei = find(xSol(nSim).ie); % == eventIndex
        if length(ei) < eventNumber, return, end    % Not enough events detected
        
        if t == xSol(nSim).xe(ei(eventNumber))
            ye      = deval(xSol(nSim), t, 1:m.nX);
            f       = m.f(t, ye, u);
            dfdx    = m.dfdx(t, ye, u);
%             dfdp    = m.dfdp(t, ye, u);
            dfidx   = C*dfdx;
            
            d2fidx2  = m.d2fdx2(t, ye, u);
            d2fidx2  = d2fidx2((1:m.nX)+(speciesIndex-1)*m.nX,:);
%             d2fidxdp = m.d2fdxdp(t, ye, u);
%             d2fidxdp = d2fidxdp((1:m.nX)+(speciesIndex-1)*m.nX,:);
%             d2fidpdx = m.d2fdpdx(t, ye, u);
%             d2fidpdx = d2fidpdx((1:m.nP)+(speciesIndex-1)*m.nP,:);
            val     = dfidx'*(1/(dfidx*f)^2*(dfidx*dfdx + f'*d2fidx2))-1/(dfidx*f)*d2fidx2;
        end
    end
    function val = d2Gdp2(t, xSol, u, nSim)
        val = zeros(m.nP, m.nP);
        if isempty(xSol(nSim).ie), return, end      % No events detected
        
        ei = find(xSol(nSim).ie); % == eventIndex
        if length(ei) < eventNumber, return, end    % Not enough events detected
        
        if t == xSol(nSim).xe(ei(eventNumber))
            ye      = deval(xSol(nSim), t, 1:m.nX);
            f       = m.f(t, ye, u);
            dfdx    = m.dfdx(t, ye, u);
            dfdp    = m.dfdp(t, ye, u);
            dfidx   = C*dfdx;
            dfidp   = C*dfdp;
            
%             d2fidx2  = m.d2fdx2(t, ye, u);
%             d2fidx2  = d2fidx2((1:m.nX)+(speciesIndex-1)*m.nX,:);
            d2fidp2  = m.d2fdp2(t, ye, u);
            d2fidp2  = d2fidp2((1:m.nP)+(speciesIndex-1)*m.nP,:);
            d2fidxdp = m.d2fdxdp(t, ye, u);
            d2fidxdp = d2fidxdp((1:m.nX)+(speciesIndex-1)*m.nX,:);
%             d2fidpdx = m.d2fdpdx(t, ye, u);
%             d2fidpdx = d2fidpdx((1:m.nP)+(speciesIndex-1)*m.nP,:);
            
            val     = dfidp'*(1/(dfidx*f)^2*(dfidx*dfdp + f'*d2fidxdp))-1/(dfidx*f)*d2fidp2;
        end
    end
    function val = d2Gdpdx(t, xSol, u, nSim)
        val = zeros(m.nX, m.nP);
        if isempty(xSol(nSim).ie), return, end      % No events detected
        
        ei = find(xSol(nSim).ie); % == eventIndex
        if length(ei) < eventNumber, return, end    % Not enough events detected
        
        if t == xSol(nSim).xe(ei(eventNumber))
            ye      = deval(xSol(nSim), t, 1:m.nX);
            f       = m.f(t, ye, u);
            dfdx    = m.dfdx(t, ye, u);
            dfdp    = m.dfdp(t, ye, u);
            dfidx   = C*dfdx;
            dfidp   = C*dfdp;
            
%             d2fidx2  = m.d2fdx2(t, ye, u);
%             d2fidx2  = d2fidx2((1:m.nX)+(speciesIndex-1)*m.nX,:);
%             d2fidp2  = m.d2fdp2(t, ye, u);
%             d2fidp2  = d2fidp2((1:m.nP)+(speciesIndex-1)*m.nP,:);
            d2fidxdp = m.d2fdxdp(t, ye, u);
            d2fidxdp = d2fidxdp((1:m.nX)+(speciesIndex-1)*m.nX,:);
            d2fidpdx = m.d2fdpdx(t, ye, u);
            d2fidpdx = d2fidpdx((1:m.nP)+(speciesIndex-1)*m.nP,:);
            
            val     = dfidx'*(1/(dfidx*f)^2*(dfidx*dfdp + f'*d2fidxdp))-1/(dfidx*f)*d2fidxdp;
        end
    end
    function val = d2Gdxdp(t, xSol, u, nSim)
        val = zeros(m.nP, m.nX);
        if isempty(xSol(nSim).ie), return, end      % No events detected
        
        ei = find(xSol(nSim).ie); % == eventIndex
        if length(ei) < eventNumber, return, end    % Not enough events detected
        
        if t == xSol(nSim).xe(ei(eventNumber))
            ye      = deval(xSol(nSim), t, 1:m.nX);
            f       = m.f(t, ye, u);
            dfdx    = m.dfdx(t, ye, u);
            dfdp    = m.dfdp(t, ye, u);
            dfidx   = C*dfdx;
            dfidp   = C*dfdp;
            
            d2fidx2  = m.d2fdx2(t, ye, u);
            d2fidx2  = d2fidx2((1:m.nX)+(speciesIndex-1)*m.nX,:);
%             d2fidp2  = m.d2fdp2(t, ye, u);
%             d2fidp2  = d2fidp2((1:m.nP)+(speciesIndex-1)*m.nP,:);
%             d2fidxdp = m.d2fdxdp(t, ye, u);
%             d2fidxdp = d2fidxdp((1:m.nX)+(speciesIndex-1)*m.nX,:);
            d2fidpdx = m.d2fdpdx(t, ye, u);
            d2fidpdx = d2fidpdx((1:m.nP)+(speciesIndex-1)*m.nP,:);
            
            val     = dfidp'*(1/(dfidx*f)^2*(dfidx*dfdx + f'*d2fidx2))-1/(dfidx*f)*d2fidpdx;
        end
    end
end