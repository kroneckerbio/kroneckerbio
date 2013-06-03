function obj = constructObjectiveTimingOutputValue(m, outputIndex, eventNumber)
% 
% 
% function obj = constructObjectiveTimingOutputValue(m, outputIndex, eventNumber)
% 
% This function constructs an objective function that returns the value of
% the timing of an event, as well as its derivatives.
% 
% Inputs
%       m               -   the model from which output is collected.
%       outputIndex     -   the index of the output whose peak was
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
%                           obj.dGdp is the derivative with respect to
%                           parameters
% 
% See also:
%   constructObjectiveOutputValue
%   constructObjectiveOutputISE
%   constructObjectiveOutputSSE
%   constructObjectiveSpeciesValue
%   constructObjectiveSpeciesISE
%   constructObjectiveSpeciesSSE
%   constructObjectiveTimingSpeciesValue
%   constructObjectiveTimingSpeciesSSE
%   constructObjectiveTimingOutputSSE

% This work is licensed under the Creative Commons Attribution-Noncommercial-No Derivative 
% Works 3.0 United States License. To view a copy of this license, visit 
% http://creativecommons.org/licenses/by-nc-nd/3.0/us/ or send a letter to Creative Commons,
% 171 Second Street, Suite 300, San Francisco, California, 94105, USA.


assert(isnumeric(outputIndex), '', 'outputIndex must have a numeric value');
assert(length(outputIndex) == 1, '', 'outputIndex must be of length 1.');

nS       = m.nP+m.nX;
nP       = m.nP;
Ix       = eye(m.nX);
Ip       = eye(m.nP);

obj.update      = @update;
obj.G           = @G;
obj.dGdx        = @dGdx;
obj.dGdp        = @dGdp;
obj.d2Gdx2      = @d2Gdx2;
obj.d2Gdp2      = @d2Gdp2;
obj.d2Gdxdp     = @d2Gdxdp;
obj.d2Gdpdx     = @d2Gdpdx;

%% update function
    function objNew = update(mNew)
        objNew = constructObjectiveTimingOutputValue(mNew, outputIndex, eventNumber);
    end
%% G
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

%% dGdx
    function val = dGdx(t, xSol, u, nSim)
        val = zeros(1, m.nX);
        if isempty(xSol(nSim).ie), return, end      % No events detected
        
        ei = find(xSol(nSim).ie); % == eventIndex
        if length(ei) < eventNumber, return, end    % Not enough events detected
        
        if t == xSol(nSim).xe(ei(eventNumber))
            ye      = deval(xSol(nSim), t, 1:m.nX);
            f       = m.f(t, ye, u);
            dfdx    = m.dfdx(t, ye, u);
            dfidx   = m.c(outputIndex,:)*dfdx;
            val     = -1/(dfidx*f)*dfidx;
        end
    end

%% dGdp
    function val = dGdp(t, xSol, u, nSim)
        val = zeros(1, m.nP);
        
        if isempty(xSol(nSim).ie), return, end      % No events detected
        
        ei = find(xSol(nSim).ie); % == eventIndex
        if length(ei) < eventNumber, return, end    % Not enough events detected
        
        if t == xSol(nSim).xe(ei(eventNumber))
            ye      = deval(xSol(nSim), t, 1:m.nX);
            f       = m.f(t, ye, u);
            dfdx    = m.dfdx(t, ye, u);
            dfidx   = m.c(outputIndex,:)*dfdx;
            dfdp    = m.dfdp(t, ye, u);
            dfidp   = m.c(outputIndex,:)*dfdp;
            val     = -1/(dfidx*f)*dfidp;
        end
    end

%% d2Gdx2
    function val = d2Gdx2(t, xSol, u, nSim)
        val = zeros(m.nX, m.nX);
        if isempty(xSol(nSim).ie), return, end      % No events detected
        
        ei = find(xSol(nSim).ie); % == eventIndex
        if length(ei) < eventNumber, return, end    % Not enough events detected
        
        if t == xSol(nSim).xe(ei(eventNumber))
            ye       = deval(xSol(nSim), t, 1:m.nX);
            f        = m.f(t, ye, u);
            dfdx     = m.dfdx(t, ye, u);
            dfidx    = m.c(outputIndex,:)*dfdx;
            dfdp     = m.dfdp(t, ye, u);
            dfidp    = m.c(outputIndex,:)*dfdp;
            d2fdp2   = m.d2fdp2(t, ye, u);
            d2fidp2  = kron(m.c(outputIndex,:),Ip)*d2fdp2;
            d2fdx2   = m.d2fdx2(t, ye, u);
            d2fidx2  = kron(m.c(outputIndex,:),Ix)*d2fdx2;
            d2fdxdp  = m.d2fdxdp(t, ye, u);
            d2fidxdp = kron(m.c(outputIndex,:),Ix)*d2fdxdp;
            d2fdpdx  = m.d2fdpdx(t, ye, u);
            d2fidpdx = kron(m.c(outputIndex,:),Ip)*d2fdpdx;
            val      = -1/(dfidx*f)*(d2fidx2  + -1/(dfidx*f)*d2fidx2*f*dfidx  + -1/(dfidx*f)*dfdx'*dfidx'*dfidx) + 1/(dfidx*f)^2*dfidx'*(f'*(-1/(dfidx*f)*d2fidx2*f*dfidx + d2fidx2)  + dfidx*(-1/(dfidx*f)*dfdx*f*dfidx + dfdx));
        end
    end

%% d2Gdp2
    function val = d2Gdp2(t, xSol, u, nSim)
        val = zeros(m.nP, m.nP);
        
        if isempty(xSol(nSim).ie), return, end      % No events detected
        
        ei = find(xSol(nSim).ie); % == eventIndex
        if length(ei) < eventNumber, return, end    % Not enough events detected
        
        if t == xSol(nSim).xe(ei(eventNumber))
            ye       = deval(xSol(nSim), t, 1:m.nX);
            f        = m.f(t, ye, u);
            dfdx     = m.dfdx(t, ye, u);
            dfidx    = m.c(outputIndex,:)*dfdx;
            dfdp     = m.dfdp(t, ye, u);
            dfidp    = m.c(outputIndex,:)*dfdp;
            d2fdp2   = m.d2fdp2(t, ye, u);
            d2fidp2  = kron(m.c(outputIndex,:),Ip)*d2fdp2;
            d2fdx2   = m.d2fdx2(t, ye, u);
            d2fidx2  = kron(m.c(outputIndex,:),Ix)*d2fdx2;
            d2fdxdp  = m.d2fdxdp(t, ye, u);
            d2fidxdp = kron(m.c(outputIndex,:),Ix)*d2fdxdp;
            d2fdpdx  = m.d2fdpdx(t, ye, u);
            d2fidpdx = kron(m.c(outputIndex,:),Ip)*d2fdpdx;
            val      = -1/(dfidx*f)*(d2fidp2  + -1/(dfidx*f)*d2fidpdx*f*dfidp + -1/(dfidx*f)*dfdp'*dfidx'*dfidp) + 1/(dfidx*f)^2*dfidp'*(f'*(-1/(dfidx*f)*d2fidx2*f*dfidp + d2fidxdp) + dfidx*(-1/(dfidx*f)*dfdx*f*dfidp + dfdp));
        end
    end

%% d2Gdxdp
    function val = d2Gdpdx(t, xSol, u, nSim)
        val = zeros(m.nX, m.nP);
        if isempty(xSol(nSim).ie), return, end      % No events detected
        
        ei = find(xSol(nSim).ie); % == eventIndex
        if length(ei) < eventNumber, return, end    % Not enough events detected
        
        if t == xSol(nSim).xe(ei(eventNumber))
            ye       = deval(xSol(nSim), t, 1:m.nX);
            f        = m.f(t, ye, u);
            dfdx     = m.dfdx(t, ye, u);
            dfidx    = m.c(outputIndex,:)*dfdx;
            dfdp     = m.dfdp(t, ye, u);
            dfidp    = m.c(outputIndex,:)*dfdp;
            d2fdp2   = m.d2fdp2(t, ye, u);
            d2fidp2  = kron(m.c(outputIndex,:),Ip)*d2fdp2;
            d2fdx2   = m.d2fdx2(t, ye, u);
            d2fidx2  = kron(m.c(outputIndex,:),Ix)*d2fdx2;
            d2fdxdp  = m.d2fdxdp(t, ye, u);
            d2fidxdp = kron(m.c(outputIndex,:),Ix)*d2fdxdp;
            d2fdpdx  = m.d2fdpdx(t, ye, u);
            d2fidpdx = kron(m.c(outputIndex,:),Ip)*d2fdpdx;
            val      = -1/(dfidx*f)*(d2fidxdp + -1/(dfidx*f)*d2fidx2*f*dfidp  + -1/(dfidx*f)*dfdx'*dfidx'*dfidp) + 1/(dfidx*f)^2*dfidx'*(f'*(-1/(dfidx*f)*d2fidx2*f*dfidp + d2fidxdp) + dfidx*(-1/(dfidx*f)*dfdx*f*dfidp + dfdp));
        end
    end

%% d2Gdpdx
    function val = d2Gdxdp(t, xSol, u, nSim)
        val = zeros(m.nP, m.nX);
        
        if isempty(xSol(nSim).ie), return, end      % No events detected
        
        ei = find(xSol(nSim).ie); % == eventIndex
        if length(ei) < eventNumber, return, end    % Not enough events detected
        
        if t == xSol(nSim).xe(ei(eventNumber))
            ye       = deval(xSol(nSim), t, 1:m.nX);
            f        = m.f(t, ye, u);
            dfdx     = m.dfdx(t, ye, u);
            dfidx    = m.c(outputIndex,:)*dfdx;
            dfdp     = m.dfdp(t, ye, u);
            dfidp    = m.c(outputIndex,:)*dfdp;
            d2fdp2   = m.d2fdp2(t, ye, u);
            d2fidp2  = kron(m.c(outputIndex,:),Ip)*d2fdp2;
            d2fdx2   = m.d2fdx2(t, ye, u);
            d2fidx2  = kron(m.c(outputIndex,:),Ix)*d2fdx2;
            d2fdxdp  = m.d2fdxdp(t, ye, u);
            d2fidxdp = kron(m.c(outputIndex,:),Ix)*d2fdxdp;
            d2fdpdx  = m.d2fdpdx(t, ye, u);
            d2fidpdx = kron(m.c(outputIndex,:),Ip)*d2fdpdx;
            val      = -1/(dfidx*f)*(d2fidpdx + -1/(dfidx*f)*d2fidpdx*f*dfidx + -1/(dfidx*f)*dfdp'*dfidx'*dfidx) + 1/(dfidx*f)^2*dfidp'*(f'*(-1/(dfidx*f)*d2fidx2*f*dfidx + d2fidx2)  + dfidx*(-1/(dfidx*f)*dfdx*f*dfidx + dfdx));
        end
    end
end
