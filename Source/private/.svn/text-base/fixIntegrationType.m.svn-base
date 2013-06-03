function [continuous complex tGet] = fixIntegrationType(con, obj)
%FIXINTEGRATIONTYPE Standardize continuous, complex, and tGet
%
%   There are two attributes of objective functions that affect how the
%   system is integrated. Continuous indicates that there is a continuous
%   aspect to the objective function. If any objective function has this
%   atribute, then the integration must integrate it. Complex indicates
%   that the objective must have the complete solution in order to compute
%   the objective (like to optimize the peak time or peak value of a
%   species). Objective functions that are not complex must have
%   discreteTimes indicated to tell when the objective function needs to be
%   evaluated.

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Constants
[nObj nCon nTop] = size(obj);

% Initialize
continuous = false(nCon,nTop);
complex = false(nCon,nTop);
tGet = cell(nCon,nTop);

% Loop through each objective to determine if any require special treatment
for iTop = 1:nTop
    for iCon = 1:nCon
        for iObj = 1:nObj
            continuous(iCon,iTop) = continuous(iCon,iTop) || obj(iObj,iCon,iTop).Continuous;
            complex(iCon,iTop) = complex(iCon,iTop) || obj(iObj,iCon,iTop).Complex;
        end
        
        % Only for non-complex integrations is tGet useful
        if ~complex(iCon,iTop)
            discreteTimes = [];
            for iObj = 1:nObj
                discreteTimes = [discreteTimes; vec(obj(iObj,iCon,iTop).DiscreteTimes)];
            end
            if continuous(iCon,iTop)
                % We need the final point as well
                discreteTimes = [discreteTimes; con(iCon).tF];
            end
            tGet{iCon,iTop} = unique(discreteTimes);
        end
    end
end