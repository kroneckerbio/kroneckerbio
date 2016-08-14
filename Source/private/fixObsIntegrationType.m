function [complex, t_get] = fixObsIntegrationType(obs)
% fixObsIntegrationType Standardize complex and tGet for observations

% (c) 2015 David R Hagen
% This work is released under the MIT license.

% Constants
[nObj, nCon, nTop] = size(obs);

% Initialize
complex = false(nCon,nTop);
t_get = cell(nCon,nTop);

% Loop through each objective to determine if any require special treatment
for iTop = 1:nTop
    for iCon = 1:nCon
        for iObj = 1:nObj
            complex(iCon,iTop) = complex(iCon,iTop) || obs(iObj,iCon,iTop).Complex;
        end
        
        % Only for non-complex integrations is tGet useful
        if ~complex(iCon,iTop)
            discrete_times = [];
            for iObj = 1:nObj
                discrete_times = [discrete_times, row(obs(iObj,iCon,iTop).t)];
            end
            t_get{iCon,iTop} = unique2013a(discrete_times);
        end
    end
end
