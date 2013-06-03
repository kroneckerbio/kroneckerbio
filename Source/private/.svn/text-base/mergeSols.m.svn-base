function xSol = mergeSols(xSol1, xSol2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Merge solution structure to preserve deval() ability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This work is licensed under the Creative Commons Attribution-Noncommercial-No Derivative 
% Works 3.0 United States License. To view a copy of this license, visit 
% http://creativecommons.org/licenses/by-nc-nd/3.0/us/ or send a letter to Creative Commons,
% 171 Second Street, Suite 300, San Francisco, California, 94105, USA.

xSol = xSol2;

if isempty(xSol1)
    return
end

% append xSol2 onto the end of xSol1
xSol.x = [xSol1.x, xSol2.x(2:end)];
xSol.y = [xSol1.y, xSol2.y(:,2:end)];

switch xSol1.solver
    case 'ode15s'
        size1 = size(xSol1.idata.dif3d, 2);
        size2 = size(xSol2.idata.dif3d, 2);
        if(size1 == size2)
            dif3d1 = xSol1.idata.dif3d;
            dif3d2 = xSol2.idata.dif3d;
        elseif(size1 < size2)
            dif3d1 = padarray(xSol1.idata.dif3d, [0, size2-size1, 0], 'post');
            dif3d2 = xSol2.idata.dif3d;
        else % size1 bigger
            dif3d1 = xSol1.idata.dif3d;
            dif3d2 = padarray(xSol2.idata.dif3d, [0, size1-size2, 0], 'post');
        end

        xSol.idata.kvec             = [xSol1.idata.kvec  xSol2.idata.kvec(2:end)];
        xSol.idata.dif3d            = cat(3, dif3d1, dif3d2(:, :, 2:end));
        xSol.idata.idxNonNegative   = [xSol1.idata.idxNonNegative xSol2.idata.idxNonNegative(2:end)];
        
    case 'ode23s'
        xSol.idata.k1 = [xSol1.idata.k1 xSol2.idata.k1(:,2:end)];
        xSol.idata.k2 = [xSol1.idata.k2 xSol2.idata.k2(:,2:end)];
    case 'ode45'
        xSol.idata.f3d = cat(3, xSol1.idata.f3d, xSol2.idata.f3d(:,:,2:end));
        xSol.idata.idxNonNegative = [xSol1.idata.idxNonNegative xSol2.idata.idxNonNegative(2:end)];
end

% if there were events in either, append these as well
if isfield(xSol1, 'xe') && isfield(xSol2, 'xe')
    xSol.xe = [xSol1.xe xSol2.xe];
    xSol.ye = [xSol1.ye xSol2.ye];
    xSol.ie = [xSol1.ie xSol2.ie];
elseif isfield(xSol1, 'xe')
    xSol.xe = xSol1.xe;
    xSol.ye = xSol1.ye;
    xSol.ie = xSol1.ie;
elseif isfield(xSol2, 'xe')
    xSol.xe = xSol2.xe;
    xSol.ye = xSol2.ye;
    xSol.ie = xSol2.ie;
end