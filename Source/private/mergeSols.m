function xSol = mergeSols(varargin)
% mergeSols Merge multiple ode solutions structures in a way that preserves the
% ability to use deval on them

% (c) 2015 David R Hagen, Joshua F Apgar, Jared E Toettcher, & Bruce Tidor
% This work is released under the MIT license.

xSol = varargin{2};

if isempty(varargin{1})
    return
end

% Append xSol2 onto the end of xSol1
xSol.x = cellfun(@(sol){sol.x}, varargin);
xSol.x = [xSol.x{:}];
xSol.y = cellfun(@(sol){sol.y}, varargin);
xSol.y = [xSol.y{:}];

switch varargin{1}.solver
    case 'ode15s'
        for i = 1:numel(varargin)-1
            xSol1 = varargin{i};
            xSol2 = varargin{i+1};
            size1 = size(xSol1.idata.dif3d, 2);
            size2 = size(xSol2.idata.dif3d, 2);
            if(size1 == size2)
                % No change
%                 dif3d1 = xSol1.idata.dif3d;
%                 dif3d2 = xSol2.idata.dif3d;
            elseif(size1 < size2)
                xSol1.idata.dif3d = padarray_post(xSol1.idata.dif3d, [0, size2-size1, 0]);
%                 dif3d2 = xSol2.idata.dif3d;
            else % size1 bigger
%                 dif3d1 = xSol1.idata.dif3d{1};
                xSol2.idata.dif3d = padarray_post(xSol2.idata.dif3d, [0, size1-size2, 0]);
            end
            varargin{i} = xSol1;
            varargin{i+1} = xSol2;
        end
        temp = cellfun(@(sol){sol.idata.kvec}, varargin); %[xSol1.idata.kvec, xSol2.idata.kvec];
        xSol.idata.kvec = [temp{:}];
        temp = cellfun(@(sol){sol.idata.dif3d}, varargin);
        xSol.idata.dif3d = cat(3, temp{:});
    case 'ode23s'
        temp = cellfun(@(sol){sol.idata.k1}, varargin);
        xSol.idata.k1 = [temp{:}];
        temp = cellfun(@(sol){sol.idata.k2}, varargin);
        xSol.idata.k2 = [temp{:}];
%         xSol.idata.k1 = [xSol1.idata.k1, xSol2.idata.k1];
%         xSol.idata.k2 = [xSol1.idata.k2, xSol2.idata.k2];
    case 'ode45'
        temp = cellfun(@(sol){sol.idata.f3d}, varargin);
        xSol.idata.f3d = cat(3, temp{:});
%         xSol.idata.f3d = cat(3, xSol1.idata.f3d, xSol2.idata.f3d);
end

% If there were events, append them
[xe, ye, ie] = deal(cell(nargin,1));
any_event = false;
for i = 1:nargin
    if isfield(varargin{i}, 'xe')
        any_event = true;
        xe{i} = varargin{i}.xe;
        ye{i} = varargin{i}.ye;
        ie{i} = varargin{i}.ie;
    end
end
if any_event
    xSol.xe = [xe{:}];
    xSol.ye = [ye{:}];
    xSol.ie = [ie{:}];
end

% if isfield(xSol1, 'xe') && isfield(xSol2, 'xe')
%     xSol.xe = [xSol1.xe, xSol2.xe];
%     xSol.ye = [xSol1.ye, xSol2.ye];
%     xSol.ie = [xSol1.ie, xSol2.ie];
% elseif isfield(xSol1, 'xe')
%     xSol.xe = xSol1.xe;
%     xSol.ye = xSol1.ye;
%     xSol.ie = xSol1.ie;
% elseif isfield(xSol2, 'xe')
%     xSol.xe = xSol2.xe;
%     xSol.ye = xSol2.ye;
%     xSol.ie = xSol2.ie;
% end
end

function B = padarray_post(A, pad_size)
% padarray_post Pad an array to a given size with zeros.
%
%   B = padarray_post(A, pad_size)
% 
% This function implements the post method of the full padarray function
% found in the Image Processing Toolbox. It pads the input array A with
% zeros to the right and bottom until the final array size is B.

n_dims = numel(pad_size);

ind   = cell(1,n_dims);
size_B = zeros(1,n_dims);
for k = 1:n_dims
    M = size(A,k);
    ind{k}   = 1:M;
    size_B(k) = M + pad_size(k);
end

B         = zeros(size_B);
B(ind{:}) = A;
end
