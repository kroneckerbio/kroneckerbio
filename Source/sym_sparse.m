function S = sym_sparse(varargin)
% Redefined sparse() function that constructs full matrices using the
% sparse syntax if input arguments are symbolic.  Note that this function
% needs to be accessible on the path (i.e., not in the private folder) for
% analytic model function handles to be able to find it.

switch nargin
    case 1
        % One input argument normally converts to sparse matrix. Here, just
        % return the input if it is symbolic, and convert to sparse
        % otherwise.
        if isa(varargin{1}, 'sym')
            S = varargin{1};
        else
            S = builtin('sparse', varargin{:});
        end
    case 2
        % Two input arguments initializes an empty sparse matrix.
        % Here, initialize with zeros.
        S = builtin('sparse', varargin{:});
    case 3
        % Three input arguments i,j,v generates a matrix
        % S(i(k),j(k)) = v(k). If v is symbolic, create a full symbolic
        % matrix. Otherwise just call sparse.
        [i,j,v] = deal(varargin{:});
        if isa(v,'sym')
            imax = max(i(:));
            jmax = max(j(:));
            S = zeros(imax, jmax, 'like', v);
            for k = 1:numel(i)
                S(i(k),j(k)) = S(i(k),j(k)) + v(k);
            end
        else
            S = builtin('sparse', varargin{:});
        end
    case 5
        % Five input arguments i,j,v,m,n generates a matrix
        % S(i(k),j(k)) = v(k), with S of size m-by-n. If v is symbolic,
        % create a full symbolic matrix. Otherwise just call sparse.
        [i,j,v,m,n] = deal(varargin{:});
        if isempty(i) || isempty(j)
            % If no nonzero elements are provided, just initialize an empty
            % S of the same type as v
            if isa(v, 'sym')
                S = zeros(m, n, 'like', v);
            else
                S = sparse(varargin{:});
            end
        elseif isa(v, 'sym')
            imax = max(i(:));
            jmax = max(j(:));
            assert(imax <= m, 'Index exceeds array bounds.')
            assert(jmax <= n, 'Index exceeds array bounds.')
            S = zeros(m, n, 'like', v);
            for k = 1:numel(i)
                S(i(k),j(k)) = S(i(k),j(k)) + v(k);
            end
        else
            S = builtin('sparse', varargin{:});
        end
    case 6
        % Six input arguments is the same as 5, but allocates space for
        % more nonzero elements, which has no use for full matrices.
        % Just use the code for five input arguments if symbolic.
        [i,j,v,m,n,~] = deal(varargin{:});
        if isa(v, 'sym')
            S = sparse(i,j,v,m,n);
        else
            S = builtin('sparse', varargin{:});
        end
    otherwise
        error('Unsupported number of input arguments.')
end

end