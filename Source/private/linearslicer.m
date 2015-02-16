function sliced = linearslicer(dim, varargin)


box = reshape(1:prod(dim), dim);
sliced = vec(box(varargin{:}));
