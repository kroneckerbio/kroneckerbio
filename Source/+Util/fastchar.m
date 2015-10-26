function fstr = fastchar(f)

origsize = size(f);

filename = which('fastchar.mu');
read(symengine, filename);
fstr = feval(symengine,'fastchar',f(:));
fstr = reshape(eval(fstr),origsize);

end