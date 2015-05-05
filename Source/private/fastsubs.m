function out = fastsubs(f, old, new)

filename = which('fastsubs.mu');
read(symengine, filename);
out = feval(symengine,'fastsubs',f,old,new);

end