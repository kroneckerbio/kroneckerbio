function read( symengine, file )
%READ Read MuPAD program file into symbolic engine

file = regexprep(file, '\\', '\\\');
feval(symengine, 'read',[' "' file '" '], 'Plain');

end

