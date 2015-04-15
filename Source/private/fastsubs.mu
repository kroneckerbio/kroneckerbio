fastsubs := proc(f, oldvars, newvars)
begin
    [nsubs,temp] := linalg::matdim(oldvars);
    for ii from 1 to nsubs do
        f := subs(f, oldvars[ii] = newvars[ii]);
    end_for;
    return(f);
end_proc