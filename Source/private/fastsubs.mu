fastsubs := proc(f, oldvars, newvars)
begin
    [numberofsubstitutions,temporaryunusedvariable] := linalg::matdim(oldvars);
    for substitutionindex from 1 to numberofsubstitutions do
        f := subs(f, oldvars[substitutionindex] = newvars[substitutionindex]);
    end_for;
    return(f);
end_proc