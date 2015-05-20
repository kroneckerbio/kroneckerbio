fastsubs := proc(f, oldvars, newvars)
begin
    if testtype(oldvars,Dom::Matrix()) then
        numberofsubstitutions := nops(oldvars);
    else
        numberofsubstitutions := 1;
    end_if;
    if numberofsubstitutions = 1 then
        f := subs(f, oldvars = newvars);
    else
        for substitutionindex from 1 to numberofsubstitutions do
            f := subs(f, oldvars[substitutionindex] = newvars[substitutionindex]);
        end_for;
    end_if;
    return(f)
end_proc