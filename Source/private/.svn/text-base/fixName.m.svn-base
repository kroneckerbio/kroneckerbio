function name = fixName(name)

if isempty(name)
    name = '';
end

assert(isempty(regexp(name, '\.|,|"', 'once')), 'KroneckerBio:fixName:ComponentName', 'Names of components of Kronecker Bio models cannot contain "." or "," or """ like "%s"', name)
