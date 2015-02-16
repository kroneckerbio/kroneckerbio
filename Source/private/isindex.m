function tf = isindex(index, length)
if islogical(index)
    if nargin < 2
        tf = true;
    else
        tf = ~any(index(length:end));
    end
elseif isnumeric(index)
    if nargin < 2
        tf = all(index >= 1 & floor(index) == index);
    else
        tf = all(index >= 1 & index <= length & floor(index) == index);
    end
else
    tf = false;
end
