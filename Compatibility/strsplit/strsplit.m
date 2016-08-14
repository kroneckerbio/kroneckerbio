function [ C ] = strsplit( str, delimiter )
%STRSPLIT Compatibility version of strsplit for MATLAB versions older than
%R2013a
%   Splits a string at a specified delimiter. Only implements the
%   functionality needed for Kronecker, which is the (str, delimiter) input
%   argument case.

delimiter = sprintf(delimiter); % Accounts for escaped characters like \n
di = strfind(str, delimiter);
ncells = length(di)+1;
C = cell(1,ncells);
for ii = 1:ncells
    if ii == 1
        starti = 1;
    else
        starti = di(ii-1)+1;
    end
    if ii == ncells
        endi = length(str);
    else
        endi = di(ii)-1;
    end
    C{ii} = str(starti:endi);
end

end

