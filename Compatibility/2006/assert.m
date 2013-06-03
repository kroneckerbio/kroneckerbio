function  assert(condition, identifier, message, varargin)
% ASSERT Tests if a condition is true and throws an error if it is false.
%    ASSERT(condition, identifier, message, ...) reports an error if
%    condition is false.  The message may be a sprintf format string with
%    variable arguments following.  This function is intended to replace
%    code such as:
%    
%    if (somethingbadhappend)
%       error("something bad happend in my code!");
%    end
%
% ARGUMENTS:
%    condition : A boolean value.  If it is false the error is reported
%    identifier: The systematic name of the error.
%    message   : The human readable message
%
% See also error, sprintf
% 
% This work is licensed under the Creative Commons Attribution-Noncommercial-No Derivative 
% Works 3.0 United States License. To view a copy of this license, visit 
% http://creativecommons.org/licenses/by-nc-nd/3.0/us/ or send a letter to Creative Commons,
% 171 Second Street, Suite 300, San Francisco, California, 94105, USA.



if(~condition)
    stack   = dbstack;
    stack   = stack(2:(end));
    
    errText = sprintf('Error using ==> %s (%s line %d)', stack(end).name, stack(end).file, stack(end).line);
    msgText = sprintf(message, varargin{:});
    
    msgStruct.message    = sprintf('%s\n%s', errText, msgText);
    msgStruct.identifier = identifier;
    msgStruct.stack      = stack;
    rethrow(msgStruct);
end