function uid = genUID
% Generate random unique ID from UUID suitable for variable name
uid = ['bt', strrep(char(java.util.UUID.randomUUID()), '-', '_')];
end