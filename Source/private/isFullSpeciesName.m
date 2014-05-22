function val = isFullSpeciesName(name)

if any(name == '.')
    val = true;
else
    val = false;
end