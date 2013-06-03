function val = isoctave()

val = (exist('OCTAVE_VERSION', 'built-in') == 5);
