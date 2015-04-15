function val = isoctave()

val = (exist('OCTAVE_VERSION', 'builtin') == 5);
