function [A] = itrimatrand (m, n, radexp, density, midbounds)
  A = imatrand(m, n, radexp, density);
  for i = 1:(n - 2)
    A((i + 2):m, i) = 0;
    A(1:i, i + 2) = 0;
  end
endfunction