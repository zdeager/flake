% mat2gray.m
% Converts the matrix M to an intensity image I
% 
% zde

function I = mat2gray (M)
  Mmin = min (min (M));
  Mmax = max (max (M));
  I = (M < Mmin) .* 0;
  I = I + (M >= Mmin & M < Mmax) .* (1 / (Mmax - Mmin) * (M - Mmin));
  I = I + (M >= Mmax);
end
