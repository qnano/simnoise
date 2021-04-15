function [TriangleFilter] = triangle_apodization(L, deltax, lambda_em, NA)
% compute triangular apodization filter with cutoff frequency = 4NA/lambda

F = -1/2/deltax:1/L:1/2/deltax-1/L;
[Fx,Fy] = meshgrid(F,F);

fnorm = sqrt(Fx.^2 + Fy.^2)*lambda_em./2./NA;
f_cutoff = 2;
TriangleFilter = 1-fnorm/f_cutoff;
TriangleFilter(TriangleFilter<0) = 0;

end

