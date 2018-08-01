% Wetted perimeter correction factor for asymmetric triangular cross
% section
%
% July 10, 2018
% Eventually, want to extend this to a general exponential cross section
% Just need to find the shear transformation that makes the shape symmetric
% Also, try out a wetted-perimeter conserved shape, possibly
%
% INPUTS
% Width, depth of asymmetric triangle
% One side slope
%
% OUTPUTS
% Ratio of wetted perimeter of symmetric triangle to the asymmetric
% triangle with the same width and depth (area-conserved)

s1 = 1/2;
s2 = 1/6;
sp = 1/4; % side slope of isosceles triangle

num = sqrt(1+(1/s1)^2) + sqrt(1+(1/s2)^2);
denom = 2*sqrt(1+(1/sp)^2);

k = num/denom; % correction factor