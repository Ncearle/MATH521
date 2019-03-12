function [msh] = meshRectangle(x,N)
%   meshRectangle creates a 2-dimensional rectangular domain given
%   the coordinates of the rectangle and the number of nodes in each
%   dimension.

xv = linspace(x(1), x(2), N(1)+1);
yv = linspace(x(3), x(4), N(2)+1);

[X1, X2] = meshgrid(xv, yv);

h1 = (x(2) - x(1))/N(1);
h2 = (x(4) - x(3))/N(2);
h = [h1, h2];

msh = struct('X1', X1, 'X2', X2, 'N', N, 'h', h);

end
