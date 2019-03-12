%% MATH 521 - HW3

close all; clear; clc;

msh = meshRectangle([0, 1, 2, 3],[20, 60]); 

% Source term
f = @(x1, x2) 40*pi^2*cos(2*pi*x1).*cos(6*pi*x2);

% Boundaries
g = @(x1, x2) cos(2*pi*x1).*cos(6*pi*x2);

[A, b] = discretisePoisson(f, g, msh);
u = vec2msh(A\b, msh); 

figure();
surf(u);
xlabel('X');
ylabel('Y');
zlabel('u');
title('u');

figure();
spy(A);
title('Sparsity of A matrix');


