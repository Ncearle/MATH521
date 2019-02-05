%% MATH 521 - HW5

clear; close all; clc;

% load('kiwi.mat');
load('video.mat');

msh = struct('P', P, 'E', E, 'T', T);

Mbar = discretiseLinearElasticity(msh);

figure(1);
pdemesh(P, E, T);
axis equal
title('Domain');
xlabel('x_1'); ylabel('x_2');

figure(2);
spy(Mbar);
title('Sparsity Pattern of the Mass Matrix');

% u = @(x1, x2) 1 - abs(x1+x2) - abs(x2-x1);   %Pyramid function
% u = u(P(1,:), P(2,:));
b = @circleb1;
c = 1;
a = 5;
f = 10;
u = assempde(b, P, E, T, c, a, f);

figure(3);
pdeplot(P, E, T, 'XYData', u, 'ZData', u)
