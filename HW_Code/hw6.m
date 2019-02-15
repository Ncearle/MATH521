%% MATH 521 - HW6

clear; close all; clc;

% load('video.mat');
% load('kiwi.mat');
load('maple.mat');
% load('pi.mat');
% load('ubc.mat');

msh = struct('P', P, 'E', E, 'T', T);

c = 1;
a = 3;

% f = @(x1,x2) sin(x1.*x2);
% f = @(x1,x2) .75./exp((x1*5).^2.*(x2*5).^2);        %Intersecting Fences
f = @(x1, x2) (sign(0.2-((x1-0.5).^2+(x2-0.5).^2)) + sign(0.2-((x1-0.5).^2/3+(x2-0.5).^2/3)))/3-1;   %Top hat
g = @(x1,x2) 0.*x1;

fh = f(msh.P(1,:),msh.P(2,:));

GammaD = @(x1,x2) true(size(x1));

[A, b, uD, Kbar, Mbar, Pf, PD] = discretiseLinearElasticity(c,a,f,g,GammaD,msh);

uh = A \ b;

ubar = Pf'*uh + PD'*uD;

% Norms
fprintf('L²-norm = %3.2f\n', sqrt(ubar'*Mbar*ubar));
fprintf('H1-norm = %3.2f\n', sqrt(ubar'*Kbar*ubar + ubar'*Mbar*ubar));
fprintf('energy norm = %3.2f\n', sqrt(c*(ubar'*Kbar*ubar) + a*(ubar'*Mbar*ubar)));

% % 
% figure(1);
% pdemesh(P, E, T);
% axis equal
% title('Domain');
% xlabel('x_1'); ylabel('x_2');
% % 
% figure(2);
% spy(Mbar);
% title('Sparsity Pattern of the Mass Matrix');
% 
figure(3);
pdeplot(P, E, T, 'XYData', ubar, 'ZData', ubar, 'mesh','on', 'colormap','hot');
title('Maple leaf; solution');
xlabel('{\itx}_1');
ylabel('{\itx}_2');
zlabel('{\itu^h}({\itx}_1,{\itx}_2)');

figure(4);
pdeplot(P, E, T, 'XYData', fh, 'ZData', fh, 'mesh', 'on', 'colormap', 'hot');
title('Maple leaf: forcing function, "Top Hat"');
xlabel('{\itx}_1');
ylabel('{\itx}_2');
zlabel('{\itf^h}({\itx}_1,{\itx}_2)');
