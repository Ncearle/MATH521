%% MATH 521 HW1
% Nicholas Earle

close all; clc;

% sample function
u = @(x1,x2) cos(2.*pi.*x1).*sin(6.*pi.*x2);

% mesh the rectange [0,1] x [2,3] with 20 / 60 subintervals in x1- / x2-direction, respectively
msh = meshRectangle([0, 1, 2, 3],[20, 60]); 

% evaluate u on msh and draw a surface plot
surf(msh.X1,msh.X2,u(msh.X1,msh.X2));

% axis labels
xlabel('x');
ylabel('y');
zlabel('u');

saveas(gcf, 'u_graph', 'svg');