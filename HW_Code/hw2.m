%% MATH 521 - HW2
% Nicholas Earle

close all; clear; clc;

run('hw1.m');

% Not sure if you wanted inline or separate functions considering they can
% be only a line of code each

% msh2vec = @(U, msh) reshape(fliplr(U'), (msh.N(1)-1)*(msh.N(2)-1), 1);
% vec2msh = @(u, msh) fliplr(reshape(u, msh.N(1)-1, msh.N(2)-1))';

U = u(msh.X1(2:end-1, 2:end-1), msh.X2(2:end-1, 2:end-1));
V = vec2msh(msh2vec(U, msh), msh);

if (U == V)
    disp('U and V are the same');
end