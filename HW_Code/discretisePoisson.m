function [A, b] = discretisePoisson(f, g, msh)

e = ones((msh.N(1)-1)*(msh.N(2)-1), 1);

L1 = spdiags([-e 2*e -e], [-1 0 1], (msh.N(1)-1)*(msh.N(2)-1), (msh.N(1)-1)*(msh.N(2)-1));
L2 = spdiags([-e 2*e -e], [-(msh.N(1)-1) 0 (msh.N(1)-1)], (msh.N(1)-1)*(msh.N(2)-1), (msh.N(1)-1)*(msh.N(2)-1));

A = 1/msh.h(1) * L1 + 1/msh.h(2) * L2;

% Source Term Matrix
s = f(msh.X1(2:end-1, 2:end-1), msh.X2(2:end-1, 2:end-1));

% Boundary Values
bE = g(msh.X1(2:end-1,1), msh.X2(2:end-1,1));
bW = g(msh.X1(2:end-1,end), msh.X2(2:end-1,end));
bN = g(msh.X1(1,2:end-1), msh.X2(1,2:end-1));
bS = g(msh.X1(end,2:end-1), msh.X2(end,2:end-1));

s(:,1) = s(:,1) + 1/msh.h(2)*bE;
s(:,end) = s(:,end) + 1/msh.h(2)*bW;
s(end,:) = s(end,:) + 1/msh.h(1)*bN;   % North boundary applied to bottom row of matrix (y = 3)
s(1,:) = s(1,:) + 1/msh.h(1)*bS;       % South boundary applied to top row of matrix (y = 2)

b = msh2vec(s, msh);