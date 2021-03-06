function [A, b, uD, Kbar, Mbar, Pf, PD] = discretiseLinearElasticity(c,a,f,g,GammaD,msh)
%DISCRETISELINEARELASTICITY Applies linear finite elements on the
%triangulation defined in the structure msh to discretise the equation
% -cΔu + au = f
% with Dirichlet boundary conditions u = g on GammaD.
%   Input:
%       c: a positive real number
%       a: a positive real number
%       f: a function handle for the inhomogeneity f(x1,x2)
%       g: a function handle for the boundary values g(x1,x2)
%       GammaD: a function handle to a function of x1 and x2, which returns
%       true, if Dirichlet conditions are imposed at the boundary point
%       (x1,x2) and false otherwise
%       msh: a mesh structure with fields P, E, T
%   Output:
%       A: a sparse array corresponding to the discrete elliptic operator
%       b: a column vector for the discretised right hand side and boundary
%       terms
%       uD: a column vector with the prescribed values of the solution on
%       the Dirichlet boundary
%       Kbar: the sparse stiffness matrix, including all boundary points
%       Mbar: the sparse mass matrix, including all boundary points
%       Pf: a sparse matrix that projects onto all free nodes
%       PD: a sparse matrix that projects onto all Dirichlet nodes

% -------------------------------------------------------------------------
% COMPUTE EXTRA MESH DATA
% -------------------------------------------------------------------------

% edge vectors
msh.D32 = msh.P(:, msh.T(3,:)) - msh.P(:, msh.T(2,:));
msh.D13 = msh.P(:, msh.T(1,:)) - msh.P(:, msh.T(3,:));
msh.D21 = msh.P(:, msh.T(2,:)) - msh.P(:, msh.T(1,:));

% row vector of triangle areas
msh.A = (msh.D13(1,:).*msh.D21(2,:) - msh.D21(1,:).*msh.D13(2,:))./2;
% = [|T_1|, |T_2|, |T_3|, ..., |T_nt|]

% number of points
msh.np = size(msh.P, 2);

% number of triangles
msh.nt = size(msh.T, 2);

% -------------------------------------------------------------------------
% ASSEMBLE THE DISCRETE SYSTEM
% -------------------------------------------------------------------------

% complete system including all boundary points
Kbar = assembleStiffness(msh);
Mbar = assembleMass(msh);
fbar = assembleLoad(f, msh);

% -------------------------------------------------------------------------
% ELIMINATE DIRICHLET BOUNDARY CONDITIONS
% -------------------------------------------------------------------------

[Pf, PD, uD] = eliminateDirichletBC(g, GammaD, msh);

% contributions from points with prescribed Dirichlet data
kD = Pf*Kbar*PD'*uD;
mD = Pf*Mbar*PD'*uD;

% reduced system for the free points only
K = Pf*Kbar*Pf';
M = Pf*Mbar*Pf';
A = c*K + a*M;
b = Pf*fbar - (c*kD + a*mD);
end

function Kbar = assembleStiffness(msh)
%ASSEMBLESTIFFNESS Assembles the complete stiffness matrix including all
%boundary points, discretised with linear finite elements

% element stiffness matrix
k11 = sum(msh.D32.^2)./(4*msh.A);
% = [k11 for T_1, k11 for T_2, k11 for T_3, ..., k11 for T_nT]
k12 = sum(msh.D32.*msh.D13)./(4*msh.A);
k13 = sum(msh.D32.*msh.D21)./(4*msh.A);
k22 = sum(msh.D13.^2)./(4*msh.A);
k23 = sum(msh.D13.*msh.D21)./(4*msh.A);
k33 = sum(msh.D21.^2)./(4*msh.A);

% the following three matrices are 3x3 block matrices, where each block is
% a row vector of length nt

i = [msh.T(1,:) msh.T(1,:) msh.T(1,:);
    msh.T(2,:) msh.T(2,:) msh.T(2,:);
    msh.T(3,:) msh.T(3,:) msh.T(3,:)]; % row indices

j = [msh.T(1,:) msh.T(2,:) msh.T(3,:);
    msh.T(1,:) msh.T(2,:) msh.T(3,:);
    msh.T(1,:) msh.T(2,:) msh.T(3,:)]; % column indices

kij = [k11 k12 k13;
    k12 k22 k23;
    k13 k23 k33]; % entries of the element stiffness matrices

% Note that with the repmat command, the previous lines could be made more
% compact and efficient, but the more verbose commands here are probably
% more illustrative.

% complete stiffness matrix
Kbar = sparse(i, j, kij, msh.np, msh.np);
end


function Mbar = assembleMass(msh)
%ASSEMBLEMASS Assembles the complete mass matrix including all boundary
%points, discretised with linear finite elements

% element mass matrix
mDiag = msh.A/6;
mOffDiag = msh.A/12;

% the following three matrices are 3x3 block matrices, where each block is
% a row vector of length nt

i = [msh.T(1,:) msh.T(1,:) msh.T(1,:);
    msh.T(2,:) msh.T(2,:) msh.T(2,:);
    msh.T(3,:) msh.T(3,:) msh.T(3,:)]; % row indices

j = [msh.T(1,:) msh.T(2,:) msh.T(3,:);
    msh.T(1,:) msh.T(2,:) msh.T(3,:);
    msh.T(1,:) msh.T(2,:) msh.T(3,:)]; % column indices

mij = [mDiag mOffDiag mOffDiag;
    mOffDiag mDiag mOffDiag;
    mOffDiag mOffDiag mDiag]; % entries of the element mass matrices

% Note that with the repmat command, the previous lines could be made more
% compact and efficient, but the more verbose commands here are probably
% more illustrative.

% complete mass matrix
Mbar = sparse(i, j, mij, msh.np, msh.np);
end

function qf = assembleLoad(f, msh)
%ASSEMBLELOAD Assembles the complete load vector including all boundary
%points, discretised with linear finite elements

% triangle midpoints
xm = (msh.P(:, msh.T(1,:)) + msh.P(:, msh.T(2,:)) + msh.P(:, msh.T(3,:)))./3;

% f(xm)
fm = f(xm(1,:)', xm(2,:)');

% quadrature operator (midpoint rule) on one element
i = [msh.T(1,:);
   msh.T(2,:);
   msh.T(3,:)]; % row indices

j = [1:msh.nt;
   1:msh.nt;
   1:msh.nt]; % column indices

qij = [msh.A./3;
   msh.A./3;
   msh.A./3]; % midpoint rule on each triangle

% full quadrature operator
Q = sparse(i, j, qij, msh.np, msh.nt);

% complete load vector
qf = Q*fm;
end

function [Pf, PD, uD] = eliminateDirichletBC(g, GammaD, msh)
%ELIMINATEDIRICHLETBC Calculates the vector uD = g(x1, x2) for all boundary
%points (x1, x2) on the Dirichlet section GammaD of the boundary. Also
%assembles the projection matrices Pf and PD.

% indices of all boundary points
indGamma = unique(msh.E(1:2,:));

% coordinates of all boundary points
Gamma = msh.P(:,indGamma);

% points on that part of the boundary where Dirichlet conditions are
% imposed
indD = indGamma(GammaD(Gamma(1,:), Gamma(2,:))); % indices of Dirichlet points
nD = length(indD); % number of Dirichlet points

% all other points of the mesh
indf = setdiff(1:msh.np, indD); % indices of free points
nf = msh.np - nD; % number of free points

% projection onto the Dirichlet points
PD = sparse(1:nD, indD, ones(1,nD), nD, msh.np);

% projection onto the free points
Pf = sparse(1:nf, indf, ones(1,nf), nf, msh.np);

% boundary values of u on the Dirichlet points
uD = g(msh.P(1,indD), msh.P(2,indD))';
end
