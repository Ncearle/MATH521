% 3.6  femcode.m
% Example Finite Element Solver, 
% See p303, Section 3.6,
% Computational Science and Engineering 
% Gilbert Strang 
% & thanks to Per-Olof Persson
%
% Solves Poisson's equation
%
%   - u_xx   - u_yy   =   1
%
% L shape domain
%
clear
clc
close all

%% Mesh the L shape domain
% p     an N by 2 matrix listing x,y coordinates of all N=mn nodes
% t     lists 3 node numbers of all triangles in T by 3 matrix
% Start the mesh with some points on the boundary and on the corners.
% We do not have to do this, but it may help to get a good mesh
h=0.5; n=ceil(1/h); h=1/n; z=h:h:1-h; N=length(z); z2=z(1:floor(N/2)); z3=z(floor(N/2)+1:end); 
top1=[z2' ones(length(z2),1)]; top2=[z3' 0.5*ones(length(z3),1)]; bottom=[z' zeros(N,1)]; left=[zeros(N,1) z']; right1=[ones(length(z2),1) z2']; right2=[0.5*ones(length(z3),1) z3'];
corners=[0 0
         0 1
         1 0
         0.5 1
         0.5 0.5
         1 0.5];
sides=[top1
       top2
       bottom
       left
       right1
       right2];
pfix=[corners
      sides];
  
disp('getting mesh...')
h=0.05; % THIS is the h to make smaller for a finer mesh
[p,t]=distmesh2d(@fd_L,@huniform,h,[0,0;1,1],pfix);

% refine if desired, 
% very simple refinement (better to use distmesh with smaller h)
nr=0; % number of refinements 
for i=1:nr
    disp('refining mesh...')
    tic; [p,t]=refine(p,t); toc
end

tol=1e-5*h; ib=abs(fd_L(p))<tol; b=1:size(p,1); b=b(ib); % boundary nodes
ix0=abs(p(:,1)-0)<tol; % y-axis 
ix0(1:2) = 0;          % y-axis, minus the corners [0 0] and [0 1]
% ibminusyaxis=ib & ~ix0; b=1:size(p,1); b=b(ibminusyaxis);

disp('Getting stiffness matrix...')
% [K,F] = assemble(p,t) % K and F for any mesh of triangles: linear phi's
N=size(p,1);T=size(t,1); % number of nodes, number of triangles
% p lists x,y coordinates of N nodes, t lists triangles by 3 node numbers
K=sparse(N,N); % zero matrix in sparse format: zeros(N) would be "dense"
F=zeros(N,1); % load vector F to hold integrals of phi's times load f(x,y)
for e=1:T  % integration over one triangular element at a time
  nodes=t(e,:); % row of t = node numbers of the 3 corners of triangle e
  Pe=[ones(3,1),p(nodes,:)]; % 3 by 3 matrix with rows=[1 xcorner ycorner] 
  Area=abs(det(Pe))/2; % area of triangle e = half of parallelogram area
  C=inv(Pe); % columns of C are coeffs in a+bx+cy to give phi=1,0,0 at nodes
  % now compute 3 by 3 Ke and 3 by 1 Fe for element e
  grad=C(2:3,:);Ke=Area*grad'*grad; % element matrix from slopes b,c in grad
  Fe=Area/3; % integral of phi over triangle is volume of pyramid: f(x,y)=1
  % multiply Fe by f at centroid for load f(x,y): one-point quadrature!
  % centroid would be mean(p(nodes,:)) = average of 3 node coordinates
  K(nodes,nodes)=K(nodes,nodes)+Ke; % add Ke to 9 entries of global K
  F(nodes)=F(nodes)+Fe; % add Fe to 3 components of load vector F
end   % all T element matrices and vectors now assembled into K and F

% [Kb,Fb] = dirichlet(K,F,b) % assembled K was singular! K*ones(N,1)=0
% Implement Dirichlet boundary conditions U(b)=0 at nodes in list b
K(b,:)=0; K(:,b)=0; F(b)=0; % put zeros in boundary rows/columns of K and F 
K(b,b)=speye(length(b),length(b)); % put I into boundary submatrix of K
Kb=K; Fb=F; % Stiffness matrix Kb (sparse format) and load vector Fb
% Solving for the vector U will produce U(b)=0 at boundary nodes
U=Kb\Fb;  % The FEM approximation is U_1 phi_1 + ... + U_N phi_N

% Plot the FEM approximation U(x,y) with values U_1 to U_N at the nodes 
%figure(1)
%trisurf(t,p(:,1),p(:,2),0*p(:,1),U,'edgecolor','k','facecolor','interp'); % would be nice to find of other ways to visualise this...
%view(2), axis equal, colorbar

figure(2)
trimesh(t,p(:,1),p(:,2),U); % hand in this plot
xlabel('x'), ylabel('y'), zlabel('u(x,y)'), az = 70; el = 60; view(az, el);