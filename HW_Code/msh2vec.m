function [u] = msh2vec(U, msh)
% msh2vec reorders an array of any size to a column vector

% u = reshape(fliplr(U'), (msh.N(1)-1)*(msh.N(2)-1), 1);    % Grid Conversion

u = reshape(U', (msh.N(1)-1)*(msh.N(2)-1), 1);  % Matrix Conversion

end 