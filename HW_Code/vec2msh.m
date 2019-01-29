function [U] = vec2msh(u, msh)
% undoes the effects of msh2vec

% U = fliplr(reshape(u, msh.N(1)-1, msh.N(2)-1))';  % Grid Conversion

U = reshape(u, msh.N(1)-1, msh.N(2)-1)';    % Matrix conversion

end


