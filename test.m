function test
import mFEM.*
% 
% Ke = [11,12,13,14; 21, 22, 23, 24; 31, 32, 33, 34; 41, 42, 43, 44];
% 
% dof = (1:4)';
% i = repmat(dof, length(dof),1)
% j = sort(i)
% 
% dof = [4,3,5,6]';
% 
% [dof(i,1),dof(j,1)]


% reshape(Ke, numel(Ke), 1)

% k = ones(4,4);
% M = Matrix(5,5);
% dof1 = [1,2,3,4];
% dof2 = [4,3,5,6];
% M.add_matrix(k,dof1);
% M.add_matrix(2*k, dof2);
% 
% m = M.init; full(m)
% 
% m2 = zeros(size(m));
% m2(dof1,dof1) = k;
% m2(dof2,dof2) = m2(dof2,dof2) + 2*k;
% m2

example6a('N',2,'Method','normal');
example6a('N',2,'Method','alt');