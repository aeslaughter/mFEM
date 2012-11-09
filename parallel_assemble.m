clear;
import mFEM.*

mesh = FEmesh();
mesh.grid('Quad4',0,1,0,1,100,100);
mesh.init();

sys = System(mesh);
sys.add_matrix('M','N''*N');
sys.assemble('M');


nel = mesh.n_elements;
% n_dof = mesh.n_dof;
% n_dof_elem = mesh.element(1).n_dof;
% nI = nel * n_dof_elem^2;

tic;
spmd
   %M = sparse(n_dof,n_dof,codistributor1d());
   elem_index = codistributed(1:nel); 
   local_index = getLocalPart(elem_index);
   
   I = []; 
   J = [];
   Mij = [];
   for e = 1:length(local_index);
       elem = mesh.element(local_index(e));
%        disp(['Lab: ', num2str(labindex),'; id: ',num2str(elem.id)]);

       [qp,w] = elem.quad.rules();
       Me = zeros(elem.n_dof);
       for i = 1:length(qp);
           for j = 1:length(qp);
               Me = Me + elem.shape(qp(i),qp(j))'*elem.shape(qp(i),qp(j))*elem.detJ(qp(i),qp(j));
           end
       end
       
       dof = elem.get_dof()';
       ii = repmat((1:length(dof))', length(dof), 1);
       jj = sort(ii);
       I = [I;dof(ii)'];
       J = [J;dof(jj)'];
       Mij = [Mij; reshape(Me,numel(Me),1)];
   end
   
   labs = 1:numlabs;
   other = labs(labs ~= labindex);
   labSend(length(I), other);
   
   data = zeros(1,numlabs);
   for i = 1:length(other)
        data(other(i)) = labReceive(other(i));
   end
    data(labindex) = length(I);

   codist = codistributor1d(1,data,[sum(data),1]);
   I_dist = codistributed.build(I, codist);
   J_dist = codistributed.build(J, codist);
   Mij_dist = codistributed.build(Mij, codist);
end

M = sparse(I_dist,J_dist,Mij_dist);
full(M);
toc;
