function example10

import mFEM.*

mesh = FEmesh();
mesh.grid('Line2',-1,1,1);
mesh.init();

sys = System(mesh);
sys.add_matrix('M1', 'N''*N');
sys.add_matrix('M2', 'B''*N');

M1 = sys.assemble('M1'); full(M1)
M2 = sys.assemble('M2'); full(M2)    
    
    
    





