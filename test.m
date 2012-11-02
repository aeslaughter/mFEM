function test
import mFEM.*



mesh = FEmesh('Space','vector');
mesh.grid('Quad4',0,1,0,1,2,2);
mesh.init();


mesh.add_boundary('right',1)
mesh

