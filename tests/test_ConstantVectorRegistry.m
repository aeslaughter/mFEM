function T = test_ConstantVectorRegistry(varargin)

T = mFEM.Test('Name','ConstantVectorRegistry',varargin{:});

mesh = mFEM.Mesh('Space','vector');
mesh.grid('Quad4',0,1,0,1,2,2);
mesh.addBoundary(1,'left');
mesh.update();

reg = mFEM.registry.ConstantVectorRegistry(mesh);
reg.add('f0', 0);
f = reg.get('f0');
T.compare(f.init(), zeros(mesh.n_dof,1), 'Scalar, zero intilization');

reg.add('f1', 11);
f = reg.get('f1');
T.compare(f.init(), 11*ones(mesh.n_dof,1), 'Scalar, scalar intilization');
 
reg.add('f2', 7*ones(mesh.n_dof,1));
f = reg.get('f2');
T.compare(f.init(), 7*ones(mesh.n_dof,1), 'Complete vector intilization');

reg.add('f3', 1, 'Tag', 1);
f = reg.get('f3');
fex = zeros(mesh.n_dof,1);
fex([1:2, 7:8, 13:14]) = 1;
T.compare(f.init(), fex, 'Partial dof, scalar input');

fin = 8*ones(6,1);
reg.add('f4', fin, 'Tag', 1);
f = reg.get('f4');
T.compare(f.init(), 8*fex, 'Partial dof, vector input');

reg.add('f4',1,'-add');
f = reg.get('f4');
T.compare(f.init(), 8*fex+1, 'Partial dof, vector input');