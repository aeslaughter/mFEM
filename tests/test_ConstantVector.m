function T = test_ConstantVector(varargin)

T = mFEM.Test('Name','ConstantVector',varargin{:});

mesh = mFEM.Mesh('Space','vector');
mesh.grid('Quad4',0,1,0,1,2,2);
mesh.addBoundary(1,'left');
mesh.update();

kern = mFEM.kernels.ConstantVector(mesh, 'f', 0);
f = kern.get();
T.compare(f.init(), zeros(mesh.n_dof,1), 'Scalar, zero intilization');

kern = mFEM.kernels.ConstantVector(mesh, 'f', 11);
f = kern.get();
T.compare(f.init(), 11*ones(mesh.n_dof,1), 'Scalar, scalar intilization');

kern = mFEM.kernels.ConstantVector(mesh, 'f', 7*ones(mesh.n_dof,1));
f = kern.get();
T.compare(f.init(), 7*ones(mesh.n_dof,1), 'Complete vector intilization');

kern = mFEM.kernels.ConstantVector(mesh, 'f', 1, 'Tag', 1);
f = kern.get();
fex = zeros(mesh.n_dof,1);
fex([1:2, 7:8, 13:14]) = 1;
T.compare(f.init(), fex, 'Partial dof, scalar input');

fin = 8*ones(6,1);
kern = mFEM.kernels.ConstantVector(mesh, 'f', fin, 'Tag', 1);
f = kern.get();
T.compare(f.init(), 8*fex, 'Partial dof, vector input');

input = [1,0.5,1,0.5,0,0,1,0.5,1,1,0,0,0,0,0,0,0,0];
kern = mFEM.kernels.ConstantVector(mesh, 'f', input);
elem = mesh.getElements(1,'-gather');
pv = kern.pointValue(elem, [0.5,0.5]);
T.compare(pv,[0;0.3125],'Point value');

pg = kern.pointGradient(elem, [0.25,0.25]);
T.compare(pg,[0,0;0,1],'Point gradient');


