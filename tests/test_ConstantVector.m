function T = test_ConstantVector

T = mFEM.Test();

mesh = mFEM.FEmesh('Space','vector');
mesh.grid(0,1,0,1,2,2,'Element','Quad4');
mesh.init();
mesh.addBoundary(1,'left');

kern = mFEM.kernels.ConstantVector(mesh, 'f', 0);
f = kern.get();
T.compare(f.init(), zeros(mesh.n_dof,1), 'Scalar, zero intilization');

kern = mFEM.kernels.ConstantVector(mesh, 'f', 11);
f = kern.get();
T.compare(f.init(), 11*ones(mesh.n_dof,1), 'Scalar, scalar intilization');

kern = mFEM.kernels.ConstantVector(mesh, 'f', 7*ones(mesh.n_dof,1));
f = kern.get();
T.compare(f.init(), 7*ones(mesh.n_dof,1), 'Complete vector intilization');

kern = mFEM.kernels.ConstantVector(mesh, 'f', 1, 'Boundary', 1);
f = kern.get();
fex = zeros(mesh.n_dof,1);
fex([1:2, 7:8, 11:12]) = 1;
T.compare(f.init(), fex, 'Partial dof, scalar input');

fin = 8*ones(6,1);
kern = mFEM.kernels.ConstantVector(mesh, 'f', fin, 'Boundary', 1);
f = kern.get();
T.compare(f.init(), 8*fex, 'Partial dof, vector input');

input = [0,0,0,0,0,0.5,0,0.5,0,1,0,1,0,0,0,0.5,0,1];
kern = mFEM.kernels.ConstantVector(mesh, 'f', input);
elem = mesh.element(1);
pv = kern.pointValue(elem, [0.25,0.25]);
T.compare(pv,[0;0.3125],'Point value');

pg = kern.pointGradient(elem, [0.25,0.25]);
T.compare(pg,[0,0;0,1],'Point gradient');


