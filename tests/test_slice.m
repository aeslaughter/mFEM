function test_slice

mesh = mFEM.FEmesh();
mesh.grid(0,1,0,1,0,1,6,6,6,'Element','Hex8');
mesh.init();

T_exact = @(x,y,z,t) exp(-t)*sin(pi*x).*sin(pi*y).*sin(pi*z);

nodes = mesh.getNodes();
x = nodes(:,1);
y = nodes(:,2);
z = nodes(:,3);
T = T_exact(x,y,z,0);


mesh.plot(T,'slice',{'x',[0.5,1],'y',[0.5,1],'z',[0,0.5]},'-new','-ShowElements')