function test
import mFEM.*



mesh = FEmesh('Space','vector');
mesh.grid('Quad4',0,10,-0.5,0.5,20,5);
mesh.init();

u = [];
x = mesh.get_nodes;
for i = 1:length(x);
    u = [u, exact(x(i,:))];
end

mesh.plot(u,'Deform',true);




function u = exact ( x )

c = 0.5;
E = 1.0e7;
nu = 0.3;
I = 2.0/3*c^3;
G = E/(2*(1+nu));
P = 100.0;
l = 10;
u(1) =-P*x(1)^2*x(2)/(2*E*I) - nu*P*x(2)^3/(6*E*I) + P*x(2)^3/(6*I*G) + (P*l^2/(2*E*I) - P*c^2/(2*I*G))*x(2);
u(2) =nu*P*x(1)*x(2)^2/(2*E*I) + P*x(1)^3/(6*E*I) - P*l^2*x(1)/(2*E*I) + P *l^3/(3*E*I);