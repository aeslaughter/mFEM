function example11b

% Set the initial conditions for the signed distance function
d_0 = @(x,y) (x-0.5).^2 + (y-0.75).^2 - 0.15^2;

% Create the mesh
mesh = mFEM.FEmesh('Element','Tri3');
mesh.grid(0,1,0,1,1,1);
mesh.init();
%mesh.plot();

% Create system
sys = mFEM.System(mesh, '-time');
dx = mesh.element(1).hmin();
sys.add_constant('dt', 0.01, 'h', dx);
sys.add_function('v',@(x,t) velocity(x,t));
% sys.add_matrix('K','N''*N - dt/2*N''*v''*B + N*h/(2*norm(v))*B''*v + dt/2*h/(2*norm(v))*B''*v*v''*B');

%sys.add_matrix('K','N''*N');
sys.add_matrix('K','-dt/2*N''*(v''*B)');

% sys.add_matrix('K','N*h/(2*norm(v))*B''*v');
% sys.add_matrix('K','dt/2*h/(2*norm(v))*B''*v*v''*B');

[x,y] = mesh.get_nodes;
sys.add_vector('phi0', d_0(x,y));
% sys.add_vector('f', 'phi0*N'' + dt/2*B''*v + phi0*h/(2*norm(v))*B''*v + dt/2*h/(2*norm(v))*B''*v*v''*grad(phi0)');
sys.add_vector('f', 'phi0*N''');
% sys.add_vector('f', 'dt/2*B''*v');
% sys.add_vector('f', 'phi0*h/(2*norm(v))*B''*v');
% sys.add_vector('f','dt/2*h/(2*norm(v))*B''*v*v''*grad(phi0)');

K = sys.assemble('K');
%ix = [1,2,9,13,4,3,10,14,6,5,11,15,8,7,12,16];
ix = [1,2,4,3];
full(K(ix,ix))
%full(K)

%f = sys.assemble('f');


%dt = sys.get('dt');
% [X,Y] = meshgrid(x,y);
% phi = d_0(x,y);
 
 return;

 
 mesh.plot(phi','-new'); hold on;
 title('t = 0');

% Z = griddata(x,y,phi,X,Y);
% [~,chandle] = contour(X,Y,Z,[0,0],'-w','LineWidth',2);
% title('t=0');
drawnow;


solver = LinearSolver(sys);

for i = dt:dt:10*dt;
    sys.time = sys.time + dt;
    phi = solver.solve();
    sys.add_vector('phi0', phi, '-overwrite');

    mesh.plot(phi); hold on;
%     Z = griddata(x,y,phi,X,Y);
%     delete(chandle);
%     [~,chandle] = contour(X,Y,Z,[0,0],'-w','LineWidth',2);   
    drawnow;
    title(['t = ',sys.time]);
end

end

function v = velocity(x,t)
   v(1,1) = cos(pi*t/8).*sin(2*pi*x(2)).*sin(pi*x(1)).^2; 
   v(2,1) = -cos(pi*t/8).*sin(2*pi*x(1)).*sin(pi*x(2)).^2; 
end

