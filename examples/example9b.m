% MAE4700/5700 HW8, Prob. 2
function example9b

% Import the mFEM library
import mFEM.*;
  
% Create a FEmesh object, add the single element, and initialize it
mesh = FEmesh('Space','vector','Element','Quad4');
mesh.grid(0,2,0,1,1,1);
mesh.init();

% Label the boundaries
mesh.add_boundary(1, 'bottom','left'); 
mesh.add_boundary(2, 'right','top');
mesh.add_boundary(3, {'x==0','y==0'});
mesh.add_boundary(4, {'x==2','y==0'});

% Create system and add matrix components
sys = System(mesh);
sys.add_constant('E', 3e11, 'v', 0.3, 'r', [1000;1000]);
sys.add_constant('D', 'E / (1-v^2) * [1, v, 0; v, 1, 0; 0, 0, (1-v)/2]');
sys.add_matrix('K', 'B''*D*B');

% Add force components
sys.add_vector('f', 'N''*r', 'Boundary', 1);
sys.add_vector('f', 'N''*-r', 'Boundary', 2);

% Assemble and solve
solver = solvers.LinearSolver(sys);
solver.add_essential_boundary({'id',3,'value',0},{'id',4,'component','y','value',0});
u = solver.solve()

% Display the displacement results
mesh.plot(u,'-Deform','-ShowNodes','Colorbar','Magnitude of Disp. (m)',...
    'Patch',{'EdgeColor','k'});
xlabel('x'); ylabel('y');

% Compute the stress and strain at the Gauss points
for e = 1:mesh.n_elements; % (include for illustration, but not needed)
    
    % Extract the current element from the mesh object
    elem = mesh.element(e);
    
    % Collect the local values of T
    d(:,1) = u(elem.get_dof());
    
    % Compute the stress and strain at the Gauss points
    k = 1;
    qp = elem.quad.rules();
    for i = 1:length(qp);
        for j = 1:length(qp);
            strain(:,k) = elem.shape_deriv(qp(i),qp(j))*d;
            stress(:,k) = sys.get('D')*strain(:,k);
            k = k + 1;
        end
    end
end    

% Display the stress and strain vectors
strain, stress
