% MAE4700/5700 HW8, Prob. 2
function example9b

% Import the mFEM library
import mFEM.*;
  
% Create a FEmesh object, add the single element, and initialize it
mesh = FEmesh('Space','vector','Element','Quad4');
mesh.grid(0,2,0,1,1,1);
mesh.init();

% Label the boundaries
mesh.addBoundary(1, 'bottom','left'); 
mesh.addBoundary(2, 'right','top');
mesh.addBoundary(3, {'x==0','y==0'});
mesh.addBoundary(4, {'x==2','y==0'});

% Create system and add matrix components
sys = System(mesh);
sys.addConstant('E', 3e11, 'v', 0.3, 'r', [1000;1000]);
sys.addConstant('D', 'E / (1-v^2) * [1, v, 0; v, 1, 0; 0, 0, (1-v)/2]');
sys.addMatrix('K', 'B''*D*B');

% Add force components
sys.addVector('f', 'N''*r', 'Boundary', 1);
sys.addVector('f', 'N''*-r', 'Boundary', 2);

% Assemble and solve
solver = solvers.LinearSolver(sys);
solver.addEssential({'boundary',3,'value',0},{'boundary',4,'component','y','value',0});
u = solver.solve()

% Display the displacement results
mesh.plot(u,'Colorbar','Magnitude of Disp. (m)','-Deform','-ShowNodes',...
    'Patch',{'EdgeColor','k'});
xlabel('x'); ylabel('y');

% Compute the stress and strain at the Gauss points
for e = 1:mesh.n_elements; % (include for illustration, but not needed)
    
    % Extract the current element from the mesh object
    elem = mesh.element(e);
    
    % Collect the local values of T
    d(:,1) = u(elem.getDof());
    
    % Compute the stress and strain at the Gauss points
    for i = 1:length(elem.qp);
        strain(:,i) = elem.shapeDeriv(elem.qp{i})*d;
        stress(:,i) = sys.get('D')*strain(:,i);
    end
end    

% Display the stress and strain vectors
strain, stress
