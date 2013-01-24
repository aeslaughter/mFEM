% Example 9.2 of Fish & Belytschko (2007)
function [stress, strain] = example4a(varargin)

% Gather options
opt.display = true;
opt = gatherUserOptions(opt, varargin{:});

% Import the mFEM library
import mFEM.*;
  
% Create a FEmesh object, add the single element, and initialize it
mesh = FEmesh('Space','vector');
mesh.addElement('Quad4',[0,1; 0,0; 2,0.5; 2,1]);
mesh.init();

% Label the boundaries
mesh.addBoundary(1, 'left');    % essential boundaries
mesh.addBoundary(2, 'top');     % distributed load (t = -20)

% Definethe constants for the problem
E = 3e7;            % modolus of elasticity
v = 0.3;            % posion's ratio
t_top = [0;-20];    % distributed load
D = E / (1-v^2) * [1, v, 0; v, 1, 0; 0, 0, (1-v)/2]; % constitutive matrix

% Create the stiffness matrix and force vector by looping over the
% elements, which in this case is a single element.
for e = 1:mesh.n_elements;
    
    % Extract the current element from the mesh object
    elem = mesh.element(e);
    
    % Define short-hand function handles for the element shape functions
    % and determinate of the Jacobian
    B = @(xi) elem.shapeDeriv(xi);
    J = @(xi) elem.detJ(xi);

    % Initialize the stiffness matrix (K) and the force vector (f), for
    % larger systems K should be sparse.
    K = zeros(elem.n_dof);
    f = zeros(elem.n_dof,1);
    
    % Loop over the quadrature points in the two dimensions to perform the
    % numeric integration
    for i = 1:length(elem.qp);
        K = K + elem.W(i)*B(elem.qp{i})'*D*B(elem.qp{i})*J(elem.qp{i});
    end
    
    % Loop through the sides of the element, if the side has the boundary
    % id of 1 (top), then add the distributed load to the force vector
    % using numeric integration via the quadrature points for element side.
    for s = 1:elem.n_sides;
        if any(elem.side(s).boundary_id == 2);
            dof = elem.getDof('side',s,'-local'); % local dof for side
            side = elem.buildSide(s);             % create side element

            for i = 1:length(side.qp);  
                f(dof) = f(dof) + side.W(i)*side.shape(side.qp{i})'*t_top*side.detJ(side.qp{i});
            end
            delete(side); % delete the side element
        end
    end          
end

% Define dof indices for the essential dofs and non-essential dofs
ess = mesh.getDof('Boundary', 1); % 5-8   
non = ~ess; % 1-4       

% Solve for the temperatures
u = zeros(size(f));         % initialize the displacement vector
u(ess) = 0;                 % apply essential boundary condtions
u(non) = K(non,non)\f(non); % solve for T on the non-essential boundaries

% Solve for the reaction fluxes
r = K*u - f;

% Display the results
if opt.display;
    u, r
end

% Compute the stress and strain at the Gauss points
% Loop through the elements
for e = 1:mesh.n_elements; % (include for illustration, but not needed)
    
    % Extract the current element from the mesh object
    elem = mesh.element(e);
    
    % Collect the local values of T
    d(:,1) = u(elem.getDof());
    
    % Compute the stress and strain at the Gauss points
    for i = 1:length(elem.qp);
        strain(:,i) = B(elem.qp{i})*d;
        stress(:,i) = D*strain(:,i);
    end
end    

% Display the stress and strain vectors
if opt.display;
    strain, stress
end
