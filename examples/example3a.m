% Example 8.2 of Fish & Belytschko (2007).
%
% Syntax:
%   example3a
%
% Description:
%   example3a solves a simple single element heat conduction problem.
function example3a
   
% Import the mFEM library
import mFEM.*;

% Create a FEmesh object, add the single element, and initialize it
mesh = FEmesh();
mesh.addElement('Quad4',[0,1; 0,0; 2,0.5; 2,1]);
mesh.init();

% Label the boundaries
mesh.addBoundary(1, 'top');     % q = 20 boundary
mesh.addBoundary(2, 'right');   % q = 0 boundary
mesh.addBoundary(3);            % essential boundaries (all others)

% Definethe constants for the problem
D = 5*eye(2);   % thermal conductivity matrix
s = 6;          % heat source (defined over entire domain)
q_top = -20;    % top boundary prescribed heat flux
T_bar = 0;      % known temperatures

% Create the stiffness matrix and force vector by looping over the
% elements, which in this case is a single element.
for e = 1:mesh.n_elements;
    
    % Extract the current element from the mesh object
    elem = mesh.element(e);
    
    % Define short-hand function handles for the element shape functions
    % and shape function derivatives
    B = @(xi) elem.shapeDeriv(xi);
    N = @(xi) elem.shape(xi);

    % Initialize the stiffness matrix (K) and the force vector (f), for
    % larger systems K should be sparse.
    K = zeros(elem.n_dof);
    f = zeros(elem.n_dof,1);
    
    % Loop over the quadrature points in the two dimensions to perform the
    % numeric integration
    for i = 1:length(elem.qp);
        f = f + elem.W(i)*s*N(elem.qp{i})'*elem.detJ(elem.qp{i});
        K = K + elem.W(i)*B(elem.qp{i})'*D*B(elem.qp{i})*elem.detJ(elem.qp{i});
    end
  
    % Loop throught the sides of the element, if the side has the boundary
    % id of 1 (top), then add the prescribed flux term to the force vector
    % using numeric integration via the quadrature points for element side.
    for s = 1:elem.n_sides;
        if any(elem.side(s).boundary_id == 1);
            dof = elem.getDof('Side',s,'-local'); % local dofs for the current side
            side = elem.buildSide(s);
            N = @(xi) side.shape(xi);
            for i = 1:length(side.qp);
                f(dof) = f(dof) + q_top*side.W(i)*N(side.qp{i})'*side.detJ(side.qp{i});              
            end
            delete(side);
        end
    end  
end

% Define dof indices for the essential dofs and non-essential dofs
ess = mesh.getDof('Boundary',3); % 4
non = ~ess;

% Solve for the temperatures
T = zeros(size(f));         % initialize the temperature vector
T(ess) = T_bar;             % apply essential boundary condtions
T(non) = K(non,non)\f(non); % solve for T on the non-essential boundaries

% Solve for the reaction fluxes
r = K*T - f;

% Display the results
T,r

% Compute the flux values at the Gauss points
% Loop through the elements
for e = 1:mesh.n_elements; % (include for illustration, but not needed)
    
    % Extract the current element from the mesh object
    elem = mesh.element(e);
    
    % Collect the local values of T
    d(:,1) = T(elem.getDof());
    
    % Compute the flux at the Gauss points
    for i = 1:length(elem.qp);
        q(:,i) = -D*B(elem.qp{i})*d;
    end
end    

% Display the flux vectors
q