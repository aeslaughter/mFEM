% Example 8.2 of Fish & Belytschko (2007).
%
% Syntax:
%   example3
%
% Description:
%   example3 solves a simple single element heat conduction problem.
function example3a
   
% Import the mFEM library
import mFEM.*;

% Create a FEmesh object, add the single element, and initialize it
mesh = FEmesh();
mesh.add_element('Quad4',[0,1; 0,0; 2,0.5; 2,1]);
mesh.initialize();

% Label the boundaries
mesh.add_boundary('top', 1);     % q = 20 boundary
mesh.add_boundary('right', 2);   % q = 0 boundary
mesh.add_boundary(3);            % essential boundaries (all others)

% Create Gauss objects for performing integration on the element and
% elements sides.
q_elem = Gauss(2,'quad');
[qp, W] = q_elem.rules();

q_face = Gauss(1,'line');
[qp_side, W_side] = q_face.rules();

% Definethe constants for the problem
D = 5*eye(2);   % thermal conductivity matrix
s = 6;          % heat source (defined over entire domain)
q_top = 20;     % top boundary prescribed heat flux
T_bar = 0;      % known temperatures

% Create the stiffness matrix and force vector by looping over the
% elements, which in this case is a single element.
for e = 1:mesh.n_elements;
    
    % Extract the current element from the mesh object
    elem = mesh.element(e);
    
    % Define short-hand function handles for the element shape functions
    % and shape function derivatives
    B = @(xi,eta) elem.shape_deriv(xi,eta);
    N = @(xi,eta) elem.shape(xi,eta);

    % Initialize the stiffness matrix (K) and the force vector (f), for
    % larger systems K should be sparse.
    K = zeros(elem.n_dof);
    f = zeros(elem.n_dof,1);
    
    % Loop over the quadrature points in the two dimensions to perform the
    % numeric integration
    for i = 1:length(qp);
        for j = 1:length(qp);
            f = f + W(j)*W(i)*s*N(qp(i),qp(j))'*elem.detJ(qp(i),qp(j));
            K = K + W(j)*W(i)*B(qp(i),qp(j))'*D*B(qp(i),qp(j))*elem.detJ(qp(i),qp(j));
        end
    end
  
    % Loop throught the sides of the element, if the side has the boundary
    % id of 1 (top), then add the prescribed flux term to the force vector
    % using numeric integration via the quadrature points for element side.
    for s = 1:elem.n_sides;
        if elem.side(s).boundary_id == 1;
            dof = elem.get_dof(s); % local dofs for the current side
            side = elem.build_side(s);
            N = @(xi) side.shape(xi);
            for i = 1:length(qp_side);
                f(dof) = f(dof) + -q_top*W_side(i)*N(qp_side(i))'*side.detJ(qp_side(i));              
            end
            delete(side);
        end
    end  
end

% Define dof indices for the essential dofs and non-essential dofs
non = mesh.get_dof(3,'ne'); % 4
ess = mesh.get_dof(3);      % 1,2,3

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
    d(:,1) = T(elem.get_dof());
    
    % Compute the flux at the Gauss points
    k = 1;
    for i = 1:length(qp);
        for j = 1:length(qp);
            q(:,k) = -D*B(qp(i),qp(j))*d;
            k = k + 1;
        end
    end
end    

% Display the flux vectors
q