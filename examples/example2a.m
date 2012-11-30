% Example 8.1 of Fish & Belytschko (2007).
%
% Syntax:
%   example2
%
% Description:
%   example2 solves a simple two element heat conduction problem.
function example2a
   
% Import the mFEM library
import mFEM.*;

% Create a FEmesh object, add the single element, and initialize it
mesh = FEmesh('Element','Tri3');
mesh.add_element([0,0; 2,0.5; 0,1]);
mesh.add_element([2,0.5; 2,1; 0,1]);
mesh.init();

% Label the boundaries
mesh.add_boundary(1, 'top');     % q = 20 boundary
mesh.add_boundary(2, 'right');   % q = 0 boundary
mesh.add_boundary(3);       % essential boundaries (all others)

% Create Gauss objects for performing integration on the element and
% elements sides.
q_elem = Gauss('Order',3,'Type','tri');
[qp, W] = q_elem.rules();

q_face = Gauss('Order', 1, 'Type', 'line');
[qp_side, W_side] = q_face.rules();

% Definethe constants for the problem
D = 5*eye(2);   % thermal conductivity matrix
b = 6;          % heat source (defined over entire domain)
q_top = 20;     % top boundary prescribed heat flux
T_bar = 0;      % known temperatures

% Create an empty sparse matrix
K = sparse(mesh.n_dof, mesh.n_dof);
f = zeros(mesh.n_dof, 1);

% Create the stiffness matrix and force vector by looping over the
% elements, which in this case is a single element.
for e = 1:mesh.n_elements;
    
    % Extract the current element from the mesh object
    elem = mesh.element(e);
    
    % Define short-hand function handles for the element shape functions
    % and shape function derivatives
    N = @(xi1,xi2) elem.shape(xi1, xi2);    
    B = @() elem.shape_deriv();

    % Initialize the stiffness matrix (K) and the force vector (f), for
    % larger systems K should be sparse.
    Ke = zeros(elem.n_dof);
    fe = zeros(elem.n_dof,1);
    
    % Loop over the quadrature points to perform the numeric integration
    qp
    for i = 1:size(qp,1);
        fe = fe + W(i)*b*N(qp(i,1),qp(i,2))'*elem.detJ();
        Ke = Ke + W(i)*B()'*D*B()*elem.detJ();
    end
    
    % Loop throught the sides of the element, if the side has the boundary
    % id of 1 (top), then add the prescribed flux term to the force vector
    % using numeric integration via the quadrature points for element side.
    for s = 1:elem.n_sides;
        if elem.side(s).boundary_id == 1;
            % Local dofs for the current side
            dof = elem.get_dof('Side',s,'-local'); 
           
            % Build the side element
            side = elem.build_side(s);
            
            % Perform Guass quadrature
            for i = 1:length(qp_side);
                fe(dof) = fe(dof) + -q_top*W_side(i)*side.shape(qp_side(i))'*side.detJ(qp_side(i));              
            end
            delete(side)
        end
    end  
    
    % Add the local stiffness and force to the global (this method is slow)
    dof = elem.get_dof();   
    K(dof, dof) = K(dof, dof) + Ke;
    f(dof) = f(dof) + fe;
end

% Define dof indices for the essential dofs and non-essential dofs
ess = mesh.get_dof('Boundary',3); % 4
non = ~ess;

% Solve for the temperatures
T = zeros(size(f));         % initialize the temperature vector
T(ess) = T_bar;             % apply essential boundary condtions
T(non) = K(non,non)\f(non); % solve for T on the non-essential boundaries

% Solve for the reaction fluxes
r = K*T - f;

% Display the results
T,r

% Compute the flux values for each element
% Loop through the elements
for e = 1:mesh.n_elements; 
    
    % Extract the current element from the mesh object
    elem = mesh.element(e);
    
    % Collect the local values of T
    d(:,1) = T(elem.get_dof());
    
    % Compute the flux at the Gauss points
    q(:,e) = -D*B()*d;

end    

% Display the flux vectors
q