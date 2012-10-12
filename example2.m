% Example 8.1 of Fish & Belytschko (2007).
%
% Syntax:
%   example2
%
% Description:
%   example2 solves a simple two element heat conduction problem.
function example2
   
% Import the mFEM library
import mFEM.*;

% Create a FEmesh object, add the single element, and initialize it
mesh = FEmesh('Tri3');
mesh.add_element([0,0; 2,0.5; 0,1]);
mesh.add_element([2,0.5; 2,1; 0,1]);
mesh.initialize();

% Label the boundaries
mesh.add_boundary_id('top', 1);     % q = 20 boundary
mesh.add_boundary_id('right', 2);   % q = 0 boundary
mesh.add_boundary_id(3);            % essential boundaries (all others)

% Create Gauss objects for performing integration on the element and
% elements sides.
q_elem = Gauss(2,'tri');
[qp, W] = q_elem.rules();

q_face = Gauss(1);
[qp_side, W_side] = q_face.rules();

% Definethe constants for the problem
D = 5*eye(2);   % thermal conductivity matrix
s = 6;          % heat source (defined over entire domain)
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
    
    % Loop over the quadrature points in the two dimensions to perform the
    % numeric integration
    for i = 1:size(qp,1);
        fe = fe + W(i)*s*N(qp(i,1),qp(i,2))'*elem.detJ();
        Ke = Ke + W(i)*B()'*D*B()*elem.detJ();
    end

    % Re-define the N short-hand function handle for use on sides
    N = @(id, beta) elem.side_shape(id, beta);
    
    % Loop throught the sides of the element, if the side has the boundary
    % id of 1 (top), then add the prescribed flux term to the force vector
    % using numeric integration via the quadrature points for the element
    % side.
    for s = 1:elem.n_sides;
        if elem.side(s).boundary_id == 1;
            for i = 1:length(qp_side);
                elem.side_detJ(s,qp_side(i))
                fe = fe + -q_top*W_side(i)*N(s,qp_side(i))'*elem.side_detJ(s,qp_side(i));              
            end
        end
    end   
    
    % Add the local stiffness and force to the global 
    % (this method is slow, see example5 for faster method)
    dof = elem.get_dof();   
    K(dof, dof) = K(dof, dof) + Ke;
    f(dof) = f(dof) + fe;
end

full(K)
f

% % Define dof indices for the essential dofs and non-essential dofs
% non = mesh.get_dof(3,'ne'); % 4
% ess = mesh.get_dof(3);      % 1,2,3
% 
% % Solve for the temperatures
% T = zeros(size(f));         % initialize the temperature vector
% T(ess) = T_bar;             % apply essential boundary condtions
% T(non) = K(non,non)\f(non); % solve for T on the non-essential boundaries
% 
% % Solve for the reaction fluxes
% r = K*T - f;
% 
% % Display the results
% T,r