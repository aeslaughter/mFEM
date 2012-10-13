% A transient heat transfer example
%
% Syntax:
%   example5
%
% Description:
%   example5 solves a simple transient heat conduction problem.
function example5
   
% Import the mFEM library
import mFEM.*;

% Create a FEmesh object, add the single element, and initialize it
mesh = FEmesh('Quad4');
mesh.grid(0,1,0,1,20,20);

% Label the boundaries
mesh.add_boundary_id(1);            % essential boundaries (all others)

% Create Gauss objects for performing integration on the element and
% elements sides.
q_elem = Gauss(2);
[qp, W] = q_elem.rules();

q_face = Gauss(1);
[qp_side, W_side] = q_face.rules();

% Definethe constants for the problem
exact = @(x,y,t) 1 + exp(-t) * sin(pi*x) * sin(pi*y);

% Thermal conducivity
k = 1 / (2 * pi * pi);   % thermal conductivity

% Initialize storage
I = zeros(mesh.n_elements * mesh.n_dim * mesh.n_dim, 1);
J = I;
K = I;
M = I;
% Ms = sparse(mesh.n_dof());
f = zeros(mesh.n_dof(),1);
tic;
% Create the stiffness matrix and force vector by looping over the
% elements, which in this case is a single element.
m = 1;
for e = 1:mesh.n_elements;

    % Extract the current element from the mesh object
    elem = mesh.element(e);
    
    % Define short-hand function handles for the element shape functions
    % and shape function derivatives
    B = @(xi,eta) elem.shape_deriv(xi,eta);
    N = @(xi,eta) elem.shape(xi,eta);

    % Initialize the local matrices and vector
    Me = zeros(elem.n_dof);     % mass matrix
    Ke = zeros(elem.n_dof);     % stiffness matrix
    fe = zeros(elem.n_dof,1);   % force vector
    
    % Loop over the quadrature points in the two dimensions to perform the
    % numeric integration
    for i = 1:length(qp);
        for j = 1:length(qp);
            Me = Me + W(j)*W(i)*N(qp(i),qp(j))'*N(qp(i),qp(j))*elem.detJ(qp(i),qp(j));
            Ke = Ke + W(j)*W(i)*B(qp(i),qp(j))'*k*B(qp(i),qp(j))*elem.detJ(qp(i),qp(j));
        end
    end
    
    dof = elem.get_dof();
%     Ms(dof,dof) = Me;
     n = numel(Me);
     r = m:(m-1)+n;
    [i,j,M(r)] = find(Me >= 0);
  
    I(r) = dof(i);
    J(r) = dof(j);
    m = m + n;
    
%     % Re-define the N short-hand function handle for use on sides
%     N = @(xi) side.shape(xi);
%     
%     % Loop throught the sides of the element, if the side has the boundary
%     % id of 1 (top), then add the prescribed flux term to the force vector
%     % using numeric integration via the quadrature points for element side.
%     for s = 1:elem.n_sides;
%         if elem.side(s).boundary_id == 1;
%             side = elem.build_side(s);
%             N = @(xi) side.shape(xi);
%             for i = 1:length(qp_side);
%                 dof = elem.side(s).dof; % local dofs for the current side
%                 f(dof) = f(dof) + -q_top*W_side(i)*N(qp_side(i))'*side.detJ();              
%             end
%             delete(side)
%         end
%     end      




end

M = sparse(I,J,M);
toc;
%full(M)

% Clean up
clear classes;


