%EXAMPLE1A Reproduces Example 5.1 of Fish & Belytschko (2007) using direct assembly
% of the stiffness matrix and force vector.
% 
% Syntax:
%   example1a
%   example1a(n)
%
% Description
%   example1a runs example exactly as done in textbook, with two elements
%   example1a(n) runs the example with n number of elements
%
% See also EXAMPLE1B

function example1a(varargin)

% Import the mFEM library
import mFEM.* elements.*;

% Parse optional input
nel = 2;
if nargin == 1 && isnumeric(varargin{1}); 
    nel = varargin{1};
end
    
% Create a FEmesh object, 2 node linear mesh from 0 to 4
mesh = FEmesh();
mesh.grid('Line2',0,4,nel); 
mesh.init();

% Label the boundaries
mesh.add_boundary(1,'left');    % T = 0 boundary (essential)    
mesh.add_boundary(2,'right');   % q = 20 boundary   

% Create Gauss objects for performing integration on the element
q_elem = Gauss(1);
[qp, W] = q_elem.rules();

% Definethe constants for the problem
k = 2;          % thermal conductivity 
A = 0.1;        % cross sectional area
b = 5;          % heat source (defined over entire domain)
q_bar = 5;      % right boundary prescribed heat flux
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
    B = @(xi) elem.shape_deriv();
    N = @(xi) elem.shape(xi);

    % Initialize the stiffness matrix (K) and the force vector (f), for
    % larger systems K should be sparse.
    Ke = zeros(elem.n_dof);
    fe = zeros(elem.n_dof,1);
    
    % Loop over the quadrature points in the two dimensions to perform the
    % numeric integration
    for i = 1:length(qp);

        % Account for the source term
        fe = fe + W(i)*N(qp(i))'*b*elem.detJ();
        
        % Build stiffness matrix
        Ke = Ke + W(i)*B(qp(i))'*k*A*B(qp(i))*elem.detJ();
    end
    
    % Re-define the N short-hand function handle for use on sides
    N = @(id) elem.shape(id);
    
    % Loop throught the sides of the element, if the side has the boundary
    % id of 2 (right), then add the prescribed flux term to force vector,
    % which for 1D elements is an evaluation at a single point so it
    % requires no integration.
    for s = 1:elem.n_sides; 
        if elem.side(s).boundary_id == 2;
            side = elem.build_side(s);
            dof = elem.get_dof('Side', s, '-local');
            fe(dof) = fe(dof) + -q_bar*A*side.shape()';  
            delete(side);
        end
    end   
    
    % Add the local stiffness and force to the global 
    % (this method is slow, see example5 for faster method)
    dof = elem.get_dof();   
    K(dof, dof) = K(dof, dof) + Ke;
    f(dof) = f(dof) + fe;
end

% Define dof indices for the essential dofs and non-essential dofs
ess = mesh.get_dof('Boundary', 1);  % 1
non = ~ess;                         % 2,3
ess
non
% Solve for the temperatures
T = zeros(size(f));         % initialize the temperature vector
T(ess) = T_bar;             % apply essential boundary condtions
T(non) = K(non,non)\f(non); % solve for T on the non-essential boundaries

% Compute the temperature gradients for the elements
% Loop through the elements
for e = 1:mesh.n_elements;
    
    % Extract the current element from the mesh object
    elem = mesh.element(e);
    
    % Collect the local values of T
    d(:,1) = T(elem.get_dof());
    
    % Compute the temperature gradient at the gauss point, store the value
    % twice for each element for creating graph, TGx is the node locations
    % used for plotting
    TG(1:2,e) = B(qp(1))*d;
    TGx(1:2,e) = elem.nodes;

end    

% Create exact solutions
x0 = 0:0.1:4;
Tex = -12.5*x0.^2 + 97.5*x0;
TGex = -25*x0 + 97.5;

% Generate figure for T and TG solutions
figure('Color','w','Name','Example 1 Results');

% Create Temperature plot
h = subplot(2,1,1);
mesh.plot(T,'ShowNodes',true); hold on;
plot(h,x0,Tex,'k-','LineWidth',1);
legend({'FEM','Exact'},'location','best');
xlabel('x (m)','interpreter','tex');
ylabel('Temperature (\circC)','interpreter','tex');

% Create TG plot
h = subplot(2,1,2);
plot(h,x0,TGex,'k-',TGx,TG,'b-o','LineWidth',1);
legend({'Exact','FEM'},'location','best');
xlabel('x (m)','interpreter','tex');
ylabel('Temp. Gradient (\circC/m)','interpreter','tex');
