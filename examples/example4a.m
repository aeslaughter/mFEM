% Example 9.2 of Fish & Belytschko (2007)
function [stress, strain] = example4a(varargin)

% Gather options
opt.display = true;
opt = gatherUserOptions(opt, varargin{:});

% Create a Mesh object, add the single element, and initialize it
mesh = mFEM.Mesh('Space','Vector');
mesh.createNode([0,1; 0,0; 2,0.5; 2,1]);
mesh.createElement('Quad4',1:4);
mesh.init();

% Label the boundaries
mesh.addBoundary('essential','left');    % essential boundaries
mesh.addBoundary('distributed','top');   % distributed load (t = -20)
mesh.update();

% Definethe constants for the problem
E = 3e7;            % modolus of elasticity
v = 0.3;            % posion's ratio
t_top = [0;-20];    % distributed load
D = E / (1-v^2) * [1, v, 0; v, 1, 0; 0, 0, (1-v)/2]; % constitutive matrix

% Create the stiffness matrix and force vector by looping over the
% elements, which in this case is a single element.
elements = mesh.getElements();
for e = 1:mesh.n_elements;
    
    % Extract the current element from the mesh object
    elem = elements(e);
    
    % Define short-hand function handles for the element shape functions
    % and determinate of the Jacobian
    W = @(i) elem.W(i);
    B = @(i) elem.shapeDeriv(elem.qp{i});
    J = @(i) elem.detJ(elem.qp{i});

    % Initialize the stiffness matrix (K) and the force vector (f), for
    % larger systems K should be sparse.
    K = zeros(elem.n_dof);
    f = zeros(elem.n_dof,1);
    
    % Loop over the quadrature points in the two dimensions to perform the
    % numeric integration
    for i = 1:length(elem.qp);
        K = K + W(i)*B(i)'*D*B(i)*J(i);
    end
    
    % Loop through the sides of the element, if the side has the boundary
    % id of 1 (top), then add the distributed load to the force vector
    % using numeric integration via the quadrature points for element side.
    [~,sid] = elem.hasTag('distributed')
    for i = 1:length(sid);
        s = sid(i);                           % current side
        dof = elem.getDof('side',s,'-local'); % local dof for side
        side = elem.buildSide(s);             % create side element
        W = @(i) side.W(i);                   % quadrature weights  
        N = @(i) side.shape(side.qp{i});      % side shape functions
        J = @(i) side.detJ(side.qp{i});       % side |J|
        for j = 1:length(side.qp);  
            f(dof) = f(dof) + W(j)*N(j)'*t_top*J(j);
        end
        delete(side); % delete the side element
    end          
end

% Define dof indices for the essential dofs and non-essential dofs
ess = mesh.getDof('Tag','essential'); % 5-8   
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
    
    % Get the current element
    elem = elements(e);
    
    % Collect the local values of T
    d(:,1) = u(elem.getDof());
    
    % Compute the stress and strain at the Gauss points
    for i = 1:length(elem.qp);
        strain(:,i) = B(i)*d;
        stress(:,i) = D*strain(:,i);
    end
end    

% Display the stress and strain vectors
if opt.display;
    strain, stress
end
