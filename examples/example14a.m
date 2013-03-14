%EXAMPLE14a Beam element
% Fish and Belytschko (2009), Ex. 10.1, p. 262
%
% This is a specialized example problem, refer that replicates the example
% problem from the text exactly. For a more rubust approach refer to
% example 14b
%
% See Also EXAMPLE14b

function example14a

% Load the mFEM library
import mFEM.*

% Create the 2-element mesh
mesh = FEmesh();
mesh.addElement('Beam',[0;8]);
mesh.addElement('Beam',[8;12]);
mesh.init();

% Add flags for the boundary elements
mesh.addBoundary(1,'left');    % fixed connection
mesh.addBoundary(2,'right');   % prescribed force and rotation

% Add flag for the distributed load
mesh.addSubdomain(10, 'x<8');    % distributed load, b

% Difine the paramters for the problem
b = -1;         % body force
P1 = -10;       % load at x = 4
P2 = 5;         % load at x = 8
EI = 10^4;      % modulus and moment of intertia
c = [-20;20];   % presecribed force and moment on boundary

% Initialize global stiffness matrix and force vector
K = Matrix(mesh);
f = zeros(mesh.n_dof,1);

% Loop through the elements
for e = 1:mesh.n_elements;
    
    % Current element
    elem = mesh.element(e);
    
    % Initilize local K and f
    Ke = zeros(elem.n_dof);    
    fe = zeros(elem.n_dof,1);    
    
    % Create functions for N and B
    B = @(xi) elem.shapeDeriv(xi);
    N = @(xi) elem.shape(xi);

    % Perform quadrature for stiffness matrix
    for i = 1:length(elem.qp);
        Ke = Ke + elem.W(i)*EI*B(elem.qp{i})'*B(elem.qp{i})*elem.detJ(elem.qp{i});
    end

    % Account for the body force, only on subdomain
    if any(elem.subdomain == 10);
        for i = 1:length(elem.qp);  
            fe = fe + elem.W(i)*N(elem.qp{i})'*b*elem.detJ(elem.qp{i});
        end
    end
    
    % Account for point loads and the boundary loads at right end
    if e == 1;
        fe = fe + N(0)'*P1;
    elseif e == 2;
        fe = fe + N(-1)'*P2;
        dof = elem.getDof('Side',2,'-local');
        fe(dof) = fe(dof) + c;
    end

    % Add the local matrix and vector to the global
    dof = elem.getDof();
    K.add(Ke, dof);
    f(dof) = f(dof) + fe;
end

% Create the stiffness matrix
K = K.init();
full(K)
f
% Solve for the unknowns
ess = mesh.getDof('Boundary' , 1);
u(ess) = 0;
u(~ess) = K(~ess,~ess)\f(~ess);
u
% Plot the results
h = subplot(2,1,1);
%mesh.plot(u,'Axes', h, 'Component', 1, 'colorbar', 'Displacement (m)');
mesh.plot(u, 'Axes', h, 'Component', 1, '-deform', 'patch', {'EdgeColor','k'});
xlabel('x (m)');
ylabel('Displacement (m)');

h = subplot(2,1,2);
% mesh.plot(u,'Axes', h, 'Component',2, 'colorbar', 'Rotation');
mesh.plot(u, 'Axes', h, 'Component', 2, '-deform', 'patch', {'EdgeColor','k'});
xlabel('x (m)');
ylabel('Rotation (rad.)');