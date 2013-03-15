%EXAMPLE13b A truss example.
% Beer and Johnston, Ex. 6.2, p. 227

function varargout = example13b(varargin)

% Define user options
opt.debug = false;
opt = gatherUserOptions(opt, varargin{:});

% Call subfunction to create the truss structure
mesh = buildTruss(); 

% Label the various boundaries
mesh.addBoundary(1,{'x==8','y==10'},{'x==24','y==10'}); % 28 kip load
mesh.addBoundary(2,{'x==40','y==10'}); % 16 kip load
mesh.addBoundary(3,{'x==0','y==0'});   % pin connections
mesh.addBoundary(4,{'x==32','y==0'});  % roller connection

% Define the system
sys = mFEM.System(mesh);
sys.addConstant('A', 1, 'E', 10^7);
sys.addMatrix('K', 'A*E/L*Ke');
sys.addConstantVector('f',-28e3,'Boundary',1,'Component','y');
sys.addConstantVector('f',16e3,'Boundary',2,'Component','x', '-add');

% Solve the system
solver = mFEM.solvers.LinearSolver(sys);
solver.addEssential({'Boundary',3,'Value',0},{'Boundary',4,'Component','y','Value',0});
[u,r] = solver.solve();

% Display the reactions
ess = mesh.getDof({'Boundary',3},{'Boundary',4,'Component','y'});
r = r(ess);

% Display Results
if ~opt.debug;
    disp(['Bx = ', num2str(r(1)/1000), ' kips']);
    disp(['By = ', num2str(r(2)/1000), ' kips']);
    disp(['J = ', num2str(r(3)/1000), ' kips']);

    % Plot results
    mesh.plot([],'-new','NodeLabels',false,'ElementLabels',false);
    mesh.plot(u,'-deform','Component',2,'Scale',3,'colorbar','Magnitude of Displacement (in)');
else
    varargout = {u,r};
end
end

function mesh = buildTruss()
    %BUILDTRUSS creates the mesh for the example truss

    mesh = mFEM.FEmesh();

    % Positions of the nodes
    x = [0,8,16,24,32,40];
    y = 10;

    % Create the first 4 boxes with cross brace
    for i = 2:length(x)-1;
        mesh.addElement('Truss',[x(i-1),y; x(i-1),0]);
        mesh.addElement('Truss',[x(i-1),0; x(i),0]);
        mesh.addElement('Truss',[x(i),y; x(i-1),y]);
        mesh.addElement('Truss',[x(i),y; x(i-1),0]);   
    end

    % Add the final triangle
    mesh.addElement('Truss',[x(end-1),y; x(end-1),0]);
    mesh.addElement('Truss',[x(end-1),0; x(end),y]);
    mesh.addElement('Truss',[x(end-1),y; x(end),y]);
    mesh.init();
end
