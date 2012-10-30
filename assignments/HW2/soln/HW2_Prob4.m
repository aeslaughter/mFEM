% Solution: Homework2, Problem 4 
function HW2_Prob4(type,nel)

% Import mFEM library
import mFEM.*

% Create mesh
mesh = FEmesh();
mesh.grid(type,0,4,nel);
mesh.init();

%Add boundary identification
mesh.add_boundary('left', 1);   % convective boundary
mesh.add_boundary('right', 2);  % q = 5 boundary  

% Define system
sys = System(mesh);
sys.add_constant('k',2,'A',0.1,'b',5,'q_bar',5,'h',0.1,'T_inf',10); 
sys.add_matrix('K', 'B''*k*A*B');
sys.add_matrix('K_h', 'N''*h*N');
sys.add_vector('f_s', 'N''*b');
sys.add_vector('f_h', '-T_inf*N''',1);
sys.add_vector('f_q', '-q_bar*A*N''', 2);

% Assembly stiffness matrix and force vector
K = sys.assemble('K') + sys.assemble('K_h');
f = sys.assemble('f_s') + sys.assemble('f_q') + sys.assemble('f_h');

% Solve for the temperature
T = K\f;

% Create Temperature plot
h = subplot(2,1,1);
mesh.plot(T,'AxesHandle',h,'ShowNodes','on'); hold on;
xlabel('x (m)','interpreter','tex');
ylabel('Temperature (\circC)','interpreter','tex');

% Compute the temperature gradients for the elements
% Compute the temperature gradients for the elements
k = 1;
for e = 1:mesh.n_elements; % loop through the elements
    
    % Extract the current element from the mesh object
    elem = mesh.element(e);
    
    % Get the Gauss points
    qp = elem.quad.rules();
    
    % Collect the local values of T
    d(:,1) = T(elem.get_dof());
    
    % Compute the temperature gradient at the gauss points
    for i = 1:length(qp);
        TG(k) = elem.shape_deriv(qp(i))*d;
        TGx(k) = elem.get_position(qp(i)); 
        k = k + 1;
    end
end 

% Create TG plot
lims = get(h,'XLim');
h = subplot(2,1,2);
plot(h,TGx,TG,'bo','LineWidth',1);
xlabel('x (m)','interpreter','tex');
ylabel('Temp. Gradient (\circC/m)','interpreter','tex');
axis auto;
set(h,'XLim',lims);


