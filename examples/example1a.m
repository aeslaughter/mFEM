function T = example1a(varargin) 
    %% Example 1a: Reproduces Example 5.1 of Fish & Belytschko (2007)
    % This example uses manually assembly of the finite element stiffness
    % matrix and force vector for a 1D problem.
    % 
    % Syntax
    %   example1a
    %   example1a('PropertyName',PropertyValue,...)
    %
    % Description
    %   example1a runs example exactly as done in textbook, with two 
    %   elements.
    %
    %   example1a('PropertyName',PropertyValue,...) allows the user to
    %   customize the behavior with the properties given below.
    %
    % EXAMPLE1A Property Descriptions
    %   N
    %       integer
    %       The number of elements, the default is 2.
    %
    %   Debug
    %       true | {false}
    %       If set to true, only the temperature is returned and the
    %       plots are not created.
    %
    % See also EXAMPLE1B EXAMPLE1C
    %
    %----------------------------------------------------------------------
    %  mFEM: A Parallel, Object-Oriented MATLAB Finite Element Library
    %  Copyright (C) 2013 Andrew E Slaughter
    % 
    %  This program is free software: you can redistribute it and/or modify
    %  it under the terms of the GNU General Public License as published by
    %  the Free Software Foundation, either version 3 of the License, or
    %  (at your option) any later version.
    % 
    %  This program is distributed in the hope that it will be useful,
    %  but WITHOUT ANY WARRANTY; without even the implied warranty of
    %  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %  GNU General Public License for more details.
    % 
    %  You should have received a copy of the GNU General Public License
    %  along with this program. If not, see <http://www.gnu.org/licenses/>.
    %
    %  Contact: Andrew E Slaughter (andrew.e.slaughter@gmail.com)
    %----------------------------------------------------------------------

    % Gather options
    opt.N = 2;
    opt.debug = false;
    opt = gatherUserOptions(opt,varargin{:});

    %% Create a Mesh Object 
    % Build mesh of 2-node linear elemnts from 0 to 4
    mesh = mFEM.Mesh();
    mesh.grid('Line2',0,4,opt.N); 

    %% Label The Boundaries
    mesh.addBoundary('essential','left');    % T = 0 boundary (essential)    
    mesh.addBoundary('flux','right');        % q = 20 boundary  
    mesh.update();

    %% Define the constants for the problem
    k = 2;          % thermal conductivity 
    A = 0.1;        % cross sectional area
    b = 5;          % heat source (defined over entire domain)
    q_bar = 5;      % right boundary prescribed heat flux
    T_bar = 0;      % known temperatures

    %% Initilize Stiffness Matrix and Force Vector
    K = sparse(mesh.n_dof, mesh.n_dof);
    f = zeros(mesh.n_dof, 1);

    %% Manual Assembly
    % Create the stiffness matrix and force vector by looping over the
    % elements, which in this case is a single element.
    elem = mesh.getElements();
    for e = 1:mesh.n_elements;

        % Define short-hand function handles for the element shape functions
        % and shape function derivatives
        B = @(xi) elem(e).shapeDeriv([]);
        N = @(xi) elem(e).shape(xi);

        % Initialize the stiffness matrix (K) and the force vector (f), for
        % larger systems K should be sparse.
        Ke = zeros(elem(e).n_dof);
        fe = zeros(elem(e).n_dof,1);

        % Loop over the quadrature points in the two dimensions to perform the
        % numeric integration
        [qp,W] = elem(e).quad.rules();
        for i = 1:length(qp);

            % Account for the source term
            fe = fe + W(i)*N(qp(i))'*b*elem(e).detJ([]);

            % Build stiffness matrix
            Ke = Ke + W(i)*B(qp(i))'*k*A*B(qp(i))*elem(e).detJ([]);
        end

        % Loop throught the sides of the element, if the side has the boundary
        % id of 2 (right), then add the prescribed flux term to force vector,
        % which for 1D elements is an evaluation at a single point so it
        % requires no integration.
        [~,sid] = elem(e).hasTag('flux');
        for s = 1:length(sid); 
            dof = elem(e).getDof('Side',sid(s),'-local');
            fe(dof) = fe(dof) + -q_bar*A;  
        end   

        % Add the local stiffness and force to the global 
        % (this method is slow, see example5 for faster method)
        dof = elem(e).getDof();
        K(dof, dof) = K(dof, dof) + Ke;
        f(dof) = f(dof) + fe;
    end

    
    K
    f
    
    return;
    
    %% Define Variables for Essential and Non-essential Degrees-of-freedom
    ess = mesh.getDof('Tag', 'essential');    % 1
    non = ~ess;                                         % 2,3

    %% Solve for the Temperatures
    T = zeros(size(f));         % initialize the temperature vector
    T(ess) = T_bar;             % apply essential boundary condtions
    T(non) = K(non,non)\f(non); % solve for T on the non-essential boundaries

    % Stop here if in debug mode
    if opt.debug; 
        return;
    end
    
    %% Compute the Temperature Gradients
    % Loop through the elements
    for e = 1:mesh.n_elements;

        % Collect the local values of T
        d(:,1) = T(elem(e).getDof());

        % Compute the temperature gradient at the gauss point, store the value
        % twice for each element for creating graph, TGx is the node locations
        % used for plotting
        TG(1:2,e) = B(qp(1))*d;
        TGx(1:2,e) = elem(e).nodes.getCoord();
    end    

    %% Generate Figure for T and TG Solutions
    % Create Exact Solutions
    x0 = 0:0.1:4;
    Tex = -12.5*x0.^2 + 97.5*x0;
    TGex = -25*x0 + 97.5;

    % Initilize the figure
    figure('Color','w','Name','Example 1 Results');

    % Create Temperature Plot
    h = subplot(2,1,1); hold on;
    plot(h,x0,Tex,'b-','LineWidth',2);
    mesh.plot(T,'-ShowNodes','Patch',{'EdgeColor','k'});
    legend({'Exact','FEM'},'location','best');
    xlabel('x (m)','interpreter','tex');
    ylabel('Temperature (\circC)','interpreter','tex');

    % Create TG Plot
    h = subplot(2,1,2);
    plot(h,x0,TGex,'b-',TGx,TG,'k-o','LineWidth',1);
    legend({'Exact','FEM'},'location','best');
    xlabel('x (m)','interpreter','tex');
    ylabel('Temp. Gradient (\circC/m)','interpreter','tex');
