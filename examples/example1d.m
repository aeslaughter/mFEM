function T = example1d(varargin)
    %EXAMPLE1D Reproduces Example 5.1 of Fish & Belytschko (2007)
    %   This example is a fully parallel, manually assembled example
    % 
    % Syntax
    %   example1d
    %   example1d('PropertyName',PropertyValue,...)
    %
    % Description
    %   example1d runs example exactly as done in textbook, with two 
    %   elements.
    %
    %   example1d('PropertyName',PropertyValue,...) allows the user to
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

tic;
    % Gather options
    opt.n = 200;
    opt.debug = true;
    opt = gatherUserOptions(opt,varargin{:});
   
    % Build mesh of 2-node linear elemnts from 0 to 4
    mesh = mFEM.Mesh();
    mesh.grid('Line2',0,4,opt.n); 

    %% Label The Boundaries
    mesh.addBoundary(1,'left');    % T = 0 boundary (essential)    
    mesh.addBoundary(2,'right');   % q = 20 boundary   
%     mesh.update();

    %% Define Constants
    k = 2;          % thermal conductivity 
    A = 0.1;        % cross sectional area
    b = 5;          % heat source (defined over entire domain)
    q_bar = 5;      % right boundary prescribed heat flux
    T_bar = 0;      % known temperatures
toc;
tic;
    % Parallel Manual Assembly
    % Extract the elements (a MATLAB Composite)
    elements = mesh.getElements();
    
    % Begin parallel block
    spmd   
        n_dof = sum([elements.n_dof]);
        K = mFEM.Matrix(n_dof);
        f = mFEM.Matrix(n_dof,1);
        
        for e = 1:length(elements);
            
            % Get the current element
            elem = elements(e);

            % Initialize the stiffness matrix (K) and the force vector (f), for
            % larger systems K should be sparse.
            Ke = zeros(elem.n_dof);
            fe = zeros(elem.n_dof,1);

            % Loop over the quadrature points in the two dimensions to perform the
            % numeric integration
            [qp,W] = elem.quad.rules();
            for i = 1:length(qp);

                % Account for the source term
                fe = fe + W(i)*elem.shape(qp(i))'*b*elem.detJ([]);

                % Build stiffness matrix
                Ke = Ke + W(i)*elem.shapeDeriv(qp(i))'*k*A*elem.shapeDeriv(qp(i))*elem.detJ([]);
            end

            % Loop throught the sides of the element, if the side has the boundary
            % id of 2 (right), then add the prescribed flux term to force vector,
            % which for 1D elements is an evaluation at a single point so it
            % requires no integration.
            [~,sid] = elem.hasTag(2);
            for s = 1:length(sid); 
                dof = elem.getDof('Side',sid(s),'-local');
                fe(dof) = fe(dof) + -q_bar*A;  
            end   
            
            % Add the local stiffness and force to the global 
            dof = elem.getDof();
            K.add(Ke,dof);
            f.add(fe,dof,1); 
        end

        ipart = buildCodistributedPartition(length(K.I));
        jpart = buildCodistributedPartition(length(K.J));       
        spart = buildCodistributedPartition(length(K.Aij));   
        fpart = buildCodistributedPartition(length(f.init())); 
        
        idist = codistributor1d(1,ipart,[sum(ipart),1]);
        jdist = codistributor1d(1,jpart,[sum(jpart),1]);  
        sdist = codistributor1d(1,jpart,[sum(spart),1]);  
        fdist = codistributor1d(1,fpart,[sum(fpart),1]);  
        
        I = codistributed.build(K.I,idist);
        J = codistributed.build(K.J,jdist);
        S = codistributed.build(K.Aij,sdist);
        f = codistributed.build(f.init(),fdist);
    end
    toc;
    tic;
    % Parallel mldivide on sparse does not yet work!
    K = gather(sparse(I,J,S));
    f = gather(f);

    %% Define Variables for Essential and Non-essential Degrees-of-freedom
    ess = mesh.getDof('Tag',1,'-gather');    % 1
    non = ~ess;                               % 2,3

    %% Solve for the Temperatures
    T = zeros(size(f));         % initialize the temperature vector
    T(ess) = T_bar;             % apply essential boundary condtions
    T(non) = K(non,non)\f(non); % solve for T on the non-essential boundaries
    toc;