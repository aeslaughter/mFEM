classdef ConstantVector < mFEM.kernels.base.Kernel;
    %CONSTANTVECTOR Vector of length mesh.n_dof
    
    properties (Access = public)%(Access = ?mFEM.registry.base.Registry)
        value;
        reserved = {'grad'};
    end
    
%     properties (Access = protected)
%         mesh;
%     end
    
    methods 
        function obj = ConstantVector(name, mesh, input, varargin)
            
            %   Component
            %       scalar | 'x' | 'y' | 'z'
            %       Returns the dof associated with the vector space or
            %       multipe dof per node elements like the Beam element.
            %
            %   Boundary
            %       scalar
            %       Extract the dofs for the boundary specified, where the
            %       scalar value is the numeric id added using the
            %       ADD_BOUNDARY method.
            %
            %   Subdomain
            %       scalar
            %       Extract the dofs for the subdomain specified, where the
            %       scalar is the numeric tag added using ADD_SUBDOMAIN.
            obj = obj@mFEM.kernels.base.Kernel(name);
            obj.testName(name);
            obj.name = name;
            obj.value = mFEM.Vector(mesh);
            
            dof = mesh.getDof(varargin{:},'-index');
            if isscalar(input);
                input = ones(length(dof),1)*input;
            end

            if ~iscolumn(input); input = input'; end
            obj.value.add(input, dof);
        end        

        function str = apply(obj, str, elem, x)
            
            % search for u and grad(u)
            % apply by calling shape or shapeDeriv -> pointValue pointGrad
            
            if ~ischar(str);
                error('Kernel:apply', 'The supplied input (str) must be a character string');
            end 
            
            expr = ['\<',obj.name,'\>'];
            repstr = mat2str(obj.valu.getLocal(dof));
            str = regexprep(str, expr, repstr); 
        end
        
        function output = pointValue(obj, elem, x)
            %POINTVALUE extract the value for a vector at a point x
            
            % Degrees of freedom
            dof = elem.getDof();

            % Extract local vector
            u = obj.value.getLocal(dof);

            % Output the value
            output = elem.shape(x)*u;
        end
        
        function M = pointGradient(obj, elem, x)
            %POINTGRADIENT extract the gradient of a vector at a point x
            
            % Degrees of freedom
            dof = elem.getDof();

            % Extract local vector
            u = obj.value.getLocal(dof);

            % Output the value
            p = elem.shapeDeriv(x)*u;
            
            % Build proper matrix output
            M = ones(elem.n_dim);
            if elem.n_dim == 1;
                M = p;
            elseif elem.n_dim == 2;
                M(1,1) = p(1);  M(1,2) = p(3);
                M(2,1) = p(3);  M(2,2) = p(2);
            elseif elem.n_dim == 3;
                M(1,1) = p(1);  M(1,2) = p(4); M(1,3) = p(6);
                M(2,1) = p(4);  M(2,2) = p(2); M(2,3) = p(5);
                M(3,1) = p(6);  M(3,2) = p(5); M(3,3) = p(3);
            end           
        end
        
        function value = get(obj)
            value = obj.value();
        end
        
        function eval(~, varargin)
            error('ConstantVector:eval', 'The eval function is not implmented for the ConstantVector kernel.');
            %value = eval(obj.value);
        end 
        
    end

end