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
            
            obj.value = mFEM.Vector(mesh);
            
            dof = mesh.getDof(varargin{:},'-index');
            if isscalar(input);
                input = ones(length(dof),1)*input;
            end

            obj.value.add(input,dof);
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
%         
%         function output = point_value(obj, name, elem, x)
%             %POINT_VALUE extract the value for a vector at a point x
%             
%             [type,idx] = obj.locate(name);
%             if ~strcmp(type,'vector');
%                 error('System:point_value', 'This method only works with vectors');
%             end
%             
%             % Degrees of freedom
%             dof = elem.get_dof();
% 
%             % Extract local vector
%             u = obj.vec(idx).vector.get_local(dof);
% 
%             % Shape functions at point x
%             N = elem.shape(x{:});
% 
%             % Output the value
%             output = N*u;
%         end
%         
%         function output = point_gradient(obj, name, elem, x)
%             %POINT_GRADIENT extract the gradient of a vector at a point x
%             
%             [type,idx] = obj.locate(name);
%             if ~strcmp(type,'vector');
%                 error('System:point_value', 'This method only works with vectors');
%             end
%             
%             % Degrees of freedom
%             dof = elem.get_dof();
% 
%             % Extract local vector
%             u = obj.vec(idx).vector.get_local(dof);
% 
%             % Shape functions at point x
%             B = elem.shape_deriv(x{:});
% 
%             % Output the value
%             output = B*u;
%         end
        function value = get(obj)
            value = obj.value();
        end
        
        function eval(~, varargin)
            error('ConstantVector:eval', 'The eval function is not implmented for the ConstantVector kernel.');
            %value = eval(obj.value);
        end  
    end
end