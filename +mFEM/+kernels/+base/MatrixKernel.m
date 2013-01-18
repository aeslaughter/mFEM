classdef MatrixKernel < mFEM.kernels.base.Kernel;

    properties
       options = struct(...
           'boundary',[],'subdomain',[]);
    end
    
    
    methods 
        function obj = MatrixKernel(name,varargin)
            obj = obj@mFEM.kernels.base.Kernel(name,[]);
            
            obj.options = gather_user_options(obj.options, varargin{:});
            
        end

%         function output = init(obj)
%             
%             output = init@mFEM.base.Kernel(obj);
%             
% %             % Insert element shape function and shape function derivatives
% %             output = regexprep(output,'\<N\>','elem.shape(qp)');
% %             output = regexprep(output,'\<B\>','elem.shape_deriv(qp)');
%         end
        

    end

end