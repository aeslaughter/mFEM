classdef MatrixKernel < mFEM.base.Kernel;

%     properties
%        const; 
%         
%     end
    
    
    methods 
        function obj = MatrixKernel(name)
            obj = obj@mFEM.base.Kernel(name);
            
            
            
            
            
            
        end

%         function output = init(obj)
%             
%             output = init@mFEM.base.Kernel(obj);
%             
% %             % Insert element shape function and shape function derivatives
% %             output = regexprep(output,'\<N\>','elem.shape(qp{:})');
% %             output = regexprep(output,'\<B\>','elem.shape_deriv(qp{:})');
%         end
        

    end

end