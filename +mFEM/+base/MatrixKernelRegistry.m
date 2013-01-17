classdef MatrixKernelRegistry < handle
    %CONSTANTKERNELREGISTRY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        matrices;% = mFEM.base.MatrixKernel.empty();
    end
    
    methods %(Access = Public)
        function obj = MatrixKernelRegistry()
        end
        
        function add(obj,name,value,varargin)
            
            obj.test_name(name);
            [idx,found] = obj.locate(name,obj.constants);

            p = mfilename('fullpath')

            
            
            if isa(value,'mFEM.base.MatrixKernel');
                obj.matrices(idx) = value;
                
            %elseif ischar
                
                
            end    
               
            
%             new_const = mFEM.base.ConstantKernel(name,value);
% 
%             if found && opt.add;    
%                 new_const.merge(obj.constants(idx));
%             end
%             
%             obj.constants(idx) = new_const;
%             obj.apply(new_const);        
        end
        
        
    end
    
end

