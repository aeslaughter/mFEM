classdef Gauss
    properties (SetAccess = private, GetAccess = public)
        order = [];
        type = 'quad';
    end
    
    methods (Access = public)
        function obj = Gauss(order, varargin)
            obj.order = order;     
            
            if nargin >= 2;
                obj.type = varargin{1};
            end    
        end
        
        function [qp, w] = rules(obj, varargin)
            
            build_cell = false;
            if nargin >= 2 && strcmpi(varargin{1},'cell');
                build_cell = varargin{1};
            end


            switch lower(obj.type)
                case {'line','quad','hex'};
                     [qp, w] = obj.rect_rules();
                    
                case 'tri';
                     [qp, w] = obj.tri_rules();
                     
                case 'tet';
                    error('Not yet supported');
                   
                otherwise
                    error('Gauss:rules', 'The %s type of quadrature is not supported.', obj.type);
            end
            
            
            if build_cell;
                [qp, w] = obj.cell_rules(qp,w);
            end
            
            
            
            
            
        end
    end
    
    methods (Access = private)
        function [qp, w] = rect_rules(obj)
            
            switch obj.order
                case 1;
                    qp = 0;
                    w = 2;
                    
                case 2;
                    qp = [-1/sqrt(3), 1/sqrt(3)]';
                    w = [1, 1]';
                    
                case 3;
                    qp = [0, -sqrt(3/5), sqrt(3/5)]';
                    w = [8/9, 5/9, 5/9]';
                    
                case 4;
                    a = sqrt((3 - 2*sqrt(6/5))/7);
                    b = sqrt((3 + 2*sqrt(6/5))/7);
                    qp = [-a, a, -b, b]';
                    
                    a = (18+sqrt(30))/36;
                    b = (18-sqrt(30))/36;
                    w = [a, a, b, b]';
                    
                case 5;
                    a = 1/3*sqrt(5 - 2*sqrt(10/7));
                    b = 1/3*sqrt(5 + 2*sqrt(10/7));
                    qp = [0, -a, a, -b, b]';
                    
                    a = (322+13*sqrt(70))/900;
                    b = (322-13*sqrt(70))/900;
                    w = [128/225, a, a, b, b]';
                    
                otherwise
                    error('Gauss:rect_rules', 'The specified order of %d is not supported.', obj.order);
            end
        end   
        
        function [qp, w] = tri_rules(obj)
            % Gauss quadrature rules for triangle domains 
            % (see Fish, p. 181, Table 7.7) 
            switch obj.order         
                case 1;
                    a  = 1/3;
                    qp = [a,a];
                    w  = [1,1];
                
                case 3;
                    a = 0.166666666;
                    b = 0.666666666;
                    qp = [a,a; b,a; a,b];
                    w = [a, a, a];
                    
                case 4;
                    a = 1/3;
                    b = 0.6;
                    c = 0.2;
                    qp = [a,a; b,c; c,b; c,c];
                    
                case 7;
                    a = 0.1012865073;
                    b = 0.7974269853;
                    c = 0.4707420641;
                    d = 0.0597158717;
                    e = 1/3;
                    qp = [a,a; b,a; a,b; c,d; c,c; d,c; e,e];
                    
                    a = 0.0629695903;
                    b = 0.0661970764;
                    c = 0.1125;
                    w = [a,a,a,b,b,b,c];

                otherwise
                    error('Gauss:rect_rules', 'The specified order of %d is not supported.', obj.order);
            end
        end
        
        function [qp, w] = cell_rules(obj,qp,w)
            
            if any(strcmpi(obj.type,{'line','tet','tri'}));
            	qp = num2cell(qp);
                
            else
                idx = 1:length(qp);
                out{1} = repmat(idx', length(idx), 1);
                out{2} = sort(repmat(idx', length(idx), 1),1);

                if strcmpi('hex',obj.type);
                   out{1} =  repmat(out{1}, n_dim, 1);
                   out{3} =  sort(repmat(out{2}, n_dim, 1),1); 
                   out{2} =  repmat(out{2}, n_dim, 1); 
                end 

                idx = cell2mat(out);                
                qp = num2cell(qp(idx));
                
                w = prod(w(idx),2);
                

            end   
            
            
            
        end
    end
end