classdef Gauss
    properties (SetAccess = private, GetAccess = public)
        order = [];
        type = 'quad';
    end
    
    methods (Access = public)
        function obj = Gauss(order, varargin)
            
            if nargin == 2 && ischar(varargin{1});
                obj.type = varargin{1};
            end
            
            obj.order = order;
        end
        
        function [qp, w] = rules(obj)
            switch obj.type
                case {'linear', 'quad', 'hex'};
                    [qp, w] = obj.quad_rules();
                    
                case 'tri';
                     [qp, w] = obj.tri_rules();
                     
                case 'tet';
                    error('Not yet supported');
                   
                otherwise
                    error('Gauss:rules', 'The %s type of quadrature is not supported.', obj.type);
            end
        end
    end
    
    methods (Access = private)
        function [qp, w] = quad_rules(obj)
            
            switch obj.order
                case 1;
                    qp = 0;
                    w = 2;
                    
                case 2;
                    qp = [-1/sqrt(3), 1/sqrt(3)];
                    w = [1, 1];
                    
                case 3;
                    qp = [0, -sqrt(3/5), sqrt(3/5)];
                    w = [8/9, 5/9, 5/9];
                    
                case 4;
                    a = sqrt((3 - 2*sqrt(6/5))/7);
                    b = sqrt((3 + 2*sqrt(6/5))/7);
                    qp = [-a, a, -b, b];
                    
                    a = (18+sqrt(30))/36;
                    b = (18-sqrt(30))/36;
                    w = [a, a, b, b];
                    
                case 5;
                    a = 1/3*sqrt(5 - 2*sqrt(10/7));
                    b = 1/3*sqrt(5 + 2*sqrt(10/7));
                    qp = [0, -a, a, -b, b];
                    
                    a = (322+13*sqrt(70))/900;
                    b = (322-13*sqrt(70))/900;
                    w = [128/225, a, a, b, b];
                    
                otherwise
                    error('Gauss:rect_rules', 'The specified order of %d is not supported.', obj.order);
            end
        end   
        
        function [qp, w] = tri_rules(obj)
            % Gauss quadrature rules for triangle domains 
            % (see Fish, p. 181, Table 7.7) 
            switch obj.order                    
                case 2;
                    a = 0.166666666;
                    b = 0.666666666;
                    qp = [a,a; b,a; a,b];
                    w = [a, a, a];
                    
                case 5;
                    a = 0.1012865073;
                    b = 0.7974269853;
                    c = 0.4707420641;
                    d = 0.0597158717;
                    e = 0.3333333333;
                    qp = [a,a; b,a; a,b; c,d; c,c; d,c; e,e];
                    
                    a = 0.0629695903;
                    b = 0.0661970764;
                    c = 0.1125;
                    w = [a,a,a,b,b,b,c];

                otherwise
                    error('Gauss:rect_rules', 'The specified order of %d is not supported.', obj.order);
            end
        end
    end
end