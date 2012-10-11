classdef Gauss
    properties (SetAccess = private, GetAccess = public)
        order = [];
    end
    
    methods
        function obj = Gauss(order)
            obj.order = order;
        end
        
        function [qp, w] = rules(obj)
            
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
            end    
        end
    end
end