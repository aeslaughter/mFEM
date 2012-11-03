classdef Gauss
    %GAUSS A class for Gauss quadrature points and weight functions.
    % Includes rectanglar and triangular quadrature rules and weight
    % functions for use in finite element calculations. An instance of this
    % class should be attached to each Element that provides the correct
    % rules for the element.
    %
    %----------------------------------------------------------------------
    % Copyright 2012 Andrew E. Slaughter
    % This software is for educational purposes only and may not be used
    % without written permession.
    %----------------------------------------------------------------------
    
    properties (SetAccess = private, GetAccess = public)
        opt = ...   % Structure of available options
            struct('order', 1, 'type', 'line');
    end
    
    methods (Access = public)
        function obj = Gauss(varargin)
            %GAUSS Class constructor
            %
            % Syntax
            %   Gauss('PropertyName',PropertyValue)
            %
            % Description
            %   Gauss('PropertyName',PropertyValue) creates an instance of 
            %   the Gauss class for the, the available properties are
            %   discussed below.
            %
            % GAUSS Property Description
            %   Order
            %       integer
            %       Specifies the number of Gauss quadrature points to use
            %       for the analysis, the default is 1. See the list of 
            %       available orders for each type below.
            %
            %   Type
            %       {'line'} | 'quad' | 'hex' | 'tri' | 'tet'
            %       Five different quadrature types are availble, the
            %       'line', 'quad', and 'hex' types are identical unless 
            %       the '-cell' option is used in the rules method, see the
            %       documentation for rules for more details.
            %
            % List of Available Gauss Quadrature Rules for Each Type
            %   
            %   Type                    Order
            %   'line', 'quad', 'hex'   1,2,3,4,5
            %   'tri'                   1,3,4,7

            % Gather/set the options
            obj.opt = gather_user_options(obj.opt, varargin{:});   
        end
        
        function [qp, w] = rules(obj, varargin)
            %RULES Returns the quadrature rules.
            %
            % Syntax
            %   [qp,w] = rules()
            %   [qp,w] = rules('-cell')
            %
            % Description
            %   [qp,w] = rules() returns the quadrature rules based on the
            %   settings established when the object was created. For the
            %   'line', 'quad', and 'hex' types it returns a vector for
            %   quadrature points (qp) and weights (w). In multidimensional
            %   cases these points are looped over multiple times, for
            %   example:
            %       [qp,w] = gauss_object.rules()
            %       for i = 1:length(qp);
            %           for j = 1:length(qp);
            %               Ke = Ke + w(i)*w(j)*B(qp(i),qp(j))'*k*B(qp(i),qp(j))*elem.detJ(qp(i),qp(j));
            %           end
            %       end
            %   For the 'tri' and 'tet' types, the qp is an array where the
            %   rows are the differenct quadrature points and the columns
            %   are contain the various components, based on the dimension.
            %   For example, the above example with the 'tri' type would be:
            %       [qp,w] = gauss_object.rules()
            %       for i = 1:size(qp,1);
            %           Ke = Ke + w(i)*B(qp(i,1),qp(i,2))'*k*B(qp(i,1),qp(i,2))*elem.detJ(qp(i,1),qp(i,2));
            %       end
            %
            %   [qp,w] = rules('-cell') this is the same as above, but
            %   outputs the qp vector or matrix in a consistent manner
            %   between the two types of output discussed above. It also
            %   outputs qp as a cell structure for dimension independant
            %   application. For example, the above examples would both
            %   work regardless of the type if this flag is used as
            %   follows:
            %       [qp,w] = gauss_object.rules(-cell)
            %       for i = 1:size(qp,1)
            %           Ke = Ke + w(i)*B(qp{i,:})'*k*B(qp{i,:})*elem.detJ(qp{i,:});
            %       end
            %   The '-cell' flag sets the 'cell' property, see the
            %   descrition below for an alternative syntax.
            %
            % RULES Property Description
            %   cell
            %       true | {false}
            %       Toggles the use of the cell type output described
            %       above. As shown in the available syntax the flag type
            %       specification is available. The following two methods
            %       for specifying this property are identical.
            %           [qp,w] = rules('-cell') or
            %           [qp,w] = rules('cell',true)
            
            build_cell = false;
            if nargin >= 2 && strcmpi(varargin{1},'cell');
                build_cell = varargin{1};
            end


            switch lower(obj.type)
                case {'line','quad','hex'};
                     [qp, w] = obj.rect_rules(obj.order);
                    
                case 'tri';
                     [qp, w] = obj.tri_rules(obj.order);
                     
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
    
    methods (Static, Access = private)
        function [qp, w] = rect_rules(order)
            %RECT_RULES Returns the linear quadrature points and weights
            %
            % Syntax
            %   [qp,w] = rect_rules(order)
            %
            % Description
            %   [qp,w] = rect_rules(order) returns the linear quadrature 
            %   points and weights for the specified order
 
            % Specify the qp and w vectors based on the specified order
            switch order
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
        
        function [qp, w] = tri_rules(order)
            %TRI_RULES Returns the triangular quadrature points and weights
            %
            % Syntax
            %   [qp,w] = tri_rules(order)
            %
            % Description
            %   [qp,w] = rect_rules(order) returns the triangular quadrature 
            %   points and weights for the specified order
            
            % Specify the qp and w vectors based on the specified order
            switch order         
                case 1; 
                    a  = 1/3;
                    qp = [a,a];
                    w  = [1,1];
                
                case 3;
                    a = 1/6;
                    b = 2/3;
                    qp = [a,a; b,a; a,b];
                    w = [a, a, a];
                    
                case 4;
                    a = 1/3;
                    b = 3/5;
                    c = 1/5;
                    qp = [a,a; b,c; c,b; c,c];
                    
                    a = -27/48;
                    b = 25/48;
                    w = [a,b,b,b];
                    
                    
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
    end
        
    methods (Access = private)
        function [qp, w] = cell_rules(obj, qp, w)
            %CELL_RULES Converts the qp and w vector to the a cell array
            %
            % Syntax
            %   [qp,w] = cell_rules(qp,w)
            %
            % Description
            %   [qp,w] = cell_rules(qp,w) converts the numeric arrays for
            %   qp and w to a consistent cell based format, see the help
            %   for the rules method for details.
            
            % These types only require conversion to a cell array
            if any(strcmpi(obj.type,{'line','tet','tri'}));
            	qp = num2cell(qp);
            
            % The 'quad' and 'hex' must be repeated to build a complete set
            % Example:
            % Given the following 2 point quadrature rules:
            %   qp = [a,b] w = [w1,w2]
            % 
            % For the 'quad' type the following is output:
            %   qp = {a,a; a,b; b,a; b,b};
            %   w = [w1*w1, w1*w2, w2*w1, w2*w2];
            %
            % For the 'hex' type the following is output:
            %   qp = {a,a,a; a,a,b; a,b,b; a,b,a; b,a,a,...}
            %   w  = {w1*w1*w1, w1*w1*w2, ...}
            %
            else
                % Creates the 2D 
                idx = 1:length(qp);
                out{1} = repmat(idx', length(idx), 1);
                out{2} = sort(repmat(idx', length(idx), 1),1);

                % Build the 3D from the 2D
                if strcmpi('hex',obj.type);
                   out{1} =  repmat(out{1}, n_dim, 1);
                   out{3} =  sort(repmat(out{2}, n_dim, 1),1); 
                   out{2} =  repmat(out{2}, n_dim, 1); 
                end 

                % Convert matrices to a cell array
                idx = cell2mat(out);                
                qp = num2cell(qp(idx));
                
                % Compute the product weights
                w = prod(w(idx),2);
            end     
        end
    end
end