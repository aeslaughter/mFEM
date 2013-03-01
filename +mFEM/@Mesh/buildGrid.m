function [pts,elms] = buildGrid(order, c, varargin)

    n = nargin-2;
    
    elms = Composite();
    spmd
        
%         E.n_elem = nx*ny;
%         E.elem_map = zeros(E.n_elem,4);
% 
%         n_nodes = length(nodes);
%         id = reshape(1:n_nodes,nx+1,ny+1);
% 
%         k = 0;
%         for j = 1:ny;
%             for i = 1:nx;
%                 k = k + 1;
%                 E.elem_map(k,[1,2,4,3]) = reshape(id(i:i+1,j:j+1),4,1);
%             end
%         end
        
        % One Dimension
        if n == 3;
            x0 = varargin{1}; x1 = varargin{2}; xn = order*varargin{3};
            pts = codistributed(x0:(x1-x0)/xn:x1);

        % Two Dimension
        elseif n == 6;
            x0 = varargin{1}; x1 = varargin{2}; xn = order*varargin{5};
            y0 = varargin{3}; y1 = varargin{4}; yn = order*varargin{6};
            x = codistributed(x0:(x1-x0)/xn:x1);
            y = codistributed(y0:(y1-y0)/yn:y1);
            [X,Y] = ndgrid(x,y);
            
            nx = length(x);
            ny = length(y);
            id = reshape(uint32(1:nx*ny),nx,ny);
            N  = false(size(id));
            k = 0;
            
            g_size = size(id);
            blk_size = [order+1,order+1];
            
            n_elem = varargin{5}*varargin{6};
            elms = zeros(n_elem, size(c,1),'uint32');
            
            
            skipped = [];
            for j = 1:order:ny-order;
                for i = 1:order:nx-order;
                    k = k + 1;
                    [ii,jj] = ndgrid(i:i+order,j:j+order);
                    g_ind = sub2ind(g_size,reshape(ii,numel(ii),1),reshape(jj,numel(jj),1));
                    ind = sub2ind(blk_size,c(:,1),c(:,2));
                    elms(k,:) = id(g_ind(ind));
                    N(g_ind(ind)) = true;
                    
                    mem = ~ismember(g_ind, g_ind(ind));
                    skipped = [skipped,g_ind(mem)];
                end
            end
            
            for s = 1:length(skipped)
                map = elms>skipped(s);
                elms(map) = elms(map) - 1;
            end
            pts = [reshape(X(N),numel(X(N)),1), reshape(Y(N),numel(Y(N)),1)];

        % Three Dimension
        elseif n == 9;
            x0 = varargin{1}; x1 = varargin{2}; xn = order*varargin{7};
            y0 = varargin{3}; y1 = varargin{4}; yn = order*varargin{8};
            z0 = varargin{5}; z1 = varargin{6}; zn = order*varargin{9};          
            x = codistributed(x0:(x1-x0)/xn:x1);
            y = codistributed(y0:(y1-y0)/yn:y1);
            z = codistributed(z0:(z1-z0)/zn:z1);
            [X,Y,Z] = ndgrid(x,y,z);
            pts = [reshape(X,numel(X),1),...
                   reshape(Y,numel(Y),1),...
                   reshape(Z,numel(Z),1)];
        else
            error('buildGrid:InputError','The number of inputs is invalid.');
        end  
    end
end

