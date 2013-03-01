function [nodes,pts] = buildGrid(varargin)
    
    n = nargin;
    spmd
        % One Dimension
        if n == 3;
            x0 = varargin{1}; x1 = varargin{2}; xn = varargin{3};
            pts = codistributed(x0:(x1-x0)/xn:x1);

        % Two Dimension
        elseif n == 6;
            x0 = varargin{1}; x1 = varargin{2}; xn = varargin{5};
            y0 = varargin{3}; y1 = varargin{4}; yn = varargin{6};
            x = codistributed(x0:(x1-x0)/xn:x1);
            y = codistributed(y0:(y1-y0)/yn:y1);
            [X,Y] = ndgrid(x,y);
            pts = [reshape(X,numel(X),1), reshape(Y,numel(Y),1)];

        % Three Dimension
        elseif n == 9;
            x0 = varargin{1}; x1 = varargin{2}; xn = varargin{7};
            y0 = varargin{3}; y1 = varargin{4}; yn = varargin{8};
            z0 = varargin{5}; z1 = varargin{6}; zn = varargin{9};          
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

        local_pts = getLocalPart(pts);
        local_pts_id = globalIndices(pts,1);
        n_nodes = size(local_pts,1);
        local_nodes = cell(n_nodes,1);
        for i = 1:n_nodes;
            local_nodes{i} = mFEM.elements.base.Node(local_pts_id(i),local_pts(i,:));
        end
        
        codist = getCodistributor(pts);
        if numlabs == 1;
            part = size(pts,1);
        else
            part = codist.Partition;
        end
        gsize = [size(pts,1),1];
        codist = codistributor1d(1, part, gsize);
        nodes = codistributed.build(local_nodes, codist);
    end
end

