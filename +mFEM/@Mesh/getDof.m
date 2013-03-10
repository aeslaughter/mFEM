function dof = getDof(obj,varargin)
    %GETDOF Returns the global degrees of freedom.
    %
    % Syntax
    %   dof = getDof()
    %   dof = getDof('PropertyName',PropertyValue,...);
    %
    % Description
    %   dof = getDof() returns the global degrees of freedom for
    %   the entire finite element mesh.
    %
    %   dof = getDof('PropertyName',PropertyValue,...) returns the
    %   global degrees of freedom for portions of the mesh
    %   depending on the properties (see descriptions below). It is
    %   possible to set multiple values. For example, the following produce
    %   the same results.
    %
    %       ess(:,1) = getDof('Tag',1,'Component','x');
    %       ess(:,2) = getDof('Tag',2,'Component','y');   
    %       ess = any(ess,2);
    %
    %       ess = getDof('Tag',[1,2],'Component',{'x','y'})
    %
    % GETDOF Property Descriptions
    %   Component
    %       vector | scalar | 'x' | 'y' | 'z' | cell string
    %       Returns the dof associated with the vector space or
    %       multipe dof per node elements like the Beam element. Note, if
    %       you mesh was constructed with nodes that have varying dofs per
    %       node this property will cause errors.
    %
    %   Tag
    %       string | cell string
    %       Extract the dofs for the tags specified, where the
    %       strings given are the tags added using the addBoundary and/or
    %       addSubdoman methods.
    %
    %   Index
    %       true | {false}
    %       Toggles the type of output, by default the GETDOF method
    %       returns the dofs as a logical array equal to the length
    %       of the number of degrees of freedom. However, if the
    %       actual numeric indices are desired set this property to
    %       true. You can also use the flag style input, the
    %       following are equivlent.
    %           getDof('index',true)
    %           getDof('-index')
    %
    %   Gather
    %       true | {false}
    %       The default behavior is to return a codistributed vector, this
    %       option will perform a gather, making the vector a serial
    %       vector.
    %
    %   Composite
    %       true | {false}
    %       Changes the output to a Composite instead of the default
    %       codistributed, the gather property does not affect the output
    %       if this option is used.
    %
    %----------------------------------------------------------------------
    %  mFEM: A Parallel, Object-Oriented MATLAB Finite Element Library
    %  Copyright (C) 2013 Andrew E Slaughter
    % 
    %  This program is free software: you can redistribute it and/or modify
    %  it under the terms of the GNU General Public License as published by
    %  the Free Software Foundation, either version 3 of the License, or
    %  (at your option) any later version.
    % 
    %  This program is distributed in the hope that it will be useful,
    %  but WITHOUT ANY WARRANTY; without even the implied warranty of
    %  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %  GNU General Public License for more details.
    % 
    %  You should have received a copy of the GNU General Public License
    %  along with this program. If not, see <http://www.gnu.org/licenses/>.
    %
    %  Contact: Andrew E Slaughter (andrew.e.slaughter@gmail.com)
    %----------------------------------------------------------------------
    
    opt.tag = {};
    opt.component = [];
    opt.composite = false;
    opt.gather = false;
    opt.index = false;
    opt = gatherUserOptions(opt,varargin{:});

    node_tag_map = obj.node_tag_map;
    nodes = obj.nodes;
%     dof_part = obj.dof_partition;
    
     tag = opt.tag;
     if ischar(tag); tag = {tag}; end
     if isnumeric(tag); tag = num2cell(tag); end
     for i = 1:length(tag);
         if isnumeric(tag{i}); tag{i} = num2str(tag{i}); end
     end
     
     cmp = opt.component;
    all_tag = obj.tag;
    
    composite_flag = opt.composite;
    index_flag = opt.index;
    
    spmd
        idx = true(length(nodes),length(tag)+1);
        map = getLocalPart(node_tag_map);
        
        for i = 1:length(tag);
            c = strcmp(tag{i},all_tag);
            idx(:,i+1) = map(:,c);
        end
                
        idx = all(idx,2);
        
        if any(idx); 
            dof(:,1) = nodes(idx).getDof('Component',cmp); 
        else
            dof(:,1) = idx; 
        end
        
        % Build index form of dofs
        if ~composite_flag;
            part = buildCodistributedPartition(length(dof));
            codist = codistributor1d(1,part,[sum(part),1]);
            dof = codistributed.build(dof,codist);
        end

        % Build logical subscript array
        if ~index_flag;
            ndof = sum([nodes.n_dof]);
            part = buildCodistributedPartition(ndof);
            codist = codistributor1d(1,part,[sum(part),1]);
            subscript = codistributed.false([sum(part),1],codist);
        end
    end
    
    % Convert from index to subscript
    if ~index_flag;
        subscript(dof) = true;
        dof = subscript;
    end
    
    if matlabpool('size') == 0 || (~opt.composite && opt.gather);
        dof = gather(dof);
    end
end
    
    
    
    
