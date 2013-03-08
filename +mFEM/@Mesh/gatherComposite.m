function out = gatherComposite(obj,varargin)
    %GATHERCOMPOSITE (protected) Returns selection of nodes or elements
    %   The gatherComposite method is a generic method for extracting data
    %   stored as a MATLAB Composite, which is the case for the node and
    %   elements objects for the mFEM.Mesh. The nodes and elements are
    %   stored as a Compsite because MATLAB currently does not support
    %   codistributed vectors of custom classes, this function allows
    %   vectors of nodes and elements to be extracted.
    %
    % Syntax
    %   no = getComposite()
    %   no = getComposite(id)
    %   no = getComposite(...,'PropertyName',PropertyValue,...)
    %
    % Description
    %   no = getComposite() returns the node or element objects for the 
    %   entire mesh.
    %
    %   no = getComposite(id) returns the objects corresponding to the
    %   global ids supplied in id.
    %
    %   no = getComposite(...,'PropertyName',PropertyValue,...) returns the 
    %   objects as above but limits the selection further based on the
    %   property parings supplied, which are detailed below.
    %
    % GETNODES Property Descriptions
    %   Name (Required)
    %       'elements' | 'nodes'
    %       Triggers which type of data to return.
    %
    %   Tag
    %       scalar | char
    %       Limits the objects returned to only those with the
    %       supplied tag, the term tag is applied to both those added with
    %       the addBoundary and addSubdomain commands.
    %
    %   Lab
    %       scalar | vector
    %       Limits the objects returned to the given processors in
    %       parallel applications
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

    % Gather the input and extract the desired objects
    [comp,map,id,all_ids,opt] = ...
        parseGatherCompositeInput(obj, varargin{:});

    % Serial case is identical to specified single lab case 
    if matlabpool('size') == 0 && isempty(opt.lab);
        opt.lab = 1;
    end
 
    % If lab is a scalar, get the values and be done
    if isscalar(opt.lab);
        out = obj.elements{opt.lab};
        if ~all_ids;
            out = out(id);
        end
        return %ends operation of this method
    end
              
    % Re-assign lab, dot ref. ont allowed in spmd and account for empty
    lab = opt.lab;
    if isempty(lab);
        lab = 1:matlabpool('size');
    end
    
    % Only perform parallel operations if paritial list is desired
    if ~all_ids;
        
        % Define variables for use in spmd
        tag = opt.tag;
        all_tag = obj.tag;
        
        spmd % begin parallel
            
            % Limit operations to the specified labs
            if any(labindex == lab);
                
                % IDs of local values
                local_id = globalIndices(map,1); 
                
                % Desired values on this processor
                idx = id >= min(local_id) & id <= max(local_id); 
                
                % Limit to defined tags
                if ~isempty(tag);
                    for i = 1:length(tag);
                        c = strcmp(tag,all_tag); 
                        idx(:,i+1) = map(:,c);
                    end
                    idx = all(idx,2);
                end
                
                % Build a new limited Composite
                comp = comp(idx); 
            end
        end
    end
    
    % Build the output vector
    out = comp{lab(1)};
    for i = 2:length(lab);
        out = [out;comp{lab(i)}];
    end
end

function [comp,map,id,all_ids,opt] = ...
    parseGatherCompositeInput(obj, varargin)
    %PARSEGATHERCOMPOSITEINPUT Prepares input for use by main function
    %   The gatherComposite method is a generic method for extracting data
    %   stored as a MATLAB Composite, which is the case for the node and
    %   elements objects for the mFEM.Mesh. This function extracts the
    %   desired information depending on the 'name' property.
    
    % No input given by user, get all nodes
    if nargin == 1;
        id = [];
        properties = {};
    
    % Ids without properties given by user    
    elseif nargin == 2 && isnumeric(varargin{1});
        id = varargin{1};
        properties = {};
        
    % Ids with properties given by the user    
    elseif nargin >= 3 && isnumeric(varargin{1});
        id = varargin{1};
        properties = varargin(2:end);
        
    % Only properties supplied    
    else
        id = [];
        properties = varargin;
    end
    
    % Gather the properties
    opt.tag = {};
    opt.lab = [];
    opt.name = '';
    opt = gatherUserOptions(opt,properties{:});
    
    % Extract the data to work with
    switch lower(opt.name)
        case {'elem','element','elements'};
            map = obj.elem_tag_map;
            comp = obj.elements;
            n = obj.n_elements;
        case {'nodes','node'};
            map = obj.node_tag_map;
            comp = obj.nodes;
            n = obj.n_nodes;
        otherwise
            error('Mesh:gatherComposite:InvalidName','The name input must be either ''elements'' or ''nodes''.');
    end
    
    % Set the ids for case empty id cases
    if isempty(id);
        id = 1:n;
    end
        
    % Determine if all of the elements are requested
    all_ids = false;
    if islogical(id) && sum(id) == n;
        all_ids = true;
    elseif length(id) == n;
        all_ids = true;
    end
end