function elem = getElement(obj, varargin)
    %GETELEMENT returns the desired elements
    %
    % Syntax
    %   elem = getElement(id)
    %   elem = getElement(id,lab)
    %
    % Description
    %   elem = getElement(id) returns an array of elements with the ids
    %   given the numeric vector id.
    %
    %   elem = getElement(id,lab) operates the same as above but limits the
    %   search to the specified processor given in lab.
    %
    %----------------------------------------------------------------------
    %  mFEM: An Object-Oriented MATLAB Finite Element Library
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
    %  along with this program.  If not, see <http://www.gnu.org/licenses/>.
    %
    %  Contact: Andrew E Slaughter (andrew.e.slaughter@gmail.com)
    %----------------------------------------------------------------------

    % Parse input
    if nargin == 1;
        id = 1:obj.n_elements;
        flag = true;
    else
        id = varargin{1};
        flag = false;
    end

    % If serial returning the elements is easy
    if matlabpool('size') == 0;
        elem = obj.elements{1};
        if ~flag;
            elem = elem(id);
        end
        return;
    end
    
    % Determine the value of the lab input
    if nargin == 3;
        lab = varargin{2};
    else
        lab = [];
    end
    
    % Begin parallel operations
    elem = gatherComposite(flag,id,obj.elem_map,obj.elements,lab);
end

function out = gatherComposite(flag,id,map,comp,lab)

    if ~flag;
        spmd
            % IDs of the elements stored on this lab
            local_id = globalIndices(map,1);

            % local ids for the desired elements
            idx = id >= min(local_id) & id <= max(local_id);

            % Build a Composite
            comp_idx = comp(idx);
        end
    else
        comp_idx = comp;
    end
    
    if isempty(lab)
        out = comp_idx{1};
        for i = 2:length(comp_idx);
            out = [out; comp_idx{i}];
        end   
    else
        out = comp_idx{lab};
    end
end


