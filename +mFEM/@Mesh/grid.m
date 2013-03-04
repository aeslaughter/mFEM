        function grid(obj, type, varargin)
            
            % Test for valid type
            d = dir(fullfile(cd,'+mFEM','+elements','*.m'));
            c = struct2cell(d);
            if ~any(strcmp([type,'.m'],c(1,:)));
                msg = 'The %s element was not found in the +mFEM/+elements directory.';
                error('Mesh:grid:InvalidType', msg, type);
            end

            % Display wait message
            if obj.options.time;
                ticID = tMessage('Generating Grid...');
            end

            [node_map, elem_map] = feval(['mFEM.elements.',type,'.buildMaps'],varargin{:});
            if ~iscodistributed(node_map)
                spmd
                    codist = codistributor1d(1, codistributor1d.unsetPartition, size(node_map));
                    node_map = codistributed(node_map, codist);
                end
            end
            if ~iscodistributed(elem_map)
                spmd
                    codist = codistributor1d(1, codistributor1d.unsetPartition, size(elem_map));
                    elem_map = codistributed(elem_map, codist);
                end
            end
       
            obj.node_map = node_map;
            obj.elem_map = elem_map;
            
            obj.nodes = obj.buildNodes(node_map);
            obj.elements = obj.buildElements(type, obj.elem_map, obj.nodes);

            % Complete message
            if obj.options.time;
                tMessage(ticID);
            end; 