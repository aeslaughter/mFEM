classdef pMatrix
    %PMATRIX A parallel matrix class
    %   Detailed explanation goes here
    
    properties
        M = Composite();
    end
    
    methods
        function obj = pMatrix(M)
            obj.M = M;
        end
        
        function K = assemble(obj,varargin)

            opt.parallel = false;
            opt = gatherUserOptions(opt,varargin{:});
            
            M = obj.M;
            spmd
                ipart = buildCodistributedPartition(length(M.I));
                jpart = buildCodistributedPartition(length(M.J));    
                spart = buildCodistributedPartition(length(M.Aij)); 

                idist = codistributor1d(1,ipart,[sum(ipart),1]);
                jdist = codistributor1d(1,jpart,[sum(jpart),1]);  
                sdist = codistributor1d(1,jpart,[sum(spart),1]);  

                I = codistributed.build(M.I,idist);
                J = codistributed.build(M.J,jdist);
                S = codistributed.build(M.Aij,sdist);
            end
            
            K = sparse(I,J,S);
            if ~opt.parallel;
                K = gather(K);
            end
        end
        
        function M = getComposite(obj)
            M = obj.M;
        end
        
        function zero(obj)
            M = obj.M;
            spmd
                M.zero();
            end
        end
    end
end

