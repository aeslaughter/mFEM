function D = transform_dof(d, n)
    %TRANSFROM_DOF Converts the dofs for vector element space
    %
    % Syntax
    %   D = transform_dof(d,n)
    %   
    % Description
    %   D = transform_dof(d,n) converts the scalar degrees of freedom
    %       for to vector based degrees of freedom. For example,
    %       inputing d = [1,3], n = 2 returns D = [1,2,5,6].
    %
    %----------------------------------------------------------------------
    % Copyright 2012 Andrew E. Slaughter
    % This software is for educational purposes only and may not be used
    % without written permession.
    %----------------------------------------------------------------------
    
    % Size of vector space dofs
    D = zeros(n*length(d),1); 

    % Loop through dimensions and build vector
    for i = 1:n;
        D(i:n:end) = d*n - (n-i);
    end 