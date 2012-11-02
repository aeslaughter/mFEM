function D = transform_dof(d, n)
    %TRANSFORM_DOF Converts the dofs for vector element space
    
    % Size of vector space dofs
    D = zeros(n*length(d),1); 

    % Loop through dimensions and build vector
    for i = 1:n;
        D(i:n:end) = d*n - (n-i);
    end 