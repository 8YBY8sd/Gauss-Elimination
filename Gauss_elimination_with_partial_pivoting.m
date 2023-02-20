function x = Gauss_elimination_with_partial_pivoting(A,b)
    n = size(A,1);
    L  = eye(n);
    
    for i = 1:n
        if A(i,i) == 0
            
            error('problem divide by zero')
            
        end
        
        max = max_entry(A,i);
        A = swap(A,i,max);
        b = swap(b,i,max);

        for j = i+1 : n 
            
            L(j,i) = A(j,i)/A(i,i);
            
            for k = i+1: n
                
                A(j,k) = A(j,k) - L(j,i)*A(i,k);
                
            end
            
            A(j,i) = 0; 
            
            b(j) = b(j) - L(j,i)*b(i);
        end
                
        
    end
    % this is the function from HW2 that solves upper triangular systems
    x = uptriangsolve(A,b);
end


