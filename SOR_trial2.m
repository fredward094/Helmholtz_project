function [ x, err, iter, flag ] = SOR_trial2(A, x, b, w, max_it, tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function computes the solution vector using SOR, if the relaxation 
%value of w=1, then this function solves using the Gauss-Seidel method, the
%inputs consist of the following
%--------------------------------------------------------------------------
%A=a square matrix representing the system of equations
%x= initial vector
%b= vector b, pertaining to the expression Ax=b
%w= relaxation value
%max_it= max number of iterations before code stops
%tol= error tolerance
%--------------------------------------------------------------------------
%The outputs are for this function are the following 
%x= solution vector
%err= the estimated error
%iter= number of iteration that occured if less than max_it
%flag= allows us to know if the system converges 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flag = 0;
    iter = 0;

    norma2_b = norm(b);
    if (norma2_b == 0.0)
        norma2_b = 1.0;
    end

    r = b - A * x;
    err = norm(r) / norma2_b;
    if (err < tol)
        return
    end

    % separate A into several matrix for SOR/Gauss-Seidel
    [ M, N, b ] = matsep(A, b, w);

    for iter = 1 : max_it
        x_1 = x;
        x = M \ (N * x + b); % adjust the aproximation
        %err = norm(x - x_1) / norm(x); % compute error
        err = norm(x_1 - x, 1); % compute error
        if (err <= tol) % check for convergence
            break
        end          
    end
    b = b / w; % vector b

    if (err > tol) % no convergence
        flag = 1;
    end   
    
end

function [ M, N, b ] = matsep(A, b, w)
%SOR process is housed in here
 b = w * b;
 M =  w * tril(A, -1) + diag(diag(A));
 N = -w * triu(A,  1) + (1.0 - w) * diag(diag(A));
end
