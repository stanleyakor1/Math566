function [iter] = iterative_method(A,b,tol,method)
x_init = zeros(length(b),1);
iter = 0;
switch method
    
    case 'jacobi'
        D = diag(A);
        while (true)
            zk = (b - A*x_init)./D;
            x = x_init+ zk;
            error = norm(zk,2)*norm(b,2);
            iter = iter + 1;
            if error < tol
                break;
            end
            x_init = x;
        end

    case 'gauss_seidel'

        n = length(b);
        x0 = zeros(length(b),1);
        x =zeros(length(b),1);

        while (true)
              for i = 1:n
                x(i) = (b(i) - A(i,1:i-1) * x(1:i-1) - A(i,i+1:n) * x0(i+1:n)) / A(i,i);
              end
            iter = iter + 1;
            if norm(x - x0)*norm(b) < tol
                break;
            end
            x0 = x;   
        end
end
end