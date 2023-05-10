function x = fast_algorithm(u,v,b,A)

        p1 = A\b;

        p2 = A\u;

        alpha = 1 - v*p2;

        x = p1 + (1/alpha).*p2*v*p1;
end
