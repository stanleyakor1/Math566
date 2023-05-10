function x = fast_pentadiag(a,b,c,d,e,f)
% STEP 1
n = length(f);
x = zeros(n,1);
B = zeros(n-2,1);
K = zeros(n-2,1);
for i = 1:n-2
    B(i) = b(i)/a(i);
    K(i) = d(i)/a(i);
    a(i+1) = a(i+1) - B(i)*c(i);
    c(i+1) = c(i+1) - B(i)*e(i);
    b(i+1) = b(i+1) - K(i)*c(i);
    a(i+2) = a(i+2) - K(i)*e(i);
    f(i+1) = f(i+1) - B(i)*f(i);
    f(i+2) = f(i+2) - K(i)*f(i);
end
a(n) = a(n) - (b(n - 1)/a(n - 1))*c(n-1);
f(n) = f(n) - (b(n - 1)/a(n - 1))*f(n - 1);
x(n) = f(n)/a(n);
x(n - 1) = (f(n - 1) - x(n)*c(n - 1))/a(n - 1);
%STEP 2
for i = n-2:-1:1
    x(i) = (f(i) - x(i+1)*c(i) - x(i+2)*e(i))/a(i);
end
end