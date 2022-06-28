clc;
clear all;
close all;

x0 = 4;
y0 = 4;
h =  0.01;
x = h:h:x0-h;
y = h:h:y0-h;
b1=0;
b2=100;

n = length(x)+2;
R = 'S';
G = numgrid(R,n);
figure(1);
spy(G);
title('A Finite Difference Grid');


n1 = sum(G(:)>0);

figure(2);
A = delsq(G);
spy(A);
title('The 5-Point Laplacian');

b = sparse(n1,1);
for i=n-2:n-2:n1
       b(i) =  b2;
end

solution = A\b;
c = reshape(solution,[n-2,n-2]);
figure(3);
[X,Y] = meshgrid(x,y);
%surf(X,Y,c);
mesh(X,Y,c);
title("Ö(x,y)");

f = 0;
for i=1:2:101
   f = f + (4*b2/(i*pi*sinh((i*pi*y0)/x0))).*sin((i*pi*X)/x0).*sinh((i*pi*Y)/x0); 
end

figure(4);
[X,Y] = meshgrid(x,y);
surf(X,Y,f);
%mesh(X,Y,f);
title("Ö(x,y)");

figure(5);
U = G;
U(G>0) = full(f(G(G>0)));
clabel(contour(U));
prism
axis square ij