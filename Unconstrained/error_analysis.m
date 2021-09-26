close all
clear
clc

x_min = -0.682327803828019*ones(2,1);
format short

f = @(x) (1/4)*x(1,:).^4 + (1/2)*x(1,:).^2 + x(1,:) ...
        + (1/4)*x(2,:).^4 + (1/2)*x(2,:).^2 + x(2,:);
gradf = @(x) [x(1,:).^3 + x(1,:) + 1;
              x(2,:).^3 + x(2,:) + 1];
Hessf = @(x) [3*x(1,:).^2 + 1, 0;
              0, 3*x(2,:).^2 + 1];

x = -3:.1:3;
y = x;
[X,Y] = meshgrid(x,y);
Z = (1/4)*X.^4 + (1/2)*X.^2 + X ...
    + (1/4)*Y.^4 + (1/2)*Y.^2 + Y;

x0 = [1;2];
alpha = 1; 
kmax = 50;
tollgrad = 1e-8;


tic;
alpha0=3;
c1 = 1e-4;
rho = 0.8;
btmax = 150;
[xknb, fknb, gradfk_normnb, knb, xseqnb, btseqnb] = ...
    newton_bcktrck(x0, f, gradf, Hessf, alpha0, ...
    kmax, tollgrad, c1, rho, btmax);
xseqnb = [x0, xseqnb];
tnb = toc;


figure(1)
contour(X, Y, Z, 'k')
hold on
grid on
xlabel("x_{1}")
ylabel("x_{2}")
plot(xseqnb(1,:), xseqnb(2,:), 'o--k')
hold on
scatter(x_min(1,1), x_min(2,1), 'ro', 'filled')

iter = length(xseqnb);
error = zeros(iter,1);
for i=1:iter
   xseqnb(:,i);
   error(i) = norm(x_min-xseqnb(:,i))/norm(x_min);
end

figure(2)
semilogy(error, 'o--k')
grid on
