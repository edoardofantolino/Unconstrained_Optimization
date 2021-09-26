close all
clear
clc

format short

f = @(x) (1/4)*x(1,:).^4 + (1/2)*x(1,:).^2 + x(1,:) ...
        + (1/4)*x(2,:).^4 + (1/2)*x(2,:).^2 + x(2,:);
gradf = @(x) [x(1,:).^3 + x(1,:) + 1;
              x(2,:).^3 + x(2,:) + 1];
Hessf = @(x) [3*x(1,:).^2 + 1, 0;
              0, 3*x(2,:).^2 + 1];

limi = 3;
x = -limi:.1:limi;
y = x;
[X,Y] = meshgrid(x,y);
Z = (1/4)*X.^4 + (1/2)*X.^2 + X ...
    + (1/4)*Y.^4 + (1/2)*Y.^2 + Y;

x0 = [2*sqrt(2); 0];
x=x0;
teta = pi/64;
G = [cos(teta), -sin(teta); sin(teta) cos(teta)];
for i=1:((2*pi)/teta)-1
    x = G*x;
    x0=[x0, x];
end
s=size(x0)

kmax=50;
tollgrad=1e-8;

[xkg, fkg, gradfk_normg, kg, xseqg] = newton([2;2], f, gradf, Hessf, 1, kmax, tollgrad);
    xseqg = [[2;2], xseqg];


for k=1:3
    for i=1:s(2)
        [xkn, fkn, gradfk_normn, kn, xseqn] = newton(x0(:,i), f, gradf, Hessf, 1, kmax, tollgrad);
        xseqn = [x0(:,i), xseqn];

        figure(1)
        axis equal
        contour(X, Y, Z, 'k')
        hold on
        grid on
        xlabel("x_{1}")
        ylabel("x_{2}")
        
        if x0(1,i) < 1e-14 && x0(1,i) > -1e-14
            zerovalue=0;
            title(['Starting point: (' num2str(zerovalue), ',', num2str(x0(2,i)) , '), num Iter:', num2str(kn)])
        else
            if x0(2,i) < 1e-14 && x0(2,i) > -1e-14
                zerovalue=0;
                title(['Starting point: (' num2str(x0(1,i)), ',', num2str(zerovalue) , '), num Iter:', num2str(kn)])
            else
                title(['Starting point: (' num2str(x0(1,i)), ',', num2str(x0(2,i)) , '), num Iter:', num2str(kn)])
            end
        end
        
        
        plot(xseqn(1,:), xseqn(2,:), 'o--k')
        hold on
        plot(xseqg(1,:), xseqg(2,:), 'o--r')
        hold off
        pause(0.05)
    end
end