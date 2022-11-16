clc
clear all
close all

global delta phi eps alpha beta gamma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta = 0.4;
phi = 0;
eps = 0;
alpha = 0;
beta = -2.6;
gamma = 0;
syms x1 x2 x3
A = [1 delta phi
     0.5 2 eps
     -0.3 0 2];
B = [0 alpha*x2 0
     beta*x1 0.3*x2 0.6*x3
     0 gamma*x2 0];
[x1_sp, x2_sp, x3_sp] = solve([A*[x1 x2 x3].' + B*[x1 x2 x3].' == [0, 0, 0]], [x1 x2 x3]);
syms y1 y2 y3 lambda
g1 = subs(A(1,:)*[y1 y2 y3].' + B(1,:)*[y1 y2 y3].',[x1 x2 x3],[y1 y2 y3]);
g2 = subs(A(2,:)*[y1 y2 y3].' + B(2,:)*[y1 y2 y3].',[x1 x2 x3],[y1 y2 y3]);
g3 = subs(A(3,:)*[y1 y2 y3].' + B(3,:)*[y1 y2 y3].',[x1 x2 x3],[y1 y2 y3]);
J = [diff(g1,y1) diff(g1,y2) diff(g1,y3)
     diff(g2,y1) diff(g2,y2) diff(g2,y3)
     diff(g3,y1) diff(g3,y2) diff(g3,y3)];
%lambda1 = solve(det(subs(J - lambda*eye(3),[y1 y2 y3],[x1_sp(1) x2_sp(1) x2_sp(1)]))==0,lambda)
%lambda2 = solve(det(subs(J - lambda*eye(3),[y1 y2 y3],[x1_sp(2) x2_sp(2) x2_sp(2)]))==0,lambda)
chareq1 = det(subs(J - lambda*eye(3),[y1 y2 y3],[x1_sp(1) x2_sp(1) x2_sp(1)]));
chareq2 = det(subs(J - lambda*eye(3),[y1 y2 y3],[x1_sp(2) x2_sp(2) x2_sp(2)]));
coeff = [[-fliplr(coeffs(chareq1))];[-fliplr(coeffs(chareq2))]];
Hur1 = [coeff(2) coeff(4) 0
       coeff(1) coeff(3) 0
       0 coeff(2) coeff(4)];
Hur2 = [coeff(6) coeff(8) 0
       coeff(5) coeff(7) 0
       0 coeff(6) coeff(8)];
for i = 1:2
    for j = 1:3
        switch i
            case 1
                M = det(Hur1(1:j,1:j))
                if M <= 0
                    disp('Неустойчивый узел 1.')
                    break
                end
            case 2
                M = det(Hur2(1:j,1:j))
                if M <= 0
                    disp('Неустойчивый узел 2.')
                    break
                end
        end
    end
end

opt = odeset('OutputSel',[1 2 3],'OutputFcn','odephas3');
X0 = [0 0 1;1 0 0;0 1 0]
hold on
for i = 1:3
    [T, X] = ode45('Phase_portrait', [0 20], X0(i,:))
    plot3(X(:,1),X(:,2),X(:,3))
end
hold off




