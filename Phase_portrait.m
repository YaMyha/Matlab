function f=Phase_portrait(t, x)
f(1) = x(1) + 0.4*x(2);
f(2) = 0.5*x(1) + 2*x(2) - 2.6*x(1)^2 + 0.3*x(2)^2 + 0.6*x(3)^2;
f(3) = -0.3*x(1) + 2*x(3);
f=f';