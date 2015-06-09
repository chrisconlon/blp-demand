function myDerivCheck(func,theta0,deps)
gradstar = finite_difference(func,theta0,deps);
[fval,gval] = func(theta0);
disp(['Derivative Check:']);
disp(num2str([gval gradstar]));
end