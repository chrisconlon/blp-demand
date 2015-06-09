function [fval,gradstar]=finite_difference(func,theta0,deps)
fval = func(theta0);
gradstar = zeros(length(theta0),length(fval));

for k=1:length(theta0)
  thetastarA = theta0;
  thetastarB = theta0;

  thetastarA(k) = thetastarA(k) + deps;
  thetastarB(k) = thetastarB(k) - deps;
  gradstar(k,:) = (func(thetastarA) - func(thetastarB))./(2*deps);
end

end