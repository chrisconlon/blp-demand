    function print_results(results)
        kk=length(results.theta);
        params =[results.theta;results.beta(1:3)];
        sval = results.SE(1:kk+3);
        disp('Results ')
        disp('---------------- ')
        disp('Parameter Estimates ')
        disp(num2str([params sval]))
        disp(['GMM Objective : ' num2str(results.fval)])
        disp(['Residual Variance : ' num2str(var(results.resid))])

    end