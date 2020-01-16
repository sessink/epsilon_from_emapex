function [output] = cost_function(kb,k_rpm,chi,noise_rpm,corrdTdzsp,dof)
% Cost function for MLE to fit spectra
%
% S.Essink, December 2019

theory = kraichnan(k_rpm, chi, kb);

a = dof ./ (theory + noise_rpm);
b = corrdTdzsp .* a';

output =  -nansum(log(a)) - nansum( (dof/2.-1)*log(b)) ...
        + nansum(b/2.) + nansum(gammaln(dof/2.) + (log(2.)*dof)/2.);

end