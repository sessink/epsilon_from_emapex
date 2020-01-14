function [deriv] = lienmmderiv(x,var)
% Very crude wrapper to take derivative 
% to replicate function
   dvar = gradient(var);
   dx = gradient(x);
   deriv = dvar./dx;
end