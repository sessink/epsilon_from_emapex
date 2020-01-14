function [intgrl] = mmintgrl(x,var)
    % very crude wrapper for a trapezoidal numerical integration 
    % to replace the custom function
    intgrl = trapz(x,var);
end