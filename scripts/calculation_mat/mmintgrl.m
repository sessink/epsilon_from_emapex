function [intgrl] = mmintgrl(x,var)
    intgrl = trapz(x,var);
end