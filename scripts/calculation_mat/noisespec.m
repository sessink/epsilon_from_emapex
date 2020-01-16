function [output] = noisespec(f_cps)
    output = 1.0e-11 * [1+(f_cps/20).^3].^2;
end