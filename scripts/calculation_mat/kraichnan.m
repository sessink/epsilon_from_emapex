function [output] = kraichnan(k_rpm,chi,kb_rpm)
% Kraichnan temperature gradient spectrum
%
% S.Essink, December 2019

% Parameters
qk = 5.27;
D = 1.4e-7;

yk = sqrt(qk) * k_rpm ./ kb_rpm;
nom = chi.* sqrt(qk) .* yk .* exp(-sqrt(6) * yk);
denom = (D * kb_rpm);

output = nom./denom;

end