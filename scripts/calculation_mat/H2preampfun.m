%
%     H2 Preamp transfer function
%
%    function H2preamp = H2preampfun(Hz);
%
     function H2preamp = H2preampfun(Hz);
     Fc1 = 339; Fc2 = 339;   % in Hz
     Gd = 0.965;
     H2_1 = (1 - (Hz.^2)/Fc1/Fc2 ).^2;
     H2_2 = (Hz/Fc1 + Hz/Fc2 + 2*pi*Hz*Gd).^2;
     H2_3 = ( 1 + (Hz/Fc1).^2).*( 1 + (Hz/Fc2).^2);
     H2preamp = H2_1 + H2_2./H2_3;
