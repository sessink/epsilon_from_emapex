%
%     H2 ADC transfer function
%
%    function H2adc = H2ADCfun(Hz);
%
     function H2adc = H2ADCfun(Hz);
     Fc5 = 120; Fc3 = 210;   % in Hz
     sinc5 = sin(pi*Hz/Fc5)./(pi*Hz/Fc5);
     sinc3 = sin(pi*Hz/Fc3)./(pi*Hz/Fc3);
     H = (sinc5.^5).*(sinc3.^3);
     H2adc = H.^2;
