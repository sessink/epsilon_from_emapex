%
%     H2 Fp07 transfer function
%
%    function H2fp07 = H2FP07fun(Hz);
%
%    Hz is frequency in Hz
%    U is velocity in m/s
%
     function H2fp07 = H2FP07fun(Hz,U);
     tau = 0.006*U.^(-1/2);
     H2fp07 = ( 1 + (2*pi*Hz*tau).^2).^(-1);
