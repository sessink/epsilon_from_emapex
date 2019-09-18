%------------------------------------------------------------------------------
%
%        BATCHELOR              08-19-92               Ren-Chieh Lien
%
%        Batchelor temperature gradient spectrum
%
%        function batchelor_jhd(krpm, epsilon,chi,nu,D,q);
%
%        reference : 
%               Oakey, N. S., "Determination of the rate of dissipation of
%               turbulent energy from simultaneous temperature and velocity 
%               shear microstructure measurements", j.p.o., 12, 256-271, 1982.
%
%------------------------------------------------------------------------------
%    function tempsp = batchelor_jhd(krpm, epsilon,chi,nu,D,q);
     function tempsp = batchelor_jhd(krpm, epsilon,chi,nu,D,q);

% keyboard

     if ~exist('q');
	 q = 3.7;
     end
     kb = (epsilon/nu/D^2)^(1/4);
     a = sqrt(2 * q) * krpm / kb;
     for i=1:length(a);
         uppera(i) = erfc(a(i)/sqrt(2))*sqrt(pi/2);
     end
     g = 2 * pi * a .* (exp(-(a.^2) / 2) - a .* uppera);
     tempsp = sqrt(q / 2) * (chi / kb / D) * g;
%    loglog(krpm,tempsp);
