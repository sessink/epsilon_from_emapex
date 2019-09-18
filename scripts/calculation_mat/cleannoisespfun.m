function Tsp = noisespfun(f,Tsp,threshold);
noisesp = 1.0e-11 * [1+(f/130).^3].^2;
noisesp = 1.0e-11 * [1+(f/20).^3].^2;
bad = find(Tsp(:)./(threshold * noisesp(:)) <= 1);
Tsp(bad) = 0;
