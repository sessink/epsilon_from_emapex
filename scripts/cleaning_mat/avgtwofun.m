
   function avgx = avgtwofun(x1,x2);

   avgx = NaN*ones(size(x1));

   g = find(~isnan(x1.*x2));
   avgx(g) = (x1(g) + x2(g))/2;

   g1 = find(~isnan(x1) & isnan(x2));
   avgx(g1) = x1(g1);
   g2 = find(isnan(x1) & ~isnan(x2));
   avgx(g2) = x2(g2);
   
   
