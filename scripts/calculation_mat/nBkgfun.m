%
%    Background
%
   function [Pr,N,dTdz,W,T]= Bkgfun(ctdfn,uxt);
   B =load(ctdfn);
   Sigma = sw_pden(B.S,B.T,B.P,0);
   B.Sigma = Sigma;
   B.N2 = -9.8/1025*lienmmderiv(-B.P(:),B.Sigma(:));
   B.dTdz = lienmmderiv(-B.P(:),B.T(:));
   B.W = lienmmderiv(B.UXT(:),-B.P(:));

   [a,j]=mysort(B.UXT);
   Pr = interp1(B.UXT(j),B.P(j),uxt);
   tmp = lienmmderiv(-B.P(:),B.T(:));
   N2 = interp1(B.UXT(j),B.N2(j),uxt);
   N = sqrt(N2);
   T = interp1(B.UXT(j),B.T(j),uxt);
   dTdz = interp1(B.UXT(j),B.dTdz(j),uxt);
   W = interp1(B.UXT(j),B.W(j),uxt);
   save test B
