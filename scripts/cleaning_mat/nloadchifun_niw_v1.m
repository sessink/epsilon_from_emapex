   function Chi = loadchifun(float);

   Float = str2mat('7779a','7781a','7783a','7786a','7787a','7788a', ...
                   '7700b','7701b','7780b','7784b','7785b','7786b');
   jf = myfindstr(Float,float);
   load_good_chi_period_niw
   load(['/n/lien/NIW/matdata/' float '_grid']);

   dtdz1 = sqrt(0.5*A.chi1./A.kT1);
   chi1 = A.chi1; kT1 = A.kT1; eps1 = A.eps1;
   smchi1 = A.smchi1; smkT1 = A.smkT1; smeps1 = A.smeps1;

   good_chi_period1 = good_chi1_period(jf,:);

   bad = find(dtdz1 <= 1.5e-3 | chi1 >= 5.0e-5 | kT1 >= 1.0e-1 );
   dtdz1(bad) = NaN; 
   chi1(bad) = NaN; kT1(bad) = NaN; eps1(bad) = NaN;
   smchi1(bad) = NaN; smkT1(bad) = NaN; smeps1(bad) = NaN;
   bad = find(A.Jday_gmt < good_chi_period1(1) | A.Jday_gmt > good_chi_period1(2));
   chi1(:,bad) = NaN; kT1(:,bad) = NaN; eps1(:,bad) = NaN;
   smchi1(:,bad) = NaN; smkT1(:,bad) = NaN; smeps1(:,bad) = NaN;

   dtdz2 = sqrt(0.5*A.chi2./A.kT2);
   chi2 = A.chi2; kT2 = A.kT2; eps2 = A.eps2;
   smchi2 = A.smchi2; smkT2 = A.smkT2; smeps2 = A.smeps2;
   good_chi_period2 = good_chi2_period(jf,:);
   bad = find(dtdz2 <= 1.5e-3 | chi2 >= 5.0e-5 | kT2 >= 1.0e-1 );
   chi2(bad) = NaN; kT2(bad) = NaN; eps2(bad) = NaN;
   smchi2(bad) = NaN; smkT2(bad) = NaN; smeps2(bad) = NaN;
   bad = find(A.Jday_gmt < good_chi_period2(1) | A.Jday_gmt > good_chi_period2(2));
   chi2(:,bad) = NaN; kT2(:,bad) = NaN; eps2(:,bad) = NaN;
   smchi2(:,bad) = NaN; smkT2(:,bad) = NaN; smeps2(:,bad) = NaN;
%
%     Comparing two sensors
%
   [nchi1,nchi2]=consistentchkfun(chi1,chi2,0.5,2);
   [nkT1,nkT2]=consistentchkfun(kT1,kT2,0.5,2);
   [neps1,neps2]=consistentchkfun(eps1,eps2,0.5,2);

   nchi=avgtwofun(nchi1,nchi2); nkT=avgtwofun(nkT1,nkT2); neps=avgtwofun(eps1,neps2);
   Chi.chi = nchi; Chi.kT = nkT; Chi.eps = neps;

   [nsmchi1,nsmchi2]=consistentchkfun(smchi1,smchi2,0.5,2);
   [nsmkT1,nsmkT2]=consistentchkfun(smkT1,smkT2,0.5,2);
   [nsmeps1,nsmeps2]=consistentchkfun(smeps1,smeps2,0.5,2);

   nsmchi=avgtwofun(nsmchi1,nsmchi2); nsmkT=avgtwofun(nsmkT1,nsmkT2); 
   nsmeps=avgtwofun(nsmeps1,nsmeps2);
   Chi.smchi = nsmchi; Chi.smkT = nsmkT; Chi.smeps = nsmeps;

   Chi.Jday_gmt = A.Jday_gmt; Chi.Pr = A.Pr;
   Chi.chi1 = nchi1; Chi.kT1 = nkT1; Chi.eps1 = neps1;
   Chi.chi2 = nchi2; Chi.kT2 = nkT2; Chi.eps2 = neps2;
   Chi.smchi1 = nsmchi1; Chi.smkT1 = nsmkT1; Chi.smeps1 = nsmeps1;
   Chi.smchi2 = nsmchi2; Chi.smkT2 = nsmkT2; Chi.smeps2 = nsmeps2;

   for i = 1:length(Chi.Jday_gmt);
       j = max(find(~isnan(Chi.kT1(:,i))));
       Chi.maxPr(i) = NaN;
       if ~isempty(j);
          Chi.maxPr(i) = Chi.Pr(j);
       end
   end
