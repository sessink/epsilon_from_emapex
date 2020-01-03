   % chi_data_dir = [niskinehomedir 'data/emapex/dec/']; 
   chi_data_dir = [niskinehomedir 'data/emapex/dec/']; 
   floatid = 
   hpid =
   kzmin = 20; kzmax = 400; plotting = 1; threshold = 4;
   plotting = 1; 
   Chi = chiprofile_fun(floatid,hpid,chi_data_dir,kzmin,kzmax,plotting,threshold);
      
   if ~isempty(Chi) ;
      npr = length(Chi.P);
      CHI.chi1(1:npr,j) = Chi.chi1(:); CHI.chi2(1:npr,j) = Chi.chi2(:);
      CHI.kT1(1:npr,j) = Chi.kT1(:); CHI.kT2(1:npr,j) = Chi.kT2(:);
      CHI.eps1(1:npr,j) = Chi.eps1(:); CHI.eps2(1:npr,j) = Chi.eps2(:);
      CHI.T(1:npr,j) = Chi.T(:);
      CHI.P(1:npr,j) = Chi.P(:);
      CHI.jday_gmt(1:npr,j) = Chi.jday_gmt(:);
      bad = find(CHI.chi1 == 0);
      CHI.chi1(bad) = NaN; CHI.chi2(bad) = NaN; CHI.kT1(bad) = NaN;
      CHI.kT2(bad) = NaN; CHI.eps1(bad) = NaN; CHI.eps2(bad) = NaN;
      CHI.P(bad) = NaN;
   else
      CHI.chi1(1,j) = NaN; CHI.chi2(1,j) = NaN;
      CHI.eps1(1,j) = NaN; CHI.eps2(1,j) = NaN;
      CHI.kT1(1,j) = NaN; CHI.kT2(1,j) = NaN;
      CHI.P(1,j) = NaN; CHI.jday_gmt(1,j) = NaN;
      CHI.T(1,j) = NaN;
   end
