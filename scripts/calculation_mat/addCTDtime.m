   clear; clf
   Float = [7781 7783 7786 7787 7788];
   niwhomedir = '/bb/lien/NIW/'
   ctd_data_dir = ['/bb/lien/NIW/NIW2016/data/emapex/dec/'];
   beghpid = 1;
   for i = 1:length(Float);
       floatid=Float(i);
       HPID = gethpidfun(floatid);
       for j = beghpid:length(HPID);
           [HPID(j) max(HPID)]
           hpid = HPID(j);
           Ctd = ctdprofile_fun(floatid,hpid,ctd_data_dir);
           npr = length(Ctd.P);
           CTD.T(1:npr,j) = Ctd.T(:); CTD.S(1:npr,j) = Ctd.S(:);
           CTD.P(1:npr,j) = Ctd.P(:); CTD.Sigma(1:npr,j) = Ctd.Sigma(:);
           CTD.N2(1:npr,j) = Ctd.N2(:);
           CTD.jday_gmt(1:npr,j) = Ctd.jday_gmt;
           bad = find(CTD.T == 0);
           CTD.T(bad) = NaN; CTD.S(bad) = NaN; CTD.Sigma(bad) = NaN;
           CTD.P(bad) = NaN; CTD.N2(bad) = NaN;
       end
       CTD
       pause
%      save(['/bb/lien/NIW/NIW2016/Cronjob/matdata/' int2str(floatid)],'CTD','-append');
       clear CTD HPID floatid
   end


