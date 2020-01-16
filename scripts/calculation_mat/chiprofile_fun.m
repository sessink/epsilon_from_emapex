   function [Chi,Fit] = chiprofile_fun(floatid,hpid,data_dir,kzmin,kzmax,plotting,threshold);
   %%    
       tmsdata_dir = data_dir;
       ctddata_dir = tmsdata_dir;
       tmsfile=[tmsdata_dir int2str(floatid) 'b/ema-' int2str(floatid) 'b-' ...
                sprintf( '%04d', hpid ) '-tms.mat'];
       disp(tmsfile)

       if ~exist(tmsfile)
           disp("doesn't exist")
           Chi = [];
       else
          A = load(tmsfile);
          for i = 1:A.nobs;
              jblock = i;
              Fit=tmsspecfun(tmsdata_dir,ctddata_dir,floatid,hpid,jblock,kzmin,kzmax,...
              plotting,threshold);
              
              Pr(i) = Fit.Pr; W(i) = Fit.W; N(i) = Fit.N;
              chi1(i) = Fit.chi1; chi2(i) = Fit.chi2;
              epsilon1(i) = Fit.epsilon1; epsilon2(i) = Fit.epsilon2;
              KT1(i) = Fit.KT1; KT2(i) = Fit.KT2;
              Chi.f_cps(:,i) = Fit.f_cps(:);
              Chi.corrTsp1_cps(:,i) = Fit.corrTsp1_cps(:);
              Chi.corrTsp2_cps(:,i) = Fit.corrTsp2_cps(:);
              Chi.H2fp07 = Fit.H2fp07;
              Chi.rawTsp1_cps = Fit.rawTsp1_cps;
              Chi.rawTsp2_cps = Fit.rawTsp2_cps;
              T(i) = Fit.T;
              Chi.jday_gmt(i) = Fit.jday_gmt; 
    %         Pr(i)
    %         disp('hit enter to continue')
        %      pause
              clf
          end
          if exist('Pr');
             Chi.P = Pr; Chi.W = W; Chi.N = N; Chi.chi1 = chi1; Chi.chi2 = chi2;
             Chi.eps1 = epsilon1; Chi.eps2 = epsilon2; Chi.kT1= KT1; Chi.kT2 = KT2;
             Chi.T = T;
             if plotting;
                x0 = 0.1; y0 = 0.7; dx = 0.2; dy = 0.25; ddx= 0.02;
                axes('position',[x0 y0 dx dy],'box','on','xscale','log','ydir','reverse');
                hold on
                plot(chi1,-Pr,chi2,-Pr);
                xlabel('\chi (C^2/s)');
                ylabel('Depth (m)'); 

                x0 = x0+dx+ddx; 
                axes('position',[x0 y0 dx dy],'box','on','xscale','log','ydir','reverse');
                hold on
                plot(KT1,-Pr,KT2,-Pr);
                xlabel('k_T (m^2/s)');
                set(gca,'yticklabel','');

                x0 = x0+dx+ddx; 
                axes('position',[x0 y0 dx dy],'box','on','xscale','log','ydir','reverse');
                hold on
                plot(epsilon1,-Pr,epsilon2,-Pr);
                xlabel('\epsilon (W/kg)');
                set(gca,'yticklabel','');
             end
          else 
            Chi = [];
          end
   end
