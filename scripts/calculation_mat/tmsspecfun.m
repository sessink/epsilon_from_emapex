%
%    tmsspecun.m
%
%        compute turbulence temperature spectrum, chi, kT, and epsilon
%
   function Fit=tmsspecfun(tmsdata_dir,ctddata_dir,floatid,hpid,jblock,kzmin,kzmax,plotting, threshold);
%
%     dT = beta * dV = beta * dADC * Vref / 2^23 / Inet;
%
   D = 1.4e-7; nu = 1.2e-6; q =3.7;
   tmsfile = [tmsdata_dir int2str(floatid) 'b/' 'ema-' int2str(floatid) 'b-' sprintf( '%04d', hpid ) '-tms.mat'];
   % changed sprintf from original code
   
   A = load(tmsfile);
   A.f = (A.flabeg +A.flaend)/2;
   A.jday_gmt = nanmean(A.uxt/86400) + datenum(1970,1,1,0,0,0);

   A.Slad1 = (A.Sla1-A.logavgoff)/A.logavgsf;
   A.Slad2 = (A.Sla2-A.logavgoff)/A.logavgsf;
   beta = 25; Vref = 4; Inet = 0.8;
   scale2 = (beta * Vref / 2^23 / Inet)^2;
   A.rawTsp1 = 10.^(A.Slad1(jblock,:)/10)*scale2;
   A.rawTsp2 = 10.^(A.Slad2(jblock,:)/10)*scale2;

   Fit.rawTsp1_cps = A.rawTsp1;
   Fit.rawTsp2_cps = A.rawTsp2;
   Fit.f_cps = A.f;
%
%    Find Vertical Velocity and Stratification and convert to wavenumber
%
   ctdfn = [ctddata_dir int2str(floatid) 'b/' 'ema-' int2str(floatid) 'b-' sprintf( '%04d', hpid ) '-ctd.mat'];
   % changed sprintf from original code
   [Pr,N,dTdz,W,T]= nBkgfun(ctdfn,A.uxt(jblock));
   N = abs(N); W = abs(W); A.k = A.f/W;
   Fit.k_cpm = A.k;
   f_rps = A.f*2*pi;
   k_rpm = f_rps/W;
   Fit.f_rps = f_rps; Fit.k_rpm = k_rpm;
   Fit.N = N; Fit.W = W; Fit.dTdz = dTdz; Fit.Pr = Pr;
   Fit.T = T; Fit.jday_gmt = A.uxt(jblock)/86400+datenum(1970,1,1,0,0,0);
%
%  Applying transfer functions
%
   H2adc = H2ADCfun(A.f);
   H2preamp = H2preampfun(A.f);
   H2fp07 = H2FP07fun(A.f,W);
   H2total = H2adc.*H2preamp.*H2fp07;
   Fit.H2total_cps = H2total;
   Fit.H2fp07 = H2fp07;

   Gamma = 0.2;

   A.corrTsp1_cps = A.rawTsp1./H2total;
   A.corrTsp1_cps= cleannoisespfun(Fit.f_cps,A.corrTsp1_cps,threshold);
   A.corrTsp1_rpm = A.corrTsp1_cps*W/2/pi;
   A.corrdTdzsp1_rpm = (k_rpm(:).^2).*A.corrTsp1_rpm(:);
   g  = find(k_rpm <= kzmax & k_rpm >= kzmin);
   if length(g) >= 3;
      chi1 = 6*D*max(mmintgrl(k_rpm(g),A.corrdTdzsp1_rpm(g)));
   else
      chi1 = NaN;
   end
   Fit.corrTsp1_rpm = A.corrTsp1_rpm;
   Fit.corrdTdzsp1_rpm = A.corrdTdzsp1_rpm;
   Fit.chi1 = chi1;
   Fit.corrTsp1_cps = A.corrTsp1_cps;

   KT1 = 0.5*chi1/(dTdz^2);
   epsilon1 = KT1*N^2/Gamma;
   Batchelorsp1 = batchelor_jhd(k_rpm, epsilon1,chi1, nu, D, q) /(2*pi);
   Fit.KT1 = KT1; Fit.epsilon1 = epsilon1; Fit.Batchelorsp1 = Batchelorsp1;

   A.corrTsp2_cps = A.rawTsp2./H2total;
   A.corrTsp2_cps = cleannoisespfun(Fit.f_cps,A.corrTsp2_cps,threshold);
   A.corrTsp2_rpm = A.corrTsp2_cps*W/2/pi;
   A.corrdTdzsp2_rpm = (k_rpm(:).^2).*A.corrTsp2_rpm(:);
   g  = find(k_rpm <= kzmax & k_rpm >= kzmin);
   if length(g) >= 3;
      chi2 = 6*D*max(mmintgrl(k_rpm(g),A.corrdTdzsp2_rpm(g)));
   else
      chi2 = NaN;
   end
   Fit.corrTsp2_rpm = A.corrTsp2_rpm;
   Fit.corrdTdzsp2_rpm = A.corrdTdzsp2_rpm;
   Fit.chi2 = chi2;
   Fit.corrTsp2_cps = A.corrTsp2_cps;

   KT2 = 0.5*chi2/(dTdz^2);
   epsilon2 = KT2*N^2/Gamma;
   Batchelorsp2 = batchelor_jhd(k_rpm, epsilon2,chi2, nu, D, q) /(2*pi);
   Fit.KT2 = KT2; Fit.epsilon2 = epsilon2; Fit.Batchelorsp2 = Batchelorsp2;

   if plotting;

      Bkgtitle1 = sprintf('Float %4.0f  hpid %4.0f block %6.0f',floatid,hpid,jblock);
      Bkgtitle2 = sprintf('Pr  = %4.0f m ',Pr);
      Bkgtitle3 = sprintf('\\partial_z T = %6.4f C/m',dTdz);
      Bkgtitle4 = sprintf('N   = %6.2e 1/s',N);
      Bkgtitletext = str2mat(Bkgtitle1,Bkgtitle2,Bkgtitle3,Bkgtitle4);

      chititle1 = sprintf('\\chi_1 = %6.2e C^2/s',chi1);
      chititle2 = sprintf('K_T_1  = %6.2e m^2/s',KT1);
      chititle3 = sprintf('\\epsilon_1  = %6.2e m^2/s^3',epsilon1);
      chititletext1 = str2mat(chititle1,chititle2,chititle3);

      chititle1 = sprintf('\\chi_2 = %6.2e C^2/s',chi2);
      chititle2 = sprintf('K_T_2  = %6.2e m^2/s',KT2);
      chititle3 = sprintf('\\epsilon_2  = %6.2e m^2/s^3',epsilon2);
      chititletext2 = str2mat(chititle1,chititle2,chititle3);

      xtick = 10.^(-1:2);
      x0 = 0.18; y0 = 0.65; dx = 0.3; dy = 0.2; ddx = 0.15; ddy = 0.25;
      axes('position',[x0 y0 dx dy],'box','on','xscale','log','yscale','log','xtick',xtick);
      hold on
      plot(A.f,H2adc,A.f,H2preamp,A.f,H2fp07,A.f,H2total,'linewidth',3);
      xlabel('Frequency (Hz)');
      ylabel('Transfer function squared H2');
%       myaxis;
      title(Bkgtitletext)
   %
   %  Raw spectrum and Corrected Spectrum

      x0 = x0+dx+ddx;
      axes('position',[x0 y0 dx dy],'box','on','xscale','log','yscale','log','xtick',xtick);
      hold on
      plot(A.f,A.rawTsp1,A.f,A.corrTsp1_cps,'linewidth',3);
      xlabel('Frequency (Hz)');
      ylabel('Raw and Corrected \Phi_T (C^2 / Hz)');
%       myaxis;
  
   %
   %  Convert from frequency to wavenumber spectrum
   %
      xtick = 10.^(1:3);
      x0 = 0.18; y0 = y0-dy-ddy;
      axes('position',[x0 y0 dx dy],'box','on','xscale','log','yscale','log','xtick',xtick);
      hold on
      plot(k_rpm,A.corrTsp1_rpm,'linewidth',3);
      xlabel('kz (m^{-1})');
      ylabel('Corrected \Phi_T (C^2 m)');
%       myaxis;
      title(chititletext1)

      x0 = x0+dx+ddx;
      axes('position',[x0 y0 dx dy],'box','on','xscale','log','yscale','log','xtick',xtick);
      hold on
      plot(k_rpm,A.corrdTdzsp1_rpm,k_rpm,A.corrdTdzsp2_rpm,'linewidth',3);
      plot(k_rpm(g),Batchelorsp1(g),'k--',k_rpm(g),Batchelorsp2(g),'r--','linewidth',3);
      xlabel('kz (m^{-1})');
      ylabel('Corrected \Phi_{\partial_z T} (C^2 m^{-1})');
%       ax = myaxis;
 %    set(gca,'ylim',[1.0e-9 Inf])
      title(chititletext2)
 %    pause
   %  saveas(gcf,'chispec_example','pdf')
   end
