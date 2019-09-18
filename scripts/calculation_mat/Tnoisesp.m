   clear; clf
   set(0,'defaulttextfontsize',12,'defaultaxesfontsize',12);
   hpid = 130:2:140;
   hpid = 132; jblock = 11;
%    
%     dT = beta * dV = beta * dADC * Vref / 2^23 / Inet;
%
% 
   W = 0.15; kzmax = 400;
   D = 1.4e-7; nu = 1.2e-6; q =3.7;
   A = load(['ema-7489a-0' int2str(hpid) '-tms.mat']); 
   A.f = (A.flabeg +A.flaend)/2;
   A.k = A.f/W; 

   A.Slad1 = (A.Sla1-A.logavgoff)/A.logavgsf; 
   A.Slad2 = (A.Sla2-A.logavgoff)/A.logavgsf; 
   beta = 25; Vref = 4; Inet = 0.8;
   scale2 = (beta * Vref / 2^23 / Inet)^2;
   A.rawTsp1 = 10.^(A.Slad1(jblock,:)/10)*scale2; A.rawTsp2 = 10.^(A.Slad2(jblock,:)/10)*scale2;
%
   ctdfn = ['ema-7489a-0' int2str(hpid) '-ctd.mat'];
   [Pr,N,dTdz,W]= Bkgfun(ctdfn,A.uxt(jblock));
   N = abs(N);
   W = abs(W);
%
%  trnasfer functina
%
   H2adc = H2ADCfun(A.f);
   H2preamp = H2preampfun(A.f);
   H2fp07 = H2FP07fun(A.f,W);
   H2total = H2adc.*H2preamp.*H2fp07;

   A.corrTsp1 = A.rawTsp1./H2total;

   f_rps = A.f*2*pi; k_rpm = f_rps/W;
   A.corrTsp1_rpm = A.corrTsp1*W/2/pi;
   A.corrdTdzsp1_rpm = (k_rpm(:).^2).*A.corrTsp1_rpm(:);
   g  = find(k_rpm <= kzmax);
   chi = 6*D*max(mmintgrl(k_rpm(g),A.corrdTdzsp1_rpm(g))); 
   KT = 0.5*chi/(dTdz^2);
   Gamma = 0.2;
   epsilon = KT*N^2/Gamma;
   Batchelorsp = batchelor_jhd(k_rpm, epsilon,chi, nu, D, q) /(2*pi);

   Bkgtitle1 = sprintf('Float 7498 hpid %4.0f block %6.0f',hpid,jblock);
   Bkgtitle2 = sprintf('Pr  = %4.0f m ',Pr);
   Bkgtitle3 = sprintf('\\partial_z T = %6.4f C/m',dTdz);
   Bkgtitle4 = sprintf('N   = %6.2e 1/s',N);
   Bkgtitletext = str2mat(Bkgtitle1,Bkgtitle2,Bkgtitle3,Bkgtitle4);

   chititle1 = sprintf('\\chi = %6.2e C^2/s',chi);
   chititle2 = sprintf('K_T  = %6.2e m^2/s',KT);
   chititle3 = sprintf('\\epsilon  = %6.2e m^2/s^3',epsilon);
   chititletext = str2mat(chititle1,chititle2,chititle3);
   
   x0 = 0.1; y0 = 0.7; dx = 0.35; dy = 0.2; ddx = 0.12; ddy = 0.1;
   axes('position',[x0 y0 dx dy],'box','on','xscale','log','yscale','log');
   hold on
   plot(A.f,H2adc,A.f,H2preamp,A.f,H2fp07,A.f,H2total,'linewidth',3);
   xlabel('Frequency (Hz)');
   ylabel('Transfer function squared H2');
   myaxis;
   title(Bkgtitletext)
%
%  Raw spectrum and Corrected Spectrum
%
   
   x0 = x0+dx+ddx;
   axes('position',[x0 y0 dx dy],'box','on','xscale','log','yscale','log');
   hold on
   plot(A.f,A.rawTsp1,A.f,A.corrTsp1,'linewidth',3);
   xlabel('Frequency (Hz)');
   ylabel('Raw and Corrected \Phi_T (C^2 / Hz)');
   myaxis;
%
%  Convert from frequency to wavenumber spectrum
%
   xtick = 10.^(1:3);
   x0 = 0.1; y0 = y0-dy-ddy;
   axes('position',[x0 y0 dx dy],'box','on','xscale','log','yscale','log','xtick',xtick);
   hold on
   plot(k_rpm,A.corrTsp1_rpm,'linewidth',3);
   xlabel('kz (m^{-1})');
   ylabel('Corrected \Phi_T (C^2 m)');
   myaxis;

   x0 = x0+dx+ddx;
   axes('position',[x0 y0 dx dy],'box','on','xscale','log','yscale','log','xtick',xtick);
   hold on
   plot(k_rpm,A.corrdTdzsp1_rpm,'linewidth',3);
   plot(k_rpm(g),Batchelorsp(g),'r-','linewidth',3);
   xlabel('kz (m^{-1})');
   ylabel('Corrected \Phi_{\partial_z T} (C^2 m^{-1})');
   ax = myaxis;

   ddy = 0.13;
   x0 = 0.1; y0 = y0-dy-ddy;
   axes('position',[x0 y0 dx dy],'box','on','xscale','log','xtick',xtick);
   hold on
   g = find(k_rpm <= kzmax);
   plot(k_rpm(g),col(k_rpm(g)).*col(A.corrdTdzsp1_rpm(g)),'linewidth',3);
   xlabel('kz (m^{-1})');
   ylabel('kz * Corrected \Phi_T (C^2 m)');
   myaxis; set(gca,'xlim',ax(1:2));
   title(chititletext)

   saveas(gcf,'chispec_example','pdf')
