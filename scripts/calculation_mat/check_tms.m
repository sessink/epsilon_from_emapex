   clear; clf
   co = get(0,'defaultaxescolororder');
   hpid = 130:2:140;
   for i = 1:length(hpid);
       A = load(['ema-7489a-0' int2str(hpid(i)) '-tms.mat']); 
       plot(A.Sla1,A.Sla2,'.','color',co(i,:))
       drawnow
       hold on
%      pause
   end
   xlabel('Sla1'); ylabel('Sla2');
   plot([20 100],[20 100],'k-')
   clf
   A.f = (A.flabeg +A.flaend)/2;
   A.k = A.f/0.15; A.Tsp1 = 10.^(A.Sla1/10); A.Tsp2 = 10.^(A.Sla2/10);
   A.k2 = ones(size(A.Sla1,1),1)*A.k(:)';
   A.dTdzsp1 = (A.k2.^2).*A.Tsp1;
   A.dTdzsp2 = (A.k2.^2).*A.Tsp2;
   loglog(A.k,A.dTdzsp1(100:210,:))


   break
   loglog(f,Sla1,f,Sla2)
   for i = 1:size(Sla1,1);
       loglog(f,Sla1(i,:),f,Sla2(i,:))
       [i size(Sla1,1)]
   end
   loglog(f,nanmean(Sla1,1),f,nanmean(Sla2,1))
    
