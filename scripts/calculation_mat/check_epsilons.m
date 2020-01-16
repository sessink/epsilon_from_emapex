figure;
semilogx(Chi.eps1,-Chi.P,'LineWidth',2)
hold on
semilogx(Chi.eps2,-Chi.P,'LineWidth',2)
semilogx(Chi.eps1_goto,-Chi.P,'--','LineWidth',2)
semilogx(Chi.eps2_goto,-Chi.P,'--','LineWidth',2)

%%
figure;
plot(Chi.kb1,-Chi.P,'LineWidth',2)
hold on
plot(Chi.kb2,-Chi.P,'LineWidth',2)