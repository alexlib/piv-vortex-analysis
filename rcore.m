
X2 = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21];

V2 = 0.75*(1-exp(-X2/4));
V2(2) = 0.31;
V2(3) = 0.35;
V2(4) = 0.41;
V2(5) = 0.45;
V2(6) = 0.48;
V2(7) = 0.539;
V2(8) = 0.55;
V2(9) = 0.57;
V2(10) = 0.62;
V2(12) = 0.66;
V2(13) = 0.710;
V2(14) = 0.732;
V2(15) = 0.745;



% X3 = [0 1.5 3.6 4.5];
% V3 = [0 3 3 0];

plot(X2,V2,'b')
% plot(X2,V2,'k',X3,V3,'b')
axis([0 21 0 1])
% legend('Velocity Profile (2m/s)','Location','Best')
% legend('Velocity Profile (2m/s)','Velocity Profile (3m/s)','Location','Best')
% set(gca, 'XTick', 1:7, 'XTickLabel', labels);
xlabel('Grid size, n')
% ylabel('\frac{\Gamma}{Vc}')
ylabel('$\displaystyle\frac{\Gamma}{Vc}$','interpreter','latex')
