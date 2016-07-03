
X2 = [0 1 3.6 4.5];
V2 = [0 2 2 0];
% X3 = [0 1.5 3.6 4.5];
% V3 = [0 3 3 0];

plot(X2,V2,'b')
% plot(X2,V2,'k',X3,V3,'b')
axis([0 5 0 3])
% legend('Velocity Profile (2m/s)','Location','Best')
% legend('Velocity Profile (2m/s)','Velocity Profile (3m/s)','Location','Best')
% set(gca, 'XTick', 1:7, 'XTickLabel', labels);
xlabel('Distance (m)')
ylabel('Velocity (m/s)')