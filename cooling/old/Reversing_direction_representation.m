data=get(gca,'Children');

Bond_out_surface=data(1).YData;
Inner_gas=data(2).YData;


Inner_gas=Inner_gas(end:-1:1);
Bond_out_surface=Bond_out_surface(end:-1:1);

Liner_length=0.650;
X_axis=linspace(0,Liner_length,length(Inner_gas));

figure;
plot(X_axis,Bond_out_surface);
hold on;
plot(X_axis,Inner_gas);
title('Temperature of Surface outside');
% ylabel('Temperature ^{\circ}C')
xlabel('Length (m)');
ylabel('Temperature ^{\circ}C');

legend('Inner Gas','Bonded outside Surface')
