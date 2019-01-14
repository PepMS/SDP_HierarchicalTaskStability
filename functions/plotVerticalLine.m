function plotVerticalLine(fig, t)
figure(fig);
hold on;
p = plot([t t], [-10000 +10000],'k-.');
set(get(get(p(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
hold off;

end
