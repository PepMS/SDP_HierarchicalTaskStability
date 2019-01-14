function plotHorizontalLine(fig, t, y);
figure(fig)
hold on
p = plot(t, y*ones(1, length(t)),'k--');
set(get(get(p(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
hold off
end