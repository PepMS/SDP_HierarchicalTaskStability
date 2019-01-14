function [hfig] = genericPlotDataLog(fig,x,y)
hfig = figure(fig.num);

% hold on;

ymin = 0.0;
ymax = 0.0;

for ii = 1:size(y,2)
    semilogy(x.data,y{ii}.data,'Color',y{ii}.color,'LineWidth',fig.linewidth,'LineStyle',y{ii}.linestyle,'DisplayName',y{ii}.name);
    if (min(y{ii}.data) < ymin)
        ymin = min(y{ii}.data);
    end
    if (max(y{ii}.data) > ymax)
        ymax = max(y{ii}.data);
    end
end

xlabel(fig.labels.x);
ylabel(fig.labels.y);

% axis([x.data(1) x.data(end) ymin ymax]);

% hold off;

genericFigureParams(hfig,fig.title,fig.fontsize,fig.position);

end