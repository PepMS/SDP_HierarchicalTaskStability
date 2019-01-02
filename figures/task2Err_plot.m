function [hfig] = task2Err_plot(nfig, t, data)
% Figure propierties
fig.num = nfig;
fig.title = 'Task Error';
fig.fontsize = 25;
fig.position =  [100, 100, 1000, 400];
fig.linewidth = 3;
fig.labels.x = '[s]';
fig.labels.y = '[rad]';

% names
x.name = '$t$';
y{1}.name = '$||\mbox{\boldmath $e$}_2||$ - Ct. gains';
y{2}.name = '$||\mbox{\boldmath $e$}_2||$ - Gain Sch.';

% colors
color = lines(size(data,1));

% line styles
ls{1} = '-';
ls{2} = '--';

% fill structures
x.data = t;
for ii=1:size(data,1)
    y{ii}.data = data(ii,:);
    y{ii}.color = color(ii,:);
    y{ii}.linestyle = ls{ii};
end

hfig = genericPlotData(fig,x,y);
end