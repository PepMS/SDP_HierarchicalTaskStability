function [hfig] = svd_plot(nfig, t, data)
% Figure propierties
fig.num = nfig;
fig.title = 'SVD';
fig.fontsize = 25;
fig.position =  [100, 100, 1000, 400];
fig.linewidth = 3;
fig.labels.x = '[s]';
fig.labels.y = '[rad]';

% names
x.name = '$t$';
y{1}.name = '$\underline{\mbox{\boldmath $\sigma$}}_{1}$';
y{2}.name = '$\underline {\mbox{\boldmath $\sigma$}}_{2}$';

% colors
color = lines(size(data,1));

% line styles
ls{1} = '-';
ls{2} = '--';
ls{3} = ':';


% fill structures
x.data = t;
for ii=1:size(data,1)
    y{ii}.data = data(ii,:);
    y{ii}.color = color(ii,:);
    y{ii}.linestyle = ls{ii};
end

hfig = genericPlotData(fig,x,y);
end