function hfig = plotData(nfig_, title, t, data, name_x, name_y, label_x, label_y)

% Figure propierties
fig.num = nfig_;
fig.fontsize = 25;
fig.position =  [100, 100, 1000, 400];
fig.linewidth = 3;

% Title
fig.title = title;

% Label
fig.labels.x = label_x;
fig.labels.y = label_y;

% Names
x.name = name_x;

% colors
color = lines(size(data,1));

% line styles
ls = {'-.', '-', '--', ':'};

% fill structures
x.data = t;
for ii=1:size(data,1)
    y{ii}.name = strcat(name_y, num2str(ii),'$');
    y{ii}.data = data(ii, :);
    y{ii}.color = color(ii,:);
    y{ii}.linestyle = ls{mod(ii,4) + 1};
end

hfig = genericPlotData(fig,x,y);
end