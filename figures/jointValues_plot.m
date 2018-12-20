function [hfig] = jointValues_plot(nfig, t, data)
% Figure propierties
    fig.num = nfig;
    fig.title = 'Joint velocity';
    fig.fontsize = 25;
    fig.position =  [100, 100, 1000, 400];
    fig.linewidth = 3;
    fig.labels.x = '[s]';
    fig.labels.y = '[rad/s]';

    % names
    x.name = '$t$';
    y{1}.name = '$\dot{\bfq}_1$';
    y{2}.name = '$e_{\dot{y}}$';
    y{3}.name = '$e_{\dot{z}}$';

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