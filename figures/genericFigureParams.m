function [hleg] = genericFigureParams(hfig,title,fontsize,pos)

    hleg = legend(gca,'show');
    set(hleg,'Interpreter','latex','FontSize',fontsize);
    
    set(hfig,'NumberTitle','off','Name',title,'Renderer','Painters','defaulttextinterpreter','latex');
    set(gca,'FontSize',fontsize,'box','off','layer','top','TickDir','out')
    set(gcf,'Color',[1,1,1])
    set(findall(gcf,'-property','FontSize'),'FontSize',fontsize)
    set(findall(gcf,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
    set(findall(gca,'-property','Interpreter'),'Interpreter','latex')
    
    if ~isempty(pos)
        set(hfig, 'Position', pos);
    end

    % set(gca,'LooseInset',get(gca,'TightInset'));

    grid ON
    % axis tight

    hax = hfig.CurrentAxes;
    xval = get(hax,'XLim');
    yval = get(hax,'YLim');
    rectangle('Position',[xval(1) yval(1) xval(2)-xval(1) yval(2)-yval(1)],'EdgeColor',[0.8 0.8 0.8]);

    hax.GridLineStyle = '--';
    
    plotPos = get(hax, 'Position');
    legendSize = get(hleg, 'Position');
    legendSize(1:2) = [];
    % legendPos = [plotPos(1)+plotPos(3)-0.005, 0.94-legendSize(2), legendSize];
    % set(hleg,'Position',legendPos);

return