function genericPrintFig(hfig, name)
    set(hfig,'Units','centimeters');
    pos = get(hfig,'Position');
    set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3)+2, pos(4)+2])
    print(hfig,name,'-dpdf','-r0')
end