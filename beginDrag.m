function beginDrag(ax)
    fig = ancestor(ax,'figure');
    pos0 = get(ax,'Position');
    pt0  = get(fig,'CurrentPoint');
    
    set(fig,'WindowButtonMotionFcn',@(s,~)drag(ax,pos0,pt0));
    set(fig,'WindowButtonUpFcn',@(s,~)stopDrag(fig));
end

function drag(ax,pos0,pt0)
    pt = get(ancestor(ax,'figure'),'CurrentPoint');
    dx = (pt(1)-pt0(1))/get(0,'ScreenPixelsPerInch');
    dy = (pt(2)-pt0(2))/get(0,'ScreenPixelsPerInch');
    set(ax,'Position',pos0 + [dx dy 0 0]);
end

function stopDrag(fig)
    set(fig,'WindowButtonMotionFcn','');
    set(fig,'WindowButtonUpFcn','');
end