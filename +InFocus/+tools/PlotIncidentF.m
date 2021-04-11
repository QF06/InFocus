function [] = PlotIncidentF(ax,intensity,name,XTick,XTickLabel,YTick,YTickLabel,map,scale)

if (map == 'fire')
  map = InFocus.tools.fire;
  imagesc(intensity,'Parent',ax);
  set(ax,'YDir','normal');
  set(ax,'XTick',XTick);
  set(ax,'XTickLabel',XTickLabel);
  set(ax,'YTick',YTick);
  set(ax,'YTickLabel',YTickLabel);
  set(ax,'Colorscale',scale);
  title(ax, ['Intensity at: ', name]);
  colormap(ax,map);
  close;
else
  imagesc(intensity,'Parent',ax);
  set(ax,'YDir','normal');
  set(ax,'XTick',XTick);
  set(ax,'XTickLabel',XTickLabel);
  set(ax,'YTick',YTick);
  set(ax,'YTickLabel',YTickLabel);
  set(ax,'Colorscale',scale);
  title(ax, ['Intensity at: ', name]);
  colormap(ax,map);
  close;
end
end

