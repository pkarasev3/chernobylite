function setFigure(fh,pos,scaleinx,scaleiny)
    if( nargin < 3 )
        scaleinx =1;
        scaleiny =1;
    elseif(nargin < 4)
        scaleiny = scaleinx;
    end
    set(fh,'Position', [pos 400*scaleinx 300*scaleiny]);
end
