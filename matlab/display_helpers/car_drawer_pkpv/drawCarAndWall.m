function sh = drawCarAndWall( d, dstar, sh_in)

  if exist('sh_in','var')
    sh = sh_in;
  else
    sh = sfigure(1);
  end
  set(sh,'Resize','off');
  clf; 
  
  carLen = 1.0; 
  hold on;
  
  plot( [-0.1 -0.1], [-0.2 1], 'k-','LineWidth',3 );
  plot( [-0.1 -0.1], [-0.2 1], 'r-.','LineWidth',1 );
  plot( [0 0], [-0.2 1], 'r-','LineWidth',1 );
  
  plot( [dstar dstar], [-0.2 0.75], 'g--','LineWidth',2 );
  
  plot( d+0.5,-0.07,'o','MarkerFaceColor',[0 0 0],'MarkerSize',12);
  
  plot( [0 d], [-0.5 -0.5],'--o','MarkerFaceColor',[0 0 0],'MarkerSize',8);
  plot( 0, -0.5,'o','MarkerFaceColor','red','MarkerSize',8);
  plot( d, -0.5,'o','MarkerFaceColor','green','MarkerSize',8);
  
  plot( d+2.0,-0.07,'o','MarkerFaceColor',[0 0 0],'MarkerSize',12);
  plot( -3:100,-0.2+0*(-3:100),'k-','LineWidth',4); 
  plot( -3:100,-0.2+0*(-3:100),'y--','LineWidth',2); 
  text( d, 0.0,'    CAR    ','EdgeColor','none',...
             'BackgroundColor','cyan','HorizontalAlignment','Left','FontSize',8);

  text( 0, 1.1,'WALL','EdgeColor','red',...
             'BackgroundColor','none','HorizontalAlignment','Center','FontSize',8);
  
  text( d*0.5, -0.75,sprintf('d=%+2.4f',d),'EdgeColor','none',...
             'BackgroundColor','none','HorizontalAlignment','Center','FontSize',8);
                        
  if d < 0 
    text( 0.5, 1, '  {\bf \color{red}FAIL!} Crashed','HorizontalAlignment','Left','FontSize',36);
  end
           
  axis([ -3, 20, -1, 3 ] );
  grid on;  
  axis off;
  hold off; 

end
