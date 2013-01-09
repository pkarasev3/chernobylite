set(0,'defaultaxesfontsize',16);  
set(0,'defaulttextfontsize',18);
set(0,'defaulttextfontname','Arial');
set(0,'defaultaxesfontweight','bold');
set(0,'defaultlinelinewidth',2);
set(0,'defaultlinemarkersize',12);
dbstop if error;

d     = 5;
dstar = 1; 
sh = drawCarAndWall( d, dstar );

sh2   = sfigure(2);
d     = -1;
dstar =  2;
sh2   = drawCarAndWall( d, dstar, sh2);
