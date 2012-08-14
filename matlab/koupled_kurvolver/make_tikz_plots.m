clear variables
set(0,'defaultaxesfontsize',16);
set(0,'defaulttextfontsize',18);
set(0,'defaulttextfontname','Arial');
set(0,'defaultaxesfontweight','bold');
set(0,'defaultlinelinewidth',2);
set(0,'defaultlinemarkersize',4);
dbstop if error;
addpath('~/source/chernobylite/matlab/util/');
addpath('~/source/chernobylite/matlab/display_helpers/');
addpath('~/source/chernobylite/matlab/LevelSetMethods/');
addpath('~/source/chernobylite/matlab/LSMlibPK/');
addpath('~/source/matlab2tikz_0.2.2/src');

load ./bridge1/all_data.mat
sfigure(1); 
plot( t_all, Dval_all, 'r-.'); hold on;
plot( t_all, Fval_all, 'b--'); hold off;

matlab2tikz( 'bridge1_DandF.tikz', 'height', '4cm', 'width', '62mm');
!cp -vu bridge1_DandF.tikz  ~/source/visioncontrol/kslice-TAC/


load ./bridge2/all_data.mat
sfigure(2); 
plot( t_all, Dval_all, 'r-.'); hold on;
plot( t_all, Fval_all, 'b--'); hold off;

matlab2tikz( 'bridge2_DandF.tikz', 'height', '4cm', 'width', '62mm');
!cp -vu bridge2_DandF.tikz  ~/source/visioncontrol/kslice-TAC/
