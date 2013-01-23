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
sfigure(1); t_all = linspace(0,1,numel(t_all));
plot( t_all, Dval_all, 'r-', 'LineWidth', 1); hold on;
plot( t_all, Fval_all, 'b--'); hold off; grid on;
legend('D(t)','F(t)'); xlabel('time (scaled)'); ylabel('functional value'); pause(0.1);
matlab2tikz( 'bridge1_DandF.tikz', 'height', '4cm', 'width', '62mm');
Ffinal=mean(Fval_all(end-5:end))
Dfinal=mean(Dval_all(end-5:end))
!cp -vu bridge1_DandF.tikz  ~/source/visioncontrol/kslice-TAC/


load ./bridge2/all_data.mat
sfigure(2); t_all = linspace(0,1,numel(t_all));
plot( t_all, Dval_all, 'r-', 'LineWidth', 1); hold on;
plot( t_all, Fval_all, 'b--'); hold off; grid on;
legend('D(t)','F(t)'); xlabel('time (scaled)'); ylabel('functional value'); pause(0.1);
matlab2tikz( 'bridge2_DandF.tikz', 'height', '4cm', 'width', '70mm');
Ffinal=mean(Fval_all(end-5:end))
Dfinal=mean(Dval_all(end-5:end))
!cp -vu bridge2_DandF.tikz  ~/source/visioncontrol/kslice-TAC/



load run_openloop_bridge_demo.mat
sfigure(3); t_all = linspace(0,1,numel(t_all));
plot( t_all, imgFunc_all, 'm.', 'LineWidth', 2);
legend('E(t)'); 
xlabel('time (scaled)'); ylabel('functional value'); grid on; pause(0.1);
matlab2tikz( 'imgFunc_OpenLoop.tikz', 'height', '3.5cm', 'width', '80mm');
!cp -vu imgFunc_OpenLoop.tikz  ~/source/visioncontrol/kslice-TAC/


% % Test to demo matlab2tikz
addpath('~/source/matlab2tikz_0.2.2/src');
figure(4); 
t = logspace(-4,1,1000);
x = sin( 1./t );
semilogx( t, x, 'g--','LineWidth',1);  legend('x(t)');
xlabel('x axis'); ylabel('y axis');
grid on;
matlab2tikz( 'aug17_tikz.tex', 'height', '3.5cm', 'width', '80mm');
!cp -vu aug17_tikz.tex  ~/source/visioncontrol/kslice-slides/overview/


