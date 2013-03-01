addpath('~/source/chernobylite/matlab/display_helpers/');
addpath('~/source/chernobylite/matlab/util/');
addpath('~/source/chernobylite/matlab/display_helpers/');
addpath('~/source/chernobylite/matlab/LSMlibPK/');                                                                            %#ok<*NOPTS>
addpath('~/source/matlab2tikz_0.2.2/src');

resnames = {'feb11C/ResultsKtrack.mat',...
             'feb11A/ResultsKtrack.mat'}

cntr_colors={'g-s','k-o'};
%close all; 
sfigure(1); clf; sfigure(2); clf; sfigure(3); clf; 
for k =1:numel(resnames)
  s = load([resnames{k}]);
  results = s.results;
  
  err_xy = sqrt( sum( (results.estm_xy - results.true_xy).^2, 2 ) );
  meanXYerr = mean( err_xy ) %#ok<NOPTS>
  
  frames     = results.nFrame_in;
  idx0       = 2:numel(frames);
  frameDelay = diff( frames );
  Ck         = cntr_colors{k}; 
  
  % % 
  sfigure(1); hold on; 
  semilogy( frames(idx0), err_xy(idx0),Ck); %title('centroid error');
  grid on; axis([0 203 2  min([128,1.01*max(err_xy(idx0) )]) ]); 
  drawnow; pause(0.001); 
  set(gca,'YTick',[2,20,40,80,120,200],'YTickLabel',[2,20,40,80,120,200],'YMinorGrid','on');
  legend('nominal system','controlled system','Location','NorthWest'); hold off;
  % % 
  sfigure(2); hold on; % title('frameskip');
  plot( frames(idx0), frameDelay, Ck ); drawnow; pause(0.001);
  grid on; set(gca,'YTick',[0 1 2 3 4 5 6],'YTickLabel',[0 1 2 3 4 5 6]);
  axis([0 203 0 1.01*max(frameDelay(:))]);
  legend('nominal system','controlled system','Location','North'); 
  hold off; 
  % %
  sfigure(3); hold on; plot( frames(idx0), results.Area(idx0),Ck ); %title('area');
  grid on; legend('nominal system','controlled system','Location','NorthWest');
  drawnow; pause(.001); set(gca,'YTick',[0 700 1100 2000 3500 4500],...
                                       'YTickLabel',[0 700 1100 2000 3500 4500]);
  axis([0 203 0 1.01*max(results.Area(idx0))]);
  hold off;
end

figure(1); drawnow; matlab2tikz('ktrack_typeIIa_centroid.tikz.tex','width','11cm','height','4cm',...
           'extraTikzpictureSettings','\tikzset{style={font=\fontsize{10}{10}\selectfont}}');
           
figure(2); drawnow; matlab2tikz('ktrack_typeIIa_frameskip.tikz.tex','width','11cm','height','4cm',...
           'extraTikzpictureSettings','\tikzset{style={font=\fontsize{10}{10}\selectfont}}');
figure(3); drawnow; matlab2tikz('ktrack_typeIIa_area.tikz.tex','width','11cm','height','4cm',...
           'extraTikzpictureSettings','\tikzset{style={font=\fontsize{10}{10}\selectfont}}');

!cp -v ./ktrack_typeIIa*.tex   ~/source/visioncontrol/thesis-pk/figs/

% figure(1); xlabel('sent frame index'); axis([0 202 0 6.5]); ylabel('frame skip');
% legend('A','B','C','D');pause(.001); drawnow; 
% matlab2tikz('ktrack_frameskip.tikz.tex','width','10cm','height','6cm');
% 
% figure(3); xlabel('sent frame index'); axis([0 202 0 0.85]); ylabel('angular displacement [deg]');
% legend('A','B','C','D');pause(.001); drawnow; 
% matlab2tikz('ktrack_angularDisplacement.tikz.tex','width','10cm','height','6cm');
% 
% figure(4); xlabel('sent frame index'); axis([0 202 4.0 200.0]); ylabel('centroid error [pixels]');
% legend('A','B','C','D');pause(.001); drawnow; 
% matlab2tikz('ktrack_centroidError.tikz.tex','width','10cm','height','6cm');

% f2=sfigure(2); ylabel('segmentation error');
% axis([0,numel(resnames)+1,10,900]); 
% grid on; drawnow; pause(0.05);
% matlab2tikz('ktrack_boxplot_1.tikz.tex','width','10cm','height','7cm');
% 
% sfigure(3); ylabel('angular displacement (deg)');
% axis([0,numel(resnames)+1,0,2]); grid on; drawnow; pause(0.05); 
% matlab2tikz('ktrack_boxplot_2.tikz.tex','width','10cm','height','7cm');

%!cp -v ./*.tex   ~/source/visioncontrol/thesis-pk/figs/
