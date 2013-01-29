addpath('~/source/chernobylite/matlab/display_helpers/');
addpath('~/source/chernobylite/matlab/util/');
addpath('~/source/chernobylite/matlab/display_helpers/');
addpath('~/source/chernobylite/matlab/LSMlibPK/');                                                                            %#ok<*NOPTS>
addpath('~/source/matlab2tikz_0.2.2/src');

resnames = {'jan29A/ResultsKtrack.mat',...
            'jan29B/ResultsKtrack.mat',...
            'jan29C/ResultsKtrack.mat',...
            'jan29D/ResultsKtrack.mat'}

cntr_colors={'k-o','b--x','r-s','g--x'};

sfigure(1); clf; sfigure(2); clf;           sfigure(3); clf;           
for k =1:numel(resnames)
  s = load([resnames{k}]);

  frames     = s.results.nFrame_in;
  idx0       = 2:numel(frames);
  frameDelay = diff( frames );
  % error  ||psi - phi||
  
  
  f2=sfigure(1); hold on; 
  plot( frames(idx0), frameDelay, cntr_colors{k},...
                  'MarkerSize',4,'MarkerFaceColor',cntr_colors{k}(1),...
                  'LineWidth',1+(mod(k,2)==0));
  hold off;
  
  
  % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ang = s.results.ang_diff(idx0) * 180 / pi;
  m_b = mean(ang)
  sd_b = sqrt(var(ang));
  
  sfigure(2); hold on; 
  plot( [k k], [min(ang),max(ang)], 'k--s', 'LineWidth',3);
  plot( [k k], [m_b-sd_b,m_b+sd_b], 'k-o', 'LineWidth',6); 
  plot( k, m_b, cntr_colors{k},'MarkerSize',11,'MarkerFaceColor',cntr_colors{k}(1) );
  hold off;
  
  f3=sfigure(3); hold on; 
  plot( frames(idx0), ang, cntr_colors{k},...
                  'MarkerSize',4,'MarkerFaceColor',cntr_colors{k}(1),...
                  'LineWidth',1+(mod(k,2)==0));
  hold off;
  
  % One more plot: (x,y) errors
  
  fprintf('');
  %boxplot([d_a d_b d_c],'notch','on',...
  %        'labels',{'Y_5','N_10','N_30'})

end

% f2=sfigure(2); ylabel('segmentation error');
% axis([0,numel(resnames)+1,10,900]); 
% grid on; drawnow; pause(0.05);
% matlab2tikz('ktrack_boxplot_1.tikz.tex','width','10cm','height','7cm');
% 
% sfigure(3); ylabel('angular displacement (deg)');
% axis([0,numel(resnames)+1,0,2]); grid on; drawnow; pause(0.05); 
% matlab2tikz('ktrack_boxplot_2.tikz.tex','width','10cm','height','7cm');

%!cp -v ./*.tex   ~/source/visioncontrol/thesis-pk/figs/
