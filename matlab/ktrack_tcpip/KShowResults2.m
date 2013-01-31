addpath('~/source/chernobylite/matlab/display_helpers/');
addpath('~/source/chernobylite/matlab/util/');
addpath('~/source/chernobylite/matlab/display_helpers/');
addpath('~/source/chernobylite/matlab/LSMlibPK/');                                                                            %#ok<*NOPTS>
addpath('~/source/matlab2tikz_0.2.2/src');

resnames = {'jan29A/ResultsKtrack.mat',...
            'jan29B/ResultsKtrack.mat',...
            'jan29C/ResultsKtrack.mat',...
            'jan29D/ResultsKtrack.mat'}

cntr_colors={'k-o','b-d','r-x','g-*'};

sfigure(1); clf; sfigure(2); clf; sfigure(3); clf; sfigure(4); clf;          
for k =1:numel(resnames)
  s = load([resnames{k}]);

  frames     = s.results.nFrame_in;
  idx0       = 2:numel(frames);
  frameDelay = diff( frames );
  % error  ||psi - phi||
  
  
  f1=sfigure(1); hold on; 
  plot( frames(idx0), frameDelay, cntr_colors{k},...
                  'MarkerSize',4,'MarkerFaceColor',cntr_colors{k}(1),...
                  'LineWidth',1+(mod(k,2)==0));
  grid on; hold off;
  
  
  % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ang = s.results.ang_diff(idx0) * 180 / pi;
  m_b = median(ang)
  sd_bp = median(ang(ang>m_b));
  sd_bn = median(ang(ang<m_b));
  
  f2=sfigure(2); hold on; 
  plot( [k k], [min(ang),max(ang)], 'k--s', 'LineWidth',3);
  plot( [k k], [sd_bn,sd_bp], 'k-o', 'LineWidth',6); 
  plot( k, m_b, cntr_colors{k},'MarkerSize',11,'MarkerFaceColor',cntr_colors{k}(1) );
  grid on; hold off;
  
  f3=sfigure(3); hold on; 
  plot( frames(idx0), ang, cntr_colors{k},...
                  'MarkerSize',4,...
                  'LineWidth',1+(mod(k,2)==0)); grid on;
  hold off;
  %median(s.results.D_ls_err)
  err_xy = sqrt( sum( (s.results.estm_xy - s.results.true_xy).^2, 2 ) );
  med_exy=median(err_xy)
  f4=sfigure(4); 
  semilogy( frames(idx0), err_xy(idx0), cntr_colors{k},...
        'MarkerSize',4,'LineWidth',1+(mod(k,2)==0)); grid on; hold on; 
  
  
  
  fprintf('');
  %boxplot([d_a d_b d_c],'notch','on',...
  %        'labels',{'Y_5','N_10','N_30'})

end

figure(1); xlabel('sent frame index'); axis([0 202 0 6.5]); ylabel('frame skip');
legend('A','B','C','D');pause(.001); drawnow; 
matlab2tikz('ktrack_frameskip.tikz.tex','width','10cm','height','6cm');

figure(3); xlabel('sent frame index'); axis([0 202 0 0.85]); ylabel('angular displacement [deg]');
legend('A','B','C','D');pause(.001); drawnow; 
matlab2tikz('ktrack_angularDisplacement.tikz.tex','width','10cm','height','6cm');

figure(4); xlabel('sent frame index'); axis([0 202 4.0 200.0]); ylabel('centroid error [pixels]');
legend('A','B','C','D');pause(.001); drawnow; 
matlab2tikz('ktrack_centroidError.tikz.tex','width','10cm','height','6cm');

% f2=sfigure(2); ylabel('segmentation error');
% axis([0,numel(resnames)+1,10,900]); 
% grid on; drawnow; pause(0.05);
% matlab2tikz('ktrack_boxplot_1.tikz.tex','width','10cm','height','7cm');
% 
% sfigure(3); ylabel('angular displacement (deg)');
% axis([0,numel(resnames)+1,0,2]); grid on; drawnow; pause(0.05); 
% matlab2tikz('ktrack_boxplot_2.tikz.tex','width','10cm','height','7cm');

%!cp -v ./*.tex   ~/source/visioncontrol/thesis-pk/figs/
