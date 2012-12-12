addpath('~/source/chernobylite/matlab/display_helpers/');
addpath('~/source/chernobylite/matlab/util/');
addpath('~/source/chernobylite/matlab/display_helpers/');
addpath('~/source/chernobylite/matlab/LSMlibPK/');                                                                            %#ok<*NOPTS>
addpath('~/source/matlab2tikz_0.2.2/src');

resnames = {'Final_Seq1/NoCompNoU_LoC_10_results.mat',...
            'Final_Seq2/NoCompNoU_LoC_10_results.mat',...
            'Final_Seq3/NoCompNoU_LoC_10_results.mat',...
            'Final_Seq1/NoCompNoU_HiC_30_results.mat',...
            'Final_Seq2/NoCompNoU_HiC_30_results.mat',...
            'Final_Seq3/NoCompNoU_HiC_30_results.mat',...
            'Final_Seq1/YesCompNoU_5_results.mat',...
            'Final_Seq2/YesCompNoU_5_results.mat',...
            'Final_Seq3/YesCompNoU_5_results.mat'}

cntr_colors={'ms','ms','ms','kd','kd','kd','ro','ro','ro'};
          
sfigure(2); clf;           
sfigure(3); clf;           
for k =1:numel(resnames)
  s = load([resnames{k}]);

  frames = s.results.nFrame_in;
  idx0   = find( frames > 20, 1);
  
  % error  ||psi - phi||
  d_a = s.results.D_ls_err(idx0:end);
  m_a = mean(d_a); sd_a = sqrt(var(d_a));
  
  f2=sfigure(2); hold on; 
  semilogy( [k k], [min(d_a),max(d_a)], 'b--s', 'LineWidth',3);
  semilogy( [k k], [m_a-sd_a,m_a+sd_a], 'b-o', 'LineWidth',6); 
  semilogy( k, m_a, cntr_colors{k},'MarkerSize',11,'MarkerFaceColor',cntr_colors{k}(1) );
  hold off;
  
  
  % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ang = s.results.ang_diff(idx0:end) * 180 / pi;
  m_b = mean(ang); sd_b = sqrt(var(ang));
  
  sfigure(3); hold on; 
  plot( [k k], [min(ang),max(ang)], 'k--s', 'LineWidth',3);
  plot( [k k], [m_b-sd_b,m_b+sd_b], 'k-o', 'LineWidth',6); 
  plot( k, m_b, cntr_colors{k},'MarkerSize',11,'MarkerFaceColor',cntr_colors{k}(1) );
  hold off;
  
  
  %boxplot([d_a d_b d_c],'notch','on',...
  %        'labels',{'Y_5','N_10','N_30'})

end

f2=sfigure(2); ylabel('segmentation error');
axis([0,numel(resnames)+1,10,900]); 
grid on; drawnow; pause(0.05);
matlab2tikz('ktrack_boxplot_1.tikz.tex','width','10cm','height','7cm');

sfigure(3); ylabel('angular displacement (deg)');
axis([0,numel(resnames)+1,0,2]); grid on; drawnow; pause(0.05); 
matlab2tikz('ktrack_boxplot_2.tikz.tex','width','10cm','height','7cm');

!cp -v ./*.tex   ~/source/visioncontrol/thesis-pk/figs/
