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

cntr_colors={'ks','ks','ks','bd','bd','bd','go','go','go'};
          
sfigure(2); clf;           
sfigure(3); clf;           
for k =1:numel(resnames)
  s = load([resnames{k}]);

  frames = s.results.nFrame_in;
  idx0   = find( frames > 20, 1);
  
  % error  ||psi - phi||
  d_a = s.results.D_ls_err(idx0:end);
  m_a = median(d_a); sd_a = sqrt(var(d_a)); mL=median(d_a(d_a<m_a)); mU=median(d_a(d_a>m_a));
  
  sfigure(2); hold on; 
  semilogy( [k k], [min(d_a),max(d_a)], '--s', 'Color',cntr_colors{k}(1),'LineWidth',2);
  semilogy( [k k], [mL,mU], '-o','Color',cntr_colors{k}(1),'LineWidth',3); 
  semilogy( k, m_a, '+','LineWidth',2, 'Color',cntr_colors{k}(1),'MarkerSize',6,'MarkerFaceColor',cntr_colors{k}(1) );
  hold off;
  
  
  % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ang = s.results.ang_diff(idx0:end) * 180 / pi;
  m_b = median(ang); sd_b = sqrt(var(ang));
  mL=median(ang(ang<m_b)); mU=median(ang(ang>m_b));
  
  sfigure(3); hold on; 
  plot( [k k], [min(ang),max(ang)], '--s','Color',cntr_colors{k}(1), 'LineWidth',2);
  plot( [k k], [mL,mU], '-o','Color',cntr_colors{k}(1), 'LineWidth',3); 
  plot( k, m_b, '+','Color',cntr_colors{k}(1),'LineWidth',2,'MarkerSize',6,'MarkerFaceColor',cntr_colors{k}(1) );
  hold off;
  
  
  %boxplot([d_a d_b d_c],'notch','on',...
  %        'labels',{'Y_5','N_10','N_30'})

end

moreOpts={'extraTikzpictureSettings','\tikzset{style={font=\fontsize{10}{10}\selectfont}}'};

f2=sfigure(2); ylabel('segmentation error');
axis([0,numel(resnames)+1,10,900]); 
grid on; drawnow; pause(0.05);
matlab2tikz('ktrack_boxplot_SegErr_tikz.tex','width','13cm','height','6cm');

sfigure(3); ylabel('angular displacement (deg)');
axis([0,numel(resnames)+1,0,2]); grid on; drawnow; pause(0.05); 
matlab2tikz('ktrack_boxplot_AngDis_tikz.tex','width','13cm','height','6cm');

%!cp -v ./*.tex   ~/source/visioncontrol/thesis-pk/figs/
