
sfigure(2); 
plot(results.D_ls_err,'r-o');
mean(results.D_ls_err)
sqrt(var(results.D_ls_err))


% plot( results.estm_xy(:,1) - results.true_xy(:,1),'m-x'); hold on;
% plot( results.estm_xy(:,2) - results.true_xy(:,2),'r-s'); hold off;
% legend('x err','y err');
