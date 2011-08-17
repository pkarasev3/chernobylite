disp('A simple multivariate gaussian (pause)')

N           = 100;
mu          = [-5 ; -5];    %(d x 1)
sigma       = [2 0; 0 1];   %(d x d)
[Z , index] = sample_mvgm(N , mu , sigma);
disp(Z)
figure(1); title('A simple multivariate gaussian (pause)');
plot(Z(1 , :) , Z(2 , :) , '+')



pause;
disp('ND slice multivariate gaussian (pause)')


N           = 100;
mu          = [-5 ; -5];    %(d x 1)
sigma       = [2 0; 0 1];   %(d x d)
p           = 1;            %(1 x 1)
L           = 4;
[Z , index] = sample_mvgm(N , mu , sigma , p , L);

%or [Z , index] = sample_mvgm(N , mu , sigma , [] , L);

disp(Z)

pause;
disp('A simple Mixture of multivariate gaussian (pause)')


N           = 100;
V           = 3;
mu          = cat(3 , [-5 ; -5] , [0 ; 0] ,[ 5 ; 5]);                %(d x 1 x M)
sigma       = cat(3 , [2 0; 0 1] , [2 -.2; -.2 2] , [1 .9; .9 1]);   %(d x d x M)
p           = cat(3 , [0.3] , [0.2]  , [0.5]);                       %(1 x 1 x M)
[Z , index] = sample_mvgm(N , mu , sigma , p);

disp(Z)

figure(2); clf; hold on; title('A simple Mixture of multivariate gaussian (pause)');

plot(Z(1 , (index==1) ) , Z(2 , (index==1) ) , 'r+')
plot(Z(1 , (index==2) ) , Z(2 , (index==2) ) , 'bo')
plot(Z(1 , (index==3) ) , Z(2 , (index==3) ) , 'gx'); hold off;



pause;
disp('ND slice Mixture of multivariate gaussian (pause)')
pause;

N           = 100;
V           = 3;
L           = 2;
mu          = cat(3 , [-5 ; -5] , [0 ; 0] ,[ 5 ; 5]);                %(d x 1 x M)
sigma       = cat(3 , [2 0; 0 1] , [2 -.2; -.2 2] , [1 .9; .9 1]);   %(d x d x M)
p           = cat(3 , [0.3] , [0.2]  , [0.5]);                       %(1 x 1 x M)
[Z , index] = sample_mvgm(N , mu , sigma , p , L);


disp(Z);


pause;
disp('ND slice Mixture of multivariate gaussian with ND gaussian parameters (pause)')
pause;

N           = 100;
V           = 3;
L           = 2;
mu          = repmat(cat(3 , [-5 ; -5] , [0 ; 0] ,[ 5 ; 5]) , [1 , 1 , 1 , V]) + 1*randn(2 , 1 ,3 , V);                %(d x 1 x M)
sigma       = repmat(cat(3 , [2 0; 0 1] , [2 -.2; -.2 2] , [1 .9; .9 1]) , [1 , 1 , 1 , V]);   %(d x d x M)
p           = repmat(cat(3 , [0.3] , [0.2]  , [0.5]) , [1 , 1 , 1 , V]);                       %(1 x 1 x M)
[Z , index] = sample_mvgm(N , mu , sigma , p , L);

disp(Z)

