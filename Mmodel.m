function Mmodel(De,K, dataout,figout1,figout2,D10,D20,D30,D40,D,pi,sd1,sd2,sd3,sd4,psd2,psd3,psd4)

% matlab function for maximum likelihood estimation of diffusion
% coefficients for mixture models

% The program computes the maximum likelihood estimation of the diffusion coefficient
% from  models with one,two,three, and four diffusion coefficients

% Input parameters
% De: estimated diffusion coefficients from SPRIA analysis
% K: number of increments used to compute each of the value in De
% dataout: string containing the name of txt file where to save the results
% figout1, figout2: strings containing the names where to save the figures
% produced
% D10, D20, D30, D40: start values for the different diffusion
% coefficients (use uninformal guesses)
% D , pi: true diffusion coefficients and their proportions, leave as [] if unknown
% sd1,sd2,sd3,sd4,psd2,psd3,psd4 are the bootstrapped standard deviation of
% the diffusion coefficients and proportions for mixture models with
% 1,2,3,4 components

% the results are saved in the file dataout and figout1-figout2

options = optimset('MaxFunEvals', 10^4, 'MaxIter', 10^4,'Display','iter-detailed', 'TolFun', 10^(-6) );


sumK = sum(K);
J = size(K,1);
temporary = gcf;
figcount = temporary.Number + 1;
figure (figcount)
             clf
             hold on
  plot(De,K,'.')
  xlabel('D_{SPRIA}')
  ylabel('K')
             hold off
  print('-depsc',figout1)
figcount = figcount +1;
numOfBins = 100;
[histFreq, histXout] = hist(De, numOfBins);
binWidth = histXout(2)-histXout(1);
h = bar(histXout, histFreq/binWidth/sum(histFreq),1,'w'); 
axis([0 min([max(De), 10]) 0 1])
hold on  
  
dc = De' * K / sumK;
sdc = dc / sqrt( sumK/2 );
const = ( K/2 - 1 )' * log( De )-sum( log( gamma( K/2 ) ) )+ K'/2 * log( K/2 );

h1=plot([dc dc], [0 1],'g');
set(h1,'LineWidth',1.5);

if ( isempty( D ) == 0 )
    
    if (numel(unique(D))== 2)
        h5 = plot([D(1) D(1)], [ 0 pi(1) ],'c-','LineWidth', 1.5);
        h5 = plot([D(2) D(2)], [ 0 pi(2) ],'c-','LineWidth', 1.5);
    end
    
    if (numel(unique(D))== 1)
        h5 = plot([D(1) D(1)], [ 0 1 ],'m-','LineWidth', 1.5);
    end
    
    if (numel(unique(D))== 3)
        h5 = plot([D(1) D(1)], [ 0 pi(1) ],'m-','LineWidth', 1.5);
        h5 = plot([D(2) D(2)], [ 0 pi(2) ],'m-','LineWidth', 1.5);
        h5 = plot([D(3) D(3)], [ 0 pi(3) ],'m-','LineWidth', 1.5);
    end
    
end

lDhat = const - ( 1 + log( dc ) ) *sum( K/2 );

fp=fopen(dataout,'a');
nw=fprintf(fp,'%s - - - - - - - - - - - - - - - - - - - - - - - - - \n',' ');
fclose(fp);
fp=fopen(dataout,'a');
nw=fprintf(fp,'%s \n',' ');
nw=fprintf(fp,['Vector of True diffusion coefficients: ',num2str(D)]);
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'%s',['    Analysis of data']);
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'    Number of particles  = %11.0f',J);
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'%s',['    Analysis with one diffusion coefficient']);
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'    Diffusion coefficient estimate  = %11.2f',dc);
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'    Diffusion coeffic. stand. err.  = %11.4f',sd1);
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'    log likelihood  = %11.2f',lDhat);
nw=fprintf(fp,'%s \n \n',' ');
fclose(fp);



% Minus Log-likelihood calculation
 
 mllik2=@(x)-(sum ( reallog( x(1) * ( K/(2*x(2)) ).^(K/2) .* De.^(K/2-1) .* exp(-(1/(2*x(2)))*K.*De)./gamma(K/2)+...
                  (1-x(1)) * ( K/(2*x(3)) ).^(K/2) .* De.^(K/2-1) .* exp(-(1/(2*x(3)))*K.*De)./gamma(K/2)  ) ) )

  [x,fval,exitflag]=fmincon(mllik2,[0.5,D10,D20],[],[],[],[],[0 0 0 ],[1 Inf Inf],[], options);
loglik = - fval;
lDhat2 = - fval;
pi1 = x(1);
D1 = x(2);
D2 = x(3);
[D_est2, ind] = sort([D1 D2]);
p_est2 = [pi1 1-pi1];
p_est2 = p_est2(ind);
chi2 = 2 *( -fval - lDhat );
p = 2 * ( 1 - chi2cdf( chi2 , 2 ) );


fp=fopen(dataout,'a');
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'%s',['    Analysis with two diffusion coefficients']);
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'    First diffusion coeff.   = %11.2f ± %11.4f',D_est2(1),sd2(1));
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'    Second diffusion coeff.  = %11.2f ± %11.4f',D_est2(2), sd2(2) );
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'    Proportion of first diff coeff =  %8.4f ± %11.4f',p_est2(1), psd2(1) );
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'    Proportion of second diff coeff = %8.4f ± %11.4f',p_est2(2), psd2(2));
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'    log likelihood  = %11.2f',loglik);
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'    log likelihood improvement = %11.2f',lDhat2-lDhat);
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'    corresponding chisquare  = %11.2f',chi2);
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'    corresponding p-value  = %7.5f',p);
nw=fprintf(fp,'%s \n \n',' ');
fclose(fp);


h2 = plot([D1 D1],[0 pi1],'r');
set(h2,'LineWidth',1.5);
h2 = plot([D2 D2],[0 1-pi1],'r');
set(h2,'LineWidth',1.5);

 mllik3=@(x)-(sum ( reallog( x(1) * ( K/(2*x(2)) ).^(K/2) .* De.^(K/2-1) .* exp(-(1/(2*x(2)))*K.*De)./gamma(K/2)+...
                  x(3) * ( K/(2*x(4)) ).^(K/2) .* De.^(K/2-1) .* exp(-(1/(2*x(4)))*K.*De)./gamma(K/2) + ...
                  (1-x(1)-x(3)) * ( K/(2*x(5)) ).^(K/2) .* De.^(K/2-1) .* exp(-(1/(2*x(5)))*K.*De)./gamma(K/2)  ) ) )

  [x,fval,exitflag]=fmincon(mllik3,[.4,D10,.3,D20,D30],[1,0,1,0,0],1,[],[],[0 0 0 0 0],[1 Inf 1 Inf Inf],[], options);
loglik = - fval;
lDhat3 = - fval;
pi1 = x(1);
pi2 = x(3);
pi3 = 1-x(1)-x(3);
D1 = x(2);
D2 = x(4);
D3 = x(5);
[D_est3, ind] = sort([D1 D2 D3]);
p_est3 = [pi1 pi2 pi3];
p_est3 = p_est3(ind);
chi2 = 2 *( -fval - lDhat2 );
p = 2 * ( 1 - chi2cdf( chi2 , 2 ) );


fp=fopen(dataout,'a');
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'%s',['    Analysis with three diffusion coefficients']);
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'    First diffusion coeff.   = %11.2f ± %11.4f',D_est3(1), sd3(1));
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'    Second diffusion coeff.  = %11.2f ± %11.4f',D_est3(2), sd3(2));
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'    Third diffusion coeff.  = %11.2f ± %11.4f',D_est3(3), sd3(3));
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'    Proportion of first diff coeff =  %8.4f ± %11.4f',p_est3(1), psd3(1));
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'    Proportion of second diff coeff = %8.4f ± %11.4f',p_est3(2), psd3(2));
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'    Proportion of third diff coeff =  %8.4f ± %11.4f',p_est3(3), psd3(3) );
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'    log likelihood  = %11.2f',loglik);
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'    log likelihood improvement = %11.2f',lDhat3-lDhat2);
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'    corresponding chisquare  = %11.2f',chi2);
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'    corresponding p-value  = %7.5f',p);
nw=fprintf(fp,'%s \n \n',' ');
fclose(fp);

h3 = plot([D1 D1],[0 pi1],'k');
set(h3,'LineWidth',1.5);
h3 = plot([D2 D2],[0 pi2],'k');
set(h3,'LineWidth',1.5);
h3 = plot([D3 D3],[0 1-pi1-pi2],'k');
set(h3,'LineWidth',1.5);

 mllik4=@(x)-(sum ( reallog( x(1) * ( K/(2*x(2)) ).^(K/2) .* De.^(K/2-1) .* exp(-(1/(2*x(2)))*K.*De)./gamma(K/2)+ ...
                  x(3) * ( K/(2*x(4)) ).^(K/2) .* De.^(K/2-1) .* exp(-(1/(2*x(4)))*K.*De)./gamma(K/2) + ...
                  x(5) * ( K/(2*x(6)) ).^(K/2) .* De.^(K/2-1) .* exp(-(1/(2*x(6)))*K.*De)./gamma(K/2) + ...
                  (1-x(1)-x(3)-x(5)) * ( K/(2*x(7)) ).^(K/2) .* De.^(K/2-1) .* exp(-(1/(2*x(7)))*K.*De)./gamma(K/2)   ) ) )

  [x,fval,exitflag]=fmincon(mllik4,[.25,D10,.25,D20,.25,D30,D40],[1,0,1,0,1,0,0],1,[],[],[0 0 0 0 0 0 0],[1 Inf 1 Inf 1 Inf Inf],[], options);
loglik = - fval;
lDhat4 = - fval;
pi1 = x(1);
pi2 = x(3);
pi3 = x(5);
pi4 = 1-x(1)-x(3)-x(5);
D1 = x(2);
D2 = x(4);
D3 = x(6);
D4 = x(7);
[D_est4, ind] = sort([D1 D2 D3 D4]);
p_est4 = [pi1 pi2 pi3 pi4];
p_est4 = p_est4(ind);
chi2 = 2 *( -fval - lDhat3 );
p = 2 * ( 1 - chi2cdf( chi2 , 2 ) );


fp=fopen(dataout,'a');
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'%s',['    Analysis with four diffusion coefficients']);
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'    First diffusion coeff.   = %11.2f ± %11.4f',D_est4(1), sd4(1));
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'    Second diffusion coeff.  = %11.2f ± %11.4f',D_est4(2), sd4(2));
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'    Third diffusion coeff.  = %11.2f ± %11.4f',D_est4(3), sd4(3));
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'    Fourth diffusion coeff.  = %11.2f ± %11.4f',D_est4(4), sd4(4));
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'    Proportion of first diff coeff =  %8.4f ± %11.4f',p_est4(1), psd4(1));
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'    Proportion of second diff coeff = %8.4f ± %11.4f', p_est4(2),  psd4(2));
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'    Proportion of third diff coeff =  %8.4f ± %11.4f',p_est4(3),  psd4(3));
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'    Proportion of fourth diff coeff =  %8.4f ± %11.4f',p_est4(4),  psd4(4) );
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'    log likelihood  = %11.2f',loglik);
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'    log likelihood improvement = %11.2f',lDhat4-lDhat3);
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'    corresponding chisquare  = %11.2f',chi2);
nw=fprintf(fp,'%s \n \n',' ');
nw=fprintf(fp,'    corresponding p-value  = %7.5f',p);
nw=fprintf(fp,'%s \n \n',' ');
fclose(fp);

h4 = plot([D1 D1],[0 pi1],'b');
set(h4,'LineWidth',1.5);
h4 = plot([D2 D2],[0 pi2],'b');
set(h4,'LineWidth',1.5);
h4 = plot([D3 D3],[0 pi3],'b');
set(h4,'LineWidth',1.5);
h4 = plot([D4 D4],[0 1-pi1-pi2-pi3],'b');
set(h4,'LineWidth',1.5);






set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [0 0 0], ...
  'YColor'      , [0 0 0], ...
  'LineWidth'   , 1.6         );

set(gca,'FontSize',10)

if ( isempty( D ) == 0 )
    legend([h h5 h1 h2 h3 h4],'D SPRIA distribution','True Model','1 comp', '2 comp', '3 comp', '4 comp','Location','NorthEast')
    title(['Vector of True diffusion coefficients: ',num2str(D)])
    print('-dpsc2',figout2,'-append')
else
    legend([h h1 h2 h3 h4],'D SPRIA distribution','1 comp', '2 comp', '3 comp', '4 comp','Location','NorthEast')
    title(['Vector of True diffusion coefficients: ',num2str(D)])
    print('-dpsc2',figout2,'-append')
end