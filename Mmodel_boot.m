function [ D_est1 , D_est2 , p_est2 , D_est3 , p_est3 , D_est4 , p_est4 ]= Mmodel_boot( De , K , D10 , D20 , D30 , D40 )


% matlab function for maximum likelihood estimation of diffusion
% coefficients for mixture models

% The program computes the maximum likelihood estimation of the diffusion coefficient
% from  models with one,two,three, and four diffusion coefficients

% Input parameters
% De: estimated diffusion coefficients from SPRIA analysis
% K: number of increments used to compute each of the value in De
% D10, D20, D30, D40: start values for the different diffusion
% coefficients (use uninformal guesses)

% Output parameters
% D_est1: is the estimated diffusion coefficient from a model with one
% diffusing component;
% D_est2: are the two estimated diffusion coefficients from a model with
% two diffusing component;
% p_est2: are the two estimated proportions of particle with diffusion coefficient equal to D_est2 for a model with
% two diffusing component;
% D_est3, p_est3, D_est4, p_est4: they are the corresponding quantities as
% D_est2, p_est2 for models with three and four diffusing components

options = optimset('MaxFunEvals', 10^4, 'MaxIter', 10^4,'Display','off', 'TolFun', 10^(-6) );

% model with one component
sumK = sum( K );
dc = De' * K / sumK;
D_est1 = dc;


% model with two components
% Minus Log-likelihood calculation

mllik2=@(x) - (sum ( reallog( x( 1 ) * ( K / ( 2 * x( 2 ) ) ).^( K / 2) .* De .^( K / 2 - 1) .* exp( - ( 1 / ( 2 * x( 2 ) ) ) * K .* De) ./ gamma( K / 2 ) + ...
    ( 1 - x( 1 ) ) * ( K /( 2 * x( 3 ) ) ) .^( K / 2) .* De.^( K / 2 - 1) .* exp( -( 1 / ( 2 * x( 3 ) ) ) * K.* De) ./ gamma( K / 2 )  ) ) );

x=fmincon( mllik2 , [.5,D10,D20] , [] , [] , [] , [] , [ 0 0 0 ] , [ 1 Inf Inf ] , [] , options) ;

pi1 = x(1);
D1 = x(2);
D2 = x(3);
[D_est2, ind] = sort( [ D1 D2 ] );
p_est2 = [ pi1 1-pi1 ];
p_est2 = p_est2( ind );


% model with three components
mllik3=@(x)-(sum ( reallog( x(1) * ( K/(2*x(2)) ).^(K/2) .* De.^(K/2-1) .* exp(-(1/(2*x(2)))*K.*De)./gamma(K/2)+...
    x(3) * ( K/(2*x(4)) ).^(K/2) .* De.^(K/2-1) .* exp(-(1/(2*x(4)))*K.*De)./gamma(K/2) + ...
    (1-x(1)-x(3)) * ( K/(2*x(5)) ).^(K/2) .* De.^(K/2-1) .* exp(-(1/(2*x(5)))*K.*De)./gamma(K/2)  ) ) );

x=fmincon( mllik3 , [ .4 , D10 , .3 , D20 , D30 ], [ 1 , 0 , 1 , 0 , 0 ] , 1 , [] , [] , [0 0 0 0 0] , [1 Inf 1 Inf Inf] , [] , options );
pi1 = x(1);
pi2 = x(3);
pi3 = 1 - x(1) - x(3);
D1 = x(2);
D2 = x(4);
D3 = x(5);
[D_est3, ind] = sort([D1 D2 D3]);
p_est3 = [pi1 pi2 pi3];
p_est3 = p_est3(ind);


% model with four components
mllik4=@(x)-(sum ( reallog( x(1) * ( K/(2*x(2)) ).^(K/2) .* De.^(K/2-1) .* exp(-(1/(2*x(2)))*K.*De)./gamma(K/2)+ ...
    x(3) * ( K/(2*x(4)) ).^(K/2) .* De.^(K/2-1) .* exp(-(1/(2*x(4)))*K.*De)./gamma(K/2) + ...
    x(5) * ( K/(2*x(6)) ).^(K/2) .* De.^(K/2-1) .* exp(-(1/(2*x(6)))*K.*De)./gamma(K/2) + ...
    (1-x(1)-x(3)-x(5)) * ( K/(2*x(7)) ).^(K/2) .* De.^(K/2-1) .* exp(-(1/(2*x(7)))*K.*De)./gamma(K/2)   ) ) );

x=fmincon(mllik4,[.25,D10,.25,D20,.25,D30,D40],[1,0,1,0,1,0,0],1,[],[],[0 0 0 0 0 0 0],[1 Inf 1 Inf 1 Inf Inf],[], options);
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