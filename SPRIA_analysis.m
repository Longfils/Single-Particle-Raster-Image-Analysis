% perform SPRIA analysis and save results in a .txt file
% run this code after "extraction.m"

S = A.Sx; % pixel size
Tp = A.Tp; % pixel dwell time
Tl = A.Tl; % line time

estimates = zeros( particleCount , 3 ); %initialize structures where to save the estimates
Pos = cell( 1 , particleCount );

F( particleCount ) = struct('cdata',[],'colormap',[]);
f = figure('visible', 'off');

for particle = 1 : particleCount
    
    [ numline, num ] = size( data{ particle } );
    
    Pos{particle} = zeros(numline, 1);
    
    for j = 1:numline
        % line ( = position ) we analyse
        V = [ 1 : num ] * S;
        lin = data{ particle }( j , : );
        Pos{ particle }( j ) = ( V * lin' ) / sum( lin ) ; %here we compute the centroid for the j-th line
        
    end
    
    clf
    imagesc( data{ particle } ) %show the particle extracted
    hold on
    plot( Pos{ particle } / S , 1 : numline , 'k' ) %show the estimated trajectory
    F( particle ) = getframe;
    
    % extract the X-axis trajectory of the particle
    XE = Pos{ particle };
    % compute the length of the trajectory
    NN = size( XE , 1 );
    % initialize the mean square displacement
    DE = 0;
    number = 0; %number of lines used to compute the mean square displacement
    for nn = 2:NN
        DE = DE+(XE(nn)-XE(nn-1))^2; %add sequentially the displacements from one line to the next
        number = number + 1;
    end
    
    % normalize the displacement by the number of steps and the time
    % interval between steps.
    DE=DE/(2*(number)*Tl); %estimated diffusion coeff. for the particle
    
    %save the number of lines used and the MSD
    if ( number > 0 )
        if ( isnan(DE)~=1 && isinf(DE)~=1 && DE > 0 )
            estimates( particle , :)=[number  DE n_im];
        end
    end
    
end

implay( F )

%%

% check if mixture is needed



D_SPRIA = estimates(:,2);
K = estimates(:,1);
str = 'run_SPRIA.txt'; %name of the file where to save the analysis


formatSpec = 'SPRIA_scatterplot'; %name of scatterplot figure

fig1 = sprintf(formatSpec);

formatSpec = 'SPRIA_fit'; %name of SPRIA fitted models figure

fig2 = sprintf(formatSpec);


% Mmodel(D_SPRIA, K , str,fig1,fig2,.5,1,3,5,7, D)
nboot = 50;
D_est1 = zeros( nboot , 1 );
D_est2 = zeros( nboot , 2 );
D_est3 = zeros( nboot , 3 );
D_est4 = zeros( nboot , 4 );
p_est2 = zeros( nboot , 2 );
p_est3 = zeros( nboot , 3 );
p_est4 = zeros( nboot , 4 );

%here we perform the bootstrap method to obtain errors for the estimates 
for i = 1:nboot
    smp = datasample( 1 : numel( D_SPRIA ) , numel( D_SPRIA ) ); %resample the data with replacement
    D_tmp = D_SPRIA( smp );
    K_tmp = K( smp );
    
    [ D_est1t ,D_est2t , p_est2t ,D_est3t , p_est3t , D_est4t , p_est4t ]= Mmodel_boot( D_tmp , K_tmp , 2 , 4 , 6 , 8 );
    D_est1( i , : ) = D_est1t;
    D_est2( i , : ) = D_est2t;
    D_est3( i , : ) = D_est3t;
    D_est4( i , : ) = D_est4t;
    p_est2( i , : ) = p_est2t;
    p_est3( i , : ) = p_est3t;
    p_est4( i , : ) = p_est4t;
end

sd1 = std(D_est1);
sd2 = std(D_est2);
sd3 = std(D_est3);
sd4 = std(D_est4);
psd2 = std(p_est2);
psd3 = std(p_est3);
psd4 = std(p_est4);

Mmodel( D_SPRIA , K , str , fig1 , fig2 , 2 , 4 , 6 , 8 , [], [] , sd1 , sd2 , sd3 , sd4 , psd2 , psd3 , psd4 );