% Supports data saved in a .mat file containing a structure "A" with at least the following fields:
% - Imgs : a MxMxN matrix containing N raster images with resolution MxM pixels
% - Sx : pixel size in the x direction
% - Sy : pixel size in the y direction
% - Tl : line dwell time used to collect the raster images
% - Tp : pixel dwell time used to collect the raster images


% load the file containing the data, for example by specifying the path
% load('C:/Users/mydata.mat')


Img = double( A.Imgs );
[ num , ~ , nimg ] = size( Img ); %get the resolution "num" of the images and how many they are "nimg" 
T = multithresh( Img , 2 );
thres = T( 1 );
thres1 = mean( T );
particleCount = 0; %initialize the number of particles "seen"
counter = 0; %counter used to go ahead until all the particles have been extracted
for n_im = 1 : nimg
    
    Img_tmp = double( Img( : , : , n_im ) ); %create a temporary copy of the image on which we can do anything
    
    if ( max( Img_tmp( : ) ) >  thres1  )  %we define a particle if the local maximum associated is at least the threshold thres1;
        particleCount = particleCount + 1; %if the maximum is higher then the some of the thresholds then there is at least one particle
        counter = 1;
    end
    
    
    %lines below will extract for each particle the signal sent from that
    %particle (run it and see the images)
    while (counter == 1)
        
        maxValue = max( Img_tmp( : ) ); %find the maximum
        [ rowsOfMaxes , colsOfMaxes ] = find( Img_tmp == maxValue ); %find the position of the maximum inside the image
        
        
        %the lines 133-190 are used to define the area around the local
        %maximum which correspond to the particle
        tempo = all( Img_tmp < thres , 2 ); %find all the columns which has all the values above the threshold
        indic = find( tempo );
        dim3 = max( indic( indic < rowsOfMaxes( 1 ) ) ) + 1; %extract the number of the first colums to the left of the max such that all the values are below the thres.
        if ( isempty( indic( indic < rowsOfMaxes( 1 ) ) ) )  %in case we are near the border we have to correct possible errors
            dim3 = + 1;
        end
        
        dim4 = min( indic( indic > rowsOfMaxes( 1 ) ) ) - 1; %extract the number of the first colums to the right of the max such that all the values are below the thres.
        if ( isempty( indic( indic > rowsOfMaxes( 1 ) ) ) ) %correction for edges
            dim4 = num;
        end
        
        
        %do the same for rows
        tempo = all( Img_tmp( dim3 : dim4 , : ) < thres , 1 );
        indic = find( tempo );
        dim1 = max( indic( indic < colsOfMaxes( 1 ) ) ) + 1;
        if ( isempty( indic( indic < colsOfMaxes( 1 ) ) ) )
            dim1 = + 1;
        end
        
        dim2 = min( indic( indic > colsOfMaxes(1) ) ) - 1;
        
        if ( isempty( indic( indic > colsOfMaxes( 1 ) ) ) )
            dim2 = num;
        end
        
        %redo for columns to avoid parallel particles
        tempo = all( Img_tmp( : , dim1 : dim2)  < thres , 2 );
        indic = find( tempo );
        dim3 = max( indic( indic < rowsOfMaxes( 1 ) ) ) + 1;
        if (isempty( indic( indic < rowsOfMaxes( 1 ) ) ) )
            dim3 = + 1;
        end
        
        dim4 = min( indic( indic > rowsOfMaxes(1) ) )-1;
        if (isempty(indic( indic > rowsOfMaxes(1) )))
            dim4 = num;
        end
        
        
        %if the particle has both dimensions bigger than one (more than one
        %line and one row)
        sizeimg = size( Img_tmp( dim3 : dim4 , dim1 : dim2 ) );
        if ( sizeimg( 1 ) > 1 && sizeimg( 2 ) > 1 )
            data{particleCount} = Img_tmp( dim3 : dim4 , dim1 : dim2 );  %save the particle in a struct        
            
        else
            particleCount = particleCount - 1;
        end
        
        Img_tmp( dim3 : dim4 , dim1 : dim2 ) = zeros( sizeimg( 1 ) , sizeimg( 2 ) ); %once the particle is saved, we delete it from the temporary image
        counter = counter - 1 ;
        if ( max( Img_tmp( : ) ) > thres1 ) %if in the image without the previous particles we have another "high" local maxima we add a particle
            particleCount = particleCount + 1 ;
            counter = counter + 1 ;
        end
        
        
    end
end

%these lines below will just create a movie with all the extracted
%particles, comment it if not needed
F( particleCount ) = struct('cdata',[],'colormap',[]);
f = figure('visible', 'off');
for j = 1 : particleCount
   imagesc( data{ j } )
   F(j) = getframe; 
end

implay( F )