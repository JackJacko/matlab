function [surfArea,totalArea,maxArea,minArea,errArea] = occlusion3d(Nref,Dref,Dmin,slope,volLen,res)
%OCCLUSION3D Maps and evaluates visual occlusion of a 2d square plane by a 
%power law distribution of spheres in a cubic volume.
%   Using 5 inputs, i.e. Nref (number of reference particles), Dref
%   (diameter of the reference particles), Dmin (minimum diameter cutoff),
%   slope (slope of the power law distribution), volLen (side length of
%   cubic volume), res (resolution of the squarelet grid) this code
%	randomly generates a population of spheres within a cubic volume and 
%   calculates the visual occlusion caused by the spheres (i.e. their 
%   projection). The reference particles are used to define the limits of 
%   the population (they represent the largest desired particles in the 
%   population and their desired number), as is the minimum diameter 
%   cutoff. The function outputs the area of the square plane, the total 
%   unoccluded area projected by the spheres, and the maximum and minimum 
%   occluded area projected by the spheres (the occluded area is calculated 
%   by summing the combined area covered by the spheres within a grid of 
%   squarelets with set resolution - the difference between max and min
%   area estimates the error).
%
%   Nref - Desired number of reference diameter particles
%   Dref - Reference diameter (microns)
%   Dmin - Minimum diameter (microns)
%   slope  - Slope of the distribution
%   volLen - Side length of the volume (microns)
%   res - Grid resolution
%
%	Jacopo Agagliate
%	University of Strathclyde
%	7 June 2016

%   2d square plane

    Pref = (Dref/Dmin)^(-slope);                                            %Probability of generating a particle with Dref diameter;
    N = ceil(Nref/Pref);                                                    %Number of particles generated - calculated so that on average Nref reference particles are being generated;
    gridLine = linspace(0,volLen,res);                                      %Node grid for total visible area estimation
    gLN = numel(gridLine);                                                  %Number of nodes in one grid row
    areaElm = (gridLine(2)-gridLine(1))^2;                                  %Area of grid area element (microns squared)
    [gridX,gridY] = meshgrid(gridLine,gridLine);                            %Grid matrices generation - these hold the coordinates of each grid point
    nodeCheckXY = zeros(gLN,gLN);                                           %This matrix keeps track of which nodes are covered
    sqrletsXY = ones(gLN-1,gLN-1);                                          %This matrix keeps track of whether the area elements of the grid are free, partially covered or fully covered
    sqrN = (gLN-1)^2;                                                       %Number of area elements in the grid

    preDvar = zeros(N,1);                                                   %Particle sizes are pre-generated and sorted so that large particles are put in first
    for i = 1:1:N
        preDvar(i) = volLen/2+1;                                            %Maximum particle diameter primer
        while preDvar(i)>volLen/2                                           %Maximum particle diameter allowed
            preDvar(i) = Dmin*rand(1)^(-1/slope);                           %Particle diameter generation
        end
    end
    Dvar = sort(preDvar,'descend');

    tic
    for i = 1:1:N
        toc
        if i == 1
            x0 = (Dvar(i)/2)+rand(1)*(volLen-Dvar(i));                      %X coord. of first particle center - constrained to be wholly contained within the volume boundaries
            y0 = (Dvar(i)/2)+rand(1)*(volLen-Dvar(i));                      %Y coord. of first particle center - constrained to be wholly contained within the volume boundaries
            z0 = (Dvar(i)/2)+rand(1)*(volLen-Dvar(i));                      %Z coord. of first particle center - constrained to be wholly contained within the volume boundaries
            X = x0;                                                         %First element of X coord. array
            Y = y0;                                                         %First element of Y coord. array
            Z = z0;                                                         %First element of Z coord. array
            D = Dvar(i);                                                    %First element of particle diameter array
            A = pi*(Dvar(i)/2)^2;                                           %First element of particle projected area array
        else
            chk = 1;                                                        %Particles are generated so that they can't share the same volume - primer
            while chk == 1                                                  %Particles are generated so that they can't share the same volume - loop
                chk = 0;
                x0 = (Dvar(i)/2)+rand(1)*(volLen-Dvar(i));                  %X coord. of particle center - constrained to be wholly contained within the volume boundaries
                y0 = (Dvar(i)/2)+rand(1)*(volLen-Dvar(i));                  %Y coord. of particle center - constrained to be wholly contained within the volume boundaries
                z0 = (Dvar(i)/2)+rand(1)*(volLen-Dvar(i));                  %Z coord. of particle center - constrained to be wholly contained within the volume boundaries
                for u = 1:1:size(X,1)
                    if sqrt(((X(u)-x0)^2)+((Y(u)-y0)^2)+((Z(u)-z0)^2))< ...
                            (Dvar(i)+D(u))/2                                %This checks the distance between the new particle center and the centers of all other particles
                        chk = 1;
                        break
                    end
                end
            end
            X = [X;x0];                                                     %X coord. of new particle is appended to X coord. array
            Y = [Y;y0];                                                     %Y coord. of new particle is appended to X coord. array
            Z = [Z;z0];                                                     %Z coord. of new particle is appended to X coord. array
            D = [D;Dvar(i)];                                                %Diameter of new particle is appended to particle diameter array
            A = [A;pi*(Dvar(i)/2)^2];                                       %Projected area of new particle is appended to particle projected area array
        end
        for c1 = 1:1:gLN
            for c2 = 1:1:gLN
                if sqrt(((x0-gridX(c1,c2))^2)+((y0-gridY(c1,c2))^2))< ...
                        Dvar(i)/2                                           %This checks which grid nodes are covered by the projecte area of the particle
                    nodeCheckXY(c1,c2) = 1;                                 %Covered nodes are marked as 1
                end
            end
        end
    end

    for cx = 1:1:gLN-1
        for cy = 1:1:gLN-1
            if nodeCheckXY(cx,cy)+nodeCheckXY(cx+1,cy)+ ...
                    nodeCheckXY(cx,cy+1)+nodeCheckXY(cx+1,cy+1) == 0        %Grid area elements are checked one by one - if sum of four pertaining nodes is 0, area is free
                    
                sqrletsXY(cx,cy) = 0;                                       %Free area elements are marked as 0
            elseif nodeCheckXY(cx,cy)+nodeCheckXY(cx+1,cy)+ ...             %Grid area elements are checked one by one - if sum of four pertaining nodes is 4, area is fully covered
                    nodeCheckXY(cx,cy+1)+nodeCheckXY(cx+1,cy+1) == 4
                sqrletsXY(cx,cy) = 2;                                       %Fully covered area elements are marked as 2
            end                                                             %In all other cases, area is partially covered - area element is left marked as 1
        end
    end

    surfArea = volLen^2;                                                    %2D square plane surface area
    totalArea = sum(A);                                                     %Total particle projected area if unoccluded
    maxArea = (sqrN-sum(sqrletsXY(:)==0))*areaElm;                          %Maximum total particle occluded projected area
    minArea = sum(sqrletsXY(:)==2)*areaElm;                                 %Minimum total particle occluded projected area
    errArea = maxArea-minArea;                                          %Total particle occluded projected area uncertainty (improve with resolution)

    %% FIGURES
    % 3D rendering
    [x,y,z] = sphere;
    figure;
    for s = 1:1:N
        surf(x*(D(s)/2)+X(s),y*(D(s)/2)+Y(s),z*(D(s)/2)+Z(s))
        hold on
        xlabel('x (\mum)')
        ylabel('y (\mum)')
        zlabel('z (\mum)')
        axis square
    end
    %%
    % 2D contour plot
    figure;
    contour(sqrletsXY)
    xlabel('x (Area element units)')
    ylabel('y (Area element units)')
end

