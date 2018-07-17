%%% Conway's Game of Life - B3/S23
%%% 24 June 2016
%%
clear all
close all
%%
N = 200;
midP = ceil(N/2);
a = (1:1:N);
X = repmat(a,N,1);
Y = repmat(a',1,N);
S = zeros(N,N);
sizes = zeros(N,N);
nElm = numel(S);
plotX = reshape(X,1,nElm);
plotY = reshape(Y,1,nElm);
plotS = zeros(1,nElm);
%%
seed = randi(2,9,9)-1;
%%
S(midP-4:midP+4,midP-4:midP+4) = seed;
transf = S;
nGen = 10000;
cmap = colormap(copper(9));
for i = 1:1:nGen
    transf = S;
    for x = 2:1:N-1
        for y = 2:1:N-1
            neighbN = 0;
            if S(x-1,y+1) == 1
            neighbN = neighbN+1;
            end
            if S(x,y+1) == 1
            neighbN = neighbN+1;
            end
            if S(x+1,y+1) == 1
            neighbN = neighbN+1;
            end
            if S(x-1,y) == 1
            neighbN = neighbN+1;
            end
            if S(x+1,y) == 1
            neighbN = neighbN+1;
            end
            if S(x-1,y-1) == 1
            neighbN = neighbN+1;
            end
            if S(x,y-1) == 1
            neighbN = neighbN+1;
            end
            if S(x+1,y-1) == 1
            neighbN = neighbN+1;
            end
            sizes(x,y) = neighbN;
            if S(x,y) == 0 && neighbN == 3
                transf(x,y) = 1;
            elseif S(x,y) == 1 && (neighbN < 2 || neighbN > 3)
                transf(x,y) = 0;
            end
        end
    end
    plotS = reshape(S,1,nElm);
    ind = find(S == 1);
    
    figure(1)
    scatter(plotX(ind),plotY(ind),10+(10*sizes(ind)),'ko')
    axis square
    box on
    axis([1 N 1 N])
    title(['Generation ',num2str(i)])
    xlabel('X')
    ylabel('Y')
    S = transf;
    pause(0.05)
end
