function circledrop(n_input,r_input,x_input)
%CIRCLEDROP Plots a user-given number of circles with user-given radius one
%by one; the circles arrange themselves as they were falling from infinite
%height and stopped upon touching y = 0 or another circle. The range of
%possible x values for the circles goes from 0 to a user-defined value.

    n = n_input;
    r = r_input;
    xArr = x_input*rand(n);
    xRec = []; 
    yRec = [];

    xInd = [];
    xSearch = [];
    ySearch = [];
    
    figure
    daspect([1 1 1])

    for i = 1:n
        xInd = find(abs(xRec-xArr(i)) <= r*2);
        xSearch = xRec(xInd);
        ySearch = yRec(xInd);
        yCalc = ySearch;

        if isempty(xSearch)
            xRec = [xRec xArr(i)];
            yRec = [yRec r];
        else
            xRec = [xRec xArr(i)];
            for j = 1:numel(xSearch)
                yCalc(j) = ySearch(j)+sqrt(-(xArr(i)-xSearch(j))^2+(2*r)^2);
            end
            yMe = max(yCalc);
            yRec = [yRec yMe];
        end
        
        hold on
        rad = 0:pi/50:2*pi;
        xCir = r * cos(rad) + xRec(i);
        yCir = r * sin(rad) + yRec(i);
        plot(xCir, yCir);
        hold off

    end
end

