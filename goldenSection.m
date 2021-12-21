function [minPoint, minValue, iterNum] = goldenSection(g_start,g_end,func)
% Golden Section Method for exact line search
%   此处显示详细说明
    epsilon=1e-4;
    ratio = sqrt(5)/2-0.5;
    intervalLen = g_end - g_start;
    middleL = g_start + (1 - ratio) *intervalLen;
    middleR = g_start +ratio*intervalLen;
    iterNum = 0;
    while(intervalLen >= epsilon)
        if(func(middleL) > func(middleR))
            g_start = middleL;
            intervalLen = g_end - g_start;
            middleL = middleR;
            middleR = g_start + ratio * intervalLen;
        else
            g_end = middleR;
            intervalLen = g_end - g_start;
            middleR = middleL;
            middleL = g_start + (1-ratio)*intervalLen;
        end
        iterNum = iterNum+ 1;
    end
        minPoint = (g_start + g_end ) /2;
        minValue = func(minPoint);
%         print("极小值点",minPoint);
%         print("极小值",minValue)
end

