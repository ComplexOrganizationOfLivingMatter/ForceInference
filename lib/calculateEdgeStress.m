function [edgeStress] = calculateEdgeStress(actualEdgeInfo, numEdge, inverse, u, v)
%CALCULATEEDGESTRESS Summary of this function goes here
%   Detailed explanation goes here
    tension = actualEdgeInfo.TensionValue(numEdge);
    if u == 1 && ~inverse
        l_u = actualEdgeInfo.vertex2_X(numEdge) - actualEdgeInfo.vertex1_X(numEdge);
    else
        l_u = actualEdgeInfo.vertex1_Y(numEdge) - actualEdgeInfo.vertex2_Y(numEdge);
    end

    if v == 2 && ~inverse
        l_v = actualEdgeInfo.vertex2_Y(numEdge) - actualEdgeInfo.vertex1_Y(numEdge);
    else
        l_v = actualEdgeInfo.vertex1_X(numEdge) - actualEdgeInfo.vertex2_X(numEdge);
    end
    
    %l_ij = [actualEdgeInfo.vertex1_Y
    %l_ij = [l_u, l_v, 1];

    edgeLength = actualEdgeInfo.EdgeLength(numEdge);
    edgeStress = tension * ( l_u * l_v ) / edgeLength;
end

