function [T] = create_TransformationMatrices(alpha)
    % Transformation Matrices
    % Invert Z axis
    T1 = [1,0,0;0,1,0;0,0,-1];
    % Rotate XZ-plane by alpha
    alpha = - alpha;
    T2 = [cosd(alpha), 0, -sind(alpha);0, 1, 0;sind(alpha),0,cosd(alpha)];
%     alpha = - alpha;
    % Multiply Transformation Matrices
    T = T1*T2;
end

