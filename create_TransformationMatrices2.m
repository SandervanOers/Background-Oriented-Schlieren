function [T] = create_TransformationMatrices2(alpha,beta)
    % Transformation Matrices
    % Invert Z axis
    T1 = [1,0,0;0,1,0;0,0,-1];
    % Rotate XZ-plane by alpha
    alpha = - alpha;
    T2 = [cosd(alpha), 0, -sind(alpha);0, 1, 0;sind(alpha),0,cosd(alpha)];
    % Rotate YZ-plane by beta
    beta = - beta;
    T3 = [1, 0, 0; 0, cosd(beta), -sind(beta); 0, sind(beta), cosd(beta)];
%     alpha = - alpha;
    % Multiply Transformation Matrices
    T = T1*T2*T3;
%     T = T1*T3*T2;
end

