function [phi_x_1_1, phi_y_1_1, phi_x_2_1, phi_y_2_1] = SnellsLaw(phi_x_1_0, phi_y_1_0, phi_x_2_0, phi_y_2_0, alpha, n1, n2, T)

    % Form Direction Vectors in the Camera Coordinate System
    i_1 = [tand(phi_x_1_0), tand(phi_y_1_0), ones(length(phi_x_1_0),1)];
    i_2 = [tand(phi_x_2_0), tand(phi_y_2_0), ones(length(phi_x_2_0),1)];

%     i_1 = [tand(phi_x_1_0),  zeros(length(phi_x_1_0),1), ones(length(phi_x_1_0),1)];
%     i_2 = [tand(phi_x_2_0), tand(phi_y_2_0), ones(length(phi_x_2_0),1)];
    % Transformation Matrices
%     % Invert Z axis
%     T1 = [1,0,0;0,1,0;0,0,-1];
%     % Rotate XZ-plane by alpha
%     alpha = - alpha;
%     T2 = [cosd(alpha), 0, -sind(alpha);0, 1, 0;sind(alpha),0,cosd(alpha)];
% %     alpha = - alpha;
%     % Multiply Transformation Matrices
%     T = T1*T2;
%     [T] = create_TransformationMatrices(alpha);
    % Form Direction Vectors in the Snell's Law Coordinate System
    I1 = T*i_1';
    I2 = T*i_2';
    % Normalize Direction Vectors
    I1n = I1./vecnorm(I1);
    I2n = I2./vecnorm(I2);
    % Normal Direction Vector in the Snell's Law Coordinate System
    n = [0,0,1];
    % Calculate Transmission Direction Vectors
    % cos theta_i 
%     cti1 = -dot(I1n,n);
%     cti2 = -dot(I2n,n);
    cti1 = I1n(3,:);
    cti2 = I2n(3,:);
    % sin^2 theta_t
    sst1 = (n1/n2)^2*(1-cti1.^2);
    sst2 = (n1/n2)^2*(1-cti2.^2);
    % transmission direction vector // Snell's Law
    t1 = n1/n2 *I1n - (n1/n2 * cti1 + sqrt(1-sst1)).*n';
    t2 = n1/n2 *I2n - (n1/n2 * cti2 + sqrt(1-sst2)).*n';
    % Transform back to Camera Coordinate System
    tt1 = inv(T)*t1;
    tt2 = inv(T)*t2;
    % Obtain Angles
    phi_x_1_1 = atand(tt1(1,:)./tt1(3,:))';
    phi_y_1_1 = atand(tt1(2,:)./tt1(3,:))';
    phi_x_2_1 = atand(tt2(1,:)./tt2(3,:))';
    phi_y_2_1 = atand(tt2(2,:)./tt2(3,:))';
    
end

