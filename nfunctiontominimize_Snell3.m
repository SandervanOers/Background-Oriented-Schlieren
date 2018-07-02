function F = nfunctiontominimize_Snell3(Vx,datax_rawnonnan, datay_rawnonnan, gridxnonnan, gridynonnan, Lengths, Distance_to_mm, focal_length, Xcameracenter, Ycameracenter, take_y_into_account)
    % Refractive Indices
    n_air = 1;
    n_glass = 1.5;

    % Distances 
%     L_camera = Lengths(1);
%     L_camera = Vx(3);
    L_glass = Lengths(2);
    L_test = Lengths(3);
    L_screen = Lengths(4);

    
    Xbarcenter = Vx(1);
    distance_ccd_lens = 1/(1/(focal_length*1e-3)-1/(Vx(2)));
    L_tot = Vx(2)*cosd(atand((Xcameracenter-Xbarcenter)*Distance_to_mm/1e3/distance_ccd_lens));
    L_camera = L_tot-2*L_glass-L_test-L_screen;

    Ybarcenter = Vx(3);
%     L_magn = Vx(3);
%     Ycameracenter = (1080/2);

    alpha = atand((Xcameracenter-Xbarcenter)*Distance_to_mm/1e3/(distance_ccd_lens));
    
    beta = atand((Ycameracenter-Ybarcenter)*Distance_to_mm/1e3/(distance_ccd_lens));
    % Distance from Center [in pixels]
    DistanceX = gridxnonnan-Xcameracenter;
    DistanceY = gridynonnan-Ycameracenter;
    % Distance = sqrt(DistanceX.^2+DistanceY.^2);

    phi_x_1 = atand(DistanceX*Distance_to_mm/1e3/(distance_ccd_lens));
    phi_y_1 = atand(DistanceY*Distance_to_mm/1e3/(distance_ccd_lens));

    phi_x_2 = atand((DistanceX+datax_rawnonnan)*Distance_to_mm/1e3/(distance_ccd_lens));
    phi_y_2 = atand((DistanceY+datay_rawnonnan)*Distance_to_mm/1e3/(distance_ccd_lens));

%     beta  = atand((Ycameracenter-Ybarcenter)*Distance_to_mm/1e3/(distance_ccd_lens));

    phi_x_1_0 = phi_x_1;
    phi_x_2_0 = phi_x_2;
    phi_y_1_0 = phi_y_1;
    phi_y_2_0 = phi_y_2;    
    
    [T] = create_TransformationMatrices2(alpha,beta);
    [phi_x_1_1, phi_y_1_1, phi_x_2_1, phi_y_2_1] = SnellsLaw(phi_x_1_0, phi_y_1_0, phi_x_2_0, phi_y_2_0, alpha, n_air, n_glass,T);
    
    
%     phi_x_1_1 = asind(n_air/n_glass*sind(alpha+phi_x_1_0))-alpha;
%     phi_x_2_1 = asind(n_air/n_glass*sind(alpha+phi_x_2_0))-alpha;
    phi_x_1_2 = phi_x_1_0; phi_y_1_2 = phi_y_1_0;
    

    A = (L_camera+L_screen)* (tand(phi_x_1_0)./(cosd(alpha)*cosd(beta)-cosd(beta)*tand(phi_x_1_0)*sind(alpha)-cosd(alpha)*tand(phi_y_1_0)*sind(beta)) ...
            -tand(phi_x_2_0)./(cosd(alpha)*cosd(beta)-cosd(beta)*tand(phi_x_2_0)*sind(alpha)-cosd(alpha)*tand(phi_y_2_0)*sind(beta))) ...
        + 2*L_glass * (tand(phi_x_1_1)./(cosd(alpha)*cosd(beta)-cosd(beta)*tand(phi_x_1_1)*sind(alpha)-cosd(alpha)*tand(phi_y_1_1)*sind(beta)) ...
            -tand(phi_x_2_1)./(cosd(alpha)*cosd(beta)-cosd(beta)*tand(phi_x_2_1)*sind(alpha)-cosd(alpha)*tand(phi_y_2_1)*sind(beta))) ...
        + L_test *tand(phi_x_1_2)./(cosd(alpha)*cosd(beta)-cosd(beta)*tand(phi_x_1_2)*sind(alpha)-cosd(alpha)*tand(phi_y_1_2)*sind(beta));
    %B = cosd(alpha)./(L_test./A+sind(alpha));
    
    n2 = zeros(length(A),1);
    stepsize = 100;
    for i=1:stepsize:length(A)
        ICC = [1.333];
        lbb = [1.3];
        ubb = [1.4];
    
        [dn] = @(n)SnellsLaw_mini3(phi_x_1_1(i), phi_y_1_1(i), phi_x_2_1(i), phi_y_2_1(i), alpha, beta, n_glass, n, A(i), T, L_test);  
        lsqOpts = optimoptions('lsqnonlin' , 'StepTolerance', 1e-12,'MaxFunctionEvaluations', 2000,'FunctionTolerance', 1e-12, 'OptimalityTolerance',1e-12, 'MaxIterations', 1e3,'display','off');
         [Vn]  = lsqnonlin(dn,ICC, lbb, ubb,lsqOpts);
%          [round(i/length(B)*100) Vn  ]
         
         n2(i) = Vn;
    end
%     [phi_x_1_2, phi_y_1_2, phi_x_2_2, phi_y_2_2] = SnellsLaw(phi_x_1_1, phi_y_1_1, phi_x_2_1, phi_y_2_1, alpha, n_glass, Vn);
    n = n2(n2>0);
%     n(n==0)=[];
%     figure;plot(n,'k+');
% 
% GridX = gridxnonnan(1:stepsize:end);
% GridY = gridynonnan(1:stepsize:end);
% tri = delaunay(GridX,GridY);
% 
%      figure
%     hold all
%     trisurf(tri, GridX, GridY, n);
%     colormap(cbrewer('seq', 'OrRd', 100));
%     shading interp
%     colorbar EastOutside
%     title('n')
%     xlabel('Horizontal Position [pixels]')
%     ylabel('Vertical Position [pixels]')    
%     axis tight
%     axis equal
%     ax=gca;
% ax.FontSize = 20;
% 	hold off
%     
%     pause(0.1);
    F = n-1.3345;
    
%     F = n-1.334;
    [Vx]

end