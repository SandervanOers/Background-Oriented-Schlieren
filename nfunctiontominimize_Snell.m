function F = nfunctiontominimize_Snell(Vx,datax_rawnonnan, datay_rawnonnan, gridxnonnan, gridynonnan, Lengths, Distance_to_mm, focal_length, Xcameracenter, Ycameracenter, take_y_into_account)
    Vx

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

%     Ybarcenter = Vx(2);
%     L_magn = Vx(3);
%     Ycameracenter = (1080/2);

    alpha = atand((Xcameracenter-Xbarcenter)*Distance_to_mm/1e3/(distance_ccd_lens));
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
    
    [T] = create_TransformationMatrices(alpha);
    [phi_x_1_1, phi_y_1_1, phi_x_2_1, phi_y_2_1] = SnellsLaw(phi_x_1_0, phi_y_1_0, phi_x_2_0, phi_y_2_0, alpha, n_air, n_glass);
    
    
%     phi_x_1_1 = asind(n_air/n_glass*sind(alpha+phi_x_1_0))-alpha;
%     phi_x_2_1 = asind(n_air/n_glass*sind(alpha+phi_x_2_0))-alpha;
    phi_x_1_2 = phi_x_1_0;
    
%     if (take_y_into_account == 1)
% 
%     phi_y_1_1 = asind(n_air/n_glass*sind(phi_y_1_0));
%     phi_y_2_1 = asind(n_air/n_glass*sind(phi_y_2_0));
%     phi_y_1_2 = phi_y_1_0;
%     phi_y_2_2 = asind(n_air/1.333*sind(phi_y_2_0));
% %     
% %     A = (L_camera+L_screen)* (tand(phi_x_1_0)./(cosd(alpha)-tand(phi_x_1_0)*sind(alpha))./cosd(phi_y_1_0)-tand(phi_x_2_0)./(cosd(alpha)-tand(phi_x_2_0)*sind(alpha))./cosd(phi_y_2_0)) ...
% %         + 2*L_glass * (tand(phi_x_1_1)./(cosd(alpha)-tand(phi_x_1_1)*sind(alpha))./cosd(phi_y_1_1)-tand(phi_x_2_1)./(cosd(alpha)-tand(phi_x_2_1)*sind(alpha))./cosd(phi_y_2_1)) ...
% %         + L_test *tand(phi_x_1_2)./(cosd(alpha)-tand(phi_x_1_2)*sind(alpha))./cosd(phi_y_1_2);
% %     B = cosd(alpha)./(L_test./A./cosd(phi_y_2_2)+sind(alpha));
% %     
% %     n = n_air*sind(alpha+phi_x_2_0)./(sind(atand(B)+alpha));
% %     
%         Asol = (L_camera+L_screen)* (tand(phi_x_1_0)./(cosd(alpha)-tand(phi_x_1_0)*sind(alpha))./cosd(phi_y_1_0)-tand(phi_x_2_0)./(cosd(alpha)-tand(phi_x_2_0)*sind(alpha))./cosd(phi_y_2_0)) ...
%         + 2*L_glass * (tand(phi_x_1_1)./(cosd(alpha)-tand(phi_x_1_1)*sind(alpha))./cosd(phi_y_1_1)-tand(phi_x_2_1)./(cosd(alpha)-tand(phi_x_2_1)*sind(alpha))./cosd(phi_y_2_1)) ...
%         + L_test *tand(phi_x_1_2)./(cosd(alpha)-tand(phi_x_1_2)*sind(alpha))./cosd(phi_y_1_2);
% for i=1:10:length(Asol)
%     myfun = @(n)nfunction_ydependece(n, alpha, phi_x_2_0(i), L_test, Asol(i), n_air, phi_y_2_0(i));
%     % myfun = @(n)nfunction_ydependece(n, alpha, phi_x_2_0, L_test, Asol, n_air, phi_y_2_0);
%     x(i) = fzero(myfun, 1.333);
%     % [round(i/length(Asol)*100) x(i)]
% end
% n = x';

%     else
    A = (L_camera+L_screen)* (tand(phi_x_1_0)./(cosd(alpha)-tand(phi_x_1_0)*sind(alpha))-tand(phi_x_2_0)./(cosd(alpha)-tand(phi_x_2_0)*sind(alpha))) ...
        + 2*L_glass * (tand(phi_x_1_1)./(cosd(alpha)-tand(phi_x_1_1)*sind(alpha))-tand(phi_x_2_1)./(cosd(alpha)-tand(phi_x_2_1)*sind(alpha))) ...
        + L_test *tand(phi_x_1_2)./(cosd(alpha)-tand(phi_x_1_2)*sind(alpha));
    B = cosd(alpha)./(L_test./A+sind(alpha));
    
    n2 = zeros(length(B),1);
    stepsize = 100;
    for i=1:stepsize:length(B)
        ICC = [1.333];
        lbb = [1.3];
        ubb = [1.4];
    
        [dn] = @(n)SnellsLaw_mini(phi_x_1_1(i), phi_y_1_1(i), phi_x_2_1(i), phi_y_2_1(i), alpha, n_glass, n, B(i), T);  
        lsqOpts = optimoptions('lsqnonlin' , 'StepTolerance', 1e-12,'MaxFunctionEvaluations', 2000,'FunctionTolerance', 1e-12, 'OptimalityTolerance',1e-12, 'MaxIterations', 1e3,'display','off');
         [Vn]  = lsqnonlin(dn,ICC, lbb, ubb,lsqOpts);
%          [round(i/length(B)*100) Vn  ]
         
         n2(i) = Vn;
    end
%     [phi_x_1_2, phi_y_1_2, phi_x_2_2, phi_y_2_2] = SnellsLaw(phi_x_1_1, phi_y_1_1, phi_x_2_1, phi_y_2_1, alpha, n_glass, Vn);
    n = n2(n2>0);
    n(n==0)=[];
    figure;plot(n,'k+');

GridX = gridxnonnan(1:stepsize:end);
GridY = gridynonnan(1:stepsize:end);
tri = delaunay(GridX,GridY);

     figure
    hold all
    trisurf(tri, GridX, GridY, n);
    colormap(cbrewer('seq', 'OrRd', 100));
    shading interp
    colorbar EastOutside
    title('n')
    xlabel('Horizontal Position [pixels]')
    ylabel('Vertical Position [pixels]')    
    axis tight
    axis equal
    ax=gca;
ax.FontSize = 20;
	hold off
    
    pause(0.1);
%     F = n-1.333;
    
    F = n-1.335;

end