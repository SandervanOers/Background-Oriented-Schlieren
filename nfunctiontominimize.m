function F = nfunctiontominimize(Vx,datax_rawnonnan, datay_rawnonnan, gridxnonnan, gridynonnan, Lengths, Distance_to_mm, focal_length, Xcameracenter, Ycameracenter, take_y_into_account)


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
    phi_x_1_1 = asind(n_air/n_glass*sind(alpha+phi_x_1_0))-alpha;
    phi_x_2_1 = asind(n_air/n_glass*sind(alpha+phi_x_2_0))-alpha;
    phi_x_1_2 = phi_x_1_0;
    
%     if (take_y_into_account == 1)
%     phi_y_1_0 = phi_y_1;
%     phi_y_2_0 = phi_y_2;
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
    
    n = n_air*sind(alpha+phi_x_2_0)./(sind(atand(B)+alpha));
%     end

    F = n-1.333;

end