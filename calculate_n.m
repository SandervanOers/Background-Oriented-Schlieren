clc;clear;close all;
% path='/home/avanoers/Documents/DIC/matrixstorage/calibration/';
% path='/home/avanoers/Documents/DIC/New Derivation/Calibration - Air/';
% path = '/home/avanoers/Documents/DIC/New Derivation/Two Layer - Air/';
% path='/home/avanoers/Documents/DIC/New Derivation/WS - Air/';
path='/home/avanoers/Documents/DIC/New Derivation/Lyon/FreshWaterLyon/';
% path='/home/avanoers/Documents/DIC/New Derivation/Lyon/TwoLayerLyon/';
d=dir([path,'*.csv']);
store_figures=0;

% xstart = 1;%38;%60;%45;%
for i=1:length(d)
%     B=table2array(importfile1([path,d(i).name]));
    C = importfile([path,d(i).name]);
%     A(:,:,i) = B(:,xstart:end);
    A(:,:,i) = C(:,:);
    c=d(i).name;
    figure
    hold on
    surf(A(:,:,i),'EdgeColor','none','LineStyle','none','FaceLighting','phong');   
    view(90,0)
    title(c(1:end-4));
    colorbar
    ax=gca;
    ax.FontSize = 30;    
    hold off
    
if (store_figures == 1)
    path1=[path,strcat(c(1:end-4),'front')];
    print(gcf,'-dpng', path1)
pause(0.1)
end

    figure
    hold on
    surf(A(:,:,i),'EdgeColor','none','LineStyle','none','FaceLighting','phong');
    title(c(1:end-4));
    colorbar
    ax=gca;
    ax.FontSize = 30;    
    axis image
    hold off

if (store_figures == 1)    
    path1=[path,c(1:end-4)];    
   print(gcf,'-dpng', path1)
pause(0.1) 
end
end

%%
pause(0.1)
% 
% xstart = 5;
% xend = 441;
% ystart = 16;
% yend = 88;
% 
% xstart = 25;
% xend = 84;
% ystart = 1;
% yend = 94;
xstart = 38;
xend = 219;
ystart = 1;
yend = 373;

    DeltaX = A(ystart:yend,xstart:xend,4);
    DeltaY = -A(ystart:yend,xstart:xend,10);
    CorrCoeff = A(ystart:yend,xstart:xend,1);
    GridX = A(ystart:yend,xstart:xend,2);
    GridY = A(ystart:yend,xstart:xend,3);

    tri = delaunay(GridX, GridY);   
    %%%
    % Noise Removal - Median Filtering
    %%%
    % Filter
    DeltaX = medfilt2(DeltaX, 'symmetric');
    DeltaY = medfilt2(DeltaY, 'symmetric');
   
    DeltaX=-DeltaX(:);
    DeltaY=DeltaY(:);%-
    GridX=GridX(:);
    GridY=GridY(:);
    CorrCoeff=CorrCoeff(:);
    
       figure
        hold all
        trisurf(tri, GridX, GridY, DeltaX);
        shading interp
        colorbar EastOutside
        title('Displacements X [in pixels]')
        xlabel('Horizontal Position [pixels]')
        ylabel('Vertical Position [pixels]')    
        axis tight
        axis equal
        ax=gca;
        ax.FontSize = 20;
        hold off

if (store_figures == 1)
    path1=[path,'DX'];    
   print(gcf,'-dpng', path1)    
end    
        
        figure
        hold all
        trisurf(tri, GridX, GridY, DeltaY);
        shading interp
        colorbar EastOutside
        title('Displacements Y [in pixels]')
        xlabel('Horizontal Position [pixels]')
        ylabel('Vertical Position [pixels]')    
        axis tight
        axis equal
        ax=gca;
        ax.FontSize = 20;
        hold off
        
if (store_figures == 1)
    path1=[path,'DY'];    
   print(gcf,'-dpng', path1)    
end          

            % Distances
    L_camera = 0.5; %[m] distance camera - tank 
    L_glass = 5.8/1000; %[m] glass thickness
    L_test = 13.8/100-2*L_glass;%19.1/100; %[m] thickness test section
    L_screen = 0/100;%12.6/100; %[m] distance tank-screen
    Distance_to_mm = 3.1e-3;
    % Camera Property
    focal_length = 16; %[mm]
    Xcameracenter = mean(GridX);
    Ycameracenter = mean(GridY);

    % Refractive Indices
    n_air = 1;
    n_glass = 1.5;
    take_y_into_account = 0;
    Lengths = [L_camera, L_glass, L_test, L_screen];

    IC = [-1e3, 0.6];
    % boundaries Cx Cy L_tot_magn
    lb = [-inf, .1];
    ub = [+inf, 10];
    dx = @(Vx)nfunctiontominimize(Vx,DeltaX(:), DeltaY(:), GridX(:), GridY(:), Lengths, Distance_to_mm, focal_length, Xcameracenter, Ycameracenter, take_y_into_account);

    lsqOpts = optimoptions('lsqnonlin', 'MaxFunctionEvaluations', 2000,'FunctionTolerance', 1e-12, 'OptimalityTolerance',1e-12, 'MaxIterations', 1e3);
     [Vx,resnorm,residual,exitflag,output]  = lsqnonlin(dx,IC, lb, ub,lsqOpts);
     Vx
% True Values
distance_ccd_lens = 1/(1/(focal_length*1e-3)-1/(Vx(2)));

Xbarcenter = Vx(1);
% Ybarcenter = Vx(2);
    L_tot = Vx(2) *cosd(atand((Xcameracenter-Xbarcenter)*Distance_to_mm/1e3/distance_ccd_lens));
    L_camera = L_tot-2*L_glass-L_test-L_screen;

    L_total = L_camera+2*L_glass+L_test+L_screen;
% Distance from Center [in pixels]
DistanceX = GridX-Xcameracenter;
DistanceY = GridY-Ycameracenter;
% Distance = sqrt(DistanceX.^2+DistanceY.^2);

    phi_x_1 = atand(DistanceX*Distance_to_mm/1e3/(distance_ccd_lens));
    phi_y_1 = atand(DistanceY*Distance_to_mm/1e3/(distance_ccd_lens));
    
    phi_x_2 = atand((DistanceX+DeltaX)*Distance_to_mm/1e3/(distance_ccd_lens));
    phi_y_2 = atand((DistanceY+DeltaY)*Distance_to_mm/1e3/(distance_ccd_lens));

    alpha = atand((Xcameracenter-Xbarcenter)*Distance_to_mm/1e3/(distance_ccd_lens));
%     beta  = atand((Ycameracenter-Ybarcenter)*Distance_to_mm/1e3/(distance_ccd_lens));

    phi_x_1_0 = phi_x_1;
    phi_x_2_0 = phi_x_2;
    phi_x_1_1 = asind(n_air/n_glass*sind(alpha+phi_x_1_0))-alpha;
    phi_x_2_1 = asind(n_air/n_glass*sind(alpha+phi_x_2_0))-alpha;
    phi_x_1_2 = phi_x_1_0;

    Asol = (L_camera+L_screen)* (tand(phi_x_1_0)./(cosd(alpha)-tand(phi_x_1_0)*sind(alpha))-tand(phi_x_2_0)./(cosd(alpha)-tand(phi_x_2_0)*sind(alpha))) ...
        + 2*L_glass * (tand(phi_x_1_1)./(cosd(alpha)-tand(phi_x_1_1)*sind(alpha))-tand(phi_x_2_1)./(cosd(alpha)-tand(phi_x_2_1)*sind(alpha))) ...
        + L_test *tand(phi_x_1_2)./(cosd(alpha)-tand(phi_x_1_2)*sind(alpha));
    B = cosd(alpha)./(L_test./Asol+sind(alpha));

    n = n_air*sind(alpha+phi_x_2_0)./(sind(atand(B)+alpha));
   
    
    figure
    hold all
    trisurf(tri, GridX, GridY, n);
    colormap(cbrewer('seq', 'OrRd', 100));
    shading interp
    colorbar EastOutside
    title('n_x')
    xlabel('Horizontal Position [pixels]')
    ylabel('Vertical Position [pixels]')    
    axis tight
    axis equal
    ax=gca;
ax.FontSize = 20;
%     caxis([1.33 1.34])
    caxis([1.3 1.4])
	hold off

if (store_figures == 1)
    print('n_field','-dpng')
end    

    figure
    hold all
    trisurf(tri, GridX, GridY, n-1.333);
    colormap(cbrewer('seq', 'OrRd', 100));
    shading interp
    colorbar EastOutside
    title('Difference from 1.333')
    xlabel('Horizontal Position [pixels]')
    ylabel('Vertical Position [pixels]')    
    axis tight
    axis equal
    ax=gca;
ax.FontSize = 20;
%     caxis([1.33 1.34])
%     caxis([1.3 1.4])
	hold off
    
    n_x_2 =  n;
use = GridX>-inf;
% use = n_x_2 > 1.33;
[uy,iay,idy] = unique(GridY(use));

% N_X = reshape(n_x_2,[grid_setup.Ny,grid_setup.Nx]);
% SA=sum(isnan(N_X),2);

znanmean = accumarray(idy,n_x_2(use),[],@nanmean);
% znanmeany = accumarray(idy,n_y(use),[],@nanmean);
Groups = idy;
Weights = CorrCoeff;
Vals = n_x_2(use);
Vals (Weights < 0) = NaN;
Weights (Weights < 0) = 0;
Weights (isnan(Vals)) = 0;

% Weights = data(:,3); Vals = data(:,2); % pick your columns here
WeightedMeanFcn = @(ii) nansum(Vals(ii).*Weights(ii))/sum(Weights(ii));
wmeans = accumarray(Groups, 1:numel(Groups), [], WeightedMeanFcn);

figure; 
hold on
plot(wmeans,linspace(min(GridY), max(GridY), (max(GridY)-min(GridY))/5+1) );
    title('horizontally averaged n_x')
    xlabel('n_x ')
    ylabel('Vertical Position [pixels]')    
    axis tight
    ax=gca;
    ax.FontSize = 20;
	hold off

if (store_figures == 1)
    print('n_field_averaged','-dpng')
end    
% 
% 
% figure; 
% hold on
% Y= linspace(min(GridY), max(GridY), (max(GridY)-min(GridY))/5+1) ;
% DistanceYY = Y-Ycameracenter;
% phi_Y = atan(DistanceYY*Distance_to_mm/1e3/(distance_ccd_lens));
% plot(phi_Y(1:end/2),wmeans(1:end/2));
%     title('horizontally averaged n_x')
%     xlabel('n_x ')
%     ylabel('Vertical angle [deg]')    
%     axis tight
%     ax=gca;
%     ax.FontSize = 20;
% 	hold off
    
    n_x_2 =  n;
use = GridY>-inf;
% use = n_x_2 > 1.33;
[uy,iay,idy] = unique(GridX(use));

% N_X = reshape(n_x_2,[grid_setup.Ny,grid_setup.Nx]);
% SA=sum(isnan(N_X),2);

znanmean = accumarray(idy,n_x_2(use),[],@nanmean);
% znanmeany = accumarray(idy,n_y(use),[],@nanmean);
Groups = idy;
Weights = CorrCoeff;
Vals = n_x_2(use);
Vals (Weights < 0) = NaN;
Weights (Weights < 0) = 0;
Weights (isnan(Vals)) = 0;

% Weights = data(:,3); Vals = data(:,2); % pick your columns here
WeightedMeanFcn = @(ii) nansum(Vals(ii).*Weights(ii))/sum(Weights(ii));
wmeans = accumarray(Groups, 1:numel(Groups), [], WeightedMeanFcn);

figure; 
hold on
plot(linspace(min(GridX), max(GridX), (max(GridX)-min(GridX))/5+1), wmeans);
    title('vertically averaged n_x')
    ylabel('n_x ')
    xlabel('Horizontal Position [pixels]')    
    axis tight
    ax=gca;
    ax.FontSize = 20;
	hold off
% % 
% % % K_ave = (mean(n_x_2)-1);
% % K_ave = 0.333007236582526;
% % conv = Distance_to_mm/1000/distance_ccd_lens*Vx(3);
% % 
% % 
% % rho = (n_x_2-1)/K_ave;
% % 
% % 
% %     figure
% %     hold all
% %     trisurf(tri, GridX, GridY, rho);
% %     colormap(cbrewer('seq', 'OrRd', 100));
% %     shading interp
% %     colorbar EastOutside
% % %     title('\rho', 'interpreter', 'latex')
% %     title('\rho')
% %     xlabel('Horizontal Position [pixels]')
% %     ylabel('Vertical Position [pixels]')    
% %     axis tight
% %     axis equal
% %     ax=gca;
% % ax.FontSize = 20;
% % %     caxis([1 1.2])
% % %     caxis([1.3 1.4])
% % 	hold off
% % if (store_figures == 1)
% %     print('rho','-dpng')
% % end        
% % 
% % MiddleDomain = (max(GridY)+min(GridY))/2;
% % CoG=sum(GridY.*rho)/sum(rho)
% % meanrho = mean(rho)