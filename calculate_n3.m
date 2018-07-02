clc;clear;close all;
% close all
% path='/home/avanoers/Documents/DIC/matrixstorage/calibration/';
% path='/home/avanoers/Documents/DIC/New Derivation/Calibration - Air/';
% path = '/home/avanoers/Documents/DIC/New Derivation/Two Layer - Air/';
% path='/home/avanoers/Documents/DIC/New Derivation/WS - Air/';
% path='/home/avanoers/Documents/DIC/New Derivation/Lyon/FreshWaterLyon/';
% path='/home/avanoers/Documents/DIC/New Derivation/Lyon/TwoLayerLyon/';
% path='/home/avanoers/Documents/DIC/Lyon/Air - Calibration/';
% path='/home/avanoers/Documents/DIC/Lyon/Air - TwoLayer/';
% path='/home/avanoers/Documents/DIC/Lyon/Air - Stratified/';
% path='/home/avanoers/Documents/DIC/Lyon/Air - Calibration (copy)/';
% % 
path = '/home/avanoers/Documents/DIC/New Derivation/New/Accurate - 5 21 3 3 1 0/Air - Strat/';

% path = '/home/avanoers/Documents/DIC/New Derivation/New/Accurate - 5 21 3 3 1 0/Air - TwoL/';

% path = '/home/avanoers/Documents/DIC/New Derivation/New/Accurate - 5 21 3 3 1 0/Air - Cal/';
% 
% path='/home/avanoers/Documents/DIC/New Derivation/New/Accurate - 5 11 3 3 1 0/Air - Cal/';
% path='/home/avanoers/Documents/DIC/New Derivation/New/Accurate - 5 11 3 3 1 0/Air - TwoL/';

% path = '/home/avanoers/Documents/DIC/New Derivation/New/Air - Cal/';

d=dir([path,'*.csv']);
store_figures=0;

xstart = 95;%40;
for i=1:length(d)
    if (store_figures==1 || (i==1||i==2||i==3||i==4||i==10||i==14||i==5 || i==6))
%         tic
%     C = importfile([path,d(i).name]);
%         toc
%         tic
%     C = importfile([path,d(i).name]);
%     toc
    C=txt2mat([path,d(i).name],1);
        
    % Filter
     if (store_figures==1)
        C = medfilt2(C, 'symmetric'); 
     end
    A(:,:,i) = C(:,xstart:end);
    c=d(i).name;
    figure
    hold on
    surf(A(:,:,i),'EdgeColor','none','LineStyle','none','FaceLighting','phong');   
    view(90,0)
    name = [strcat(c(1:end-4),' front')];
    title(name);
    display([num2str(i), ' ', c(1:end-4)]);
    colorbar
%     caxis([mean(C(:))-std(C(:)) mean(C(:))+std(C(:))])
    ax=gca;
    ax.FontSize = 30;    
    hold off
    
if (store_figures == 1)
    path1=[path,strcat(c(1:end-4),'front')];
    print(gcf,'-dpng', path1)
% pause(0.1)
end

    figure
    hold on
    surf(A(:,:,i),'EdgeColor','none','LineStyle','none','FaceLighting','phong');
    title(c(1:end-4));
    colorbar
    
%     caxis([mean(C(:))-std(C(:)) mean(C(:))+std(C(:))])
    ax=gca;
    ax.FontSize = 30;    
    axis image
    hold off

if (store_figures == 1)    
    path1=[path,c(1:end-4)];    
   print(gcf,'-dpng', path1)
% pause(0.1) 
end
    end
end

if (store_figures == 1) 
    close all
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
% % yend = 94;
% xstart = 38;
% xend = 219;
% ystart = 1;
% yend = 373;

% ystart = 1;
% yend = 373;
% xstart = 1;
% xend = 178;%142
% xstart = 1;
% xend = 238;
% ystart = 1;
% yend = 378;
% xstart = 25;
% xend = 198;
% ystart = 15;
% yend = 367;
% xstart =1;
% xend=124;
% ystart=1;
% yend=189;
xstart=1;%95;
xend=size(A,2); %243
ystart=1;
yend=size(A,1);

    DeltaX = A(ystart:yend,xstart:xend,4);
    DeltaY = A(ystart:yend,xstart:xend,10);
    Ux = A(ystart:yend,xstart:xend,5);
    Uxx = A(ystart:yend,xstart:xend,6);
    Vy = A(ystart:yend,xstart:xend,14);
    CorrCoeff = A(ystart:yend,xstart:xend,1);
    GridX = A(ystart:yend,xstart:xend,2);
    GridY = A(ystart:yend,xstart:xend,3);

%  
        % Filter
    DeltaX = medfilt2(DeltaX, 'symmetric');
    DeltaY = medfilt2(DeltaY, 'symmetric');
    
I = [];%abs(Ux)>0.65 | Ux > 0   | CorrCoeff < 0.9;%[];%(abs(Uxx)>1e-3 | abs(Ux)>0.1) || CorrCoeff < 0.9); %[];%

    GridX(I)=[];
    GridY(I)=[];
    DeltaX(I)=[];
    DeltaY(I)=[];
    Vy(I)=[];
    Ux(I)=[];
    CorrCoeff(I)=[];
    Uxx(I)=[];

    
    
    %     
% % % cutoff = -2;%0.95;%-2;   
%     GridX(CorrCoeff<cutoff)=[];
%     GridY(CorrCoeff<cutoff)=[];
%     DeltaX(CorrCoeff<cutoff)=[];
%     DeltaY(CorrCoeff<cutoff)=[];
%     Vy(CorrCoeff<cutoff)=[];
%     Ux(CorrCoeff<cutoff)=[];
%     CorrCoeff(CorrCoeff<cutoff)=[];
%     
%     mUx = mean(Ux(:));
%     sUx = std(Ux(:));
 
%  GridX(Ux>0) = [];
%  GridY(Ux>0) = [];
%  DeltaX(Ux>0) = [];
%  DeltaY(Ux>0) = [];
%  Vy(Ux>0) = [];
%  CorrCoeff(Ux>0) = [];
%  Ux(Ux>0) = [];
%  
%  GridX(Ux>mUx+sUx)=[];
%  GridY(Ux>mUx+sUx)=[];
%  DeltaX(Ux>mUx+sUx)=[];
%  DeltaY(Ux>mUx+sUx)=[];
%  Vy(Ux>mUx+sUx)=[];
%  CorrCoeff(Ux>mUx+sUx)=[];
%  Ux(Ux>mUx+sUx)=[];
%  
%  GridX(Ux<mUx-sUx)=[];
%  GridY(Ux<mUx-sUx)=[];
%  DeltaX(Ux<mUx-sUx)=[];
%  DeltaY(Ux<mUx-sUx)=[];
%  Vy(Ux<mUx-sUx)=[];
%  CorrCoeff(Ux<mUx-sUx)=[];
%  Ux(Ux<mUx-sUx)=[];
 
    tri = delaunay(GridX, GridY);  
    %%%
    % Noise Removal - Median Filtering
%     %%%
%     % Filter
%     DeltaX = medfilt2(DeltaX, 'symmetric');
%     DeltaY = medfilt2(DeltaY, 'symmetric');
%     

    DeltaX=-DeltaX(:);
    DeltaY=DeltaY(:);%-
    GridX=GridX(:);
    GridY=GridY(:);
    CorrCoeff=CorrCoeff(:);
    Vy = Vy(:);
    Ux = Ux(:);
    Uxx = Uxx(:);
      
        figure
        hold all
        plot3(GridX, GridY, CorrCoeff, 'k+');
        title('Correlation Coefficient [-1 1]')
        xlabel('Horizontal Position [pixels]')
        ylabel('Vertical Position [pixels]') 
        axis tight
        box on
        ax=gca;
        ax.FontSize = 20;
        hold off         
        
        
        figure
        hold all
         trisurf(tri,GridX, GridY, CorrCoeff);
        shading interp
        colorbar EastOutside
        title('Correlation Coefficient [-1 1]')
        xlabel('Horizontal Position [pixels]')
        ylabel('Vertical Position [pixels]') 
        axis tight
%         box on
        ax=gca;
        ax.FontSize = 20;
        hold off         
        
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

  
        figure
        hold all
        trisurf(tri, GridX, GridY, Vy);
        shading interp
        colorbar EastOutside
        title('V_y [-]')
        xlabel('Horizontal Position [pixels]')
        ylabel('Vertical Position [pixels]')    
        axis tight
        axis equal
        ax=gca;
        ax.FontSize = 20;
        hold off        
        
        figure
        hold all
        trisurf(tri, GridX, GridY, Ux);
        shading interp
        colorbar EastOutside
        title('U_x [-]')
        xlabel('Horizontal Position [pixels]')
        ylabel('Vertical Position [pixels]')    
        axis tight
        axis equal
        ax=gca;
        ax.FontSize = 20;
        hold off
        
        figure
        hold all
        trisurf(tri, GridX, GridY, Uxx);
        shading interp
        colorbar EastOutside
        title('U_{xx} [1/pixel]')
        xlabel('Horizontal Position [pixels]')
        ylabel('Vertical Position [pixels]')    
        axis tight
        axis equal
        ax=gca;
        ax.FontSize = 20;
        hold off        
        
        
p=[GridX,GridY];


figure;
hold all
 tricontour(p,tri, DeltaY,15); 
        shading interp
        colorbar EastOutside
        title('Displacements Y [in pixels]')
        xlabel('Horizontal Position [pixels]')
        ylabel('Vertical Position [pixels]')    
        axis tight
        axis equal
        box on
        ax=gca;
        ax.FontSize = 20;
        hold off
figure
hold all
tricontour(p,tri, DeltaY,15); 
        shading interp
        colorbar EastOutside
        title('Displacements Y [in pixels]')
        xlabel('Horizontal Position [pixels]')
        ylabel('Vertical Position [pixels]')    
        axis tight
        axis equal
        box on
        ax=gca;
        ax.FontSize = 20;
        hold off        
        
% 
% figure;
% hold all
% tricontour(p,tri, Vy,15); 
%         shading interp
%         colorbar EastOutside
%         title('V_y [-]')
%         xlabel('Horizontal Position [pixels]')
%         ylabel('Vertical Position [pixels]')    
%         axis tight
%         axis equal
%         box on
%         ax=gca;
%         ax.FontSize = 20;
%         hold off        
% %   
% figure;
% hold all
% [c,h] = tricontour(p,tri, DeltaX,15); 
%         shading interp
%         colorbar EastOutside
%         title('Displacements X [pixels]')
%         xlabel('Horizontal Position [pixels]')
%         ylabel('Vertical Position [pixels]')    
%         axis tight
%         axis equal
%         box on
%         ax=gca;
%         ax.FontSize = 20;
%         hold off     
% if (store_figures == 1)
%     path1=[path,'DX_contour'];    
%    print(gcf,'-dpng', path1)    
% end          

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
    n_air = 1.0003;
    n_glass = 1.5;
    take_y_into_account = 0;
    Lengths = [L_camera, L_glass, L_test, L_screen];
% 
% %         IC = [-1e3, 0.6];
%         IC = [-932, 0.565];
%         % boundaries Cx Cy L_tot_magn
%         lb = [-inf, .1];
%         ub = [+inf, 10];
pause(0.1)
        %%
                IC = [-1011 0.6186, 1232];
%                 IC = [-1031, 0.621, 1007];
        % boundaries Cx Cy L_tot_magn
        lb = [-inf, .1, -inf];
        ub = [+inf, 10, +inf];
        dx = @(Vx)nfunctiontominimize_Snell3(Vx,DeltaX(:), DeltaY(:), GridX(:), GridY(:), Lengths, Distance_to_mm, focal_length, Xcameracenter, Ycameracenter, take_y_into_account);

        lsqOpts = optimoptions('lsqnonlin' , 'StepTolerance', 1e-12,'MaxFunctionEvaluations', 2000,'FunctionTolerance', 1e-12, 'OptimalityTolerance',1e-12, 'MaxIterations', 1e3);
%          [Vx,resnorm,residual,exitflag,output]  = lsqnonlin(dx,IC, lb, ub,lsqOpts);
%          Vx
%          Vx = 1e3*[-0.990555438769635   0.000614569872955   1.257769046731914];
%         Vx = 1e3*[  -1.016924129058061   0.000618343679763   1.223315233588966];
% Vx = 1e3*[-1.010969792776444   0.000619657574061   1.223353830030936];
% 1.335
% Vx = 1e3*[-1.012220715916666   0.000619310034272   1.246942638222305]; %1.3345

%     Vx = 1e3*[-1.063886184076514   0.000629441701773   1.2228633r53282208];

Vx = 1e3*[-0.984896945941551   0.000614170743748   1.234471138949774];
        % True Values
        distance_ccd_lens = 1/(1/(focal_length*1e-3)-1/(Vx(2)));

        Xbarcenter = Vx(1);
        Ybarcenter = Vx(3); %Ycameracenter;%
            L_tot = Vx(2)*cosd(atand((Xcameracenter-Xbarcenter)*Distance_to_mm/1e3/distance_ccd_lens));


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
    beta  = atand((Ycameracenter-Ybarcenter)*Distance_to_mm/1e3/(distance_ccd_lens));

    phi_x_1_0 = phi_x_1;
    phi_x_2_0 = phi_x_2;
    phi_y_1_0 = phi_y_1;
    phi_y_2_0 = phi_y_2;    
%     phi_x_1_1 = asind(n_air/n_glass*sind(alpha+phi_x_1_0))-alpha;
%     phi_x_2_1 = asind(n_air/n_glass*sind(alpha+phi_x_2_0))-alpha;
    phi_x_1_2 = phi_x_1_0;
    phi_y_1_2 = phi_y_1_0;
%     phi_x_2_2 = asind(n_air/1.333*sind(alpha+phi_x_2_0))-alpha;
  
%         [T] = create_TransformationMatrices(alpha);
        
        [T] = create_TransformationMatrices2(alpha,beta);
    [phi_x_1_1_new, phi_y_1_1_new, phi_x_2_1_new, phi_y_2_1_new] = SnellsLaw(phi_x_1_0, phi_y_1_0, phi_x_2_0, phi_y_2_0, alpha, n_air, n_glass, T);

 
%     [phi_x_1_2_new, phi_y_1_2_new, phi_x_2_2_new, phi_y_2_2_new1] = SnellsLaw(phi_x_1_1_new, phi_y_1_1_new, phi_x_2_1_new, phi_y_2_1_new, alpha, n_glass, 1.333);
    
    phi_x_1_1 = phi_x_1_1_new;
    phi_y_1_1 = phi_y_1_1_new;
    phi_x_2_1 = phi_x_2_1_new;
    phi_y_2_1 = phi_y_2_1_new;

    

    A = (L_camera+L_screen)* (tand(phi_x_1_0)./(cosd(alpha)*cosd(beta)-cosd(beta)*tand(phi_x_1_0)*sind(alpha)-cosd(alpha)*tand(phi_y_1_0)*sind(beta)) ...
            -tand(phi_x_2_0)./(cosd(alpha)*cosd(beta)-cosd(beta)*tand(phi_x_2_0)*sind(alpha)-cosd(alpha)*tand(phi_y_2_0)*sind(beta))) ...
        + 2*L_glass * (tand(phi_x_1_1)./(cosd(alpha)*cosd(beta)-cosd(beta)*tand(phi_x_1_1)*sind(alpha)-cosd(alpha)*tand(phi_y_1_1)*sind(beta)) ...
            -tand(phi_x_2_1)./(cosd(alpha)*cosd(beta)-cosd(beta)*tand(phi_x_2_1)*sind(alpha)-cosd(alpha)*tand(phi_y_2_1)*sind(beta))) ...
        + L_test *tand(phi_x_1_2)./(cosd(alpha)*cosd(beta)-cosd(beta)*tand(phi_x_1_2)*sind(alpha)-cosd(alpha)*tand(phi_y_1_2)*sind(beta));
    %%
    %     profile on
    tic;
%         [T] = create_TransformationMatrices2(alpha,beta);
        
       lsqOpts = optimoptions('lsqnonlin' , 'StepTolerance', 1e-12,'MaxFunctionEvaluations', 2000,'FunctionTolerance', 1e-12, 'OptimalityTolerance',1e-12, 'MaxIterations', 1e3,'display','off');
    flag_halfway_done = 0;
    stepsize = 100;
%     n2 = zeros(floor(length(B)/stepsize),1);
    n2 = zeros(length(DeltaX),1);
        ICC = [1.334];
        lbb = [1.3];
        ubb = [1.4];
    for i=1:stepsize:length(A)
        [dn] = @(n)SnellsLaw_mini3(phi_x_1_1(i), phi_y_1_1(i), phi_x_2_1(i), phi_y_2_1(i), alpha, beta, n_glass, n, A(i), T, L_test);  
        [Vn]  = lsqnonlin(dn,ICC, lbb, ubb,lsqOpts);
         if(flag_halfway_done == 0 && i/length(A)-0.5 > 0) 
            display('Halfway Done');
            flag_halfway_done = 1;
         end
         n2(i) = Vn;
    end
    display('Done');
    toc;

GridX = GridX(1:stepsize:end);
GridY = GridY(1:stepsize:end);
% tri1=delaunay(GridX(1:1000:end), GridY(1:1000:end));
tri = delaunay(GridX,GridY);

    
    %%
%     profile off
%     profile viewer
    figure;plot(n2(n2>0),'k+');
In=isoutlier(n2, 'ThresholdFactor', 8); %4
t=1:size(n2);
figure;hold all;plot(t,n2); plot(t(In),n2(In),'k+');hold off;

    n = n2';
    n(n==0)=[];
    
if (store_figures == 1)
    print('n_field','-dpng')
end

if (stepsize == 1)
        figure
        hold all
        plot3(GridX(In), GridY(In), CorrCoeff(In), 'k+');
        title('Correlation Coefficient [-1 1]')
        xlabel('Horizontal Position [pixels]')
        ylabel('Vertical Position [pixels]') 
        axis tight
        box on
        ax=gca;
        ax.FontSize = 20;
        axis([min(GridX) max(GridX) min(GridY) max(GridY)])
        hold off         
        
        
        figure
        hold all
        plot(t,CorrCoeff);
        plot(t(In),CorrCoeff(In), 'k+');
        title('C');
        hold off
        figure
        hold all
        plot(t,Ux);
        plot(t(In),Ux(In), 'k+');
        title('U_x');
        hold off
        figure
        hold all
        plot(t,Uxx);
        plot(t(In),Uxx(In), 'k+');
        title('U_{xx}');
        hold off    
        figure
        hold all
        plot(t,Vy);
        plot(t(In),Vy(In), 'k+');
        title('V_{y}');
        hold off       
end
%%
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
%     caxis([1.33 1.34])
%     caxis([mean(n)-3*std(n) mean(n)+3*std(n)])
	hold off

if (store_figures == 1)
    print('n','-dpng')
end       



     figure
    hold all
    trisurf(tri, GridX, GridY, n-1.333);
%     colormap(cbrewer('seq', 'OrRd', 100));
    shading interp
    colorbar EastOutside
    title('diff from 1.333')
    xlabel('Horizontal Position [pixels]')
    ylabel('Vertical Position [pixels]')    
    axis tight
    axis equal
    ax=gca;
ax.FontSize = 20;
%     caxis([1.33 1.34])
%     caxis([mean(n)-3*std(n) mean(n)+3*std(n)])
	hold off

%        
%     n2 = medfilt2(n,'symmetric');
%     figure
%     hold all
%     trisurf(tri, GridX, GridY, n2);
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
% %     caxis([1.33 1.34])
% %     caxis([mean(n)-3*std(n) mean(n)+3*std(n)])
% 	hold off
% %      
%     
%     
%     figure
%     hold all
%     trisurf(tri, GridX, GridY, n-1.333);
%     colormap(cbrewer('seq', 'OrRd', 100));
%     shading interp
%     colorbar EastOutside
%     title('Difference from 1.333')
%     xlabel('Horizontal Position [pixels]')
%     ylabel('Vertical Position [pixels]')    
%     axis tight
%     axis equal
%     ax=gca;
% ax.FontSize = 20;
% %     caxis([1.33 1.34])
% %     caxis([1.3 1.4])
% 	hold off
    
    n_x_2 =  n';
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
cutoff = 0;
Vals (Weights < cutoff) = NaN;
Weights (Weights < cutoff) = 0;
Weights (isnan(Vals)) = 0;

% Weights = data(:,3); Vals = data(:,2); % pick your columns here
WeightedMeanFcn = @(ii) nansum(Vals(ii).*Weights(ii))/sum(Weights(ii));
wmeans = accumarray(Groups, 1:numel(Groups), [], WeightedMeanFcn);

figure; 
hold all
plot(wmeans,linspace(min(GridY), max(GridY), length(wmeans) )); 
    title('horizontally averaged n_x')
    xlabel('n_x ')
    ylabel('Vertical Position [pixels]')    
    axis tight
    ax=gca;
    ax.FontSize = 20;
	hold off

if (store_figures == 1)
    print('n_hor_averaged','-dpng')
end    


    n_x_2 =  n';
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
cutoff = 0;
Vals (Weights < cutoff) = NaN;
Weights (Weights < cutoff) = 0;
Weights (isnan(Vals)) = 0;

% Weights = data(:,3); Vals = data(:,2); % pick your columns here
WeightedMeanFcn = @(ii) nansum(Vals(ii).*Weights(ii))/sum(Weights(ii));
wmeans_y = accumarray(Groups, 1:numel(Groups), [], WeightedMeanFcn);

figure; 
hold on
plot(linspace(min(GridX), max(GridX), length(wmeans_y)),wmeans_y); %(max(GridY)-min(GridY))/3+1)
    title('vertically averaged n_x')
    ylabel('n_x ')
    xlabel('Horizontal Position [pixels]')    
    axis tight
    ax=gca;
    ax.FontSize = 20;
	hold off

if (store_figures == 1)
    print('n_ver_averaged','-dpng')
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
%     
% %     n_x_2 =  n;
% use = GridY>-inf;
% % use = n_x_2 > 1.33;
% [uy,iay,idy] = unique(GridX(use));
% 
% % N_X = reshape(n_x_2,[grid_setup.Ny,grid_setup.Nx]);
% % SA=sum(isnan(N_X),2);
% 
% znanmean = accumarray(idy,n_x_2(use),[],@nanmean);
% % znanmeany = accumarray(idy,n_y(use),[],@nanmean);
% Groups = idy;
% Weights = CorrCoeff;
% Vals = n_x_2(use);
% cutoff = 0.95;
% Vals (Weights < cutoff) = NaN;
% Weights (Weights < cutoff) = 0;
% Weights (isnan(Vals)) = 0;
% 
% % Weights = data(:,3); Vals = data(:,2); % pick your columns here
% WeightedMeanFcn = @(ii) nansum(Vals(ii).*Weights(ii))/sum(Weights(ii));
% wmeans = accumarray(Groups, 1:numel(Groups), [], WeightedMeanFcn);
% 
% figure; 
% hold on
% plot(linspace(min(GridX), max(GridX), (max(GridX)-min(GridX))/3+1), wmeans);
%     title('vertically averaged n_x')
%     ylabel('n_x ')
%     xlabel('Horizontal Position [pixels]')    
%     axis tight
%     ax=gca;
%     ax.FontSize = 20;
% 	hold off
%     
% 
% if (store_figures == 1)
%     print('n_ver_averaged','-dpng')
% end        
% 
% K_ave = (mean(n)-1)/0.9998;
% K_ave =   0.333066367557091;
% % K_ave = 0.333;
% % % 
% % % 
% rho = (n-1)/K_ave;
% % % 
% % % 
%     figure
%     hold all
%     trisurf(tri, GridX, GridY, rho);
%     colormap(cbrewer('seq', 'OrRd', 100));
%     shading interp
%     colorbar EastOutside
% %     title('\rho', 'interpreter', 'latex')
%     title('\rho')
%     xlabel('Horizontal Position [pixels]')
%     ylabel('Vertical Position [pixels]')    
%     axis tight
%     axis equal
%     ax=gca;
% ax.FontSize = 20;
% %     caxis([1 1.2])
% %     caxis([1.3 1.4])
% 	hold off
% % % if (store_figures == 1)
% % %     print('rho','-dpng')
% % % end        
% % % 
% % % MiddleDomain = (max(GridY)+min(GridY))/2;
% % % CoG=sum(GridY.*rho)/sum(rho)
% % % meanrho = mean(rho)
% % 
% % K_ave = 0.327;
% % K_ave=0.333;
% rhomeans = (wmeans-1)/K_ave;
% % n=n';
% T = 16;
% n_temp_corrected = n +1.3531e-4*T+5.1e-8*T^2;
% K_ave = (mean(n)-1)/0.9998;
% K_ave =   0.333066367557091;
% K_LL = (mean(n)^2-1)/(mean(n)^2+2)/.9998;
% K_LL = 0.205736490652667;
% 
% K_LL = 0.205736490652667;%(mean(n)^2-1)/(mean(n)^2+2)/.9998;
% K_LL_T = 0.206957019085785;% (mean(n_temp_corrected)^2-1)/(mean(n_temp_corrected)^2+2)/.9998;
% K_GL = 0.333066367557091;% (mean(n)-1)/0.9998
% K_GL_T = 0.335244819247428; % (mean(n_temp_corrected)-1)/0.9998
% 
% rho_GL = (wmeans-1)/K_GL;
% rho_LL = (wmeans.^2-1)./(wmeans.^2+2)/K_LL;
% % 
% % figure; 
% % hold on
% LL = L_tot-L_screen-L_glass-L_test/2;
% LM = LL/cosd(alpha);
% conv = Distance_to_mm/1000/distance_ccd_lens*LM;
% Y = linspace(max(GridY), min(GridY), length(wmeans)); %(max(GridY)-min(GridY))/3+1
% Y=Y*conv;
% % plot(rhomeans-1, Y*100);
% % box on
% % grid on
% % axis([0 0.07 0 25])
% % xticks([0 0.01 0.02 0.03 0.04 0.05 0.06 0.07])
% %     title('horizontally averaged rho')
% %     xlabel('\Delta\rho ')
% %     ylabel('Vertical Position [cm]')    
% % %     axis tight
% %     ax=gca;
% %     ax.FontSize = 20;
% % 	hold off
% % %  
% % n_LL = (wmeans.^2-1)./(wmeans.^2+2);
% % figure; 
% hold on
LL = L_tot-L_screen-L_glass-L_test/2;
LM = LL/cosd(alpha)/cosd(beta);
conv = Distance_to_mm/1000/distance_ccd_lens*LM;
% % Y = linspace(max(GridY), min(GridY), length(wmeans)); %(max(GridY)-min(GridY))/3+1
% % Y=Y*conv;
% % plot(n_LL, Y*100);
% % box on
% % grid on
% % % axis([0 0.07 0 25])
% % % xticks([0 0.01 0.02 0.03 0.04 0.05 0.06 0.07])
% %     title('n_LL')
% %     xlabel('\Delta\rho ')
% %     ylabel('Vertical Position [cm]')    
% % %     axis tight
% %     ax=gca;
% %     ax.FontSize = 20;
% % 	hold off
%     
% %     
% %     figure; 
% % hold on
% % LL = L_tot-L_screen-L_glass-L_test/2;
% % LM = LL/cosd(alpha);
% % conv = Distance_to_mm/1000/distance_ccd_lens*LM;
% % Y = linspace(max(GridY), min(GridY), length(wmeans)); %(max(GridY)-min(GridY))/3+1
% % Y=Y*conv;
% % plot(rhomeans-1, Y*100);
% % plot(rho_profiel-1, Y*100, 'b');
% % box on
% % grid on
% % axis([0 0.07 0 25])
% % xticks([0 0.01 0.02 0.03 0.04 0.05 0.06 0.07])
% %     title('horizontally averaged rho')
% %     xlabel('\Delta\rho ')
% %     ylabel('Vertical Position [cm]')    
% % %     axis tight
% %     ax=gca;
% %     ax.FontSize = 20;
% % 	hold off
% % %     
% % var(n)
% % 
% % T = 17;
% % wmeans_temp = wmeans+1.351e-4*T+5.1e-8*T^2;  %-1.3373
% %     figure; 
% % hold on
% % plot(wmeans_temp, Y*100);
% % box on
% % grid on
% %     title('horizontally averaged rho')
% %     xlabel('n_temp ')
% %     ylabel('Vertical Position [cm]')    
% % %     axis tight
% %     ax=gca;
% %     ax.FontSize = 20;
% % 	hold off
% 
%     
% % T = 16;
% wmeans_temp = wmeans+1.351e-4*T+5.1e-8*T^2;  %-1.3373
% % K_LL = (1.333^2-1)/(1.333^2+2);
% rho_LL_T = (wmeans_temp.^2-1)./(wmeans_temp.^2+2)/K_LL_T;
% 
% %     K_ave=0.333;
% rho_GL_T = (wmeans_temp-1)/K_GL_T;
% % 
% % rho_LL = (wmeans.^2-1)./(wmeans.^2+2)/K_LL;
% % rho_means_LL = (wmeans_temp.^2-1)./(wmeans_temp.^2+2)/K_LL;
%     figure; 
% hold on
% plot(rhomeans-1, Y*100, 'b');
% plot(rho_GL_T-1, Y*100, 'r');
% plot(rho_LL-1, Y*100, 'g');
% plot(rho_LL_T-1, Y*100, 'y');
% legend('GL', 'GL Temp', 'LL', 'LL Temp');
% box on
% grid on
%     title('horizontally averaged rho')
%     xlabel('\Delta\rho ')
%     ylabel('Vertical Position [cm]')    
% %     axis tight
%     ax=gca;
%     ax.FontSize = 20;
%     hold off
%     
%     
% p=[GridX,GridY];
% 
% 
% figure;
% hold all
% [c,h] = tricontour(p,tri, n',25); 
%         shading interp
%         colorbar EastOutside
%         title('n')
%         xlabel('Horizontal Position [pixels]')
%         ylabel('Vertical Position [pixels]')  
%         axis tight
%         axis equal
%         box on
%         ax=gca;
%         ax.FontSize = 20;
%         hold off
%         
% 
%         
% % Tan

% 
% % c1 = (1.7682e-3/5.8e-6/2)+0.5*sqrt((1.7682e-3/5.8e-6)^2-4*lhs/5.8e-6);
% 
% c2 = (1.7682e-3/5.8e-6/2)-0.5*sqrt((1.7682e-3/5.8e-6)^2-4*lhs/5.8e-6);
% % figure;plot(c1)
% figure;plot(c2)
% 
% % figure;plot(c1*10+1000)
% figure;plot(c2*10+1000)

%%
rhomeans = interp1(T(:,2), T(:,1),wmeans);
figure; 
hold all
plot(rhomeans,linspace(min(GridY), max(GridY), length(wmeans) )); 
    title('horizontally averaged rho')
    xlabel('rho')
    ylabel('Vertical Position [pixels]')    
    axis tight
    ax=gca;
    ax.FontSize = 20;
	hold off
    
Ypix =     linspace(min(GridY), max(GridY), length(wmeans)) ;
Yimag = conv*Ypix;
Yimag = flipud(Yimag);

Temp = 16;
wmeansT = wmeans+1.3531e-4*Temp+5.1e-8*Temp^2;
Table2;
rhomeans = interp1(T(:,2), T(:,1), wmeans);%, 'spline', 'extrap');

rhomeansT = interp1(T(:,2), T(:,1),wmeansT);

figure;
plot(flipud(rhomeans),Yimag-0.22)


figure;
plot(Yimag-0.22, flipud(rhomeans))


load 5_step.mat;
figure; 
hold all
% plot(y.Drho,-y.z1)
plot((flipud(rhomeans)-1)*1000,Yimag-0.22); 
% plot((flipud(rhomeansT)-1)*1000,Yimag-0.22); 
    title('horizontally averaged rho')
    xlabel('\Delta\rho')
    ylabel('Vertical Position [cm]')    
    axis tight
    ax=gca;
%     xlim([0 70])
    ax.FontSize = 20;
	hold off    
var(n)