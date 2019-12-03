clc; clear;close all;

path='C:\cygwin64\home\sande\T2\20-6-2018\SmallTank\Side\air-freshtest\';
d=dir([path,'Sfile.csv']);
C=txt2mat([path,d.name],1, 5);

[V,I] = min(C(:,5))

a = C(:,1);
b = C(:,2);
c = C(:,3);
d = C(:,4);
s = C(:,5);

[a(I) b(I) c(I) d(I) s(I)]
avalue = a(I);
bvalue = b(I);
cvalue = c(I);
dvalue = d(I);
svalue = s(I);
%%
% cvalue = 0.52;
% dvalue = -1.2;

[I] = find(C(:,3) == cvalue & C(:,4) == dvalue);

A = a(I);
B = b(I);
S = s(I);

        tri12 = delaunay(A, B);
        figure
        hold all
         trisurf(tri12, A, B, S);
        shading interp
        colorbar EastOutside
        title(['c = ', num2str(cvalue),' d = ', num2str(dvalue)])
        xlabel('a')
        ylabel('b') 
        axis tight
%         box on
        ax=gca;
        ax.FontSize = 20;
        hold off 
        
%%
% bvalue = 0.04;
% dvalue = -1.2;
[I] = find(C(:,2) == bvalue & C(:,4) == dvalue);

A = a(I);
B = c(I);
S = s(I);

        tri12 = delaunay(A, B);
        figure
        hold all
         trisurf(tri12, A, B, (S));
        shading interp
        colorbar EastOutside
        title(['b = ', num2str(bvalue),' d = ', num2str(dvalue)])
        xlabel('a')
        ylabel('c') 
        axis tight
%         box on
        ax=gca;
        ax.FontSize = 20;
        hold off            
        %%
%         bvalue = 0.04;
%         cvalue = 0.52;
[I] = find(C(:,2) == bvalue & C(:,3) == cvalue);

A = a(I);
B = d(I);
S = s(I);

        tri12 = delaunay(A, B);
        figure
        hold all
         trisurf(tri12, A, B, (S));
        shading interp
        colorbar EastOutside
        title(['b = ', num2str(bvalue),' c = ', num2str(cvalue)])
        xlabel('a')
        ylabel('d') 
        axis tight
%         box on
        ax=gca;
        ax.FontSize = 20;
        hold off      
        %%
%         avalue = 0.04;
%         dvalue = -1.2;
[I] = find(C(:,1) == avalue & C(:,4) == dvalue);

A = b(I);
B = c(I);
S = s(I);

        tri12 = delaunay(A, B);
        figure
        hold all
         trisurf(tri12, A, B, (S));
        shading interp
        colorbar EastOutside
        title(['a = ', num2str(avalue),' d = ', num2str(dvalue)])
        xlabel('b')
        ylabel('c') 
        axis tight
%         box on
        ax=gca;
        ax.FontSize = 20;
        hold off    
                  
                %%
%                 avalue = 0.04;
%                 cvalue = 0.44;
[I] = find(C(:,1) == avalue & C(:,3) == cvalue);

A = b(I);
B = d(I);
S = s(I);

        tri12 = delaunay(A, B);
        figure
        hold all
         trisurf(tri12, A, B, (S));
        shading interp
        colorbar EastOutside
        title(['a = ', num2str(avalue),' c = ', num2str(cvalue)])
        xlabel('b')
        ylabel('d') 
        axis tight
%         box on
        ax=gca;
        ax.FontSize = 20;
        hold off    
         
%%
% avalue = 0.04;
% bvalue = 0.04;
[I] = find(C(:,1) == avalue & C(:,2) == bvalue);

A = c(I);
B = d(I);
S = s(I);

        tri12 = delaunay(A, B);
        figure
        hold all
         trisurf(tri12, A, B, (S));
        shading interp
        colorbar EastOutside
        title(['a = ', num2str(avalue),' b = ', num2str(bvalue)])
        xlabel('c')
        ylabel('d') 
        axis tight
%         box on
        ax=gca;
        ax.FontSize = 20;
        hold off    
        
