close all

%% Set the variables and environment

i=1;


% choose the good functional group (:,j)
capture_percent = mesozooplankton.capture_efficacity_percent(:,i)';

% Set the axes
ESD = ESDmatrix(:,1);
ESD = ESD';
x = [min(ESD)  max(ESD)];
y = [0.0001  max(sum(community.CW_community_total,2))];

%% figure

figure

opacity = 0.4;
im = image(x,y,capture_percent, ...
   "CDataMapping","scaled",...           
   "AlphaData", opacity);                  % Set Transparency as through the alphaData property 
set(gca,'Ydir','normal')
xlim([min(ESD) , max(ESD)])
ylim([min(sum(community.CW_community_total,2)) max(sum(community.CW_community_total,2))])
ax = gca;
ax.FontSize = 12; 

hold on

cb = colorbar;
cb.Title.String = "capture efficiency";
cb.Title.FontSize = 14;
cb.FontSize = 12; 

plot(ESD,sum(community.CW_community_total,2) , LineWidth=2, Color=[0,0,0])
plot(ESD,mesozooplankton.CW_size_class(:,i), LineWidth=2,  Color=[1.00,0.41,0.16])

xlabel('Equivalent spherical diameter (ESD)','FontSize', 14)
ylabel('Total carbon mass (gC.m^{-2})','FontSize',14)  

%191
xlim([6*10.^-1 , 2*10.^4])
ylim([10.^-6 , 2*10.^0])