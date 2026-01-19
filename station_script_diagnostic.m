%clear 
close all 
clearvars -except metaplanktonC metaplankton Mergedbase base            

%% This script call the function : function_diagnostic_tool_map
% Outputs of this function

% abundance_community : abundance of the community
% biovolume_community : biovolue of the community
% abundance_j_size_class : abundance of each functional group 
% biovolume_j_size_class : biovolume of each functional group 
% WW_j_size_class : wet mass of each functional group 
% CW_j_size_class : carbone mass of each functional group 
% clearance_group_station : clearance rate of higher trophic levels 
% respiration_group_station : respiration rate of higher trophic levels 
% growth_group_station : growth rate of higher trophic levels 
% capture_efficacity_percent : capture efficacity of each functional group for each size class
% nbr_prey_consumed_size_abundance : number of prey items consumed per day by each functional group (in terms abundance)
% nbr_prey_consumed_size_abundance : number of prey items consumed per day by each functional group (in terms biovolume)
% nbr_prey_consumed_size_carbon : number of prey items consumed per day by each functional group (in terms carbon)
% after_grazing_abundance : abundance after grazing for each functional group
% after_grazing_biovolume : biovolume after grazing for each functional group

%% Inputs for the function
station_number = 168;

load metaplanktonC.mat
metaplanktonC = metaplankton;
load metaplankton.mat

X=metaplanktonC(1).X;
X1=metaplanktonC(1).X1;
X=exp(X);
ESDmatrix=metaplanktonC(station_number).merged.ESDmatrix;
ESD=ESDmatrix(:,1);

ordre=metaplanktonC(station_number).merged.ordre;
namepft2=metaplanktonC(station_number).WP2.namepft;
namepft2 = namepft2';

table_bv =metaplankton(station_number).merged4devicesmin.datapft; %NBSS of all functional groups at station 178
table_c = metaplanktonC(station_number).merged4devicesmin.datapft; %NBSS of all functional groups at station 178

%table_bv = metaplankton(station_number).merged2netsUVPH20.datapft; %NBSS of all functional groups at station 178%
%table_c = metaplanktonC(station_number).merged2netsUVPH20.datapft; %NBSS of all functional groups at station 178

table_bv = array2table(table_bv);
table_c = array2table(table_c);

table_bv = renamevars(table_bv, {'table_bv1','table_bv2','table_bv3','table_bv4','table_bv5','table_bv6','table_bv7','table_bv8','table_bv9','table_bv10','table_bv11','table_bv12','table_bv13','table_bv14','table_bv15','table_bv16','table_bv17','table_bv18','table_bv19','table_bv20','table_bv21','table_bv22','table_bv23','table_bv24','table_bv25','table_bv26','table_bv27','table_bv28'},namepft2); 
table_c = renamevars(table_c, {'table_c1','table_c2','table_c3','table_c4','table_c5','table_c6','table_c7','table_c8','table_c9','table_c10','table_c11','table_c12','table_c13','table_c14','table_c15','table_c16','table_c17','table_c18','table_c19','table_c20','table_c21','table_c22','table_c23','table_c24','table_c25','table_c26','table_c27','table_c28'},namepft2); 

X1 = array2table(X1);
X = array2table(X);
ESD = array2table(ESD);

table_bv = [X1,X,ESD, table_bv];
table_c  = [X1,X,ESD, table_c];

sizespectre_bv = renamevars(table_bv,  {'X1','X','ESD'}, {'largeurBV','BvMedian','ESD_micron_AllUnitsInMm3_Mm3_m3'});
sizespectre_c = renamevars(table_c,  {'X1','X','ESD'}, {'largeurBV','BvMedian','ESD_micron_AllUnitsInMm3_Mm3_m3'});

% carbon percent 
copepod = 9.48; chaetognatha = 3.56; gelatinous_carnivorous = 0.18; hydrozoa = 0.384141406; scyphozoa = 0.491428569; ctenophora = 0.205286; 
gelatinous_filter_feeders = 0.56; appendicularia = 0.07; thaliacea = 1.081628667; large_crustacean = 9.64; pteropod = 5.52;

carbon_percent = [appendicularia;gelatinous_carnivorous;scyphozoa;thaliacea;ctenophora;copepod;chaetognatha;hydrozoa;gelatinous_filter_feeders;pteropod;large_crustacean];

% temperature of the station
temperature = metaplanktonC(station_number).hydro.Temp_10m; %from the hydro table en °C at 10m

%% function

[community,mesozooplankton] = function_diagnostic_tool_station(sizespectre_bv, sizespectre_c, carbon_percent, temperature);
