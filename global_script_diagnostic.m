
clearvars -except metaplanktonC metaplankton
close all

%% This script call the function : function_diagnostic_tool_map
% Outputs of this function

% abundance_community : abundance of each functional group at each station 
% biovolume_community : biovolume of each functional group at each station
% WW_community : wet mass of each functional group at each station
% CW_community : carbone mass of each functional group at each station
% respiration rate of the entire community for each station
% clearance_group : clearance rate of higher trophic levels at each station
% respiration_group : respiration rate of higher trophic levels at each station
% growth_group_station : growth rate of higher trophic levels at each station
% sum_nbr_prey_consumed_carbon : number of prey items consumed per day by each functional group and station (in terms carbon)
% percentage_prey_consumed : percentage of the community consumed per day


%% Inputs for the function : 

%Data set

% load metaplanktonC.mat
% metaplanktonC = metaplankton;
% load metaplankton.mat

% carbon percent: 
copepod = 9.48; chaetognatha = 3.56; gelatinous_carnivorous = 0.18; hydrozoa = 0.384141406; scyphozoa = 0.491428569; ctenophora = 0.205286;
gelatinous_filter_feeders = 0.56; appendicularia = 0.07; thaliacea = 1.081628667; large_crustracean = 9.64; pteropod = 5.52;

carbon_percent = [appendicularia;gelatinous_carnivorous;scyphozoa;thaliacea;ctenophora;copepod;chaetognatha;hydrozoa;gelatinous_filter_feeders;pteropod;large_crustracean];

%% function

[diagnostic_global] = function_diagnostic_global(carbon_percent,metaplankton,metaplanktonC);
