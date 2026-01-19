
function [diagnostic_global] = global_function_diagnostic(carbon_percent,metaplankton,metaplanktonC)

%% 1. Initialization

%% 1.1. Matrix preparation

group_names = {'Appendicularia' 'gelatinous_carnivorous' 'Scyphozoa' 'Thaliacea' 'Ctenophora' 'all_copepoda' 'chaetognatha' 'Hydrozoa' 'gelatinous_filter_feeders' 'pteropoda' 'large_crustracean'};

num_steps = length(group_names); 
descriptive_global.abundance_community = zeros(size(1, 1), 28);
descriptive_global.biovolume_community = zeros(size(1, 1), 28);
descriptive_global.WW_community = zeros(size(1, 1), 28);
descriptive_global.CW_community = zeros(size(1, 1), 28);
functional_global.respiration_community = zeros(size(1, 1));
feeding_global.clearance_group = zeros(size(1, 1), num_steps);
respiration_global.respiration_group = zeros(size(1, 1), num_steps);
growth_global.growth_group = zeros(size(1, 1), num_steps);
feeding_global.prey_consumed_per_meso_carbon = zeros(size(1, 1));

%% 1.2. Data set preparation 
  
     for s=1:210
   if size(metaplanktonC(s).merged2netsUVPH20,2)>0
   %if size(metaplankton(s).merged4devicesmin,2)>0

ln_bvmedian=metaplanktonC(1).X;
widthbv=metaplanktonC(1).X1;%mm3
bvmedian=exp(ln_bvmedian);%mm3
ESDmatrix=metaplanktonC(178).merged.ESDmatrix;%µm
ESD=ESDmatrix(:,1);%µm

namepft=metaplanktonC(178).WP2.namepft;
namepft = namepft';

sizespectre_bv = metaplankton(s).merged2netsUVPH20.datapft; %mm3/mm3/m3
sizespectre_c = metaplanktonC(s).merged2netsUVPH20.datapft; %mg/mm3/m3

%sizespectre_bv = metaplankton(s).merged4devicesmin.datapft; %mm3/mm3/m3
%sizespectre_c = metaplanktonC(s).merged4devicesmin.datapft; %mg/mm3/m3

sizespectre_bv = array2table(sizespectre_bv);
sizespectre_c = array2table(sizespectre_c);

sizespectre_bv = renamevars(sizespectre_bv, {'sizespectre_bv1','sizespectre_bv2','sizespectre_bv3','sizespectre_bv4','sizespectre_bv5','sizespectre_bv6','sizespectre_bv7','sizespectre_bv8','sizespectre_bv9','sizespectre_bv10','sizespectre_bv11','sizespectre_bv12','sizespectre_bv13','sizespectre_bv14','sizespectre_bv15','sizespectre_bv16','sizespectre_bv17','sizespectre_bv18','sizespectre_bv19','sizespectre_bv20','sizespectre_bv21','sizespectre_bv22','sizespectre_bv23','sizespectre_bv24','sizespectre_bv25','sizespectre_bv26','sizespectre_bv27','sizespectre_bv28'},namepft); 
sizespectre_c = renamevars(sizespectre_c, {'sizespectre_c1','sizespectre_c2','sizespectre_c3','sizespectre_c4','sizespectre_c5','sizespectre_c6','sizespectre_c7','sizespectre_c8','sizespectre_c9','sizespectre_c10','sizespectre_c11','sizespectre_c12','sizespectre_c13','sizespectre_c14','sizespectre_c15','sizespectre_c16','sizespectre_c17','sizespectre_c18','sizespectre_c19','sizespectre_c20','sizespectre_c21','sizespectre_c22','sizespectre_c23','sizespectre_c24','sizespectre_c25','sizespectre_c26','sizespectre_c27','sizespectre_c28'},namepft); 

widthbv = array2table(widthbv);%mm3
bvmedian = array2table(bvmedian);%mm3
ESD = array2table(ESD);%µm

sizespectre_bv = [widthbv,bvmedian,ESD, sizespectre_bv];
sizespectre_c  = [widthbv,bvmedian,ESD, sizespectre_c];

temperature = metaplanktonC(s).hydro.Temp_10m; %from the hydro table en °C at 10m

group_meso = sizespectre_bv(:,["Appendicularia" "gelatinous_carnivorous" "Scyphozoa" "Thaliacea" "Ctenophora" "all_copepoda" "chaetognatha" "Hydrozoa" "gelatinous_filter_feeders" "pteropoda" "large_crustracean"]);
group_meso_carbon = sizespectre_c(:,["Appendicularia" "gelatinous_carnivorous" "Scyphozoa" "Thaliacea" "Ctenophora" "all_copepoda" "chaetognatha" "Hydrozoa" "gelatinous_filter_feeders" "pteropoda" "large_crustracean"]);


factor = sizespectre_bv.widthbv./sizespectre_bv.bvmedian;
factor_abundance = factor(1,1);%mm3/mm3


%% 2. Station S - Entire community

%% 2.1. Descriptive analysis

abundance_total = table2array(sizespectre_bv(:,4:31)) .* factor_abundance; %nb/m3
biovolume_total = abundance_total.* sizespectre_bv.bvmedian;%mm3/m3
% WW :
CW_station = table2array(sizespectre_c(:,4:31)) .* sizespectre_c.widthbv .* 10.^-3;%gC/m3
CW_station(isinf(CW_station)) = 0;

% CW :  
%Carbon percent of all community 
appendicularia = 0.07; chrysophyceae = 2.29; cryptophyta = 2.29; ctenophora = 0.205286; Dictyochophyceae = 2.29; Hydrozoa = 0.384141406; Phaeocystis = 2.29;
Scyphozoa = 0.491428569; Thaliacea = 1.081628667; copepod = 9.3930945454545455; chaetognatha = 3.56; ciliates = 2.29; cyanobacteria = 2.29; diatoms = 2.29; dinoflagelates = 2.29; 
flagellates = 2.29; gelatinous_carnivorous = 0.18; gelatinous_filter_feeders = 0.56; large_crustracean = 9.64; nano = 2.29; other_large = 9.3095; other_small = 2.29;
other_unidentified = 5.3095; pico = 2.29; plastic = 0; pteropoda = 5.52; rhizaria = 8.9; small_grazers = 9.5687; unicellular_hole = 2.29;

%carbon_all_b = [appendicularia,chrysophyceae,cryptophyta,ctenophora,Dictyochophyceae,Hydrozoa,Phaeocystis,Scyphozoa,Thaliacea,copepod,chaetognatha,ciliates,cyanobacteria,diatoms,dinoflagelates,flagellates,gelatinous_carnivorous,gelatinous_filter_feeders,large_crustracean,nano,other_large,other_small,other_unidentified,pico,plastic,pteropoda,rhizaria,small_grazers,unicellular_hole];
carbon_all_b = [appendicularia,chrysophyceae,cryptophyta,ctenophora,Dictyochophyceae,Hydrozoa,Phaeocystis,Scyphozoa,Thaliacea,copepod,chaetognatha,ciliates,cyanobacteria,diatoms,dinoflagelates,flagellates,gelatinous_carnivorous,gelatinous_filter_feeders,large_crustracean,nano,other_large,other_small,other_unidentified,pico,plastic,pteropoda,rhizaria,small_grazers];

appendicularia_a = 1; chrysophyceae_a  = 0.8809; cryptophyta_a  = 0.8809; ctenophora_a  = 1; Dictyochophyceae_a  = 0.8809; Hydrozoa_a  = 1; Phaeocystis_a  = 0.8809;
Scyphozoa_a  = 1; Thaliacea_a  = 1; copepod_a  = 1; chaetognatha_a  = 1; ciliates_a  = 0.8809; cyanobacteria_a  = 0.8809; diatoms_a  = 0.8809; dinoflagelates_a  = 0.8809;
flagellates_a  = 0.8809; gelatinous_carnivorous_a  = 1; gelatinous_filter_feeders_a  = 1; large_crustacean_a  = 1; nano_a  = 0.8809; other_large_a  = 1; other_small_a = 0.8809;
other_unidentified_a  = 1; pico_a  = 0.8809; plastic_a  = 1; pteropoda_a  = 1; rhizaria_a = 1; small_grazers_a  = 1; unicellular_hole_a = 0.8809;

%carbon_all_a = [appendicularia_a ,chrysophyceae_a ,cryptophyta_a ,ctenophora_a ,Dictyochophyceae_a ,Hydrozoa_a ,Phaeocystis_a ,Scyphozoa_a ,Thaliacea_a ,copepod_a ,chaetognatha_a ,ciliates_a ,cyanobacteria_a ,diatoms_a ,dinoflagelates_a ,flagellates_a ,gelatinous_carnivorous_a ,gelatinous_filter_feeders_a ,large_crustacean_a ,nano_a ,other_large_a,other_small_a ,other_unidentified_a ,pico_a ,plastic_a ,pteropoda_a ,rhizaria_a ,small_grazers_a,unicellular_hole_a ];
carbon_all_a = [appendicularia_a ,chrysophyceae_a ,cryptophyta_a ,ctenophora_a ,Dictyochophyceae_a ,Hydrozoa_a ,Phaeocystis_a ,Scyphozoa_a ,Thaliacea_a ,copepod_a ,chaetognatha_a ,ciliates_a ,cyanobacteria_a ,diatoms_a ,dinoflagelates_a ,flagellates_a ,gelatinous_carnivorous_a ,gelatinous_filter_feeders_a ,large_crustacean_a ,nano_a ,other_large_a,other_small_a ,other_unidentified_a ,pico_a ,plastic_a ,pteropoda_a ,rhizaria_a ,small_grazers_a];

WW_station = CW_station.^carbon_all_a ./ (carbon_all_b .*100);%gC/m3
WW_station(isinf(WW_station)) = 0;

CW_station_sum = sum(CW_station,2);%gC/m3

% Storage :

descriptive_global.abundance_community(s,:) = sum(abundance_total.*200,1); %nb/m2 
descriptive_global.biovolume_community(s,:) = sum(biovolume_total.*200,1);%mm3/m2
descriptive_global.WW_community(s,:) = sum(WW_station.*200,1);%g/m2
descriptive_global.CW_community(s,:) = sum(CW_station.*200,1);%g/m2

%% 2.2. Functional analysis : Respiration rate of the entire community

respiration_ind_community = CW_station.^0.87 .* 2.63 .* (10.^-0.008).^carbon_all_b; %mmO2.d-1.m3 at 15°C

respiration_ind_com_temp = 10.^(log10(respiration_ind_community)-log10(2.8)*(15-temperature)./10);%mmO2.d-1.m3 at temp station

% Storage :

functional_global.respiration_community(s,:)= sum(respiration_ind_com_temp.*200,"all");%mmO2.d-1.m3 at temp station

%% 4. Taxonomic group by taxonomic group

 for j = 1 : length(group_names)
   
%% 3.1. Descriptive analysis 

biovolume_class = sizespectre_bv(:,"bvmedian");%mm3
size_class = sizespectre_bv(:,"ESD");%µm
group_NBSS = group_meso(:,j); %in mm3.mm-3.m-3

group_physio = [biovolume_class, size_class, group_NBSS]; % Data with only group j
group_physio = renamevars(group_physio, {'bvmedian', 'ESD', group_names{j}}, {'biovolume', 'size_class', 'NBSS'}); %change the name of the column

group_physio.group_abundance = group_physio.NBSS .* factor_abundance; % Abundance 

if j == 1 % Appendicularia (special case : different percentage of carbon for respiration and growth)

%% 3.2. Functional analysis : physiological rates (clearance rate (L/d), respiration rate (mmolO2/d), mass-specific growth rate (1/d)

carbon_body_appendiculaire = 4.65; 
group_physio.CW_group = table2array(group_meso_carbon(:,j)) .* sizespectre_c.widthbv .* 10.^-3;%g/m3
group_physio.WW_group = group_physio.CW_group ./ carbon_body_appendiculaire .*100;%g/m3


% Best equations from Lemoine et al., 2025 - clearance specific to the size of an individual .* abundance (already taken into account in CW_group or WW_group)

group_physio.clearance_group = 13.46.*group_physio.WW_group.^0.72;%L.d-1.m-3 at 15°C

group_physio.respiration_group = group_physio.CW_group.^0.87 .* 2.63 .* (10.^-0.008).^carbon_body_appendiculaire;%mm02.d-1.m-3 at 15°C

group_physio.growth_group = 0.01.*group_physio.CW_group.^-0.2;%d-1.m-3 at 15°C

group_physio.growth_group(isinf(group_physio.growth_group)) = 0;

% Corrects the temperature to match that of the station (without correction, rates are calculated for a temperature of 15°C) 
group_physio.clearance_group_temp = 10.^(log10(group_physio.clearance_group)-log10(2.8)*(15-temperature)./10); % at temp station
group_physio.respiration_group_temp = 10.^(log10(group_physio.respiration_group)-log10(2.8)*(15-temperature)./10); % at temp station
group_physio.growth_group_temp = 10.^(log10(group_physio.growth_group)-log10(2.8)*(15-temperature)./10); % at temp station

% Storage :
feeding_global.clearance_group(s,j) = sum(group_physio.clearance_group_temp,1);%L.d-1.m-3
respiration_global.respiration_group(s,j) = sum(group_physio.respiration_group_temp,1);%mm02.d-1.m-3
growth_global.growth_group(s,j) = sum(group_physio.growth_group_temp,1);%d-1.m-3
feeding_global.clearance_group_m2(s,j) = sum(group_physio.clearance_group_temp.* 200,1);%L.d-1.m-2

%% 3.3. Predation pressure by group j

if sum(group_physio.group_abundance,1) ~= 0 % Only work when predators are present

% Predator/prey ratio matrix (x(column) = predator size / y(row) = prey size)
length_ESD = length(sizespectre_bv.ESD);
pred_prey_ratio_matrix = zeros(length_ESD,length_ESD);
for i = 1:length_ESD
    pred_prey_ratio_matrix(:,i) = sizespectre_bv.ESD(i) ./ sizespectre_bv.ESD;
    log_pred_prey_ratio_matrix = log10(pred_prey_ratio_matrix);
end

lower_borne_normal =  0.7131475 - 0.06293929 .* carbon_percent(j); %µmpredator.µmprey-1
upper_borne_normal = 4.504562 -  0.2293139   .* carbon_percent(j); %µmpredator.µmprey-1

range_cover_normal = upper_borne_normal - lower_borne_normal; 
strip_normal = range_cover_normal ./ 0.0251; %fit to the size class of the data

% Adjustment of the lower borne to match the pred/prey ratio matrix
correction_ratio = reshape(log_pred_prey_ratio_matrix, [], 1);
[~, idx] = min(abs(correction_ratio(:) - lower_borne_normal));
lower_borne = correction_ratio(idx);

normal_strip = zeros(length(round(strip_normal)));
normal_strip(1,1) = lower_borne;
for i= 2 : round(strip_normal)
    normal_strip(i,1) = normal_strip(i-1,1)+ 0.0251;
end

% Adjustment to match the pred/prey ratio matrix
length_normal = length(normal_strip);
for i = 1 : length_normal
     [~, index] = min(abs(correction_ratio(:) - normal_strip(i)));
     normal_strip(i,2) = correction_ratio(index); %size class same than ESD
end

n = length_normal;
x = linspace(-3, 3, n); % Ajustez la plage pour obtenir une distribution plus large
mu = 0; % Moyenne
sigma = 1; % Écart-type
normal_strip(:,3) = normpdf(x, mu, sigma);
normal_strip(:,3) = normal_strip(:,3) / max(normal_strip(:,3));
normal_strip(:,3) = normal_strip(:,3) + (1 - max(normal_strip(:,3))) / 2;


% Keep only the values within the normal distribution -> the other values become 0 because they cannot be eaten by the predator.
log_pred_prey_ratio_matrix(log_pred_prey_ratio_matrix < min(normal_strip(:,2)) | log_pred_prey_ratio_matrix > max(normal_strip(:,2))) = log10(0); 

% Replaces the pred/prey ratio with capture efficiency (probability)
x = normal_strip(:,2); % pred/prey ratio
y = normal_strip(:,3); % associeted capture efficiency
[~, index] = ismembertol(log_pred_prey_ratio_matrix, x); % find and replace with an index the values of x inside the log_pred_ratio
valid_indices = index > 0;

% Extract corresponding y values using valid indices
capture_efficiency = zeros(size(log_pred_prey_ratio_matrix)); % creation of the new matrice to store the probabilties
capture_efficiency(valid_indices) = y(index(valid_indices)); % remplace the index by the Y values (=percent)

prey_catchable = zeros(size(log_pred_prey_ratio_matrix));
prey_catchable(capture_efficiency > 0) = 1;

% Matrix with prey captchable by group j (x(column) = predator y(row) = size class of prey)
pred_size = zeros(length_ESD,length_ESD);
for i = 1:length_ESD
     pred_size(:,i) = sizespectre_bv.ESD(i) ./ pred_prey_ratio_matrix(:,i);
     inf_idx = isinf(pred_size);
     pred_size(inf_idx) = 0;
end

clearance_group_j = group_physio.clearance_group_temp';
clearance_group_m3 = clearance_group_j ./1000; % L.d-1 -> m3.d-1

x = sizespectre_bv.ESD; % extract ESD values
y = CW_station_sum;
[~, index] = ismembertol(pred_size, x); % find ESD values inside the data and change them with an index
valid_indices = index > 0;
% Extract corresponding y values using valid indices
prey_available_carbon = zeros(size(pred_size));
prey_available_carbon(valid_indices) = y(index(valid_indices));%change the index with the corresponding abundance values

consumed_carbon = prey_available_carbon .* clearance_group_m3 .* prey_catchable;%gC/d-1/m3

prey_consumed_carbon = sum(consumed_carbon,1);%gC/d-1/m3


prey_consumed_carbon_total = sum(prey_consumed_carbon,2).*200;%gC/d-1/m2

feeding_global.prey_consumed_per_meso_carbon(s,j) = prey_consumed_carbon_total; %gC/d-1/m2

feeding_global.specific_grazing_rate(s,j)= mean(prey_consumed_carbon' ./ (group_physio.CW_group.*200),'omitnan');

%% Detritus - fecal pellet

assimilation = 0.6;

production_global.detritus_m2(s,j) = (1-assimilation) .* prey_consumed_carbon_total;%gC/m2


%% respiration carbon
molecular_weight_carbon = 12.011;

% Respiration carbon

respiration_mol_O2_group = group_physio.respiration_group_temp .* 10.^-3; %mol02.d-1.m3

respiration_carbon_group = respiration_mol_O2_group .* molecular_weight_carbon .* 0.97;%gC.d-1.m3

respiration_carbon_m2 = sum(respiration_carbon_group .* 200);%gC.d-1.m2

respiration_global.respiration_carbon_m2(s,j) = respiration_carbon_m2;%gC.d-1.m2

respiration_global.respiration_specific(s,j) = respiration_carbon_m2 ./ sum(group_physio.CW_group .* 200);%gC/gC.d-1.m2

else


feeding_global.prey_consumed_per_meso_carbon(s,j) = 0;
feeding_global.specific_grazing_rate(s,j)=0;
production_global.detritus_m2(s,j) = 0;

respiration_global.respiration_carbon_m2(s,j) = 0;
respiration_global.respiration_specific(s,j) = 0;

end

%% For the other group :

else

group_physio.CW_group = table2array(group_meso_carbon(:,j)) .* sizespectre_c.widthbv .* 10.^-3;%gC.m3
group_physio.WW_group = group_physio.CW_group ./ carbon_percent(j) .*100;%gC.m3


%% 3.2. Functional analysis : physiological rates (clearance rate (L/d), respiration rate (mmolO2/d), mass-specific growth rate (1/d)

% Best equations from Lemoine et Lombard, 2024 - clearance specific to the size of an individual

group_physio.clearance_group = 13.46.*group_physio.WW_group.^0.72;%L.d-1.m3 at 15°C

group_physio.respiration_group = group_physio.CW_group.^0.87 .* 2.63 .* (10.^-0.008).^carbon_percent(j);%mmol02.d-1.m3 at 15°C
group_physio.growth_group = 0.01.*group_physio.CW_group.^-0.2;%d-1.m3 at 15°C
group_physio.growth_group(isinf(group_physio.growth_group)) = 0;

% Corrects the temperature to match that of the station (without correction, rates are calculated for a temperature of 15°C) 
group_physio.clearance_group_temp = 10.^(log10(group_physio.clearance_group)-log10(2.8)*(15-temperature)./10);%L.d-1.m3 at temp station
group_physio.respiration_group_temp = 10.^(log10(group_physio.respiration_group)-log10(2.8)*(15-temperature)./10);%mmol02.d-1.m3 at temp station
group_physio.growth_group_temp = 10.^(log10(group_physio.growth_group)-log10(2.8)*(15-temperature)./10);%d-1.m3 at temp station

% Storage :
feeding_global.clearance_group(s,j) = sum(group_physio.clearance_group_temp,1);%L.d-1.m3
respiration_global.respiration_group(s,j) = sum(group_physio.respiration_group_temp,1);%mmol02.d-1.m3
growth_global.growth_group(s,j) = sum(group_physio.growth_group_temp,1);%d-1.m3
feeding_global.clearance_group_m2(s,j) = sum(group_physio.clearance_group_temp .* 200,1);%L.d-1.m2

%% 3.3. Predation pressure by group j

if sum(group_physio.group_abundance,1) ~= 0 % Only work when predators are present

% Predator/prey ratio matrix (x(column) = predator size / y(row) = prey size)
length_ESD = length(sizespectre_bv.ESD);
pred_prey_ratio_matrix = zeros(length_ESD,length_ESD);
for i = 1:length_ESD
    pred_prey_ratio_matrix(:,i) = sizespectre_bv.ESD(i) ./ sizespectre_bv.ESD;
    log_pred_prey_ratio_matrix = log10(pred_prey_ratio_matrix);
end

lower_borne_normal =  0.7131475 - 0.06293929  .* carbon_percent(j);%µmpredator_µmprey-1
upper_borne_normal = 4.504562 -  0.2293139   .* carbon_percent(j);%µmpredator_µmprey-1

range_cover_normal = upper_borne_normal - lower_borne_normal; 
strip_normal = range_cover_normal / 0.0251; %fit to the size class of the data

% Adjustment of the lower borne to match the pred/prey ratio matrix
correction_ratio = reshape(log_pred_prey_ratio_matrix, [], 1);
[~, idx] = min(abs(correction_ratio(:) - lower_borne_normal));
lower_borne = correction_ratio(idx);

normal_strip = zeros(length(round(strip_normal)));
normal_strip(1,1) = lower_borne;
for i= 2 : round(strip_normal)
    normal_strip(i,1) = normal_strip(i-1,1)+ 0.0251;
end

% Adjustment to match the pred/prey ratio matrix
length_normal = length(normal_strip);
for i = 1 : length_normal
     [~, index] = min(abs(correction_ratio(:) - normal_strip(i)));
     normal_strip(i,2) = correction_ratio(index); %size class same than ESD
end
n = length_normal;
x = linspace(-3, 3, n); % Ajustez la plage pour obtenir une distribution plus large
mu = 0; % Moyenne
sigma = 1; % Écart-type
normal_strip(:,3) = normpdf(x, mu, sigma);
normal_strip(:,3) = normal_strip(:,3) / max(normal_strip(:,3));
normal_strip(:,3) = normal_strip(:,3) + (1 - max(normal_strip(:,3))) / 2;

% Keep only the values within the normal distribution -> the other values become 0 because they cannot be eaten by the predator.
log_pred_prey_ratio_matrix(log_pred_prey_ratio_matrix < min(normal_strip(:,2)) | log_pred_prey_ratio_matrix > max(normal_strip(:,2))) = log10(0); 

% Replaces the pred/prey ratio with capture efficiency (probability)
x = normal_strip(:,2); % pred/prey ratio
y = normal_strip(:,3); % associeted capture efficiency
[~, index] = ismembertol(log_pred_prey_ratio_matrix, x); % find and replace with an index the values of x inside the log_pred_ratio
valid_indices = index > 0;
% Extract corresponding y values using valid indices
capture_efficiency = zeros(size(log_pred_prey_ratio_matrix)); % creation of the new matrice to store the probabilties
capture_efficiency(valid_indices) = y(index(valid_indices)); % remplace the index by the Y values (=percent)

prey_catchable = zeros(size(log_pred_prey_ratio_matrix));
prey_catchable(capture_efficiency > 0) = 1;

% Matrix with prey captchable by group j (x(column) = predator y(row) = size class of prey)
pred_size = zeros(length_ESD,length_ESD);
for i = 1:length_ESD
     pred_size(:,i) = sizespectre_bv.ESD(i) ./ pred_prey_ratio_matrix(:,i);
     inf_idx = isinf(pred_size);
     pred_size(inf_idx) = 0;
end

% Clearance rate of group j according to the size
clearance_group_j = group_physio.clearance_group_temp';
clearance_group_m3 = clearance_group_j ./1000; % L.d-1 -> m3.d-1

x = sizespectre_bv.ESD; % extract ESD values
y = CW_station_sum;
[~, index] = ismembertol(pred_size, x); % find ESD values inside the data and change them with an index
valid_indices = index > 0;
% Extract corresponding y values using valid indices
prey_available_carbon = zeros(size(pred_size));
prey_available_carbon(valid_indices) = y(index(valid_indices));%change the index with the corresponding abundance values

consumed_carbon = prey_available_carbon .* clearance_group_m3 .* prey_catchable;%gC/d-1/m3

prey_consumed_carbon = sum(consumed_carbon,1);%gC/d-1/m3

prey_consumed_carbon_total = sum(prey_consumed_carbon,2).*200;%gC/d-1/m2

feeding_global.prey_consumed_per_meso_carbon(s,j) = prey_consumed_carbon_total; %gC/d-1/m2

feeding_global.specific_grazing_rate(s,j)= mean(prey_consumed_carbon' ./ (group_physio.CW_group.*200),'omitnan');%gC/gC/d-1/m2

%% Detritus - fecal pellet

assimilation = 0.6;

production_global.detritus_m2(s,j) = (1-assimilation) .* prey_consumed_carbon_total;%gC/d-1/m2

%% Predicted growth

molecular_weight_carbon = 12.011;

% Respiration carbon

respiration_mol_O2_group = group_physio.respiration_group_temp .* 10.^-3; %mol02.d-1.m3

respiration_carbon_group = respiration_mol_O2_group .* molecular_weight_carbon .* 0.97;%gC.d-1.m3

respiration_carbon_m2 = sum(respiration_carbon_group .* 200);%gC.d-1.m3

respiration_global.respiration_carbon_m2(s,j) = respiration_carbon_m2;%gC.d-1.m3

respiration_global.respiration_specific(s,j) = respiration_carbon_m2 ./ sum(group_physio.CW_group .* 200);%gC/gC.d-1.m3

% RQ = 0.97

else

feeding_global.prey_consumed_per_meso_carbon(s,j) = 0;
feeding_global.specific_grazing_rate(s,j)=0;
production_global.detritus_m2(s,j) = 0;

respiration_global.respiration_carbon_m2(s,j) = 0;
respiration_global.respiration_specific(s,j) = 0;

end

end

 end

indice = 129;

sum_nbr_prey_before = sum(y(1:indice),1).*200;

sum_nbr_prey_consumed = sum(feeding_global.prey_consumed_per_meso_carbon(s,:),2); 
sum_nbr_prey_consumed(isinf(sum_nbr_prey_consumed)) = 0;

feeding_global.percentage_prey_consumed(s,:) = (sum_nbr_prey_consumed./sum_nbr_prey_before).*100;

   end

    end

row_abundance_community = sum(descriptive_global.abundance_community, 2);
nonzero_abundance_community = row_abundance_community ~= 0;
descriptive_global.abundance_community = descriptive_global.abundance_community(nonzero_abundance_community, :);

row_biovolume_community = sum(descriptive_global.biovolume_community, 2);
nonzero_biovolume_community = row_biovolume_community ~= 0;
descriptive_global.biovolume_community = descriptive_global.biovolume_community(nonzero_biovolume_community, :);

row_WW = sum(descriptive_global.WW_community, 2);
nonzero_WW = row_WW ~= 0;
descriptive_global.WW_community = descriptive_global.WW_community(nonzero_WW, :);

row_CW = sum(descriptive_global.CW_community, 2);
nonzero_CW = row_CW ~= 0;
descriptive_global.CW_community = descriptive_global.CW_community(nonzero_CW, :);

row_respiration_community = sum(functional_global.respiration_community, 2);
nonzero_respiration_community = row_respiration_community ~= 0;
functional_global.respiration_community = functional_global.respiration_community(nonzero_respiration_community, :);

row_clearance_group = sum(feeding_global.clearance_group, 2);
nonzero_clearance_group = row_clearance_group ~= 0;
feeding_global.clearance_group = feeding_global.clearance_group(nonzero_clearance_group, :);

row_respiration_group = sum(respiration_global.respiration_group, 2);
nonzero_respiration_group = row_respiration_group ~= 0;
respiration_global.respiration_group = respiration_global.respiration_group(nonzero_respiration_group, :);

row_growth_group = sum(growth_global.growth_group, 2);
nonzero_growth_group = row_growth_group ~= 0;
growth_global.growth_group = growth_global.growth_group(nonzero_growth_group, :);

row_sum_nbr_prey_consumed_j_carbon = sum(feeding_global.prey_consumed_per_meso_carbon, 2);
nonzero_sum_nbr_prey_consumed_j_carbon = row_sum_nbr_prey_consumed_j_carbon ~= 0;
feeding_global.prey_consumed_per_meso_carbon = feeding_global.prey_consumed_per_meso_carbon(nonzero_sum_nbr_prey_consumed_j_carbon, :);

row_percentage_prey_consumed = sum(feeding_global.percentage_prey_consumed, 2);
nonzero_percentage_prey_consumed = row_percentage_prey_consumed ~= 0;
feeding_global.percentage_prey_consumed = feeding_global.percentage_prey_consumed(nonzero_percentage_prey_consumed, :);

row_detritus = sum(production_global.detritus_m2, 2);
nonzero_detritus = row_detritus ~= 0;
production_global.detritus_m2 = production_global.detritus_m2(nonzero_detritus, :);

row_respiration_specific = sum(respiration_global.respiration_specific, 2);
nonzero_respiration_specific = row_respiration_specific ~= 0;
respiration_global.respiration_specific = respiration_global.respiration_specific(nonzero_respiration_specific, :);

row_respiration_carbon_m2 = sum(respiration_global.respiration_carbon_m2, 2);
nonzero_respiration_carbon_m2 = row_respiration_carbon_m2 ~= 0;
respiration_global.respiration_carbon_m2 = respiration_global.respiration_carbon_m2(nonzero_respiration_carbon_m2, :);

row_clearance_group_m2 = sum(feeding_global.clearance_group_m2, 2);
nonzero_clearance_group_m2 = row_clearance_group_m2 ~= 0;
feeding_global.clearance_group_m2 = feeding_global.clearance_group_m2(nonzero_clearance_group_m2, :);

row_specific_grazing_rate = sum(feeding_global.specific_grazing_rate, 2);
nonzero_specific_grazing_rate = row_specific_grazing_rate ~= 0;
feeding_global.specific_grazing_rate = feeding_global.specific_grazing_rate(nonzero_specific_grazing_rate, :);

    community_global.descriptive = descriptive_global;
    community_global.functional = functional_global;

    mesozooplankton_global.feeding = feeding_global;
    mesozooplankton_global.respiration = respiration_global;
    mesozooplankton_global.growth = growth_global;
    mesozooplankton_global.production = production_global;

    diagnostic_global.community = community_global;
    diagnostic_global.mesozooplankton = mesozooplankton_global;