
function [community,mesozooplankton] = station_function_diagnostic(sizespectre_bv, sizespectre_c, carbon_percent, temperature)

%% 1. Initialization

%% 1.1. Matrix preparation

group_meso = sizespectre_bv(:,["Appendicularia" "gelatinous_carnivorous" "Scyphozoa" "Thaliacea" "Ctenophora" "all_copepoda" "chaetognatha" "Hydrozoa" "gelatinous_filter_feeders" "pteropoda" "large_crustracean"]);
group_meso_carbon = sizespectre_c(:,["Appendicularia" "gelatinous_carnivorous" "Scyphozoa" "Thaliacea" "Ctenophora" "all_copepoda" "chaetognatha" "Hydrozoa" "gelatinous_filter_feeders" "pteropoda" "large_crustracean"]);
group_names = {'Appendicularia' 'gelatinous_carnivorous' 'Scyphozoa' 'Thaliacea' 'Ctenophora' 'all_copepoda' 'chaetognatha' 'Hydrozoa' 'gelatinous_filter_feeders' 'pteropoda' 'large_crustracean'};
predator_names = ["appendicularia" "gelatinous_carnivorous" "scyphozoa" "thaliacea" "ctenophora" "copepod" "chaetognatha" "hydrozoa" "gelatinous_filter_feeders" "pteropod" "large_crustracean"];

num_steps = length(group_names); 
clearance_group_station = zeros(size(1, 1), num_steps);
respiration_group_station = zeros(size(1, 1), num_steps);
growth_group_station = zeros(size(1, 1), num_steps);

factor = sizespectre_bv.largeurBV./sizespectre_bv.BvMedian;%mm3/mm3
factor_abundance = factor(1,1);

%% 2. Global community
abundance_station = table2array(sizespectre_bv(:,4:31)) .* factor_abundance;%nbr/mm3
biovolume_station = abundance_station .* sizespectre_bv.BvMedian;%mm3/mm3

% WW :
CW_station = table2array(sizespectre_c(:,4:31)) .* sizespectre_c.largeurBV .* 10.^-3;%gC/m3
CW_station(isinf(CW_station)) = 0;

% CW :  
%Carbon percent of all community 
appendicularia = 0.07; chrysophyceae = 2.29; cryptophyta = 2.29; ctenophora = 0.205286; Dictyochophyceae = 2.29; Hydrozoa = 0.384141406; Phaeocystis = 2.29;
Scyphozoa = 0.491428569; Thaliacea = 1.081628667; copepod = 9.3930945454545455; chaetognatha = 3.56; ciliates = 2.29; cyanobacteria = 2.29; diatoms = 2.29; dinoflagelates = 2.29; 
flagellates = 2.29; gelatinous_carnivorous = 0.18; gelatinous_filter_feeders = 0.56; large_crustracean = 9.64; nano = 2.29; other_large = 9.3095; other_small = 2.29;
other_unidentified = 5.3095; pico = 2.29; plastic = 0; pteropoda = 5.52; rhizaria = 8.9; small_grazers = 9.5687; 

%carbon_all_b = [appendicularia,chrysophyceae,cryptophyta,ctenophora,Dictyochophyceae,Hydrozoa,Phaeocystis,Scyphozoa,Thaliacea,copepod,chaetognatha,ciliates,cyanobacteria,diatoms,dinoflagelates,flagellates,gelatinous_carnivorous,gelatinous_filter_feeders,large_crustracean,nano,other_large,other_small,other_unidentified,pico,plastic,pteropoda,rhizaria,small_grazers,unicellular_hole];
carbon_all_b = [appendicularia,chrysophyceae,cryptophyta,ctenophora,Dictyochophyceae,Hydrozoa,Phaeocystis,Scyphozoa,Thaliacea,copepod,chaetognatha,ciliates,cyanobacteria,diatoms,dinoflagelates,flagellates,gelatinous_carnivorous,gelatinous_filter_feeders,large_crustracean,nano,other_large,other_small,other_unidentified,pico,plastic,pteropoda,rhizaria,small_grazers];

appendicularia_a = 1; chrysophyceae_a  = 0.8809; cryptophyta_a  = 0.8809; ctenophora_a  = 1; Dictyochophyceae_a  = 0.8809; Hydrozoa_a  = 1; Phaeocystis_a  = 0.8809;
Scyphozoa_a  = 1; Thaliacea_a  = 1; copepod_a  = 1; chaetognatha_a  = 1; ciliates_a  = 0.8809; cyanobacteria_a  = 0.8809; diatoms_a  = 0.8809; dinoflagelates_a  = 0.8809;
flagellates_a  = 0.8809; gelatinous_carnivorous_a  = 1; gelatinous_filter_feeders_a  = 1; large_crustacean_a  = 1; nano_a  = 0.8809; other_large_a  = 1; other_small_a = 0.8809;
other_unidentified_a  = 1; pico_a  = 0.8809; plastic_a  = 0; pteropoda_a  = 1; rhizaria_a = 1; small_grazers_a  = 1;

%carbon_all_a = [appendicularia_a ,chrysophyceae_a ,cryptophyta_a ,ctenophora_a ,Dictyochophyceae_a ,Hydrozoa_a ,Phaeocystis_a ,Scyphozoa_a ,Thaliacea_a ,copepod_a ,chaetognatha_a ,ciliates_a ,cyanobacteria_a ,diatoms_a ,dinoflagelates_a ,flagellates_a ,gelatinous_carnivorous_a ,gelatinous_filter_feeders_a ,large_crustacean_a ,nano_a ,other_large_a,other_small_a ,other_unidentified_a ,pico_a ,plastic_a ,pteropoda_a ,rhizaria_a ,small_grazers_a,unicellular_hole_a ];
carbon_all_a = [appendicularia_a ,chrysophyceae_a ,cryptophyta_a ,ctenophora_a ,Dictyochophyceae_a ,Hydrozoa_a ,Phaeocystis_a ,Scyphozoa_a ,Thaliacea_a ,copepod_a ,chaetognatha_a ,ciliates_a ,cyanobacteria_a ,diatoms_a ,dinoflagelates_a ,flagellates_a ,gelatinous_carnivorous_a ,gelatinous_filter_feeders_a ,large_crustacean_a ,nano_a ,other_large_a,other_small_a ,other_unidentified_a ,pico_a ,plastic_a ,pteropoda_a ,rhizaria_a ,small_grazers_a];

WW_station = CW_station.^carbon_all_a ./ (carbon_all_b .*100);%gC/m3
WW_station(isinf(WW_station)) = 0;

CW_station_sum = sum(CW_station,2);%gC/m3

% Storage
community.abundance_community_total = abundance_station;
community.biovolume_community_total = biovolume_station;
community.CW_community_total = CW_station;
community.WW_community_total = WW_station;

%% 3. Taxonomic group by taxonomic group

for j = 1 : length(group_names)

if j == 1 % Appendicularia (special case : different percentage of carbon for respiration and growth)
 %% 3.1. Descriptive analysis 

biovolume_class = sizespectre_bv(:,"BvMedian");%mm3
size_class = sizespectre_bv(:,"ESD_micron_AllUnitsInMm3_Mm3_m3");%µm
group_NBSS = group_meso(:,j); %in mm3.mm-3.m-3
group_physio = [biovolume_class, size_class, group_NBSS]; % Data with only group j
group_physio = renamevars(group_physio, {'BvMedian', 'ESD_micron_AllUnitsInMm3_Mm3_m3', group_names{j}}, {'biovolume', 'size_class', 'NBSS'}); %change the name of the column

group_physio.group_abundance = group_physio.NBSS .* factor_abundance; % Abundance nbr/mm3
group_physio.group_biovolume = group_physio.group_abundance .* group_physio.biovolume; % Abundance nbr/mm3

carbon_body_appendiculaire = 4.65; 

group_physio.CW_group = table2array(group_meso_carbon(:,j)) .* sizespectre_c.largeurBV .* 10.^-3;%g/m3
group_physio.WW_group = group_physio.CW_group ./ carbon_body_appendiculaire .*100;%g/m3

%% 3.2. Functional analysis : physiological rates (clearance rate (L/d), respiration rate (mmolO2/d), mass-specific growth rate (1/d)

group_physio.clearance_group = 13.46.*group_physio.WW_group.^0.72;%L.d-1.m-3 at 15°C

group_physio.respiration_group = group_physio.CW_group.^0.87 .* 2.63 .* (10.^-0.008).^carbon_body_appendiculaire;%mm02.d-1.m-3 at 15°C

group_physio.growth_group = 0.01.*group_physio.CW_group.^-0.2;%d-1.m-3 at 15°C

group_physio.growth_group(isinf(group_physio.growth_group)) = 0;

% Corrects the temperature to match that of the station (without correction, rates are calculated for a temperature of 15°C) 
group_physio.clearance_group_temp = 10.^(log10(group_physio.clearance_group)-log10(2.8)*(15-temperature)./10); % at temp station
group_physio.respiration_group_temp = 10.^(log10(group_physio.respiration_group)-log10(2.8)*(15-temperature)./10); % at temp station
group_physio.growth_group_temp = 10.^(log10(group_physio.growth_group)-log10(2.8)*(15-temperature)./10); % at temp station

% Storage :
mesozooplankton.clearance_group_station(:,j) = sum(group_physio.clearance_group_temp,1);
mesozooplankton.respiration_group_station(:,j) = sum(group_physio.respiration_group_temp,1);
mesozooplankton.growth_group_station(:,j) = sum(group_physio.growth_group_temp,1);
mesozooplankton.CW_size_class(:,j) = group_physio.CW_group;
mesozooplankton.WW_size_class(:,j) = group_physio.WW_group;
mesozooplankton.abundance(:,j) = group_physio.group_abundance;
mesozooplankton.biovolume(:,j) = group_physio.group_biovolume;

%% 3.3. Predation pressure by group j

if sum(group_physio.group_abundance ,1) ~= 0 %work only when predator are present

length_ESD = length(sizespectre_bv.ESD_micron_AllUnitsInMm3_Mm3_m3);
pred_prey_ratio_matrix = zeros(length_ESD,length_ESD);
for i = 1:length_ESD
    pred_prey_ratio_matrix(:,i) = sizespectre_bv.ESD_micron_AllUnitsInMm3_Mm3_m3(i) ./ sizespectre_bv.ESD_micron_AllUnitsInMm3_Mm3_m3;
    log_pred_prey_ratio_matrix = log10(pred_prey_ratio_matrix);
end

lower_borne_normal =  0.7131475 - 0.06293929 .* carbon_percent(j);%µmpredator.µmprey-1
upper_borne_normal = 4.504562 -  0.2293139   .* carbon_percent(j);%µmpredator.µmprey-1

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
     pred_size(:,i) = sizespectre_bv.ESD_micron_AllUnitsInMm3_Mm3_m3(i) ./ pred_prey_ratio_matrix(:,i);
     inf_idx = isinf(pred_size);
     pred_size(inf_idx) = 0;
end

clearance_group_j = group_physio.clearance_group_temp';
clearance_group_m3 = clearance_group_j ./1000; % L.d-1 -> m3.d-1

x = sizespectre_bv.ESD_micron_AllUnitsInMm3_Mm3_m3; % extract ESD values
y = CW_station_sum;
[~, index] = ismembertol(pred_size, x); % find ESD values inside the data and change them with an index
valid_indices = index > 0;
% Extract corresponding y values using valid indices
prey_available_carbon = zeros(size(pred_size));
prey_available_carbon(valid_indices) = y(index(valid_indices));%change the index with the corresponding abundance values

consumed_carbon = prey_available_carbon .* clearance_group_m3 .* prey_catchable;%gC/d-1/m3

prey_consumed_carbon = sum(consumed_carbon,1);%gC/d-1/m3

mesozooplankton.nbr_prey_consumed_size_carbon(:,j) = prey_consumed_carbon';
mesozooplankton.capture_efficacity_percent(:,j) = sum(capture_efficiency,2)/length(capture_efficiency(1,:));

end

else
 %% 3.1. Descriptive analysis 

biovolume_class = sizespectre_bv(:,"BvMedian");%mm3
size_class = sizespectre_bv(:,"ESD_micron_AllUnitsInMm3_Mm3_m3");%µm
group_NBSS = group_meso(:,j); %in mm3.mm-3.m-3
group_physio = [biovolume_class, size_class, group_NBSS]; % Data with only group j
group_physio = renamevars(group_physio, {'BvMedian', 'ESD_micron_AllUnitsInMm3_Mm3_m3', group_names{j}}, {'biovolume', 'size_class', 'NBSS'}); %change the name of the column

group_physio.group_abundance = group_physio.NBSS .* factor_abundance; % Abundance nbr/mm3
group_physio.group_biovolume = group_physio.group_abundance .* group_physio.biovolume; % Abundance nbr/mm3

group_physio.CW_group = table2array(group_meso_carbon(:,j)) .* sizespectre_c.largeurBV .* 10.^-3;%g/m3
group_physio.WW_group = group_physio.CW_group ./ carbon_percent(j) .*100;%g/m3

%% 3.2. Functional analysis : physiological rates (clearance rate (L/d), respiration rate (mmolO2/d), mass-specific growth rate (1/d)

group_physio.clearance_group = 13.46.*group_physio.WW_group.^0.72;%L.d-1.m-3 at 15°C

group_physio.respiration_group = group_physio.CW_group.^0.87 .* 2.63 .* (10.^-0.008).^carbon_percent(j);%mm02.d-1.m-3 at 15°C

group_physio.growth_group = 0.01.*group_physio.CW_group.^-0.2;%d-1.m-3 at 15°C

group_physio.growth_group(isinf(group_physio.growth_group)) = 0;

% Corrects the temperature to match that of the station (without correction, rates are calculated for a temperature of 15°C) 
group_physio.clearance_group_temp = 10.^(log10(group_physio.clearance_group)-log10(2.8)*(15-temperature)./10); % at temp station
group_physio.respiration_group_temp = 10.^(log10(group_physio.respiration_group)-log10(2.8)*(15-temperature)./10); % at temp station
group_physio.growth_group_temp = 10.^(log10(group_physio.growth_group)-log10(2.8)*(15-temperature)./10); % at temp station

% Storage :
mesozooplankton.clearance_group_station(:,j) = sum(group_physio.clearance_group_temp,1);
mesozooplankton.respiration_group_station(:,j) = sum(group_physio.respiration_group_temp,1);
mesozooplankton.growth_group_station(:,j) = sum(group_physio.growth_group_temp,1);
mesozooplankton.CW_size_class(:,j) = group_physio.CW_group;
mesozooplankton.WW_size_class(:,j) = group_physio.WW_group;
mesozooplankton.abundance(:,j) = group_physio.group_abundance;
mesozooplankton.biovolume(:,j) = group_physio.group_biovolume;

%% 3.3. Predation pressure by group j

if sum(group_physio.group_abundance ,1) ~= 0 %work only when predator are present

length_ESD = length(sizespectre_bv.ESD_micron_AllUnitsInMm3_Mm3_m3);
pred_prey_ratio_matrix = zeros(length_ESD,length_ESD);
for i = 1:length_ESD
    pred_prey_ratio_matrix(:,i) = sizespectre_bv.ESD_micron_AllUnitsInMm3_Mm3_m3(i) ./ sizespectre_bv.ESD_micron_AllUnitsInMm3_Mm3_m3;
    log_pred_prey_ratio_matrix = log10(pred_prey_ratio_matrix);
end

lower_borne_normal =  0.7131475 - 0.06293929 .* carbon_percent(j);%µmpredator.µmprey-1
upper_borne_normal = 4.504562 -  0.2293139   .* carbon_percent(j);%µmpredator.µmprey-1

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
     pred_size(:,i) = sizespectre_bv.ESD_micron_AllUnitsInMm3_Mm3_m3(i) ./ pred_prey_ratio_matrix(:,i);
     inf_idx = isinf(pred_size);
     pred_size(inf_idx) = 0;
end

clearance_group_j = group_physio.clearance_group_temp';
clearance_group_m3 = clearance_group_j ./1000; % L.d-1 -> m3.d-1

x = sizespectre_bv.ESD_micron_AllUnitsInMm3_Mm3_m3; % extract ESD values
y = CW_station_sum;
[~, index] = ismembertol(pred_size, x); % find ESD values inside the data and change them with an index
valid_indices = index > 0;
% Extract corresponding y values using valid indices
prey_available_carbon = zeros(size(pred_size));
prey_available_carbon(valid_indices) = y(index(valid_indices));%change the index with the corresponding abundance values

consumed_carbon = prey_available_carbon .* clearance_group_m3 .* prey_catchable;%gC/d-1/m3

prey_consumed_carbon = sum(consumed_carbon,1);%gC/d-1/m3

mesozooplankton.nbr_prey_consumed_size_carbon(:,j) = prey_consumed_carbon';

column_size_class_pred = find(table2array(group_NBSS(:,1)) > 0);
capture_efficiency_total = sum(capture_efficiency(:, column_size_class_pred),2);
capture_efficiency_rel = capture_efficiency_total / max(capture_efficiency_total);

mesozooplankton.capture_efficacity_percent(:,j) = capture_efficiency_rel;

end

end

end