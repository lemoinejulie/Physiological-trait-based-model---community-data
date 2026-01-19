%%%%%%%%%%%%%%%%%%%
%% Initialisation%%
%%%%%%%%%%%%%%%%%%%

X1=metaplanktonC(1).X1;
latWP2=[];
lonWP2=[];
abund=[];
site=[];
for i=1:210
    if size(metaplankton(i).merged2netsUVPH20,2)>0
        site=[site i];
        abund=[abund; metaplanktonC(i).merged2netsUVPH20.donut];
        %abund = [abund;metaplankton_pft(i).merged4devicesmin.donut];
        latWP2=[latWP2; metaplanktonC(i).Lat];
        lonWP2=[lonWP2; metaplanktonC(i).Lon];
    end
end

ocean_area_km2 = [ ...
    0 0 522000 2604000 6816000 10301000 12006000 13388000 14693000 15833000 ...
    16483000 15782000 15438000 15450000 16147000 17211000 16898000 16792000 ...
    17387000 16628000 16553000 14981000 13354000 11747000 10806000 10029000 ...
    8411000 6612000 5529000 5399000 3123000 2456000 4414000 3742000 2545000 979000 ]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot polynôme in function of latitude%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
values = sum(diagnostic_global.mesozooplankton.production.detritus_m2,2);

degree =4; % Degré du polynôme (vous pouvez ajuster ceci selon votre besoin)

coefficients = polyfit(latWP2,log10(values), degree); % Ajustement polynômial

latitude_fit = linspace(min(latWP2), max(latWP2), 100); % Valeurs de respiration pour la courbe
values_fit = polyval(coefficients, latitude_fit); % Calcul des valeurs de latitude correspondantes

hold on; % Pour superposer les graphiques

% Tracer la courbe ajustée en rouge
plot(latitude_fit, values_fit, 'r');

% Premier scatter : ronds pleins (bleu par défaut ou autre couleur)
scatter(latWP2, log10(values), 'o', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'red');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot polynôme in function of latitude comparaison predation pressure%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lat = [55; 55;52;54.30;37;43;43;43;43;43];
percentage = [5;2;20;8;50;41;49;38;32;28];
percent_log = (percentage);
lat_2 = [-40;-23.75;27.0833;27.0833;-40;27.0833;27.0833;-40;33.3333;18.3333;-40;-0.8333;33.3333;33.3333;-23.75;-0.8333;-23.75;-40;-0.8333;-0.8333;-23.75;9.5833;9.5833;9.5833;9.5833;14.1667;14.1667;33.3333;33.3333;33.3333;33.3333];
percentage_2 = [12.28;9.30;5.71;6.94;41.71;5.87;2.25;44.71;2.59;9.83;52.83;73.73;18.37;6.66;4.15;21.11;9.30;12.62;37.84;26.36;9.17;19.97;38.91;13.72;28.66;19.42;10.24;100;41.71;12.45;1.73];
percent_log2 = (percentage_2);

values = diagnostic_global.mesozooplankton.feeding.percentage_prey_consumed;
values(values > 100) = 100;

degree =4; % Degré du polynôme (vous pouvez ajuster ceci selon votre besoin)

coefficients = polyfit(latWP2,(values), degree); % Ajustement polynômial

latitude_fit = linspace(min(latWP2), max(latWP2), 100); % Valeurs de respiration pour la courbe
values_fit = polyval(coefficients, latitude_fit); % Calcul des valeurs de latitude correspondantes

hold on; % Pour superposer les graphiques

% Tracer la courbe ajustée en rouge
plot(latitude_fit, values_fit, 'r');

% Premier scatter : ronds pleins (bleu par défaut ou autre couleur)
scatter(latWP2, (values), 'o', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'red');

% Deuxième scatter : ronds creux (bord bleu, intérieur blanc)
scatter(lat, (percent_log), 'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'blue');
scatter(lat_2,(percent_log2), 'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'blue');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot polynôme in function of latitude comparaison respiration%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
lat = [-69.0149;-64.7164;-63.7612;-61.3731;-60.8955;-58.9851;-58.5075;-58.5075;-40.3582;-28.8955;-23.1642;-21.7313;-2.1493;-1.194;0.7164;2.1493;3.582;6.9254;9.791;12.6567;13.1343;17.4328;18.8657;23.6418;24.597;24.597;27.4627;28.8955;29.3731;30.806;37.4925;39.403;40.3582;43.2239;44.6567;51.8209;54.6866;69.9701;71.403;76.1791];
data = [0.0157;0.0039;0.0157;0.0404;0.0365;0.0143;0.0183;0.0078;0.1565;0.1943;0.1487;0.1748;0.2165;0.2009;0.283;0.167;0.1865;0.1552;0.2283;0.1696;0.2048;0.1852;0.1409;0.0926;0.1865;0.1043;0.1996;0.1096;0.2074;0.1487;0.1396;0.0796;0.0861;0.0678;0.0926;0.0483;0.0248;0.013;0.0209;0.0222];

close all

b =  diagnostic_global.mesozooplankton.respiration.respiration_specific(:,6);

degree =4; % Degré du polynôme (vous pouvez ajuster ceci selon votre besoin)


coefficients = polyfit(latWP2,log10(b), degree); % Ajustement polynômial

latitude_fit = linspace(min(latWP2), max(latWP2), 100); % Valeurs de respiration pour la courbe
respiration_fit = polyval(coefficients, latitude_fit); % Calcul des valeurs de latitude correspondantes

% Tracer la courbe ajustée

hold on; % Pour superposer les graphiques

% Tracer la courbe ajustée en rouge
plot(latitude_fit, respiration_fit, 'r');

% Premier scatter : ronds pleins (bleu par défaut ou autre couleur)
scatter(latWP2, log10(b), 'o', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'red');

% Deuxième scatter : ronds creux (bord bleu, intérieur blanc)
scatter(lat, log10(data), 'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'blue');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Global budget latitudinal band%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all 
group_names = {'Appendicularia' 'gelatinous_carnivorous' 'Scyphozoa' 'Thaliacea' 'Ctenophora' 'all_copepoda' 'chaetognatha' 'Hydrozoa' 'gelatinous_filter_feeders' 'pteropoda' 'large_crustracean'};
%[1:5, 8:9]
%[2,3,8]
variable_sum = sum(diagnostic_global.mesozooplankton.feeding.prey_consumed_per_meso_carbon(:,6),2);
%z = zscore(variable_sum);
%threshold = 1;
%mask_no_outlier = abs(z) < threshold;
%variable_sum_clean = variable_sum(mask_no_outlier);
%latWP2_clean = latWP2(mask_no_outlier); 

threshold = prctile(variable_sum, 95); % par exemple garder les 98% les plus bas
mask = variable_sum < threshold;
variable_sum_clean = variable_sum(mask);
latWP2_clean = latWP2(mask);

latitude_fit = [-90:5:-5, 5:5:90];

degree =4; % Degré du polynôme

logvariable_sum_clean=variable_sum_clean;
logvariable_sum_clean(logvariable_sum_clean == 0) = 1;
ft = fittype(sprintf('poly%d', degree));
[curve, gof] = fit(latWP2_clean(:), log10(logvariable_sum_clean(:)), ft, 'Robust', 'Bisquare');

respiration_fit_log = feval(curve, latitude_fit);
respiration_fit = 10.^respiration_fit_log;

max_val = max(variable_sum_clean);
respiration_fit = min(respiration_fit, max_val);


bilan_band = respiration_fit .* ocean_area_km2*10^6 .* 365 .* 10^-15;
global_bilan = sum(bilan_band,'omitnan')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Global budget incertitudes%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%y_fit_log = polyval(coefficients, latWP2_clean);
y_fit_log = feval(curve, latWP2_clean);

residuals = log10(logvariable_sum_clean) - y_fit_log;
stderr_log = std(residuals); % en base log10
bilan_upper = global_bilan * 10.^stderr_log;
bilan_lower = global_bilan / 10.^stderr_log;

incertitude = (bilan_upper - bilan_lower) / 2

%%%%%%%%%%%%%%%%%%%%%%
%% Bar plot latitude%%
%%%%%%%%%%%%%%%%%%%%%%
close all

% Noms des groupes fonctionnels, réorganisés
group_names = ["appendicularia" "gelatinous_carnivorous" "scyphozoa" "thaliacea" "ctenophora" "hydrozoa" "gelatinous_filter_feeders" "chaetognatha" "pteropod" "copepod" "large_crustracean"];

% Réorganiser les indices des groupes pour respecter l'ordre souhaité
new_order = [1, 2, 3, 4, 5, 8, 9, 7, 10, 6, 11];

name = repelem(group_names', 63);
%abondance_community_station(:,[1 4 6 8 9 10 11 17 18 19 26]
variable_per_group=diagnostic_global.mesozooplankton.respiration.respiration_carbon_m2;

variable_per_group = variable_per_group(:, new_order);


variable_per_group(variable_per_group == 0) = NaN;
variable_reshap = reshape(variable_per_group, [], 1);

% Multiplier le vecteur pour obtenir 11 colonnes
latitude = repmat(latWP2, 1, 10);
latitude = reshape(latitude, [], 1);

num_bins = 25;
edges = linspace(min(latitude), max(latitude), num_bins+1);
lat_categories = discretize(latitude, edges, 'IncludedEdge', 'right');
repname = repelem(1:10, 63)';

% Initialiser les matrices
group_contributions = zeros(num_bins, length(group_names));
total_per_bin = zeros(num_bins, 1);

% Calculer les valeurs totales et contributions des groupes
for i = 1:num_bins
    indices = (lat_categories == i);
    if any(indices)
        total_per_bin(i) = mean(variable_reshap(indices),'omitnan');
        for j=1:10
         indices_name = (repname == j);
        group_contributions(i,j) = mean(variable_reshap(indices & indices_name),'omitnan');
        end
    end
end

% Normaliser les contributions
total_per_bin(total_per_bin == 0) = NaN;
group_contributions = group_contributions ./ sum(group_contributions,2,'omitnan');

% Couleurs pastel adaptées aux catégories
color = [
    173 216 230; % Bleu pastel clair (Appendicularia)
    176 196 222; % Bleu gris pastel (Gelatinous carnivorous)
    135 206 250; % Bleu acier pastel (Scyphozoa)
    135 206 235; % Bleu ciel pastel (Thaliacea)
    175 238 238; % Turquoise pastel (Ctenophora)
    200 200 255; % Bleu lavande pastel (Hydrozoa)
    191 239 255; % Bleu clair pastel (Gelatinous filter feeders)
    255 239 160; % Jaune pastel (Chaetognatha)
    255 250 160; % Jaune clair pastel (Pteropod)
    255 160 160; % Rouge pastel (Copepod)
    255 182 193; % Rose pastel (Large crustacean)
];

% Convertir de 0-255 à 0-1 pour MATLAB
color = color / 255;
% Création du barplot empilé
figure; hold on;
h = bar(edges(1:end-1) + diff(edges)/2, group_contributions, 'stacked');
for i = 1:length(h)
    h(i).FaceColor = color(i, :); % Appliquer la couleur personnalisée
end

xlabel('Latitude');
legend(group_names, 'Location', 'bestoutside');
hold off;

%%%%%%%%%%%%%%%%%%%%%
%% Bar plot Global %%
%%%%%%%%%%%%%%%%%%%%%

close all

% Noms des groupes fonctionnels, réorganisés
group_names = ["appendicularia" "gelatinous_carnivorous" "scyphozoa" "thaliacea" "ctenophora" "hydrozoa" "gelatinous_filter_feeders" "chaetognatha" "pteropod" "copepod" "large_crustracean"];

% Réorganiser les indices des groupes pour respecter l'ordre souhaité
new_order = [1, 2, 3, 4, 5, 8, 9, 7, 10, 6, 11];

name = repelem(group_names', 63);
%abondance_community_station(:,[1 4 6 8 9 10 11 17 18 19 26]
variable_per_group=diagnostic_global.mesozooplankton.production.detritus_m2;

variable_per_group = variable_per_group(:, new_order);


% Vérifier la taille de la matrice
if size(variable_per_group,1) ~= 63 || size(variable_per_group,2) ~= 11
    error('La matrice doit être de taille 63x11.');
end

% Calcul des proportions par station (normalisation par ligne)
station_proportions = variable_per_group ./ sum(variable_per_group, 2);

% Calcul de la moyenne des proportions sur toutes les stations
mean_proportions = mean(station_proportions, 1);

% Noms des groupes fonctionnels
group_labels = strcat("G", string(1:11));

% Définition des couleurs personnalisées (converties en RGB)
% Couleurs pastel adaptées aux catégories
color = [
    173 216 230; % Bleu pastel clair (Appendicularia)
    176 196 222; % Bleu gris pastel (Gelatinous carnivorous)
    135 206 250; % Bleu acier pastel (Scyphozoa)
    135 206 235; % Bleu ciel pastel (Thaliacea)
    175 238 238; % Turquoise pastel (Ctenophora)
    200 200 255; % Bleu lavande pastel (Hydrozoa)
    191 239 255; % Bleu clair pastel (Gelatinous filter feeders)
    255 239 160; % Jaune pastel (Chaetognatha)
    255 250 160; % Jaune clair pastel (Pteropod)
    255 160 160; % Rouge pastel (Copepod)
    255 182 193; % Rose pastel (Large crustacean)
];


% Convertir de 0-255 à 0-1 pour MATLAB
color = color / 255;

% --- Barplot Empilé avec la moyenne des proportions ---
figure;
b = bar(1, mean_proportions, 'stacked'); % Barplot empilé sur une seule barre
for i = 1:length(b)
    b(i).FaceColor = color(i, :); % Appliquer la couleur personnalisée
end

% Ajouter les labels et titre
ylabel('Fraction de contribution moyenne');
ylim([0 1]); % Normalisation de l'axe Y entre 0 et 1
set(gca, 'XTick', []); % Supprimer les ticks sur l'axe X
set(gca, 'FontSize', 12);
legend(group_names, 'Location', 'bestoutside');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bar plot Global Polynome %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

% Noms des groupes fonctionnels, réorganisés
group_names = ["appendicularia" "gelatinous_carnivorous" "scyphozoa" "thaliacea" "ctenophora" "hydrozoa" "gelatinous_filter_feeders" "chaetognatha" "pteropod" "copepod" "large_crustracean"];

% Réorganiser les indices des groupes pour respecter l'ordre souhaité
new_order = [1, 2, 3, 4, 5, 8, 9, 7, 10, 6, 11];

name = repelem(group_names', 63);
%abondance_community_station(:,[1 4 6 8 9 10 11 17 18 19 26]
variable_per_group=diagnostic_global.mesozooplankton.feeding.prey_consumed_per_meso_carbon;

variable_per_group = variable_per_group(:, new_order);

degree = 4; % degré du polynôme
latitude_fit = linspace(min(latWP2), max(latWP2), 100); % points pour la courbe

% Initialisation
coefficients = zeros(degree + 1, length(group_names));
variable_fit_log = zeros(length(latitude_fit), length(group_names));

for i = 1:length(group_names)
    y = log10(variable_per_group(:, i));           % transformation log10
    y(isinf(y)) = NaN;               % remplace -Inf par NaN
    valid_idx = ~isnan(y);           % indices valides
    
    if sum(valid_idx) > 40
        coefficients(:, i) = polyfit(latWP2(valid_idx), y(valid_idx), degree); % ajustement
        variable_fit_log(:, i) = polyval(coefficients(:, i), latitude_fit);    % évaluation
    else
        coefficients(:, i) = NaN;                    % pas d'ajustement
        variable_fit_log(:, i) = NaN(size(latitude_fit)); % remplissage avec NaN
    end
end

group_names = ["appendicularia" "gelatinous_carnivorous" "scyphozoa" "thaliacea" "ctenophora" "hydrozoa" "gelatinous_filter_feeders" "chaetognatha" "pteropod" "copepod" "large_crustracean"];

% Calcul des proportions par station (normalisation par ligne)
station_proportions = 10.^(variable_fit_log) ./ sum(10.^(variable_fit_log), 2,'omitnan');

% Calcul de la moyenne des proportions sur toutes les stations
mean_proportions = mean(station_proportions, 1);
mean_proportions(isnan(mean_proportions)) = 0;


% Noms des groupes fonctionnels
group_labels = strcat("G", string(1:11));


% Définition des couleurs personnalisées (converties en RGB)
% Couleurs pastel adaptées aux catégories
color = [
    173 216 230; % Bleu pastel clair (Appendicularia)
    176 196 222; % Bleu gris pastel (Gelatinous carnivorous)
    135 206 250; % Bleu acier pastel (Scyphozoa)
    135 206 235; % Bleu ciel pastel (Thaliacea)
    175 238 238; % Turquoise pastel (Ctenophora)
    200 200 255; % Bleu lavande pastel (Hydrozoa)
    191 239 255; % Bleu clair pastel (Gelatinous filter feeders)
    255 239 160; % Jaune pastel (Chaetognatha)
    255 250 160; % Jaune clair pastel (Pteropod)
    255 160 160; % Rouge pastel (Copepod)
    255 182 193; % Rose pastel (Large crustacean)
];

% Convertir de 0-255 à 0-1 pour MATLAB
color = color / 255;

% --- Barplot Empilé avec la moyenne des proportions ---
figure;
b = bar(1, mean_proportions, 'stacked'); % Barplot empilé sur une seule barre
for i = 1:length(b)
    b(i).FaceColor = color(i, :); % Appliquer la couleur personnalisée
end

% Ajouter les labels et titre
ylabel('Fraction de contribution moyenne');
ylim([0 1]); % Normalisation de l'axe Y entre 0 et 1
set(gca, 'XTick', []); % Supprimer les ticks sur l'axe X
set(gca, 'FontSize', 12);
legend(group_names, 'Location', 'bestoutside');