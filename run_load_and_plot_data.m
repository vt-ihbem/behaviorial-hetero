function run_load_and_plot_data()
% Data on relative risk by age group for COVID-19 
% infection, hospitalization, death

clearvars; close all; clc

data = readtable('COVID-19_risk_by_age.xlsx');
data = removevars(data,'Ref_CDC_2023_RiskForCOVID_19Infection_Hospitalization_AndDeathByAgeGroup_AvailableFrom_Https___archive_cdc_gov_www_cdc_gov_coronavirus_2019_ncov_covid_data_investigations_discovery_hospitalization_death_by_age_html');
data.AgeGroup1 = data.AgeGroup;
for j = 1:length(data.AgeGroup1)
    data.AgeGroup1{j} = data.AgeGroup{j}(1:end-10);
end

f1 = figure(1);
plot(1:9,data.InfectionRiskByAge,'.-',3,data.InfectionRiskByAge(3),'*','linewidth',2,'markersize',20)
set(gca,'fontsize',16,'xtick',1:9,'xticklabel',data.AgeGroup1)
ylim([0 1.2])
xlim([0.5 9.5])
ylabel('Relative ratio compared to 18-29 age group')
xlabel('Age group')
title('Infection risk by age','fontsize',20)
xtickangle(90)
set(f1,'Position',[10 10 800 600])

f2 = figure;
plot(1:9,data.HospitalizationRiskByAge,'.-',3,data.InfectionRiskByAge(3),'*','linewidth',2,'markersize',20)
set(gca,'fontsize',16,'xtick',1:9,'xticklabel',data.AgeGroup1)
ylim([0 16])
xlim([0.5 9.5])
ylabel('Relative ratio compared to 18-29 age group')
xlabel('Age group')
title('Hospitalization risk by age','fontsize',20)
xtickangle(90)
set(f2,'Position',[10 10 800 600])

f3 = figure;
plot(1:9,data.DeathRiskByAge,'.-',3,data.InfectionRiskByAge(3),'*','linewidth',2,'markersize',20)
set(gca,'fontsize',16,'xtick',1:9,'xticklabel',data.AgeGroup1)
ylim([0 400])
xlim([0.5 9.5])
ylabel('Relative ratio compared to 18-29 age group')
xlabel('Age group')
title('Death risk by age','fontsize',20)
xtickangle(90)
set(f3,'Position',[10 10 800 600])

outdir = 'outputs_data';
if ~exist(outdir,'dir'), mkdir(outdir); end

saveas(f1,fullfile(outdir,'fig1_infectionrisk.png'))
savefig(f1,fullfile(outdir,'fig1_infectionrisk.fig'))
saveas(f2,fullfile(outdir,'fig1_hospitalizationrisk.png'))
savefig(f2,fullfile(outdir,'fig1_hospitalizationrisk.fig'))
saveas(f3,fullfile(outdir,'fig1_deathrisk.png'))
savefig(f3,fullfile(outdir,'fig1_deathrisk.fig'))

close all

fprintf('\nPlotting finished. Publication figures saved in %s\n', outdir);

