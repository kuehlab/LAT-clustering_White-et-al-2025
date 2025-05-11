clearvars

load('array_stats_array_bkgnd.mat');
load('multiseg_params.mat')
%%
ratio_names = {'1-3', '1-6', '1-12', '1-24', '1-48', '1-96', '1-192', '0'};
time_names = {'5', '10',};
for r = 1:2
    for c = 1:8
        writetable(struct2table(data{r,c}),['csvs/r' ratio_names{c} '_' time_names{r} 'min.csv'])
    end
end
