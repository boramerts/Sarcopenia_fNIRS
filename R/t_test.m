clear; clc;
data = readtable('Ttest.xlsx');

control = table();
sarco = table();

c = 1; s = 1;
for i = 1:length(data.Yas)
    if data.Sarkopenik(i) == 0
        control(c,:) = data(i,:);
        c = c+1;
    else
        sarco(s,:) = data(i,:);
        s = s+1;
    end
end

control.Sarkopenik = [];
sarco.Sarkopenik = [];

[h,p,ci,stats] = ttest2(table2array(control),table2array(sarco))