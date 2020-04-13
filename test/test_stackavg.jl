using Dates

geolist = [Date(2019, 1, 1), Date(2019, 1, 3), Date(2019, 1, 6), Date(2019, 1, 10)]
intlist = [(geolist[i], geolist[j]) for i = 1:length(geolist)-1 for j = i+1:length(geolist)]
