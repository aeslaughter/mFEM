function work

T1 = struct('name',{'3','2','4'},'type',{'subdomain','boundary','boundary'});
T2 = struct('name',{'1','2',},'type',{'subdomain','boundary'});

tic;
[Lia(:,1),~] = ismember({T1.name},{T2.name});
[Lia(:,2),Locb] = ismember({T1.type},{T2.type});


Lia = all(Lia,2);
idx = Locb(Lia)