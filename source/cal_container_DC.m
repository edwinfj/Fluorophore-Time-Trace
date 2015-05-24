%%%%%%%%%%%%%
% You FORGOT TO DISCARD BLEACHED MOLECULES!!!!!!
% You FORGOT TO DISCARD BLEACHED MOLECULES!!!!!!
% You FORGOT TO DISCARD BLEACHED MOLECULES!!!!!!
% You FORGOT TO DISCARD BLEACHED MOLECULES!!!!!!
%%%%%%%%%%%%




%%
% container of duty-cycles of individual molecules
% all the source data are saved in rec_events
% results are saved in container_DC
% 7 columns: 0~100s, 100~200s, etc.
% rows: duty-cycles of individual molecules in the period
% ???????molecule????????100s??
[row,col]=size(rec_events);
mol_num=col/6;
container_DC=zeros(mol_num,7); % ???????
% ???????????100s??on-state duration??????????????column??
for i=1:mol_num %???????
    for j=2:row %????????????
        on_t=rec_events(j,i*6-5); %emission-burst????
        if (on_t>0 && on_t<700)
            container_DC(i,ceil(on_t/100))=container_DC(i,ceil(on_t/100))+rec_events(j,i*6-3); %???on-state duration????????????
        else
             break; %??emission-burst???????????
        end
    end
end
container_DC=container_DC/100;
sum_container_DC=sum(container_DC(:,4:7),2)/4;

save container_DC container_DC sum_container_DC

%% calculate distribution of duty cycle in each 100s period
% save data in dc_dist; each rwo represents a 100s period

dc_scl=[0.0001:0.003:0.35];% -input
dc_dist=zeros(7,length(dc_scl));
for i=1:7
    dc_dist(i,:)=hist(container_DC(:,i),dc_scl);
end
%%
for i=1:1% -input
    dc_fit(i,:)=1./dc_fit(i,:);
end
save dc_dist dc_scl dc_dist dc_fit