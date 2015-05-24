% plot frequency of on/off life time
ton_scl=[0.05:0.05:10];
%toff_scl=[0.001:0.05:300];
nton=zeros(1,length(ton_scl));%register number of different on-lifetime
%ntoff=zeros(1,length(toff_scl));%register number of different off-lifetime
[row, col]=size(rec_events);
for i=1:col/6 %extracted value
    for j=2:row
        if(rec_events(j,6*i-3)==0)
            break;
        end
    end
    lngth=j-1; %extracted value
    nton=nton+hist(rec_events(2:lngth,6*i-3),ton_scl);%on-lifetime
%    ntoff=ntoff+hist(rec_events(3:lngth,6*i-5)-rec_events(2:(lngth-1),6*i-4),toff_scl);%off-lifetime
end
%lnnton=log(nton);lnton_scl=log(ton_scl);
%lnntoff=log(ntoff);lntoff_scl=log(toff_scl);
save ton_dist nton ton_scl
figure
bar(ton_scl,nton);title('Distribution of ¦Óon (647nm-50%, 405nm-0%)');
xlabel('second');ylabel('frequency');
%figure
%bar(ton_scl,lnnton);
%figure
%bar(lnton_scl,lnnton);
%figure
%bar(toff_scl,ntoff);title('Distribution of ¦Óoff (647nm-50%, 405nm-0%)');
%xlabel('second');ylabel('frequency');
%figure
%bar(toff_scl,lnntoff);
%figure
%bar(lntoff_scl,lnntoff);
%******************
%not useful, the time resolution is ~50ms, can not resolve lifetime
%the distribution is not fitted only by power-law
%some wried fitting even can do better than power-law