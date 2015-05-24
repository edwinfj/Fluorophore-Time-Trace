%molecule time trace analysis
%antibody-alexa647/alexa488
%raw data: recording time- time,
%          background(counts/px)- bkg, 
%          molecule time traces(counts/25px)- raw_trace
%exposure time:50ms, EMGain 70, Readout Rate 10MHz,recording time:1500s
% all manul inputs are labeled as "-parameter"

mol_num=68.0;%molecule number -parameter
area=25;%area(pixels) selected per molecule -parameter
exp_time=0.03;%exposure time(s) -parameter
ppc=0.1932;%photons/count  -parameter
off_thr=2;%set the threshold(on_thr*std) of a switch-off event -parameter
delta_thr=2.7;%set the threshold(on_thr*std) of a off state -parameter
on_thr=3.5;%set the threshold(on_thr*std) of a on state -parameter
stt_fr=1301;%set the first frame for analysis -parameter
%%
end_fr=length(time);%get the total frame number
fr_num=end_fr-stt_fr+1;%set the frame number for analysis

trace_N=zeros(fr_num,mol_num);%molecule time traces(photons/area)
for i=1:mol_num
    trace_N(:,i)=(raw_trace(stt_fr:end_fr,i*2)-bkg(stt_fr:end_fr)*area)*ppc;% -parameter
end
std_N=std(trace_N,1);%std(photons/area) of each molecule's time trace
%%
%trace_deltaN=zeros(fr_num,mol_num);%time trace of photon change
%for i=1:fr_num
%    if ~(i-1)
%        trace_deltaN(i,:)=traceN(i,:);
%    else
%        trace_deltaN(i,:)=traceN(i,:)-traceN(i-1,:);
%    end
%end
%only use photon number to judge on or off is a safer method


%%
%analyze switching events
rec_events=zeros(201,6*mol_num);%hope the 200 switching-events/molecule is enough -parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rec_events records switching events for each molecule
%for each molecule:
%the first row is the number of switching events
%since the second row, the columns are:
%switch-on time t_on, switch-off time t_off, on-state duration ¦Ó
%photons from all on-frames, on-frame number, photons/swiching event
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:mol_num %traverse through all molecules
    switch_num=0;%for each molecule, initialize switching number
    event=0;%flag of finding switching event
    on_fr=0;%number of on-frames during each switching event
    tt_ph=0.0;%total photons from all frames of a switching event
    temp1=6*(j-1);
    for i=1:fr_num      
        if ~event %haven't enterred a switching event
            if (i~=fr_num)&&(trace_N(i,j)>=std_N(j)*on_thr)
                %switch-on threshould and not last frame
                %enter a switching event
                switch_num=switch_num+1;
                event=1;
                rec_events(switch_num+1,1+temp1)=time(stt_fr+i-1);%record t_on(/s)
                on_fr=1;%initialize on-frame number
                tt_ph=trace_N(i,j);%initialize total photons
            end
        else %in a switching event
            if trace_N(i,j)<std_N(j)*off_thr || (i==fr_num) || trace_N(i-1,j)-trace_N(i,j)>=std_N(j)*delta_thr
                %exit from a switching event
                %or last frame, discard the last frame
                event=0;
                rec_events(switch_num+1,2+temp1)=time(stt_fr+i-1);%record t_off(/s)
                rec_events(switch_num+1,3+temp1)=time(stt_fr+i-1)-rec_events(switch_num+1,1+temp1);%¦Ó(s) of this switching event
                temp2=rec_events(switch_num+1,3+temp1);
                rec_events(switch_num+1,4+temp1)=tt_ph;%record on_frame photons of this switching event
                rec_events(switch_num+1,5+temp1)=on_fr;%record on-frame number
                rec_events(switch_num+1,6+temp1)=tt_ph/exp_time*temp2/on_fr;%record photons/switching-event
                on_fr=0;%empty on-frame number for next switching event
                tt_ph=0;%empty total photons for next switching event
            else
                on_fr=on_fr+1;
                tt_ph=tt_ph+trace_N(i,j);
            end
        end
    end
    rec_events(1,1+temp1)=switch_num;
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%give the mean photons/switching-cycle, 
%mean photons/frame during switching,
%mean switching-cycles/molecule,
%and histograms of photons/switching-cycle
%   *****and photons/frame******
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tt_switch_events=sum(rec_events(1,:));
%total switching cycles of all molecules
tt_switch_photons=0.0;
%total photons during switching cycles of all molecules
tt_on_fr=0.0;
%total on-frames of all molecules
tt_fr_photons=0.0;
%total photons during on-frames of all molecules
tt_on_t=0.0;
%total switching cycle duration of all molecules
positive_mol_num=0;%molcules that have switching events
hist_switch_photons=zeros(1,20);%20bins -parameter
%hist_fr_photons=zeros(1,20);%20 bins -parameter
bins1=500:4000:76500;%bin vector -parameter
%bins2=250:500:9750;
for i=1:mol_num
    tt_switch_photons=tt_switch_photons+sum(rec_events(:,6*i));
    tt_on_fr=tt_on_fr+sum(rec_events(:,6*i-1));
    tt_fr_photons=tt_fr_photons+sum(rec_events(:,6*i-2));
    tt_on_t=tt_on_t+sum(rec_events(:,6*i-3));
    index=rec_events(1,6*i-5);%switching events
    if index
        positive_mol_num=positive_mol_num+1;
        hist_switch_photons=hist_switch_photons+hist(rec_events(2:(index+1),i*6),bins1);
        %hist_fr_photons=hist_fr_photons+hist(rec_events(2:(index+1),i*6-2)./rec_events(2:(index+1),i*6-1),bins2);
    end
end
figure                        %****** draw histogram of photons/switching
bar(bins1,hist_switch_photons);%******
title('Histogram of Photons/Switch-Cycle');
xlabel('photons/switching-cycle'),ylabel('Frequency');
axis tight;
figure                        %****** draw histogram of photons/frame
%bar(bins2,hist_fr_photons);%******
%title('Histogram of Photons/Frame During Switch Cycles');
%xlabel('photons/frame'),ylabel('Frequency');
%axis tight;

avg_fr_photons=tt_fr_photons/tt_on_fr;
%mean photons/frame for this fluorophore under current condition
avg_switch_photons=tt_switch_photons/tt_switch_events;
%mean photons/switching-cycle for this fluorophore under current condition
avg_on_t=tt_on_t/tt_switch_events;
%mean duration(s) of one switching cycle
avg_switch_cycle=tt_switch_events/positive_mol_num;
%average switching cycles/molecule during this recording time
%under current condition for this fluorophore


%%
%begin calculating duty cycle and survival fraction
periods=15;% -parameter
matrix1=zeros(periods+1,mol_num);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate duty cycle and survival fraction for each 100s period
%first row is the last period for each molecule
%second row is the total switch-on time during 0-100s,
%third row is the total switch-on time during 100-200s, etc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:mol_num
    switch_num=rec_events(1,6*j-5);%get switching cycles for this molecule
    if switch_num
        matrix1(1,j)=floor(rec_events(switch_num+1,6*j-5)/100)+1;
        %last period
        for i=1:switch_num
            index=floor(rec_events(i+1,6*j-5)/100)+2;
            matrix1(index,j)=matrix1(index,j)+rec_events(i+1,6*j-3);
        end
    else
        matrix1(1,j)=0;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matrix2=zeros(periods,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           duty-cycle  |    survival-fraction
%0~100s         ?                    ?
%100~200s       ?                    ?
%...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:periods
    matrix2(i,1)=sum(matrix1(i+1,:));
end
%total switching-on time during each period
for j=1:mol_num
    index=matrix1(1,j);
    if index
        matrix2(1:index,2)=matrix2(1:index,2)+1;
    end
end
%survival molecule number during each period
matrix2(2:periods,1)=matrix2(2:periods,1)./matrix2(2:periods,2)/100.0;
matrix2(1,1)=matrix2(1,1)/matrix2(1,2)/(100.0-time(stt_fr));
%mean duty cycle per 100s per molecule
matrix2(:,2)=matrix2(:,2)/positive_mol_num;
%survival fraction per 100s

figure   %********** draw graph of duty cycle and survival fraction
last_p=7; % last period to show -parameter
x=1:last_p;
[AX,H1,H2]=plotyy(x,matrix2(1:last_p,1),x,matrix2(1:last_p,2));
title('Mean Duty Cycle And Survival Fraction Per 100s');
xlabel('time/100s')
set(get(AX(1),'Ylabel'),'String','Duty Cycle');
set(get(AX(2),'Ylabel'),'String','Survival Fraction');



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%save data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear i j temp1 temp2 on_fr event switch_num index tt_ph AX x H1 H2
save Results % rename file -parameter
        