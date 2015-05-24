%%
%plot the time trace of ith molecule
ith_mol=2;%-parameter
thsh=zeros(1,fr_num);
flag=0;
std_val=std_N(ith_mol);
for i=1:fr_num
    if flag
        if trace_N(i,ith_mol)<off_thr*std_val || (trace_N(i-1,ith_mol)-trace_N(i,ith_mol)>=delta_thr*std_val) || i==fr_num
            flag=0;
        else
            thsh(i)=on_thr*std_val;
        end
    else
        if (i~=fr_num)&&(trace_N(i,ith_mol)>=std_val*on_thr)
        flag=1;
        thsh(i)=on_thr*std_val;
        end
    end
end

figure
plot(time(stt_fr:end_fr),trace_N(:,ith_mol),time(stt_fr:end_fr),thsh);
xlabel('time/s'),ylabel('photons'),title('One Time Trace of a Alexa488 Molecule Under 100% 488nm and 10% 405nm Excitation');
clear thsh ith_mol flag std_val i
%%
%test cell
ith_mol=1;
hold on;
thsh=zeros(fr_num,1)+2*std_N(ith_mol);
plot(time,thsh);hold off;
clear ith_mol;