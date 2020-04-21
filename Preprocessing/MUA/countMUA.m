function countMUA (animaldir, prefix, day, epoch)

figure
spikes = zeros(32,1);
for i = 1:32
    load(sprintf('%s/%smua%02d-%d-%02d.mat',animaldir,prefix,day,epoch,i))
    spikes(i) = length(mua{day}{epoch}{i}.spiketimes);
end

plot(spikes)
xlim([1 32])
xlabel('Channel')
ylabel('Total Spikes')
title(sprintf('%s MUA %02d-%02d',prefix,day,epoch));

end