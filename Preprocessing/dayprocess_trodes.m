% Extract EEG, ripple-filtered trace, ripples, dentate spikes, MUA, 
% position tracking, and linear track performance from raw data

animals = {'Adobo_LT', 'Anise_LT','Baharat_LT', 'Basil_LT', 'Caraway_LT',...
    'Cardamom_LT', 'Chicory_LT', 'Chives_LT', 'Cinnamon_LT', 'Coriander_LT',...
    'Cumin_LT', 'Dill_LT', 'Fenugreek_LT', 'GaramMasala_LT', 'Ginger_LT',...
    'Jerk_LT', 'Mace_LT', 'Mustard_LT','Nutmeg_LT', 'OldBay_LT', 'Oregano_LT',...
    'Paprika_LT', 'Parsley_LT', 'Pepper_LT', 'Provence_LT', 'Pumpkin_LT',...
    'Rosemary_LT', 'Saffron_LT', 'Sage_LT', 'Salt_LT', 'Sumac_LT', 'Tarragon_LT', 'Thyme_LT', 'Vanilla_LT'};
multtracks = 1;
track = 'LT';
trackepoch = 1;
datadir = 'E:\CRCNS\Cohort 2\';
animalbasedir = 'E:\CRCNS\Cohort 2\';


for a = 1:length(animals)
    prefix = animals{a};
    animaldir = [animalbasedir '/' prefix];
    
    taskname = sprintf('%s/%stask.mat',animaldir,prefix);
    load(taskname);
    for d = 1:length(task)
        for e = 1:length(task{d})
            if(~contains(lower(task{d}{e}.descript), '*fail*') && ...
                    ~contains(lower(task{d}{e}.descript), '*not run*'))
                extract_eeg(datadir, animaldir, prefix, d, e, 'trackepoch', trackepoch, 'multtracks', multtracks);
                ej_rereference(animaldir, prefix, d, e);
                ej_rippledayprocess(animaldir, prefix, d, e);
                
                ag_extractdspikes(animaldir, prefix, d, .015, 5)
            	plotdspikes
                
                extract_spikes(datadir, animaldir, prefix, d, e, 'multtracks', 1, 'trackepoch', 1);
                countMUA(animaldir, prefix, d, e);
            end
        end
        ag_extractripples(animaldir, prefix, d, .015, 3);
        extract_pos(datadir, animaldir, prefix, d, 'trackepoch', trackepoch, 'track', track, 'multtracks', multtracks)
    end
end

if strcmpi(track, 'LT')
    calcLinearTrackPerformance(datadir, animaldir, prefix, d, 'trackepoch', trackepoch, 'multtracks', multtracks);
end