animals = {'Bones', 'Odo', 'Sulu', 'Worf', 'Beverly', 'Chakotay', 'Nerys',...
    'SevenOfNine', 'Rain', 'OBrien', 'Picard', 'Riker', 'Kes', 'Neelix', 'Quark',...
    'Sato', 'Garrett', 'Guinan', 'Keeler', 'Dax', 'TPol', 'Doctor', 'Tuvok', 'Bashir',...     
	'Scotty'};

animalbasedir = '\\hub.gladstone.internal\huanglab-lfp\Emily\DREADDs\Preprocessed Data';

for a = animals
    animaldir = sprintf('%s/%s/',animalbasedir,a{1});
    prefix = a{1};
    taskname = sprintf('%s/%stask.mat',animaldir,prefix);
    load(taskname);
    
    for d = 1:length(task)
        if(~contains(lower(task{d}{1}.descript), 'fail'))
            eaj_extractdspikes(animaldir, prefix, d)
        end
    end
end