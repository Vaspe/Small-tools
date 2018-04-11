function Data = logsout2data2( logsin )

    NumOutp = size(get(logsin),1);
    Data    = struct;
    Data.Channels = {};
    for i = 1:NumOutp
        try
        iSublog  = logsin.get(i);
        iTS      = iSublog.Values;
        iTSnames = fieldnames(iTS);
        Erfield  = iSublog.Name;        
        
        for ii =1: size(iTSnames,1)
            iData1   = timeseries2data(iTS.(iTSnames{ii}));
            inames  = fieldnames(iData1);
            Data.(inames{2}).(inames{1}) = iData1.(inames{1}) ;
            Data.(inames{2}).(inames{2}) = iData1.(inames{2}) ;
            
            if size(Data.Channels,1)==0
                Data.Channels{1} = inames{2};
            else
                Data.Channels{(end+1)} = inames{2};
            end
        end
        clearvars iSublog iTS iTSnames iTSnames iData1 inames

    catch ME
         ErField_int = get(logsin);
         Erfield     = ErField_int{i,1};
         disp (['Conversion of logsout to data for field ' Erfield  ' failed'])   
%          rethrow(ME)   
        end
    end

end


function Data  = timeseries2data(TimeSeriesIn)

time = TimeSeriesIn.Time;
name = TimeSeriesIn.Name;

Data.time = time;
eval(['Data.' name '= TimeSeriesIn.Data;']);


end