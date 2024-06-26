function p_myfunc_toDat(timeseries, filename)

    x1 = timeseries(:,1);
    x2 = timeseries(:,2);

    matrix = [x1;x2].';
    fid = fopen(filename', 'w');
    optionString = "%f\n";
    for i = 1:length(x1)-1
        optionString = '%f ' + optionString;
    end
    fprintf(fid, optionString, matrix.');
    fclose(fid);