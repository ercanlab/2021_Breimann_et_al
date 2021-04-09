function timestepsArray = timestepsArray(dim3_pre, timePre)
timeAr = [];
timePre1=timePre;
    for i = 1:dim3_pre
        timeAr = [timeAr, timePre];
        timePre = timePre + timePre1;
    end    
    timestepsArray = timeAr;
end
