function [] = printArray(T,sep,openBrace,closeBrace,stream)

fprintf(stream,openBrace);
for i = 1:length(T)
    fprintf(stream,num2str(T(i)));
    if i<length(T)
        fprintf(stream,sep);
    end
end
fprintf(stream,closeBrace);


end

