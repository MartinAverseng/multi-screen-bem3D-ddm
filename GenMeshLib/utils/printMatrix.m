function [] = printMatrix(T,sep,openBrace,closeBrace,stream)

fprintf(stream,openBrace);
for i = 1:size(T,1)
    printArray(T(i,:),sep,openBrace,closeBrace,stream);
    if i<size(T,1)
        fprintf(stream,sep);
    end
end
fprintf(stream,closeBrace);

end

