function [EmWaveLength, ExcWaveLength, D] = getImage3DFromFichier(chemin)

fid = fopen(chemin,'r');
trouve=0;
k = 1;
i=1;
j=1;

while feof(fid) == 0        %Test for end-of-file
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    matches = findstr(tline, 'DATA'); %Find string within another, longer string
    num = length(matches);
    if num > 0
        trouve=1;
    end
    
    if (trouve==1) && (num==0)
        k = k+1;
        if (strcmp(tline(1:2),'Ex'))
            ind      = strfind(tline,' ');
            ExcWaveLength(i) = str2double(tline(ind(2)+1:length(tline)));
            i = i+1;
        else
            ind       = strfind(tline,' ') ; % la recherche des indice de ' '
            
            EmWaveLength(j)   = str2double(tline(1:ind-1));
            D(j)      = str2double(tline(ind+1:length(tline)));
            j=j+1;                 
        end
    end
end

Nex = length(ExcWaveLength) ; 
Nem = length(EmWaveLength)/Nex;
EmWaveLength = EmWaveLength(1:Nem);
D = reshape(D,Nem,Nex);
D = D';

fclose(fid);