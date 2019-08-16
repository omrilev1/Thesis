function [finalH,sortIndices] = parse_alist(filepath,filetype,toConvertSys)
% This function parse a list files and generate the resulting parity check
% matrix. It check if the matrix has girth4, and it also convert the parity check matrix to systematic -
% thus allow to use Matlab's LDPC Encoder/Decoder.
% The alist file is build as follows :
% The The first line contains the size of the matrix as columns and rows.
% The second line contains the maximum number of ones per column followed by the maximum number of ones per row.
% Then follows a line that contains, for every column of the matrix, the number of ones in this column.
% The next line contains, for every row of the matrix, the number of zeros in that row.
% The following lines contain, for each column, the positions of the ones in that column. A 0 indicates an unused entry.
% The following lines redundantly contain for every row the positions of the non-zero entries in that row.

% The most comprehensive source for LDPC codes is mackay's site, which is
% called "Encyclopedia of sparse codes". It can be found here : 
% http://www.inference.org.uk/mackay/codes/data.html#l20

% if filetype == 'alist', then we parse the file according to the above
% format. if filetype == 'mat', then we just load the matrix, check the girth and transform to systematic form

% if toConvertSys = 1, then we try to change the parity matrix H such that
% it'll fit to matlab LDPC decoder
if strcmp(filetype,'alist')
    fid = fopen(filepath);
    currLine = 0;
    lineCNT = 1;
    matrixColCNT = 1;
    while (currLine ~= -1)
        currLine = fgetl(fid);
        currDouble = str2double(split(currLine));
        if currLine == -1
            disp('I finished the file');
            fclose(fid);
            break
        end
        
        switch lineCNT
            
            case 1
                N = currDouble(1);
                K = N - currDouble(2);
                H = sparse(N-K,N);
            case 2
                d_v = currDouble(1);
                d_c = currDouble(2);
            case 3
                varNodeDEG = zeros(1,length(currDouble));
                for i=1:length(currDouble)
                    varNodeDEG(i) = currDouble(i);
                end
            case 4
                chkNodeDEG = zeros(1,length(currDouble));
                for i=1:length(currDouble)
                    chkNodeDEG(i) = currDouble(i);
                end
            otherwise
                nanIdx = find(isnan(currDouble));
                
                if length(nanIdx) == 0
                    firstIdx = 1;
                else
                    firstIdx = nanIdx(1) + 1;
                    if (firstIdx == length(currDouble)+1)
                        firstIdx = 1;
                    end
                end
                for idx = firstIdx : firstIdx + varNodeDEG(matrixColCNT)-1
                    H(currDouble(idx),matrixColCNT) = 1;
                end
                matrixColCNT = matrixColCNT + 1;
                if matrixColCNT > N
                    fclose(fid);
                    break
                end
        end
        lineCNT = lineCNT + 1;
    end
else
    tempMat = load(filepath);
    H = tempMat.H;
    N = size(H,2);
    K = N - size(H,1);
end
%% Check the girth4 of the matrices
% girth = zeros(1,N-K);
% We can test girth 4 by using the matrix O
% Graph with cycle of 4 will have elemnts such as : 
%        - x -    ... - x - ...
% H =             ...
%        - x -    ... - x - ...
% and as so the rows inner prooduct will be greater than 2 - the matrix O
% is the wors inner product, and we check if it has elements greater than 2
% (we avoid the main diagonal which is simply the row degree)

% O=H*H';
% for i=1:(N-K)
%     O(i,i)=0;
% end
% for i=1:(N-K)
%     girth(i)=max(O(i,:));
% end
% girth4=max(girth);
% 
% if girth4<2
%     good_nonSys_H = H;
% else
%     disp('Fuck My Life');
%     good_nonSys_H = [];
% end

%% Convert H to systematic, so we'll be able to use Matlab's LDPC

if toConvertSys
    % Permute columns of a binary Matrix until the rightmost square matrix is
    % invertible over GF(2)
    
    % matrix dimensions:
    [~, n] = size(good_nonSys_H);
    
    % Initialization
    HInvertible = good_nonSys_H;
    PermuteIndex = 1:n;
    flag = true;
    counter = 0;
    
    % Initial Report
    disp('Creating a ParityCheck matrix which is suitable for MATLAB COMM Tollbox')
    
    % Permute columns
    while flag
        
        % Check if the rightmost square matrix is invertible over GF(2)
        try
            
            EncoderObject = comm.LDPCEncoder(sparse(HInvertible));
            % Check if new matrix works
            fprintf(['ParityCheck Matrix become suitable for Matlab LDPC Encoder ',...
                'after ',num2str(counter),' permutations!\n'])
            flag = false;           % Break the loop
            
        catch
            
            % Choose different columns for the rightmost part of matrix
            counter = counter+1;    %Permutation Counter
            PermuteIndex = randperm(n);
            HInvertible = H(:,PermuteIndex);
            
            if (mod(counter,100) == 0)
                disp(strcat('counter = ',num2str(counter)))
            end
            
        end
        
    end
    finalH = sparse(HInvertible);
    sortIndices = PermuteIndex;
else
    finalH = sparse(H);
    sortIndices = 1:size(H,2);
end
end

