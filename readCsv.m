function data = readCsv(file, varargin)
% FORMAT data = readCsv(file, varargin)
% Reads excel files that have been stripped of commas. Discards the top row
% as suspected headers.
%
% Inputs:
% file:			CSV file to read (e.g. file.csv). String.
% varargin{1}:	A cell array of excel columns to use (e.g. {'A' 'CX'}).
% varargin{2}:	Corresponding names for the loaded columns (e.g. {'column A'
%				'column B'}). Cell array.
% Outputs:
% data:			A structure with header and column components. 
%
% Additional Notes: This function works well with the writeCsv.m function.
% The writeCsv.m function writes out .csv files based upon the
% information in the data structure.
%
%
% readCsv.m is the proprietary property of The Regents of the
% University of California (“The Regents.”)
% Copyright © 2012-14 The Regents of the University of California, Davis
% campus. All Rights Reserved.  
% Redistribution and use in source and binary forms, with or without
% modification, are permitted by nonprofit, research institutions for
% research use only, provided that the following conditions are met:
%   - Redistributions of source code must retain the above copyright
%     notice, this list of conditions and the following disclaimer. 
%   - Redistributions in binary form must reproduce the above copyright
%     notice, this list of conditions and the following disclaimer in the
%     documentation and/or other materials provided with the distribution. 
%   - The name of The Regents may not be used to endorse or promote
%     products derived from this software without specific prior written
%     permission. 
% The end-user understands that the program was developed for research
% purposes and is advised not to rely exclusively on the program for any
% reason.
% THE SOFTWARE PROVIDED IS ON AN "AS IS" BASIS, AND THE REGENTS HAVE NO 
% OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
% MODIFICATIONS. THE REGENTS SPECIFICALLY DISCLAIM ANY EXPRESS OR IMPLIED
% WARRANTIES, INCLUDING BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
% SPECIAL, INCIDENTAL, EXEMPLARY OR CONSEQUENTIAL DAMAGES, INCLUDING BUT
% NOT LIMITED TO  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES, LOSS OF USE,
% DATA OR PROFITS, OR BUSINESS INTERRUPTION, HOWEVER CAUSED AND UNDER ANY
% THEORY OF LIABILITY WHETHER IN CONTRACT, STRICT LIABILITY OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
% THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF ADVISED OF THE POSSIBILITY
% OF SUCH DAMAGE. 
% If you do not agree to these terms, do not download or use the software.
% This license may be modified only in a writing signed by authorized
% signatory of both parties.
% For commercial license information please contact copyright@ucdavis.edu.
%
%
% Created by Benjamin Geib 141010

data = {};
if exist(file, 'file')
	fid = fopen(file); show = 0;
else
	fprintf(['ERROR: ' file ' DNE\n']);
	error('Error: Correct pathways\n');
end

% Read in the data
count = 1;
while 1
    line = fgetl(fid);
    if ~ischar(line)
        break
    end
    split{count} = strread(line, '%s', 'delimiter', ',');
    count = count + 1;
end

file_data = {};
header = split{1};
c = 1;

% Allow for dication of columns
if length(varargin) == 2
    file_data = let2num(varargin{1});
else
    file_data = 1:size(header, 1);
end

% Organize and save as a structure
for n = file_data
    % if headers are defined
    if length(varargin) == 2
        data{c}.header = varargin{2}{c}; 
        % First is header, so start with 2nd structure
        for i = 2:length(split)
            if ~isempty(split{i})
                try
                    data{c}.col(i - 1) = split{i}(n);
                catch err
                    data{c}.col(i - 1) = {''};
                end
            end
        end
        c = c + 1;
    % if headers aren't defined
    else
        data{n}.header = header(n);
       % First is header, so start with 2nd structure
        for i = 2:length(split)
            try
                data{n}.col(i - 1) = split{i}(n);
            catch err
                 data{c}.col(i - 1) = {''};
            end
                
        end
        c = c + 1;
    end
 
end

% Print out the first 5 examples to make sure looks okay
if show == 1
    for i = 1:length(data)
        fprintf([char(data{i}.header) '\t\t' ]);
    end
    fprintf('\n');
    for i = 1:5
        for j = 1:length(data)
            fprintf([char(data{j}.col(i)) '\t\t' ]);
        end
        fprintf('\n');
    end
end

fclose(fid);

end

%=========================================================================%
% Internal Function: let2num
%=========================================================================%
function num_array = let2num(array)
% FORMAT num_array = let2num(array)
% Function converts letter columns to numbers, e.g. AA = 27

for i = 1:length(array)
    sum = 0;
    N = length(array{i});
    for j = 1:length(array{i})
        switch array{i}(j)
            case 'A', L = 1;
            case 'B', L = 2;
            case 'C', L = 3;
            case 'D', L = 4;
            case 'E', L = 5;
            case 'F', L = 6;
            case 'G', L = 7;
            case 'H', L = 8;
            case 'I', L = 9;
            case 'J', L = 10;
            case 'K', L = 11;
            case 'L', L = 12;
            case 'M', L = 13;
            case 'N', L = 14;
            case 'O', L = 15;
            case 'P', L = 16;
            case 'Q', L = 17;
            case 'R', L = 18;
            case 'S', L = 19;
            case 'T', L = 20;
            case 'U', L = 21;
            case 'V', L = 22;
            case 'W', L = 23;
            case 'X', L = 24;
            case 'Y', L = 25;
            case 'Z', L = 26;
        end
        % Add the base 26
        if j == N
            sum = sum + L;
        else
            sum = sum + (L * 26 * (N - j));
        end
    end
    num_array(i) = sum;
    clear sum N;
end
end
