function [output] = expungeText(input)
%EXPUNGETEXT Removes white space and punctuation from a text.
%   The functions takes an input string and produces an output string
%   devoid of punctuation and white space. The text is left ideally suited
%   for cipher encoding.
%
%	Jacopo Agagliate
%	University of Strathclyde
%	29 April 2016

    alphabet = 'abcdefghijklmnopqrstuvwxyz';
    c_alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
    output = input;
    expIndex = [];

    for a = 1:1:size(input,2)
        pos = find(alphabet == input(a),1);
        if isempty(pos)
            pos = find(c_alphabet == input(a),1);
            if isempty(pos)
                expIndex = [expIndex;a];
            end
        end
    end
    output(expIndex) = [];
end

