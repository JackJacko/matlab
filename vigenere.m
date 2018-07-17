function [output] = vigenere(input,key,endec)
%VIGENERE Encrypts/decrypts a string using a Vigenère cipher.
%   The functions takes an input string and an input key and applies the
%   latter to the former to produce an encrypted/decrypted output string.
%   Endec value must be 1 to code, 0 to decode.
%
%	Jacopo Agagliate
%	University of Strathclyde
%	29 April 2016

    alphabet = 'abcdefghijklmnopqrstuvwxyz';
    nuInput = lower(expungeText(input));
    output = nuInput;

    keylength = size(key,2);
    plainlength = size(nuInput,2);
    reps = mod(plainlength,keylength)+1;

    extKey = repmat(key,1,reps);
    extKey(plainlength+1:end) = [];
    extKey = lower(extKey);

    for a = 1:1:size(nuInput,2)
        posA = find(alphabet == nuInput(a));
        posB = find(alphabet == extKey(a))-1;
        if endec == 1
            nupos = rem(posA+posB,size(alphabet,2));
        elseif endec == 0
            nupos = rem(posA-posB,size(alphabet,2));
        else
            error('Error. Endec input must be 1 (code) or 0 (decode).')
        end
        if nupos == 0
            nupos = size(alphabet,2);
        elseif nupos < 0
            nupos = size(alphabet,2) + nupos;
        end
        output(a) = alphabet(nupos);
    end

end

