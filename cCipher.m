function [output] = cCipher(input,shift,endec)
%CCIPHER Encrypts/decrypts a string using a Caesar cipher.
%   The functions takes an input string and an input shift and applies the
%   latter to the former to produce an encrypted/decrypted output string.
%   Non-letter characters are ignored. Endec value must be 1 to code, 0 to
%   decode.
%
%	Jacopo Agagliate
%	University of Strathclyde
%	29 April 2016

    alphabet = 'abcdefghijklmnopqrstuvwxyz';
    nuInput = lower(expungeText(input));
    output = nuInput;

    if rem(shift,1) ~= 0
        error('Error. Shift input must be an integer.')
    end

    for a = 1:1:size(nuInput,2)
        pos = find(alphabet == nuInput(a));
        if endec == 1
            nupos = rem(pos+shift,size(alphabet,2));
        elseif endec == 0
            nupos = rem(pos-shift,size(alphabet,2));
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

