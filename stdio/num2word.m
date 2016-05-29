%NUM2WORD Convert an integer N to an English word for that cardinal number. 
% If the optional parameter O is one, returns the ordinal number.

function s = num2word(n, o)
    cardword = {'zero', 'one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine', 'ten'};
    ordword = {'zeroth', 'first', 'second', 'third', 'fourth', 'fifth', 'sixth', 'seventh', 'eighth', 'ninth', 'tenth'};

    if ~exist('o', 'var'); o = []; end
    if isempty(o); o = 0; end
    
    n = n + 1;
    
    if n < 1 || n > length(cardword)
        s = 'out of range';
        return;
    end
    
    % Return either an ordinal or a cardinal
    if o
        s = ordword{n};
    else
        s = cardword{n};
    end
    
end