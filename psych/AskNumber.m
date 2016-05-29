function number = AskNumber(a, b, msg, scrInfo)
% number = AskNumber(a, b, msg, scrInfo)
% 
% Show a text screen with text MSG on screen SCRINFO.
% Wait for a response via keyboard, a single number from A to B.
% Return the response as an integer. Any non-numerical keypress is ignored.

    numberKeys = [KbName('0') KbName('1') KbName('2') KbName('3') KbName('4') KbName('5') KbName('6') KbName('7') KbName('8') KbName('9')];
    responseKeys = [(48+a):(48+b) numberKeys((1+a):(1+b))];
    [temp, tempAnswerInfo] = DisplayText(scrInfo, responseKeys, msg);
        
    % Decode the answer
    number = find(tempAnswerInfo.keyCode, 1);
    if number >= KbName('0') && number <= KbName('9')
        number = number - KbName('0');
    else
        number = number - 48;
    end
    
    clear temp tempAnswerInfo numberKeys responseKeys;
end