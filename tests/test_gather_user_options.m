function T = test_gather_user_options
%TEST_GATHER_USER_OPTIONS Tests the Registry class

% Call the Test class for this file
T = mFEM.Test(mfilename('fullfile'));

% Variable input
options.option1 = true;
options.option2 = false;
var = {'-option2','option1',false,'unknownOption',52};

% Evaluate the results
warning('off','gather_user_options:unknown')
[opt{1}, unknown] = gather_user_options(options,var{:}); 
T.compare(opt{1}.option1,false, 'Standard input');

[~, msgid] = lastwarn;
T.compare(msgid,'gather_user_options:unknown', 'Unknown option warning');

T.compare(unknown{2},52,'Unknown option collected');

opt{2} = gather_user_options(options,'-option2');
T.compare(opt{1}.option2,true, 'Flag input, multiple inputs');

gather_user_options(opt,var{:},{'-disablewarn'}); 
T.compare(opt{2}.option2,true, 'Flag input, single input');



    
    
      
