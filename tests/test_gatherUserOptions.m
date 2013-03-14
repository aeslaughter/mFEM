function T = test_gatherUserOptions
%TEST_GATHERUSEROPTIONS Tests the Registry class

% Call the Test class for this file
T = mFEM.Test(mfilename('fullfile'));

% Variable input
options.option1 = true;
options.option2 = false;
var = {'-option2','option1',false,'unknownOption',52};

% Evaluate the results
[opt{1}, unknown] = gatherUserOptions(options,var{:}); 
T.compare(opt{1}.option1,false, 'Standard input');

T.compare(unknown{2},52,'Unknown option collected');

warning('off','gatherUserOptions:UnknownProperty')
opt{1} = gatherUserOptions(options,var{:}); 
[~, msgid] = lastwarn;
T.compare(msgid,'gatherUserOptions:UnknownProperty', 'Unknown option warning');

opt{2} = gatherUserOptions(options,'-option2');
T.compare(opt{1}.option2,true, 'Flag input, multiple inputs');

gatherUserOptions(opt,var{:},'GatherUserOptions',{'-disablewarn'}); 
T.compare(opt{2}.option2,true, 'Flag input, single input, disable warning');



    
    
      
