function test_gather_user_options
%TEST_REGISTRY Tests the Registry class

% Call the Test class for this file
test = Test(mfilename('fullfile'));

% Run the code that is being tested
opt = test.func(@code_to_test);

% Evaluate the results
test.compare(opt{1}.option1,false, 'Standard input');
test.compare(opt{1}.option2,true, 'Flag input, multiple inputs');
test.compare(opt{2}.option2,true, 'Flag input, single input');

end

function opt = code_to_test

    var = {'-option2','option1',false,'unknownOption',52};
    options.option1 = true;
    options.option2 = false;
    opt{1} = gather_user_options(options,var{:}); 
    opt{2} = gather_user_options(options,'-option2');
    gather_user_options(opt,var{:},{'-disablewarn'});   
end