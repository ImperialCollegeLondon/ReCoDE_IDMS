## Testing 

### What is testing? 

Testing is the process of executing code in order to find errors and debug them. You probably already do testing, informally, to check that your code or functions are working as expected. However, it is less likely that you write functions and then test them formally on several different input and expected output values to ensure they are robust. This is formal testing, and is useful because a function that works initially with a specific set of input values may not work when even a small input value changes. 

In this project, the functions are mostly standalone, however when developing more complex code with functions that rely on the output of other functions, *unit testing* become increasingly useful to debug code. Unit testing checks small parts of the code one at a time to check they are behaving as expected. If they are not, then it is easier to debug why than when an error appears as part of a large chain of dependent functions.

In order to facilitate testing, this project splits the code into short functions with a clear and singular purpose. Therefore, by writing code that is easy to test, we improve the quality of our code, by writing more manageable functions. 


### Test that 

In order to run our unit testing, this project uses the [testthat package](https://r-pkgs.org/testing-basics.html), which tries to make testing easy and, possibly, fun. Each test takes a known input (i.e., arguments provided to a function) and checks that the output is expected (i.e., a pre-defined value, object etc. returned by the function). Once you have written all the tests for all the functions in your project, testthat will run them all and return an output which states how many tests pass, fail or throw warnings. As you build up your code and continue to test regularly, if a function changes or a new input throws an error then testthat will flag it and it will be more straightforward to debug! 

### Testing functions

In the R file *test-functions.R* there is code to test 4 functions: *simulate_data_single_var*, *simulate_data_multi_var*, *draw_init_values* and *compare_param_est*. 

For instance, the first two functions are used to simulate infectious disease transmission data by solving a set of ordinary differential equations (ODE).One of the tests we want to run is to check that the outputs are always positive, as we can't have negative numbers of people! If the outputs are negative then this could mean that our ODE equations aren't balanced. Fortunately, test.that will warn us quickly so we can address the problem before it has a knock on affect to downstream analysis. 

As well as testing regularly, the act of writing the tests helps us to improve the functions. For instance, when writing the tests for *simulate_data_single_var*, I realised that if any of the rate parameters were negative, this could also lead to the outputs being negative. Therefore I added checks at the start of the function, to ensure that all parameter values are valid. For instance: *if(n_pop <0 ) stop("n_pop cannot be negative")*. 

Once you have written tests for all the functions, the testthat package makes it easy to run all the tests in one go: simply execute *test_dir(".")*, in the console below. As long as the code to test is saved in the current working directory and has the prefix test, the testthat package will find and run all your tests. 

**Challenge: Only some of the functions here have test written for them, in order to demonstrate the testthat package. However, you would ideally write tests for all functions. Try writing some tests for additional functions (locations in the *R* folder). You can either do this by adding to the script *test-functions.R*, or write a new script with the prefix test in its name. Once you have written the tests, you can execute them as before using *test_dir(".")*. You may find that you want to edit the functions you are testing too, to ensure they are robust!**



 


