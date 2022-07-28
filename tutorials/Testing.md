## A note on Testing 

### What is testing? 

Testing is the process of executing code in order to find errors and debug them. You probably already do testing, informally, to check that your code or functions are working as expected. However, it is less likely that you write functions and then test them formally on several different input and expected output values to ensure they are robust. This is formal testing, and is useful because a function that works initially with a specific set of input values may not work when even a small input value changes. 

In this project, the functions are mostly standalone, however when developing more complex code with functions that rely on the output of other functions, *unit testing* become increasingly useful to debug code. Unit testing checks small parts of the code one at a time to check they are behaving as expected. If they are not, then it is easier to debug why than when an error appears as part of a large chain of dependent functions.

In order to facilitate testing, this project splits the code into short functions with a clear and singular purpose. Therefore, by writing code that is easy to test, we improve the quality of our code, by writing more manageable functions. 


### Test that 

In order to run our unit testing, this project uses the [test that package](https://r-pkgs.org/testing-basics.html), which tries to make testing easy and, possibly, fun. Each test takes a known input (i.e., arguments provided to a function) and checks that the output is expected (i.e., a pre-defined value, object etc. returned by the function). Once you have written all the tests for all the functions in your project, test that will run them all and return an output which states how many tests pass, fail or throw warnings. As you build up your code and continue to test regularly, if a function changes or a new input throws an error then test that will flag it and it will be more straightforward to debug! 




