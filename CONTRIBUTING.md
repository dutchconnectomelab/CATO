Contributing Guidelines
===========================

Thank you very much for your interest in contributing to CATO!

Questions about CATO
---------------------------
If you have questions you can:

* If you encounter a problem or have a feature suggestion, then you can open an issue on the issues page of this GitHub repository. This page contains a list with the outstanding issues (bugs and feature requests) that are in development.

* If you have a very specific question that you would like to discuss, feel free to contact [Martijn van den Heuvel](http://www.dutchconnectomelab.nl).

How to contribute to CATO
---------------------------
If you like to contribute code, you can open a new GitHub pull request with the patch. We can then work together to merge this pull request into CATO. You can also always send us a message to discuss your ideas.

Style guide
---------------------------
For MATLAB code, the repository follows the [MATLAB Programming Style Guidelines by Richard Johnson](https://www.ee.columbia.edu/~marios/matlab/MatlabStyle1p5.pdf).

For shell code, the repository follows the [Google Shell Style Guide](https://google.github.io/styleguide/shellguide.html).

Code changes should always be validated before being merged into the development branch. CATO uses a MATLAB test suite initiated using `run_testSuite.m`. If your code changes are not covered by tests, then create a test function to validate the changes.
