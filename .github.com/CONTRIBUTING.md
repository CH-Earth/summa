# Contributing

SUMMA is Open Source software. This means that the code is made available for free, but also that development, maintenance and support are intended to be community efforts. Our rationale for moving SUMMA model development to an open source model is that we want:

- to encourage other researchers and developers to contribute to SUMMA development, and
- to facilitate transparent development and use of the model.

## Support

There is no official support for the SUMMA model, other than the SUMMA documentation, the SUMMA source code archive and the description of the model in the literature. Any additional support relies on volunteer efforts by the SUMMA development community. The following resources are available:

- [SUMMA web site](https://www.ral.ucar.edu/projects/summa): General background, SUMMA resources, and test data sets.
- [SUMMA Source code repository](https://github.com/NCAR/SUMMA) : Source code distribution, coordination of model development, bug fixes, and releases.

We expect that the user comes prepared with some understanding of the model and scientific computing. As such, these items are specifically not supported by the SUMMA development community:

- Building and running the SUMMA model on platforms other than LINUX, UNIX, and OSX.
- Using LINUX, UNIX, or OSX operating systems.
- Development of project specific features.
- Configuring individual model applications.

## Submitting Issues
#### Submitting Bug Reports

If you think you have found a bug in SUMMA, please check whether an issue has been filed on [SUMMA's Github page](https://github.com/NCAR/SUMMA/issues). If not, go please go ahead and create an issue and include the following information in your bug report:

- Version of SUMMA that you are using (e.g. SUMMA 1.0 - even better if you can provide the specific tag or commit)
- Name and version of the fortran compiler you are using
- Operating system
- A description of relevant model settings
- A summary of the bug or error message you are getting

If you can provide more information that is great. If you know how to run the model in a debugger, you may be able to pinpoint where the problem occurs.

#### Proposing New Features

SUMMA is under active development.  If you would like to propose a new feature, driver, or extension to SUMMA, please file an issue on [SUMMA's Github page](https://github.com/NCAR/SUMMA/issues). Also, because SUMMA is an open source model with no official support for general-purpose development, be prepared to contribute to the implementation of your feature request. Features that are only of interest to you are unlikely to be implemented in the main source code repo (although you are of course free to modify the code in any way you see fit).

## Contributing to SUMMA
#### Git Workflow
We have developed some documentation to help you get started if you are new to Git but want to contribute to VIC:

- [Working with Git](https://github.com/NCAR/summa/blob/master/docs/howto/git_howto.md)
- [Git Workflow](https://github.com/NCAR/summa/blob/master/docs/howto/summa_git_workflow.md)

#### Coding Conventions
We have some simple [https://github.com/NCAR/summa/blob/master/docs/howto/summa_coding_conventions.md](coding conventions) that we would like anyone who contributes code to SUMMA to follow.
