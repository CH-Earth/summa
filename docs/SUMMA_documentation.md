# Navigating the SUMMA Documentation

The SUMMA documentation remains a work in progress. It can be navigated online on [summa.readthedocs.io](http://summa.readthedocs.io/) using the navigation panel to the left.

If you are new to SUMMA, start with the section on _Installation_ and make sure that you can run the SUMMA test suite. After that, it will depend on what you want to do. The _Development_ section is mostly of interest if you want to contribute to the SUMMA source code. Model users will want to read through the _Configuration_ and _Input/Output_ sections to understand how to configure SUMMA for there own applications.


## Contributing to SUMMA Documentation
SUMMA documentation is hosted on [summa.readthedocs.io](http://summa.readthedocs.io/), written in [Markdown](https://daringfireball.net/projects/markdown/syntax) and built using [MkDocs](http://www.mkdocs.org/). If you want to contribute to the documentation, you can do so by forking the [SUMMA repository](https://www.github.com/NCAR/summa), creating a branch for your changes and editing the documentation files in the `docs` directory in the SUMMA repo. You may need to add `mkdocs.yml` in the top level SUMMA directory to ensure that your page shows up in the right place.

You need to install [MkDocs](http://www.mkdocs.org/) locally, so that you can make sure that your edits show up correctly before you make a pull request to the NCAR SUMMA repo. If you are new to Git, please review the [SUMMA and Git](development/SUMMA_and_git.md) and [SUMMA Git Workflow](development/SUMMA_git_workflow.md) pages. We will only accept documentation updates via a pull request.
