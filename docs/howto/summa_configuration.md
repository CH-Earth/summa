# SUMMA Configuration

## Introduction

SUMMA configuration is  performed via a large number of configuration files. The format of these files is likely to change (and simplified) as we build out SUMMA's infrastructure. The following may be somewhat difficult to follow without a working example. We strongly recommend that you get the [test applications](http://ral.ucar.edu/projects/summa) to help you get started.

In the following we will assume that your SUMMA source installation is in the directory `<localInstallation>` on your  machine (this is not actually the name of the directory, but can be whatever you call) and that the test applications are installed in the directory `<localApplications>`. Note that these two directories can be the same, in which case everything would be in the same directory tree.

## SUMMA directory tree

After building SUMMA and installing the test applications, you will have the following directories (there may be some additional ones)

* `<localInstallation>/bin`

    SUMMA executable

* `<localInstallation>/build`

    SUMMA source code

* `<localInstallation>/docs/howto`

    Some documentation, including this file

* `<localApplications>/input`

    Model forcing files. There are directories for a number of different locations.

* `<localApplications>/output`

    Model output files.

* `<localApplications>/scripts`

    A number of plotting and analysis scripts.

* `<localApplications>/settings`

    Model settings and parameter files. There are directories for a number of different test cases.

* `<localApplications>/verification`

    Some analysis and plotting scripts.

## Running the model

1. The model executable `summa.exe` is in the `<localInstallation>/bin` directory (which you can add to your path to make it easier to run SUMMA). To run SUMMA, type

        > <localInstallation>/bin/summa.exe

    You should get output that looks something like

        124430.755
        1st command-line argument missing, expect text string defining the output file suffix

1. To actually run the model, you need to provide two command-line arguments, for example

        > summa.exe <case name> <configuration file>

    where `<case_name>` is simply a string that will be used to identify your simulation. The string will be used to modify the model output files. It's easiest if you start `<case_name>` with an underscore (`_`) and then follow with some text. The `<configuration file>` is the main input file to SUMMA and includes the paths to all the other input files.

## SUMMA input files

The main SUMMA configuration file (the one that you specify on the command-line) simply provides a listing of other SUMMA input files. The following provides a very brief overview of the main SUMMA configuration file. We are not going into great detail, because this will all change as we build out SUMMA. Dedicated users who want to develop their own SUMMA applications are referred to the comments and documentation in the input files for the test applications as well as to the SUMMA source code. You can trace the location of all the files by starting with the main configuration file described below.

In the following we will use the first test case that is used for figure 1 in [Clark et al., 2015](http://dx.doi.org/10.1002/2015WR017200) to discuss briefly the model setup. This test case is described in `<localApplications>/settings/README` as "Figure 1: Radiation transmission through an Aspen stand, Reynolds Mountain East"`".

A few rules apply to all the configuration files:

 * `!` signifies the start of a comment. Any text from `!` till the end of the line is discarded (this can be an entire line)
 * strings (such as file paths) need to be enclosed in single quotes, i.e. `'`.
 * the order of the entries matters. SUMMA currently does not support "free-form" configuration files (for example, using dictionaries with key-entry pairs). This means that configuration options need to be specified in the order indicated in the sample files.

### Main configuration file

The main configuration file for this test case is `<localApplications>/settings/wrrPaperTestCases/figure01/summaFileManager_riparianAspenBeersLaw.txt`. The file has the following content (note that we assume in the following that comments have been stripped):

1. `SUMMA_FILE_MANAGER_V1.0` Simply repeat this string. It is used to identify the version of the file manager with which this configuration file is compatible.

1. base path for model settings. For the test applications this is `<localApplications>/settings`.

1. base path for input files.

1. base path for output files.

1. file with all the model decisions, that is, the various model decisions that determine how summa will be configured for a particular model applications (`M_DECISIONS`).

The paths that follow are briefly explained in the configuration file for the test case. The `meta` files provide metadata for model parameters. Often these values will remain the same for different model simulations (for example, we use a fixed set of metadata files for the test applications). The next set of files (beginning with `LOCAL_ATTRIBUTES`) provide location- and application-specific information. The final line provides a prefix for the model output files that can be used to identify your simulation.

## SUMMA output files

For now, the user is referred to the IDL plotting files in `<localApplications>/verification` to get a description of the SUMMA output files.