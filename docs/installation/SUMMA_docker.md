# SUMMA using Docker

If you are not interested in compiling SUMMA locally on your machine and you are not planning to contribute to SUMMA development, then you may be interested in using a version of SUMMA that runs in a [Docker container](https://www.docker.com/what-docker).

To run SUMMA using [Docker](https://www.docker.com) you will need to do the following:

  1. Install Docker for your platform. Note that the [Docker Community Edition (Docker CE)](https://www.docker.com/community-edition) will work just fine (it's what we use) and it's free. There are Docker CE version for Mac, Windows, and a number of \*nix distributions.

  1. Once you have Docker installed, you can get the Docker containers for different versions of SUMMA from [Docker Hub](https://hub.docker.com). Currently the Docker containers are distributed as part of the [Docker Hub account of bartnijssen](https://hub.docker.com/r/bartnijssen/summa/) and you can find out which version are being built by looking at the [Build Settings](https://hub.docker.com/r/bartnijssen/summa/~/settings/automated-builds/). For example, the Docker Tag Name `latest` reflects the SUMMA version on the master branch in the SUMMA reposistory and `develop` reflects the SUMMA version on the develop branch. There may be a few others as well.

  1. To download the appropriate container, make sure that Docker is running and open a terminal. Then type

    `docker pull bartnijssen/summa:<Docker Tag Name>`

    where `<Docker Tag Name>` is equal to `latest`, `develop` or whichever tag you are interested in. This will download the container to your local machine.

  1. You are now ready to run SUMMA. No further installs are required. You can run SUMMA by typing

    `docker run bartnijssen/summa:<Docker Tag Name>`

  1. You can give SUMMA command-line arguments directly after the name of the docker container you are running. For example, to see which version of SUMMA you are running provide the `-v` command-line option. If you use `develop` as the `<Docker Tag Name>`, this would look something like this (the exact message will depend on the SUMMA version you are running)

    `docker run bartnijssen/summa:develop -v`

         SUMMA - Structure for Unifying Multiple Modeling Alternatives
                               Version: v2.0.0
                   Build Time: Fri Sep  1 22:23:26 UTC 2017
                        Git Branch: develop-0-gb9515d7
              Git Hash: b9515d7d7f44ae7043eac37caa66e7bc9f946a04

  1. To do actual SUMMA runs, you need to provide a bit more information. The main thing is that you will need to provide the path of the master file (`-m`) so that SUMMA knows where to read and write its input and output. The other part is that you need to set up a mapping between your local file paths and the path that SUMMA has access to within the Docker container. This mapping can be set up on the Docker command-line by using the `-v` option to `docker run`. See `docker run --help` for more details.

    The [SUMMA test cases](SUMMA_test_cases.md) include install and run scripts in the top level directory to configure and run the test cases using SUMMA on Docker. We recommend that you start there to understand how to make an actual SUMMA run using Docker.
