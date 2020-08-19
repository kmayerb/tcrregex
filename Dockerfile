FROM continuumio/miniconda3:4.8.2

RUN apt-get update && apt-get install -y procps && apt-get install -y nano && apt-get -y install gcc && apt-get -y install unzip && apt-get -y install curl && apt-get -y install wget

RUN pip install git+https://github.com/kmayerb/tcrregex.git@0.1.0

RUN pip install pytest
RUN pip install ipython

RUN python -c "import tcrregex as td; td.install_test_files.install_test_files()"
RUN python -c "import tcrregex as td; td.setup_db.install_all_next_gen()"
