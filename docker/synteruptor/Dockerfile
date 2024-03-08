FROM ubuntu:22.04
LABEL maintainer="mat.barba@gmail.com"
LABEL version="1.0.1"
LABEL description="Set of scripts to build a Synteruptor database."

# Disable Prompt During Packages Installation
ARG DEBIAN_FRONTEND=noninteractive

RUN apt update
RUN apt install -y git parallel
RUN apt install -y bioperl libstatistics-basic-perl ncbi-blast+ sqlite3

# Install synteruptor
RUN git clone https://github.com/i2bc/synteruptor.git && cd synteruptor && git checkout main
RUN chmod +x synteruptor/*.sh /synteruptor/*.pl
ENV PATH="$PATH:/synteruptor"

ENTRYPOINT ["run_gbk.sh"]
