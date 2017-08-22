FROM dceoy/r-tidyverse

ADD http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip /tmp/fastqc.zip
ADD https://github.com/BenLangmead/bowtie2/releases/download/v2.3.2/bowtie2-2.3.2-linux-x86_64.zip /tmp/bowtie2.zip
ADD https://github.com/deweylab/RSEM/archive/master.tar.gz /usr/local/src
ADD https://sourceforge.net/projects/prinseq/files/latest/download?source=files /usr/local/src
ADD https://github.com/samtools/htslib/archive/master.tar.gz /usr/local/src
ADD https://github.com/samtools/samtools/archive/master.tar.gz /usr/local/src

RUN set -e \
      && apt-get -y update \
      && apt-get -y upgrade \
      && apt-get -y install autoconf libtbb-dev perl pigz \
      && apt-get clean

RUN set -e \
      && unzip -d /usr/local/src /tmp/fastqc.zip \
      && chmod +x /usr/local/src/FastQC/fastqc \
      && ln -s /usr/local/src/FastQC/fastqc /usr/local/bin/fastqc

RUN set -e \
      && unzip -d /usr/local/src /tmp/bowtie2.zip \
      && mv /usr/local/src/bowtie2-* /usr/local/src/bowtie2 \
      && find /usr/local/src/bowtie2 -type f -maxdepth 1 -executable | xargs -i ln -s {} /usr/local/bin

RUN set -e \
      && cd /usr/local/src/RSEM-master \
      && make \
      && make install

RUN set -e \
      && mv /usr/local/src/prinseq-* /usr/local/src/prinseq \
      && chmod +x /usr/local/src/prinseq/*.pl \
      && find /usr/local/src/prinseq -type f -maxdepth 1 -executable | xargs -i ln -s {} /usr/local/bin

RUN set -e \
      && cd /usr/local/src/htslib-master \
      && autoheader \
      && autoconf \
      && ./configure \
      && make \
      && make install

RUN set -e \
      && cd /usr/local/src/samtools-master \
      && autoconf \
      && ./configure \
      && make \
      && make install

ENTRYPOINT ["bash", "-c"]