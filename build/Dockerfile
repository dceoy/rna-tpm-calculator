FROM dceoy/r-tidyverse

ADD http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip /tmp/fastqc.zip
ADD https://github.com/BenLangmead/bowtie2/archive/master.tar.gz /tmp/bowtie2.tar.gz
ADD https://github.com/deweylab/RSEM/archive/master.tar.gz /tmp/rsem.tar.gz
ADD https://sourceforge.net/projects/prinseq/files/latest/download?source=files /tmp/prinseq.tar.gz
ADD https://github.com/samtools/htslib/archive/master.tar.gz /tmp/htslib.tar.gz
ADD https://github.com/samtools/samtools/archive/master.tar.gz /tmp/samtools.tar.gz

RUN set -e \
      && apt-get -y update \
      && apt-get -y upgrade \
      && apt-get -y install libtbb-dev perl pigz \
      && apt-get clean

RUN set -e \
      && unzip -d /usr/local/src /tmp/fastqc.zip \
      && chmod +x /usr/local/src/FastQC/fastqc \
      && ln -s /usr/local/src/FastQC/fastqc /usr/local/bin/fastqc

RUN set -e \
      && tar xvf /tmp/bowtie2.tar.gz -C /usr/local/src \
      && cd /usr/local/src/bowtie2-master \
      && make \
      && make install

RUN set -e \
      && tar xvf /tmp/rsem.tar.gz -C /usr/local/src \
      && cd /usr/local/src/RSEM-master \
      && make \
      && make install

RUN set -e \
      && tar xvf /tmp/prinseq.tar.gz -C /usr/local/src \
      && mv /usr/local/src/prinseq-* /usr/local/src/prinseq \
      && chmod +x /usr/local/src/prinseq/*.pl \
      && find /usr/local/src/prinseq -type f -maxdepth 1 -executable | xargs -i ln -s {} /usr/local/bin

RUN set -e \
      && tar xvf /tmp/htslib.tar.gz -C /usr/local/src \
      && cd /usr/local/src/htslib-master \
      && autoheader \
      && autoconf \
      && ./configure \
      && make \
      && make install

RUN set -e \
      && tar xvf /tmp/samtools.tar.gz -C /usr/local/src \
      && cd /usr/local/src/samtools-master \
      && autoconf \
      && ./configure \
      && make \
      && make install

ENTRYPOINT ["bash", "-c"]
