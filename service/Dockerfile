FROM python:3.9-slim-bullseye

RUN apt update && apt-get install --no-install-recommends -y \
    ca-certificates \
    wget \
    bzip2 \
    unzip \
    build-essential \
    libcurl4-openssl-dev \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev

ENV BCFTOOLS_VERSION="1.19"
RUN wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
    && tar xjf bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
    && rm bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
    && cd bcftools-${BCFTOOLS_VERSION} \
    && ./configure \
    && cd plugins \
    && wget https://software.broadinstitute.org/software/score/score_${BCFTOOLS_VERSION}-dev.zip \
    && unzip score_${BCFTOOLS_VERSION}-dev.zip \
    && rm score_${BCFTOOLS_VERSION}-dev.zip \
    && cd .. \
    && make \
    && make install \
    && make clean

RUN wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver -O /usr/bin/liftOver \
    && chmod +x /usr/bin/liftOver

# downloaded from
#  https://hgdownload.soe.ucsc.edu/goldenPath/hs1/
#  https://hgdownload.soe.ucsc.edu/goldenPath/hg19/
#  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/
COPY hg38ToHs1.over.chain.gz /
COPY hs1ToHg38.over.chain.gz /
COPY hg19ToHg38.over.chain.gz /
COPY hg38ToHg19.over.chain.gz /

COPY requirements.txt /
RUN python3 -m pip install --upgrade -r /requirements.txt

COPY server.py /

ENV PORT=8080

CMD exec gunicorn --preload --bind :$PORT --workers 1 --threads 3 --timeout 0 server:app
