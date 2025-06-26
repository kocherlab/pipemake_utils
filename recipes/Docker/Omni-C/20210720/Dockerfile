FROM condaforge/mambaforge
ENV PATH /opt/conda/envs/Omni-C/bin:$PATH
ENV LD_LIBRARY_PATH /opt/conda/envs/Omni-C/lib:$LD_LIBRARY_PATH
RUN apt-get update -y && apt-get install unzip liblz4-dev build-essential zlib1g-dev -y
RUN mamba create -n Omni-C -c conda-forge bioconda::htslib bioconda::pysam=0.23.0 conda-forge::tabulate=0.9.0 bioconda::pairtools=1.1.3 bioconda::bwa=0.7.19 bioconda::samtools conda-forge::numpy conda-forge::scipy conda-forge::matplotlib conda-forge::pandas bioconda::pyBigWig bioconda::py2bit bioconda::deeptools python
RUN echo "source activate Omni-C" > ~/.bashrc
RUN wget https://github.com/smithlabcode/preseq/releases/download/v3.1.2/preseq-3.1.2.tar.gz
RUN tar xvf preseq-3.1.2.tar.gz
RUN cd preseq-3.1.2/ && ./configure --enable-hts CPPFLAGS='-I /opt/conda/envs/Omni-C/include/' LDFLAGS='-L/opt/conda/envs/Omni-C/lib' && make && make install
RUN wget https://github.com/dovetail-genomics/Omni-C/archive/refs/heads/main.zip
RUN unzip main.zip
RUN mv Omni-C-main/get_qc.py /opt/conda/envs/Omni-C/bin/get-qc
RUN chmod +x /opt/conda/envs/Omni-C/bin/get-qc