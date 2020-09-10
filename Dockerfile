FROM nfcore/base:1.9

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/quaisar_coregenome/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name quaisar_coregenome > quaisar_coregenome.yml
RUN ktUpdateTaxonomy.sh
