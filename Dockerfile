FROM rocker/rstudio

RUN Rscript -e 'install.packages("BiocManager")' -e 'BiocManager::install("pcaMethods", ask=FALSE)'
RUN Rscript -e 'install.packages("remotes")'
RUN Rscript -e 'remotes::install_version("fda", version="2.4.8.1")'
RUN Rscript -e 'remotes::install_github("ondrej-klima/growthfd", upgrage="never")'

EXPOSE 8787
CMD ["/init"]