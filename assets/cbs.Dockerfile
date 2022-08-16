FROM r-base:latest

#RUN apt update -qq
#RUN apt install --no-install-recommends -y r-cran-rjava

RUN R -e "install.packages('BiocManager',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('DNAcopy')"
RUN R -e "BiocManager::install('aCGH')"
RUN R -e "install.packages('optparse',dependencies=TRUE, repos='http://cran.rstudio.com/')"

COPY CBS_MergeLevels.R CBS_MergeLevels.R
ENTRYPOINT ["Rscript", "CBS_MergeLevels.R"]