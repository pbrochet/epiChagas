
FROM ubuntu:18.04
MAINTAINER Pauline Brochet (pauline.brochet@univ-amu.fr)

#############################################################
#                Utilitaire (wget, zip)                     #
#############################################################

RUN apt-get update -y --fix-missing
RUN apt-get install -y wget
RUN apt-get install -y zip
RUN apt-get install -y git
RUN apt-get install -y build-essential
RUN apt-get install -y zlib1g-dev
RUN apt-get install -y make
RUN apt-get install -y python3.7
RUN apt-get install -y python3-pip
RUN apt-get install -y python-dev

#############################################################
#                         pigz                              #
#############################################################

RUN apt-get update -y && apt-get install -y pigz

#############################################################
#                        fastQC                             #
#############################################################

RUN apt-get install -y perl
RUN apt-get install -y xauth
RUN apt update && apt install -y fastqc

#############################################################
#                    Trimmomatic                            #
#############################################################

RUN apt update && apt install -y openjdk-11-jre
RUN apt-get install -y trimmomatic
# RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
# RUN unzip Trimmomatic-0.39.zip
# RUN mv Trimmomatic-0.39 trimmomatic
# RUN rm Trimmomatic-0.39.zip
# RUN chmod +x trimmomatic/trimmomatic-0.39.jar
# RUN ln -s trimmomatic/trimmomatic-0.39.jar Trimmomatic

#############################################################
#                           STAR                            #
#############################################################

RUN apt-get update -y --fix-missing
RUN apt-get install -y rna-star


#############################################################
#                 featureCounts (subread)                   #
#############################################################

RUN apt-get install -y subread


#############################################################
#                         multiQC                           #
#############################################################

RUN git clone https://github.com/ewels/MultiQC.git
RUN export LC_ALL=C.UTF-8
RUN export LANG=C.UTF-8
RUN cd MultiQC && python3.7 setup.py install


#############################################################
#                       R installation                      #
#############################################################


ARG BUILD_DATE
ARG R_VERSION
ENV BUILD_DATE ${BUILD_DATE:-}
ENV R_VERSION ${R_VERSION:-3.6.2}
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8


## Dependencies

  # Dependencies
RUN apt-get update \
  && apt-get -y install \
  python2.7-dev \
  python-pip \
  build-essential \
  pkg-config \
  libffi-dev \
  libmagickwand-dev \
  python-cffi

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    bash-completion \
    ca-certificates \
    file \
    fonts-texgyre \
    g++ \
    gfortran \
    gsfonts \
    libbz2-1.0 \
    libcurl3 \
    # libicu55 \
    libturbojpeg \
    libopenblas-dev \
    libpangocairo-1.0-0 \ 
    libpcre3 \
    # libpng12-0 \
    libtiff5 \ 
    liblzma5 \
    locales \
    make \
    unzip \
    zip \
    libssl-dev \
    zlib1g

  # Configure default locale
RUN  echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
  && locale-gen en_US.utf8 \
  && /usr/sbin/update-locale LANG=en_US.UTF-8

  # Dependencies
RUN BUILDDEPS="curl \
    default-jdk \
    libbz2-dev \
    libcairo2-dev \ 
    libcurl4-openssl-dev \
    libpango1.0-dev \
    libjpeg-dev \
    libicu-dev \
    libpcre3-dev \
    libpng-dev \
    libreadline-dev \
    libtiff5-dev \
    liblzma-dev \ 
    libx11-dev \
    libxt-dev \
    perl \
    texinfo \
    texlive-extra-utils \
    texlive-fonts-recommended \
    texlive-fonts-extra \
    texlive-latex-recommended \
    x11proto-core-dev \
    xauth \
    xfonts-base \
    xvfb \
    zlib1g-dev" \
  && apt-get update \
  && apt-get install -y --no-install-recommends $BUILDDEPS # installe les dépendances listées


## Install R

RUN cd /tmp/ \
# Download and extract source code
  && curl -O https://cran.r-project.org/src/base/R-3/R-${R_VERSION}.tar.gz \
  && tar -xf R-${R_VERSION}.tar.gz \
  && cd R-${R_VERSION} \ 
# Set compiler flags
  && R_PAPERSIZE=letter \
    R_BATCHSAVE="--no-save --no-restore" \
    R_BROWSER=xdg-open \
    PAGER=/usr/bin/pager \
    PERL=/usr/bin/perl \
    R_UNZIPCMD=/usr/bin/unzip \
    R_ZIPCMD=/usr/bin/zip \
    R_PRINTCMD=/usr/bin/lpr \
    LIBnn=lib \
    AWK=/usr/bin/awk \
    CFLAGS="-g -O2 -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g" \
    CXXFLAGS="-g -O2 -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g" \
# Configure options
  ./configure --enable-R-shlib \
               --enable-memory-profiling \
               --with-readline \
               --with-blas="-lopenblas" \
               --disable-nls \
               --without-recommended-packages \
# Build and install
  && make \
  && make install
  
  # Add a default CRAN mirror
RUN echo "options(repos = c(CRAN = 'https://cran.rstudio.com/'), download.file.method = 'libcurl')" >> /usr/local/lib/R/etc/Rprofile.site

  # Add a library directory (for user-installed packages)
RUN mkdir -p /usr/local/lib/R/site-library \
  && chown root:staff /usr/local/lib/R/site-library \
  && chmod g+wx /usr/local/lib/R/site-library
  
  # Fix library path
RUN echo "R_LIBS_USER='/usr/local/lib/R/site-library'" >> /usr/local/lib/R/etc/Renviron \
  && echo "R_LIBS=\${R_LIBS-'/usr/local/lib/R/site-library:/usr/local/lib/R/library:/usr/lib/R/library'}" >> /usr/local/lib/R/etc/Renviron

  # Install packages from Microsoft R Application Network # certains packages sont sur cran microsoft 
RUN [ -z "$BUILD_DATE" ] && BUILD_DATE=$(date -I --date='TZ="America/Los_Angeles"') || true \
  && MRAN=https://mran.microsoft.com/snapshot/${BUILD_DATE} \
  && echo MRAN=$MRAN >> /etc/environment \
  && export MRAN=$MRAN \
  && echo "options(repos = c(CRAN='$MRAN'), download.file.method = 'libcurl')" >> /usr/local/lib/R/etc/Rprofile.site

  # Use littler installation scripts
RUN Rscript -e "install.packages(c('littler', 'docopt'), repo='$MRAN')" \
  && ln -s /usr/local/lib/R/site-library/littler/examples/install2.r /usr/local/bin/install2.r \
  && ln -s /usr/local/lib/R/site-library/littler/examples/installGithub.r /usr/local/bin/installGithub.r \
  && ln -s /usr/local/lib/R/site-library/littler/bin/r /usr/local/bin/r


#############################################################
#                 R package installation                    #
#############################################################

# Dependencies
RUN apt-get -y install libxml2
RUN apt-get -y install libssl1.0.0 libssl-dev
RUN apt install -y libcurl4-openssl-dev
RUN apt-get -y install pandoc

# R package
RUN R -e 'install.packages("http://cran.r-project.org/src/contrib/Archive/foreign/foreign_0.8-76.tar.gz", type="source")'
RUN R -e 'install.packages("FactoMineR", repos = "http://cran.rstudio.com/")'
RUN R -e 'install.packages("factoextra", repos = "http://cran.rstudio.com/")'
RUN R -e 'install.packages("ggplot2", repos = "http://cran.rstudio.com/")'
RUN R -e 'install.packages("gplots", repos = "http://cran.rstudio.com/")'
RUN R -e 'install.packages("gProfileR", repos = "http://cran.rstudio.com/")'
RUN R -e 'install.packages("dplyr", repos = "http://cran.rstudio.com/")'
RUN R -e 'install.packages("reshape2", repos = "http://cran.rstudio.com/")'
RUN R -e 'install.packages("http://cran.r-project.org/src/contrib/Archive/XML/XML_3.99-0.3.tar.gz", type="source")'

# R bioconductor packages
RUN R -e 'install.packages("BiocManager", repos = "http://cran.rstudio.com/")'
RUN R -e 'BiocManager::install("biomaRt")'
RUN R -e 'BiocManager::install("DESeq2")'
RUN R -e 'BiocManager::install("EnhancedVolcano")'
RUN R -e 'BiocManager::install("pathview")'
RUN R -e 'BiocManager::install("gage")'
RUN R -e 'BiocManager::install("gageData")'
RUN R -e 'BiocManager::install("apeglm")'



#############################################################
#                Qualimap (R dependencies)                  #
#############################################################

RUN mkdir -p /usr/lib/jvm && cd /usr/lib/jvm && wget https://download.java.net/java/GA/jdk11/9/GPL/openjdk-11.0.2_linux-x64_bin.tar.gz && tar xvf openjdk-11.0.2_linux-x64_bin.tar.gz
RUN update-alternatives --install /usr/bin/java java /usr/lib/jvm/jdk-11.0.2/bin/java 1 && update-alternatives --install /usr/bin/javac javac /usr/lib/jvm/jdk-11.0.2/bin/javac 1

# R dependencies
RUN R -e 'install.packages("optparse", repos = "http://cran.rstudio.com/")'
RUN R -e 'BiocManager::install("NOISeq")'
RUN R -e 'BiocManager::install("Repitools")'
RUN R -e 'BiocManager::install("Rsamtools")'
RUN R -e 'BiocManager::install("rtracklayer")'

RUN curl -OL https://bitbucket.org/kokonech/qualimap/downloads/qualimap-build-31-08-20.tar.gz && \
    tar -xf qualimap-build-31-08-20.tar.gz &&\
    cd qualimap-build-31-08-20 && Rscript scripts/installDependencies.r

#############################################################
#                     Clean installation                    #
#############################################################

# # Clean R install
# RUN cd / \
#   && rm -rf /tmp/* \
#   && apt-get remove --purge -y $BUILDDEPS \
#   && apt-get autoremove -y \
#   && apt-get autoclean -y \
#   && rm -rf /var/lib/apt/lists/*

# RUN apt-get clean \
#   && apt-get remove -y --purge build-essential

# RUN export DISPLAY=:0.0
# RUN apt-get install -y gnumeric

RUN easy_install xlsx2csv
RUN R -e 'BiocManager::install("pcaMethods")'
RUN R -e 'install.packages("ADAPTS", repos = "http://cran.rstudio.com/")'
