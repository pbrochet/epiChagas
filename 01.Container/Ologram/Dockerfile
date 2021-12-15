
FROM ubuntu:18.04
MAINTAINER Pauline Brochet (pauline.brochet@univ-amu.fr)


#############################################################
#                Utilitaire (wget, zip, git)                #
#############################################################

RUN apt-get update -y
RUN apt-get update -y --fix-missing
RUN apt-get upgrade -y

# #############################################################
# #                       R installation                      #
# #############################################################



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
  
#   # Add a default CRAN mirror
# RUN echo "options(repos = c(CRAN = 'https://cran.rstudio.com/'), download.file.method = 'libcurl')" >> /usr/local/lib/R/etc/Rprofile.site

#   # Add a library directory (for user-installed packages)
# RUN mkdir -p /usr/local/lib/R/site-library \
#   && chown root:staff /usr/local/lib/R/site-library \
#   && chmod g+wx /usr/local/lib/R/site-library
  
#   # Fix library path
# RUN echo "R_LIBS_USER='/usr/local/lib/R/site-library'" >> /usr/local/lib/R/etc/Renviron \
#   && echo "R_LIBS=\${R_LIBS-'/usr/local/lib/R/site-library:/usr/local/lib/R/library:/usr/lib/R/library'}" >> /usr/local/lib/R/etc/Renviron

#   # Install packages from Microsoft R Application Network # certains packages sont sur cran microsoft 
# RUN [ -z "$BUILD_DATE" ] && BUILD_DATE=$(date -I --date='TZ="America/Los_Angeles"') || true \
#   && MRAN=https://mran.microsoft.com/snapshot/${BUILD_DATE} \
#   && echo MRAN=$MRAN >> /etc/environment \
#   && export MRAN=$MRAN \
#   && echo "options(repos = c(CRAN='$MRAN'), download.file.method = 'libcurl')" >> /usr/local/lib/R/etc/Rprofile.site



#############################################################
#             Installing dependency for Ologram             #
#############################################################

RUN apt-get install -y wget
RUN apt-get install -y git
RUN apt-get install -y build-essential
RUN apt-get install -y make
RUN apt-get install -y python3.8
RUN apt-get install -y python3.8-dev
RUN apt-get install -y python3-pip
RUN apt-get install -y libssl-dev
RUN apt-get install -y libz-dev
RUN apt-get install -y bedtools
RUN apt-get install -y graphviz
RUN ln -sfn /usr/bin/python3.8 /usr/bin/python3 && ln -sfn /usr/bin/python3 /usr/bin/python && ln -sfn /usr/bin/pip3 /usr/bin/pip

RUN pip3 install numpy
RUN pip3 install Cython
RUN pip3 install pandas
RUN python -m pip install -U pip
RUN pip3 install --upgrade setuptools
RUN pip3 install ez_setup
RUN pip3 install statsmodels

#############################################################
#                        Ologram                            #
#############################################################

RUN git clone https://github.com/dputhier/pygtftk/
RUN cd pygtftk && git checkout develop
RUN cd pygtftk && pip3 install -r requirements.txt 
RUN export LC_ALL=C.UTF-8
RUN export LANG=C.UTF-8
RUN export LC_CTYPE=en_US.UTF-8
RUN echo "alias python=python3" >> ~/.bashrc
RUN echo "export LC_ALL=C.UTF-8" >> ~/.bashrc
RUN echo "LANG=C.UTF-8" >> ~/.bashrc
RUN echo "LC_CTYPE=en_US.UTF-8" >> ~/.bashrc
RUN cd pygtftk && python3.8 setup.py install
RUN cd pygtftk && make install