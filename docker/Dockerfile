# use the ubuntu base image, see https://hub.docker.com/_/ubuntu/
FROM ubuntu:16.04

MAINTAINER Yulin CHE yche@cse.ust.hk

# install oh-my-zsh
RUN apt-get update && apt-get install -y \
    zsh \
    git-core && rm -rf /var/lib/apt/lists/*
RUN git clone https://github.com/robbyrussell/oh-my-zsh.git ~/.oh-my-zsh \
    && cp ~/.oh-my-zsh/templates/zshrc.zsh-template ~/.zshrc \
    && chsh -s /bin/zsh

# install required packages
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    g++ \
    git \
    libboost-all-dev \
    python \
    python-dev \
    python-pip \
    && apt-get clean

# set environment
ENV BOOST_ROOT /usr

