#FROM ubuntu:22.04
FROM nvcr.io/nvidia/cuda:12.1.0-devel-ubuntu22.04 

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install cron -y && \
    apt-get install openssh-server -y && \
    apt-get install net-tools -y && \
    apt-get install wget -y && \
    apt-get install vim -y && \
    apt-get install curl -y && \
    apt-get install python3-pip -y && \
    apt-get install unzip -y && \
    apt-get install unixodbc unixodbc-dev -y && \
	apt-get install libgl1 -y && \
    apt-get install git -y && \
    apt install mpich --reinstall -y

# set image time zone
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends tzdata
RUN TZ=Asia/Taipei && \
    ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && \
    echo $TZ > /etc/timezone && \
    dpkg-reconfigure -f noninteractive tzdata

# set for workdir
WORKDIR /root
RUN mkdir /root/code/
COPY requirements.txt /root/requirements.txt
RUN pip install --upgrade pip
RUN pip install -r requirements.txt
