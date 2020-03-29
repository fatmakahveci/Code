FROM python:3.7

MAINTAINER Fatma Kahveci fatmaba@gmail.com

ENV PYTHONUNBUFFERED 1

RUN pip install docopt==0.6.2 intervaltree==3.0.2 numpy==1.18.2 pysam==0.15.4 sortedcontainers==2.1.0

COPY ./estimateSampleComposition.py /opt/
