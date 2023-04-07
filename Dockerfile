FROM alpine:3.16
COPY ./ /src

RUN pip3 install -e . 