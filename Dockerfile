FROM debian

COPY ./NASP /usr/local/bin/ 

ENTRYPOINT ./NASP
