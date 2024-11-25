FROM gcc:latest
COPY . /usr/src/multik
WORKDIR /usr/src/multik
RUN ./compile.sh
CMD ["./myapp"]