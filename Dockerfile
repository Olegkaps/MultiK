FROM gcc:latest
FROM python:latest
COPY . /usr/src/multik
WORKDIR /usr/src/multik

RUN pip install pytest
RUN mkdir exe
RUN mkdir test/exe
RUN sed -i 's/\r$//' compile.sh */compile.sh test/test.sh
RUN bash compile.sh

CMD ["bash", "./pipeline.sh"]