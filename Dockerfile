FROM gcc:14.1
FROM python:3.11

COPY . /usr/src/multik

WORKDIR /usr/src/multik

RUN pip install pytest
RUN mkdir results
RUN sed -i 's/\r$//' compile.sh */compile.sh test/test.sh
RUN bash compile.sh

CMD ["bash", "./pipeline.sh"]