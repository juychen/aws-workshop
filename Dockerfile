FROM python:3.7.5-slim
ARG KB_VERSION=0.24.4
RUN apt-get update && \
	apt-get install --no-install-recommends -y curl dpkg-dev gnupg lsb-release procps
ENV PATH="/usr/local/bin:${PATH}"
RUN pip install --upgrade pip && pip install kb-python