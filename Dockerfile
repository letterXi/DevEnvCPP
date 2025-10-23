FROM ubuntu:24.04

RUN apt-get update  && \
    apt-get install -y gcc-14  g++-14 cmake gdb && \
    update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-14 100 && \
    update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-14 100 && \
    apt-get install -y git && \
    apt-get install -y cmake-format clang-format && \
    apt-get install -y sudo rsync

RUN mkdir /app

RUN usermod -l user ubuntu && \
    groupmod -n user ubuntu && \
    usermod -d /home/user -m user && \
    echo "user ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers



RUN if ["$HOST_UID" != "1000"]; then usermod -u $HOST_UID user; fi && \
    if ["$HOST_GID" != "1000"]; then \
        if ! getent group $HOST_GID > /dev/null; then \
            groupmod -g $HOST_GID user; \
        else \
            usermod -g $HOST_GID user; \
        fi \
    fi

USER user

RUN echo 'export CC=/usr/bin/gcc-14' >> /home/user/.bashrc  && \
    echo 'export CXX=/usr/bin/g++-14' >> /home/user/.bashrc

CMD /app/scripts/init/init_container.sh