#Docker Guidance

##Install Docker and Start Service
start in fedora:

```zsh
sudo systemctl status docker
sudo systemctl start docker
```

If you want Docker to start at boot, you should also:

```zsh
sudo systemctl enable docker
```

##Check Whether It Starts

In terminal, type `sudo docker run hello-world`.

If result is as follows, then it is okay.

```zsh
Hello from Docker!
This message shows that your installation appears to be working correctly.

To generate this message, Docker took the following steps:
 1. The Docker client contacted the Docker daemon.
 2. The Docker daemon pulled the "hello-world" image from the Docker Hub.
 3. The Docker daemon created a new container from that image which runs the
    executable that produces the output you are currently reading.
 4. The Docker daemon streamed that output to the Docker client, which sent it
    to your terminal.

To try something more ambitious, you can run an Ubuntu container with:
 $ docker run -it ubuntu bash

Share images, automate workflows, and more with a free Docker ID:
 https://cloud.docker.com/

For more examples and ideas, visit:
 https://docs.docker.com/engine/userguide/
```

##Choice1, Simply Pull From Docker Hub

do as follows, you will get docker image from my docker hub repository

```zsh
sudo docker pull yche/yche-dev-env
```

##Choice2, Build Docker Image From Scratch 
###Build
build image with script [build_img.sh](build_img.sh)

```zsh
./build_img.sh
```

###Push 2 Remote Hub
give docker tag, where `10d386ca7c19` should be replaced with your docker image id

```zsh
sudo docker tag 10d386ca7c19 yche/yche-dev-env:latest
```

login docker hub

```zsh
sudo docker login
```

push, where `yche/yche-dev-env` should be replaced with your docker hub repository

```zsh
sudo docker push yche/yche-dev-env
```

##Usage of Docker Image

first load docker image

```zsh
python run_docker.py
```

second enter into workspace

```zsh
cd /opt/
```

