#TianR.
#Mar 30, 2015
#Mar 31, 2015, Tue.
#April 1, 2015, Wed.

#install docker in MaxOS 10.6+
################################
#download boot2docker, a linux virtual machine under MacOS
https://github.com/boot2docker/osx-installer/releases/download/v1.5.0/Boot2Docker-1.5.0.pkg

boot2docker init

boot2docker start

boot2docker status

#set the environment variables
$(boot2docker shellinit)  ##very important

#sudo service docker start,  in linux


docker version

docker run -d -P --name web nginx

echo 'export DOCKER_HOST="tcp://$(boot2docker ip 2>/dev/null):2375";' >> ~/.bashrc

source ~/.bashrc



###########################
#1 Hello World

https://www.docker.com/whatisdocker/

try it

two videos




#################################
#2 customized ubuntu OS
docker login
docker search ubuntu
docker pull ubuntu

docker images

docker run -t -i ubuntu

#within unbuntu apt-get install vim
#exit

docker ps

#container ID
docker commit 00088ea4bc5e
docker tag dc2beaa95ca9 ubuntu:14.04_rt_vim

#from image file to tar file
docker save ubuntu:14.04_rt_vim > /tmp/ubuntu:14.04_rt_vim.tar

#delete container docker rm 
#delete image docker rmi

#load tar file as docker image
docker load < /tmp/ubuntu\:14.04_rt_vim.tar 
docker images
docker run -i -t 
docker ps



##################################
#3 Dockfile and mongoDB
docker build --rm=true --tag="myyy/repo" .
docker images
docker run --name mongo_instance1 -d myyy/repo
docker logs mongo_instance1



#################################
#4 To mount the directory in host to a Container
docker run -v /Users/tianr/0projects/Learn/00Docker/testdir/:/tmp d0955f21bf24 ls /tmp

docker run -i -t -v /Users/tianr/Desktop/CancerSeqMachine_docker_2015Mar24-May24/RawData/Exome_Trio_Corpas/:/tmp bb4ee6fb9845 bowtie -qS /opt/local/software/bowtie-1.1.1/indexes/e_coli /tmp/test.fq 

docker run -i -t -v /Users/tianr/Desktop/CancerSeqMachine_docker_2015Mar24-May24/RawData/Exome_Trio_Corpas/:/tmp bb4ee6fb9845 bowtie -qS /tmp/e_coli.ebwt/e_coli /tmp/test.fq
