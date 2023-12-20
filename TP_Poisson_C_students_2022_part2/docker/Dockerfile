
FROM ubuntu:18.04

# name the container


RUN apt-get update 
RUN apt-get install -y gcc
RUN apt-get install -y libc-dev
RUN apt-get install -y libblas-dev 
RUN apt-get install -y liblapack-dev
RUN apt-get install -y liblapacke-dev
#RUN apt-get install -y libatlas-base-dev
RUN apt-get install -y vim
RUN apt-get install -y curl
RUN apt-get install -y make

## Add folders that will contains code before and after installation
RUN mkdir /home/tp-cn2
RUN mkdir -p /home/tp-cn2/shared_folder
VOLUME ["home/tp-cn2/shared_folder"]

## add the TP_Poisson folder to the container
ADD TP_Poisson_C_students_2022_part2 /home/tp-cn2/TP_Poisson

#make the working directory
WORKDIR /home/tp-cn2/TP_Poisson

#create a .mk file to compile the code
#mkdir machine name .mk
#copy the ambre.mk to buildkitsandbox.mk
#RUN echo "CC=gcc" > buildkitsandbox.mk
RUN cat new.mk > buildkitsandbox.mk



#test the environment
#RUN make testenv

#run ls to see the content of the folder
RUN ls -l /home/tp-cn2/TP_Poisson


CMD [ "/bin/bash" ]




