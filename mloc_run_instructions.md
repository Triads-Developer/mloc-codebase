# FortranOptimization
Optimization of code which performs multiple event relocation on a cluster of earthquakes, using the Hypocentroidal Decomposition (HD) algorithm


# How to run

Connect to the VPN and SSH into the server:


    ssh triads@mantle.wustl.edu

CD into the current working directory


    cd /RAID/users/triads/mloc_working/test1

Run the program providing a *.in file as input

    /RAID/users/triads/mloc_working/mloc_f < run01.in


note, the code will produce a bunch of files prepended with run01. If you want to run it again with that same .in file, you'll have to move the run01* files into another folder.

# How to run on the RIS

Connect to the VPN

Copy the project folder into your Shared Folder (or Home folder)

ssh into an environment

    ssh <wustl username>@compute1-client-1.ris.wustl.edu

Mount your shared volume so the jobs can see it (and the code) - change g.porter to whatever your username is

    export LSF_DOCKER_VOLUMES='/storage1/fs1/artsci/Active/g.porter:/storage1/fs1/artsci/Active/g.porter'

At this point, you have some flexibility. In this case, we are going to be using an ```interactive``` job with a image focused on ```R```. It doesn't matter because we just want an Ubuntu-based image. Ubuntu will come with GCC and GCC has Gfortran (which is what we'll use to compile the code). Again, substitude your username for g.porter for your home directory. This also assumes that you are part of A&S which gives the ability us to use the group ```compute-artsci``` and the queue ```artsci-interactive```

    PATH=/home/g.porter:$PATH LSF_DOCKER_PORTS='8081:8787' bsub -Is -G compute-artsci -q artsci-interactive -M 25G -R 'select[port8081=1]' -a 'docker(rocker/verse:4.0.2)' /bin/bash

Navigate to the previously mounted directory    

    cd /storage1/fs1/artsci/Active/g.porter

If you want to compile with gfortran, the code will be in ```mloc_par_test/mloc_gfortran```

Make a backup folder, and move all the ```*.o``` files and ```mloc_g``` into it. This is safer than just removing stuff

    mkdir backup
    mv *.o backup/
    mv mloc_g backup/

Next, build it by running ```make```
