# OmegaPlus-GPU
Author: Reinout Corts

First release: 13/12/2021

Last update: 13/12/2021

Version: 3.0.0

## About
This is a forked repo of OmegaPlus: a scalable tool for rapid detection of selective sweeps in whole-genome datasets (https://github.com/alachins/omegaplus.)

This forked repo adds GPU-acceleration functionality to the original OmegaPlus tool.

## Prerequisites
To be able to compile the GPU-accelerated source code the following packages should be installed:

* [make](https://www.gnu.org/software/make/)
* OpenCL

Installing OpenCL differs per Linux distribution and GPU manufacturer. The following pages provide some info on the process:


* [NVIDIA](https://developer.nvidia.com/opencl)
* [AMD](https://rocmdocs.amd.com/en/latest/Installation_Guide/Installation_new.html)

OpenCL headers also need to be installed for linking purposes. In most repositories the package is named ```opencl-headers```:

* [Ubuntu](https://packages.ubuntu.com/search?keywords=opencl-headers)
* [Debian](https://packages.debian.org/search?keywords=opencl-headers)
* [Arch](https://www.archlinux.org/packages/extra/any/opencl-headers)
* [Fedora](https://fedora.pkgs.org/30/fedora-armhfp/opencl-headers-2.2-4.20190205git49f07d3.fc30.noarch.rpm.html)
* [Centos 7](https://centos.pkgs.org/7/epel-x86_64/opencl-headers-2.2-1.20180306gite986688.el7.noarch.rpm.html)

OpenCL can be tested using the ```clinfo``` command which can be installed using the package manager of your distribution.

## Download and Compile
The following commands can be used to download and compile the source code of the GPU accelerated OmegaPlus version (note that ```git``` is needed for this).

```
$ cd ~/Documents/
$ git clone https://github.com/MrKzn/omegaplus.git
$ cd omegaplus
$ make -f Makefile.GPU.gcc
```

The executable named ```OmegaPlus-GPU``` is placed in the main folder

## Test Run
To verify the correct installation of the prerequisites and GPU-accelerated OmegaPlus itself, a test run can be done with the following command:

```
$ ./OmegaPlus-GPU -name TEST -grid 1000 -length 100000 -minwin 1000 -maxwin 20000 -seed 1332 -binary -input examples/ms.out
```

This command is going to execute the GPU-accelerated version of OmegaPlus that computes Omega at 1,000 equidistant grid positions with the left and right sub-region windows set to a maximum of 20,000 SNPs and minimum of 1,000 SNPs and the alignment length set to 100,000.

## GPU version exeptions
The GPU-accelerated OmegaPlus version does not support DNA data, only binary data. This can be forced using the ```-binary``` command line option or just by using binary data.

The ```-b``` ```borderTol``` approximation command line option is not implemented in the GPU-accelerated OmegaPlus version.

