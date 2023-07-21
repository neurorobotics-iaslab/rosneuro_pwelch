# ROS-Neuro pwelch package
The package provides a class and a node which compute the power spectral density (PSD) throught pwelch algorithm. Additionally, in the implementation of the pwelch algorithm with a node, a ringbuffer is asked as parameter.

# Parameters
The parameters needed to configure properly Pwelch are:
- ```psd_wlength``` (```int```, default: ```256```) 
- ```psd_novl``` (```int```, default: ```128```) 
- ```psd_dolog``` (```int```, default: ```1```) 
- ```sampling_freq``` (```int```, default: ```512```) 
- ```wtype``` (```int```, default: ```2```) 

## Usage
The package required as input a dynamicmatrix (which is a Eigen::matrixXd) in order to apply the pwlech algorithm, it must be *[samples x channels]*. Additionally, it returns a matrix with *[nfreqs x channels]*.
