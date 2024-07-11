#include "ros/ros.h"
#include "std_msgs/String.h"
#include "rosneuro_pwelch/Pwelch.hpp"
#include <rosneuro_filters/rosneuro_filters_utilities.hpp>
#include "rosneuro_buffers_ringbuffer/RingBuffer.hpp"

int main(int argc, char **argv){
    ros::init(argc, argv, "test_pwelch");
    
    std::string datapath;
	
	rosneuro::Pwelch<double>* pwelch = new rosneuro::Pwelch<double>();

	if(ros::param::get("~datapath", datapath) == false) {
		ROS_ERROR("Cannot find 'datapath' parameter");
		return 0;
	}

	int framesize = 32;
	int samplingFreq = 512;
	int wlength = 256;
	int wtype = 2; // hamming
	int novl = 128;
	int sampleRate = 512;
	int dolog = 1;

	if(pwelch->configure(wlength, wtype, novl, samplingFreq, dolog) == false){
        ROS_ERROR("FFTW plan not created");
		throw std::runtime_error("FFTW plan not created");
    }

	rosneuro::Buffer<double>* buffer = new rosneuro::RingBuffer<double>();
	buffer->configure("RingBufferCfg");
	//buffer->set(samplingFreq, framesize);

    const std::string fileinput    = datapath + "/test/rawdata.csv";
	const std::string fileout      = datapath + "/test/psd_window.csv";   
	const std::string fileFeatures = datapath + "/test/features_rosneuro.csv"; 
    
    // Load input data
	rosneuro::DynamicMatrix<double> input = readCSV<double>(fileinput); // [samples x channels]
	rosneuro::DynamicMatrix<double>  output(39345, framesize);
	//std::vector<double> allFeatures;
	//Eigen::Ref<Eigen::VectorXd> allFeatures;

	int nsamples = input.rows();
	int nchannels = input.cols();
	
	// Allocate frame data (for simulating online loop) 	
	rosneuro::DynamicMatrix<double> framedata = rosneuro::DynamicMatrix<double>::Zero(framesize, nchannels);

	// Apply the psd
	ROS_INFO("Applying psd");
	auto count = 0;
	for(auto i = 0; i<nsamples; i = i+framesize) {

		framedata = input.middleRows(i, framesize); // [samples x channels]

		// filling the ring bar
		buffer->add(framedata); 

		// check if full, if yes go on
		if(buffer->isfull() == false){
			ROS_INFO("[%s] Buffer is not full", buffer->name().c_str());
			continue;
		}

		// matrix of size [frequency x nchannels]
		rosneuro::DynamicMatrix<double> c_output = pwelch->apply(buffer->get());
		output.middleRows(count*129, 129) = c_output;
		count++;
}

	// Writing the filtered data
	// output is [(frequency x nwindows) x nchannels]
	writeCSV<double>(fileout, output);

	//writeCSV<double>(fileFeatures, allFeatures);

	ros::shutdown();
	

	return 0;
}