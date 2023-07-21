#ifndef ROSNEURO_PWELCH_PWELCHNODE_H
#define ROSNEURO_PWELCH_PWELCHNODE_H

#include <ros/ros.h>
#include <rosneuro_msgs/NeuroFrame.h>
#include <rosneuro_pwelch/Pwelch.h>
#include <rosneuro_buffers_ringbuffer/RingBuffer.h>

namespace rosneuro{

template <typename T>
class PwelchNode{
    public:
        PwelchNode(void);
        ~PwelchNode(void);

        bool configure(void);
        void run(void);

	private:
		void on_received_neurodata(const rosneuro_msgs::NeuroFrame& msg);

	private:
		ros::NodeHandle nh_;
		ros::NodeHandle p_nh_;
		ros::Subscriber sub_;
		ros::Publisher  pub_;
        unsigned int                    psd_wlength_;
        unsigned int                    psd_novl_;
	    unsigned int                    psd_dolog_; 
        unsigned int                    sampling_freq_;
        unsigned int                    wtype_; 

		Pwelch<T> * pwelch_;
		Buffer<T> * buffer_;
		rosneuro_msgs::NeuroFrame neuromsg_;
};

template <typename T>
PwelchNode<T>::PwelchNode(void) : p_nh_("~") {
    this->sub_ = this->nh_.subscribe("/neurodata_filtered", 100, &PwelchNode<T>::on_received_neurodata, this);
    this->pub_ = this->nh_.advertise<rosneuro_msgs::NeuroFrame>("/neurodata_psd", 1); // may be is better with p_nh_
}
template <typename T>
PwelchNode<T>::~PwelchNode(void){}

template <typename T>
bool PwelchNode<T>::configure(void){
    this->p_nh_.param<int>("psd_wlength", (int&) this->psd_wlength_, 256);
	this->p_nh_.param<int>("psd_novl", (int&) this->psd_novl_, 128);
	this->p_nh_.param<int>("psd_dolog", (int&) this->psd_dolog_, 1); 
    this->p_nh_.param<int>("sampling_freq", (int&) this->sampling_freq_, 512);
    this->p_nh_.param<int>("wtype", (int&) this->wtype_, 2); // 2 is hamming

	this->buffer_ = new rosneuro::RingBuffer<double>();
	if(buffer_->configure("RingBufferCfg") == false){ 
		ROS_ERROR("[%s] RingBuffer configuration failed", this->buffer_->name().c_str());
        return false;
	}
	ROS_INFO("[RingBuffer] RingBuffer has been configured");

	this->pwelch_ = new rosneuro::Pwelch<double>();
    if(this->pwelch_->configure(this->psd_wlength_, this->wtype_, this->psd_novl_, this->sampling_freq_, this->psd_dolog_) == false){
        ROS_ERROR("[Pwelch] Pwelch configuration failed");
        return false;
    }

    ROS_INFO("[Pwelch] Pwelch has been configured");
    return true;
}   

template <typename T>
void PwelchNode<T>::run(void){

	ros::Rate r(512);

	while(ros::ok()) {
		ros::spinOnce();
		r.sleep();
	}
}

template <typename T>
void PwelchNode<T>::on_received_neurodata(const rosneuro_msgs::NeuroFrame& msg) {

	// Getting information from the message
	float* ptr_in;
	rosneuro::DynamicMatrix<float> data_in;
	rosneuro::DynamicMatrix<double> data_in_T;
	rosneuro::DynamicMatrix<double> data_out_T;
	rosneuro::DynamicMatrix<double> data_out;

	unsigned int nsamples   = msg.eeg.info.nsamples;
	unsigned int nchannels  = msg.eeg.info.nchannels;

	// Getting pointer to the input data message
	ptr_in = const_cast<float*>(msg.eeg.data.data());

	// Re-map raw pointer to eigen matrix.
	// Remember input data is stored as [channels x samples]
	data_in = Eigen::Map<rosneuro::DynamicMatrix<float>>(ptr_in, nchannels, nsamples);

    // Transpose input data to [samples x channels]
	data_in_T = data_in.transpose().cast<double>();
	
	// Add data to the buffer
	this->buffer_->add(data_in_T);

	// Perform pwelch only if buffer full
	if(this->buffer_->isfull()){
		
		// Take data from the buffer [nsamples x nchannels]
		rosneuro::DynamicMatrix<double> buffer_data = this->buffer_->get();
		int nchannels = buffer_data.cols();

		// Apply pwelch [frequency x nchannels]
    	data_out_T = this->pwelch_->apply(buffer_data);
		unsigned int nfreqs = data_out_T.rows();

		// Transpose again data out to [nchannels x frequency]
		data_out = data_out_T.transpose();
       
	   	// Re-map eigen matrix to std::vector
	   	std::vector<float> neuroeeg(nfreqs * nchannels);
	   	Eigen::Map<rosneuro::DynamicMatrix<float>>(neuroeeg.data(), nchannels, nfreqs) = data_out.cast<float>();
	   	
	   	// Setting message
	   	this->neuromsg_.header.stamp = ros::Time::now();
		this->neuromsg_.eeg.data = neuroeeg;
		this->neuromsg_.sr  = msg.sr;
		this->neuromsg_.exg = msg.exg;
		this->neuromsg_.tri = msg.tri;
		this->neuromsg_.eeg.info = msg.eeg.info;
    
	   	// Publishing the message
	   	this->pub_.publish(this->neuromsg_);

	}else{
	   	ROS_INFO("[%s] Buffer is not full", this->buffer_->name().c_str());
    }
	std::cout << "seq added: " << msg.header.seq << std::endl;
}

}

#endif