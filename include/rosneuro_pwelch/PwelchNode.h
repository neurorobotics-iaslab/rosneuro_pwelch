#ifndef ROSNEURO_PWELCH_PWELCHNODE_H
#define ROSNEURO_PWELCH_PWELCHNODE_H

#include <ros/ros.h>
#include <rosneuro_msgs/NeuroFrame.h>
#include <rosneuro_pwelch/Pwelch.h>
#include <rosneuro_buffers_ringbuffer/RingBuffer.h>
#include <gtest/gtest_prod.h>

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

            ros::NodeHandle nh_;
            ros::NodeHandle p_nh_;
            ros::Subscriber sub_;
            ros::Publisher  pub_;
            unsigned int psd_wlength_;
            unsigned int psd_novl_;
            unsigned int psd_dolog_;
            unsigned int sampling_freq_;
            unsigned int wtype_;

            Pwelch<T> * pwelch_;
            Buffer<T> * buffer_;
            rosneuro_msgs::NeuroFrame neuromsg_;

            FRIEND_TEST(TestPwelchNodeSuite, Constructor);
            FRIEND_TEST(TestPwelchNodeSuite, Configure);
            FRIEND_TEST(TestPwelchNodeSuite, WrongConfigure);
            FRIEND_TEST(TestPwelchNodeSuite, OnReceivedNeuroData);
    };

    template <typename T>
    PwelchNode<T>::PwelchNode(void) : p_nh_("~") {
        this->sub_ = this->nh_.subscribe("/neurodata_filtered", 100, &PwelchNode<T>::on_received_neurodata, this);
        this->pub_ = this->nh_.advertise<rosneuro_msgs::NeuroFrame>("/neurodata_psd", 1);
    }

    template <typename T>
    PwelchNode<T>::~PwelchNode(void){
    }

    template <typename T>
    bool PwelchNode<T>::configure(void){
        this->p_nh_.param<int>("psd_wlength", (int&) this->psd_wlength_, 256);
        this->p_nh_.param<int>("psd_novl", (int&) this->psd_novl_, 128);
        this->p_nh_.param<int>("psd_dolog", (int&) this->psd_dolog_, 1);
        this->p_nh_.param<int>("sampling_freq", (int&) this->sampling_freq_, 512);
        this->p_nh_.param<int>("wtype", (int&) this->wtype_, 2);

        this->buffer_ = new RingBuffer<double>();
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
        unsigned int nsamples   = msg.eeg.info.nsamples;
        unsigned int nchannels  = msg.eeg.info.nchannels;
        float* ptr_in = const_cast<float*>(msg.eeg.data.data());

        DynamicMatrix<float> data_in = Eigen::Map<rosneuro::DynamicMatrix<float>>(ptr_in, nchannels, nsamples);
        DynamicMatrix<double> data_in_T = data_in.transpose().cast<double>();
        this->buffer_->add(data_in_T);

        if(this->buffer_->isfull()){
            DynamicMatrix<double> buffer_data = this->buffer_->get();
            int nchannels = buffer_data.cols();

            DynamicMatrix<double> data_out_T = this->pwelch_->apply(buffer_data);
            unsigned int nfreqs = data_out_T.rows();

            DynamicMatrix<double> data_out = data_out_T.transpose();
            std::vector<float> neuroeeg(nfreqs * nchannels);
            Eigen::Map<rosneuro::DynamicMatrix<float>>(neuroeeg.data(), nchannels, nfreqs) = data_out.cast<float>();

            this->neuromsg_.header.stamp = ros::Time::now();
            this->neuromsg_.eeg.data = neuroeeg;
            this->neuromsg_.sr  = msg.sr;
            this->neuromsg_.exg = msg.exg;
            this->neuromsg_.tri = msg.tri;
            this->neuromsg_.eeg.info = msg.eeg.info;

            this->pub_.publish(this->neuromsg_);

        }else{
            ROS_INFO("[%s] Buffer is not full", this->buffer_->name().c_str());
        }
        std::cout << "seq added: " << msg.header.seq << std::endl;
    }
}

#endif