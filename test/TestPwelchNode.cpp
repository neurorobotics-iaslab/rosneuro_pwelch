#include "rosneuro_pwelch/PwelchNode.h"
#include <gtest/gtest.h>

namespace rosneuro {

class TestPwelchNodeSuite : public ::testing::Test {
    public:
        TestPwelchNodeSuite() {}
        ~TestPwelchNodeSuite() {}

        void SetUp() override {
            pwelch_node = new PwelchNode<double>();
        }

        void TearDown() override {
            ros::param::del("RingBufferCfg");
            delete pwelch_node;
        }

        void config(){
            ros::NodeHandle nh;
            XmlRpc::XmlRpcValue params, validConfig;
            validConfig["name"] = "ringbuffer";
            validConfig["type"] = "RingBufferFloat";
            params["size"] = 256;
            validConfig["params"] = params;
            nh.setParam("RingBufferCfg", validConfig);
        }

        PwelchNode<double>* pwelch_node;
};

TEST_F(TestPwelchNodeSuite, Constructor) {
    ASSERT_TRUE(pwelch_node->sub_.getTopic() == "/neurodata_filtered");
    ASSERT_TRUE(pwelch_node->pub_.getTopic() == "/neurodata_psd");
}

TEST_F(TestPwelchNodeSuite, Configure) {
    config();
    ASSERT_TRUE(pwelch_node->configure());
    ASSERT_TRUE(pwelch_node->psd_wlength_ == 256);
    ASSERT_TRUE(pwelch_node->psd_novl_ == 128);
    ASSERT_TRUE(pwelch_node->psd_dolog_ == 1);
    ASSERT_TRUE(pwelch_node->sampling_freq_ == 512);
    ASSERT_TRUE(pwelch_node->wtype_ == 2);
}

TEST_F(TestPwelchNodeSuite, WrongConfigure) {
    ASSERT_FALSE(pwelch_node->configure());
}

TEST_F(TestPwelchNodeSuite, OnReceivedNeuroData) {
    config();
    pwelch_node->configure();

    rosneuro_msgs::NeuroFrame msg;
    msg.header.stamp = ros::Time::now();
    msg.header.frame_id = "neurodata";
    msg.eeg.info.nsamples = 256;
    msg.eeg.info.nchannels = 256;
    msg.eeg.data = std::vector<float>(256*256, 2);
    pwelch_node->on_received_neurodata(msg);

    ASSERT_EQ(pwelch_node->neuromsg_.eeg.info.nsamples, 256);
    ASSERT_EQ(pwelch_node->neuromsg_.eeg.info.nchannels, 256);
    ASSERT_EQ(pwelch_node->neuromsg_.eeg.data.size(), 33024);
    ASSERT_EQ(pwelch_node->neuromsg_.exg, msg.exg);
    ASSERT_EQ(pwelch_node->neuromsg_.tri, msg.tri);
}

}

int main(int argc, char **argv) {
    ros::console::set_logger_level(ROSCONSOLE_DEFAULT_NAME, ros::console::levels::Fatal);
    ros::init(argc, argv, "test_pwelch_node");
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}